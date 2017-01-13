from scipy.special import ellipk
import math
from tabulate import tabulate
import numpy as np
import cmath


class CPWResonator:
    """Coplanar wave-guide resonator class"""

    def __init__(self, length, conductorWidth, gapWidth, conductorThickness, resonatorType,
                 conductorMaterial='Niobium', substrateMaterial='Silicon',
                 temperature=4, couplingCapacitance=0E0, loadImpedance=50, loadBoundaryCondition='Short', modes=[1]):
        # Supplied parameters
        self.length = length
        self.conductorWidth = conductorWidth
        self.gapWidth = gapWidth
        self.conductorThickness = conductorThickness
        self.resonatorType = resonatorType
        self.conductor = self.conductorProperties(conductorMaterial)
        self.substrate = self.substrateProperties(substrateMaterial)
        self.temperature = temperature
        self.couplingCapacitance = couplingCapacitance
        self.loadImpedance = loadImpedance
        self.loadBoundaryCondition = loadBoundaryCondition
        self.modes = modes

        # Calculated parameters
        self.effectivePermittivity = self.effectivePermittivity(self.substrate.relativePermittivity)
        self.capacitancePerUnitLength = self.capacitancePerUnitLength(self.conductorWidth, self.gapWidth, self.effectivePermittivity)
        self.geometricInducatancePerUnitLength = self.geometricInducatancePerUnitLength(self.conductorWidth, self.gapWidth)
        self.kineticInductancePerUnitLength = self.kineticInductancePerUnitLength(self.conductorWidth, self.gapWidth, self.conductorThickness,
            self.temperature, self.conductor.criticalTemperature, self.conductor.londonPenetrationDepthZero)
        self.totalInductancePerUnitLength = self.geometricInducatancePerUnitLength + self.kineticInductancePerUnitLength
        self.characteristicImpedance = self.characteristicImpedance(self.totalInductancePerUnitLength, self.capacitancePerUnitLength)
        self.uncoupledResonantFrequency = self.uncoupledResonantFrequency(self.totalInductancePerUnitLength, self.capacitancePerUnitLength, self.length,
            self.resonatorType, self.modes)
        self.coupledResonantFrequency = self.coupledResonantFrequency(self.totalInductancePerUnitLength, self.capacitancePerUnitLength, self.length,
            self.resonatorType, self.couplingCapacitance, self.loadImpedance, self.modes)
        self.inputImpedance = self.inputImpedance(self.length, self.characteristicImpedance, self.loadBoundaryCondition,
            self.totalInductancePerUnitLength, self.capacitancePerUnitLength, self.uncoupledResonantFrequency)
        self.internalQualityFactor = self.internalQualityFactor(self.length, self.resonatorType, self.modes, self.conductor)
        self.externalQualityFactor = self.externalQualityFactor(self.resonatorType, self.modes, self.uncoupledResonantFrequency,
            self.couplingCapacitance, self.characteristicImpedance)
        self.externalQualityFactor2 = self.externalQualityFactor2(self.resonatorType, self.modes, self.uncoupledResonantFrequency,
            self.capacitancePerUnitLength, self.length, self.couplingCapacitance)
        self.loadedQualityFactor = self.loadedQualityFactor(self.internalQualityFactor, self.externalQualityFactor)

    def effectivePermittivity(self, relativePermittivity):
        return (1 + relativePermittivity)/2

    def capacitancePerUnitLength(self, conductorWidth, gapWidth, effectivePermittivity):
        # Permittivity of free space
        freeSpacePermittivity = 8.85418782E-12

        # Complete elliptic integral of the first kind
        k = conductorWidth / (conductorWidth + 2 * gapWidth)
        k2 = math.sqrt(1 - math.pow(k,2))

        # Total CPW capacitance p.u.l.
        return 4 * freeSpacePermittivity * (effectivePermittivity + 0) * (ellipk(k) / ellipk(k2))

    def geometricInducatancePerUnitLength(self, conductorWidth, gapWidth):
        # Permeability of freespace
        freeSpacePermeability = 1.25663706E-6

        # Complete elliptic integral of the first kind
        k = conductorWidth / (conductorWidth + 2 * gapWidth)
        k2 = math.sqrt(1 - math.pow(k,2))

        # Total conductor geometric inductance p.u.l.
        return (freeSpacePermeability / 4) * (ellipk(k2) / ellipk(k))

    def kineticInductancePerUnitLength(self, conductorWidth, gapWidth, conductorThickness,
                          temperature, criticalTemperature, londonPenetrationDepthZero):
        # Permeability of freespace
        freeSpacePermeability = 1.25663706E-6

        # Complete elliptic integral of the first kind
        k = conductorWidth / (conductorWidth + 2 * gapWidth)
        K = ellipk(k)

        # Penetration depth at temperature T
        londonPenetrationDepthT = self.londonPenetrationDepthT(temperature, criticalTemperature, londonPenetrationDepthZero)

        # Geometrical factor
        geometricFactor = (1 / (2 * math.pow(k, 2) * math.pow(K, 2))) * (
            - math.log(conductorThickness / (4 * conductorWidth)) + ((2 * (conductorWidth + gapWidth))
            / (conductorWidth + 2 * gapWidth)) * math.log(gapWidth / (conductorWidth + gapWidth)) - (
            conductorWidth / (conductorWidth + 2 * gapWidth)) * math.log(conductorThickness
            / (4 * (conductorWidth + 2 * gapWidth))))

        # Kinetic Inductance p.u.l.
        return freeSpacePermeability * (math.pow(londonPenetrationDepthT, 2) / (conductorWidth * conductorThickness)) * geometricFactor


    def londonPenetrationDepthT(self, temperature, criticalTemperature, londonPenetrationDepthZero):
        return londonPenetrationDepthZero / math.sqrt(1 - math.pow((temperature / criticalTemperature), 4))

    def characteristicImpedance(self, inductance, capacitance, resistance=0, conductance=0, frequency=1):
        return cmath.sqrt(
            (resistance + 1j*2*math.pi*frequency*inductance ) /
            (conductance + 1j*2*math.pi*frequency*capacitance))

    def inputImpedance(self, length, characteristicImpedance, loadBoundaryCondition, inductancePerUnitLength, capacitancePerUnitLength, frequency,
                       resistancePerUnitLength=0, conductancePerUnitLength=0):
        gamma = np.sqrt(
            (resistancePerUnitLength + 1j*2*math.pi*frequency*inductancePerUnitLength ) *
            (conductancePerUnitLength + 1j*2*math.pi*frequency*capacitancePerUnitLength))

        if loadBoundaryCondition == 'Short':
            return characteristicImpedance * np.tanh(gamma*length)
        elif loadBoundaryCondition == 'Open':
            return characteristicImpedance / np.tanh(gamma * length)
        else:
            print('Error: Load boundary condition no valid!')
            return -1

    def uncoupledResonantFrequency(self, totalInductancePerUnitLength, capacitancePerUnitLength, length, resonatorType, modes):
        m = self.getModeFactor(resonatorType, modes)

        return 1 / (math.sqrt(totalInductancePerUnitLength*capacitancePerUnitLength) * m * length)

    def coupledResonantFrequency(self, totalInductancePerUnitLength, capacitancePerUnitLength, length, resonatorType, couplingCapacitance, loadImpedane, modes):
        m = self.getModeFactor(resonatorType, modes)

        # Pre-coupled
        uncoupledResonantFrequency = 1 / (math.sqrt(totalInductancePerUnitLength*capacitancePerUnitLength) * m * length)

        # Post-coupled
        effectiveCouplingCapacitance = self.effectiveCouplingCapacitance(couplingCapacitance, uncoupledResonantFrequency, loadImpedane)

        return 1 / (np.sqrt((totalInductancePerUnitLength*length) * ((capacitancePerUnitLength*length) + effectiveCouplingCapacitance)) * m)

    def effectiveCouplingCapacitance(self, couplingCapacitance, frequency, loadImpedane):
        return couplingCapacitance / \
               (1 + np.power(frequency * couplingCapacitance * loadImpedane * math.pi, 2))

    def internalQualityFactor(self, length, resonatorType, modes, conductor):
        m = self.getModeFactor(resonatorType, modes)
        return (1/m) * (math.pi/(conductor.amplitudeAttenuationPerUnitLength*length))

    def externalQualityFactor(self, resonatorType, modes, uncoupledResonantFrequency, couplingCapacitance, characteristicImpedance):
        m = self.getModeFactor(resonatorType, modes)
        q_in = 2 * math.pi * uncoupledResonantFrequency * couplingCapacitance * characteristicImpedance
        return (1/m) * (math.pi/(q_in**2))

    def externalQualityFactor2(self, resonatorType, modes, uncoupledResonantFrequency, capacitancePerUnitLength, length, couplingCapacitance, loadResistance=50):
        omega_n = 2 * math.pi * uncoupledResonantFrequency
        r_star = (1+(omega_n*couplingCapacitance*loadResistance)**2) / ((omega_n*couplingCapacitance)**2 * loadResistance)
        C = (capacitancePerUnitLength * length)/2
        return  omega_n * r_star * C


    def loadedQualityFactor(self, internalQualityFactor, externalQualityFactor):
        return 1/((1/internalQualityFactor) + (1/externalQualityFactor))

    def getModeFactor(self, resonatorType, modes):
        if resonatorType == 'half':
            m = 4.0 / (2.0 * np.array(modes))
        elif resonatorType == 'quarter':
            m = 4.0 / ((2.0 * np.array(modes)) - 1)
        else:
            print('Error: Incorrect resonator type provided!')
            return -1
        return m

    def conductorProperties(self, material):
        return{
            'Niobium': Conductor(material=material, criticalTemperature=9.2, londonPenetrationDepthZero=33.3E-9,
                resistancePerUnitLength=0, conductancePerUnitLength=0, amplitudeAttenuationPerUnitLength=0.01),
            'Niobium Nitride': Conductor(material=material, criticalTemperature=16.2, londonPenetrationDepthZero=40E-9,
                resistancePerUnitLength=0, conductancePerUnitLength=0, amplitudeAttenuationPerUnitLength=0.01),
        }[material]

    def substrateProperties(self, material):
        return{
            'Silicon': Substrate(material=material, relativePermittivity=11.9),
            'Sapphire': Substrate(material=material, relativePermittivity=10.2),
        }[material]


class Conductor:
    """CPW conductor material class"""

    def __init__(self, material, criticalTemperature, londonPenetrationDepthZero, resistancePerUnitLength,
                 conductancePerUnitLength, amplitudeAttenuationPerUnitLength):
        self.material = material
        self.criticalTemperature = criticalTemperature
        self.londonPenetrationDepthZero = londonPenetrationDepthZero
        self.resistancePerUnitLength = resistancePerUnitLength
        self.conductancePerUnitLength = conductancePerUnitLength
        self.amplitudeAttenuationPerUnitLength = amplitudeAttenuationPerUnitLength


class Substrate:
    """Substrate material class"""

    def __init__(self, material, relativePermittivity):
        self.material = material
        self.relativePermittivity = relativePermittivity

if __name__ == '__main__':
    mCPW = CPWResonator(length=7186E-6, conductorWidth=20E-6, gapWidth=10E-6, conductorThickness=100E-9,
                    resonatorType='quarter', conductorMaterial='Niobium Nitride', substrateMaterial='Silicon',
                    temperature=4, couplingCapacitance=1E-15, loadBoundaryCondition='Short', modes=[1, 2, 3])
print(tabulate(
    [['Effective permittivity', mCPW.effectivePermittivity, ''],
     ['Substrate capacitance', mCPW.capacitancePerUnitLength * math.pow(10, 12), 'pF/m'],
     ['Geometric Inductance', mCPW.geometricInducatancePerUnitLength, 'H/m'],
     ['Kinetic Inductance', mCPW.kineticInductancePerUnitLength, 'H/m'],
     ['Characteristic Impedance', mCPW.characteristicImpedance, 'Ohms'],
     ['Input Impedance', mCPW.inputImpedance, 'Ohms'],
     ['Resonant frequency (Uncoupled)', mCPW.uncoupledResonantFrequency/math.pow(10,9), 'Ghz'],
     ['Resonant frequency (Coupled)', mCPW.coupledResonantFrequency / math.pow(10, 9), 'Ghz'],
     ['Internal Quality factor', mCPW.internalQualityFactor, ''],
     ['External Quality factor', mCPW.externalQualityFactor, ''],
     ['Loaded Quality factor', mCPW.loadedQualityFactor, '']],
    headers=['Property', 'Value', 'Units'], floatfmt=".2f"))