from scipy.special import ellipk
import math
import cmath
from tabulate import tabulate

class CPWResonator:
    'Coplanar wave-guide resonator class'

    def __init__(self, length, conductorWidth, gapWidth, conductorThickness, resonatorType,
                 conductorMaterial='Niobium', substrateMaterial='Silicon',
                 temperature=4, couplingCapacitance=0E0, loadImpedance=50, loadBoundaryCondition='Short'):
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

        # Calculated parameters
        self.effectivePermittivity = self.effectivePermittivity(self.substrate.relativePermittivity)
        self.substrateCapacitance = self.substrateCapacitance(self.conductorWidth, self.gapWidth, self.effectivePermittivity)
        self.geometricInducatance = self.geometricInducatance(self.conductorWidth, self.gapWidth)
        self.kineticInductance = self.kineticInductance(self.conductorWidth, self.gapWidth, self.conductorThickness,
            self.temperature, self.conductor.criticalTemperature, self.conductor.londonPenetrationDepthZero)
        self.totalInductance = self.geometricInducatance + self.kineticInductance
        self.characteristicImpedance = self.characteristicImpedance(self.totalInductance, self.substrateCapacitance)
        self.resonantFrequency = self.resonantFrequency(self.totalInductance, self.substrateCapacitance, self.length,
            self.resonatorType)
        self.coupledResonantFrequency = self.coupledResonantFrequency(self.totalInductance, self.substrateCapacitance, self.length,
            self.resonatorType, self.couplingCapacitance, self.loadImpedance)
        self.inputImpedance = self.inputImpedance(self.length, self.characteristicImpedance, self.loadBoundaryCondition,
            self.totalInductance, self.substrateCapacitance, self.resonantFrequency)

    def effectivePermittivity(self, relativePermittivity):
        # Filling factor, q
        # Simplification for large heights, H, H1
        q = 0.5

        return 1 + q * (relativePermittivity - 1)

    def substrateCapacitance(self, conductorWidth, gapWidth, effectivePermittivity):
        # Permittivity of free space
        freeSpacePermittivity = 8.85418782E-12

        # Complete elliptic integral of the first kind
        k = conductorWidth / (conductorWidth + 2 * gapWidth)
        k2 = math.sqrt(1 - math.pow(k,2))

        # Total substrate capacitance
        return 4 * freeSpacePermittivity * effectivePermittivity * ellipk(k) / ellipk(k2)

    def geometricInducatance(self, conductorWidth, gapWidth):
        # Permeability of freespace
        freeSpacePermeability = 1.25663706E-6

        # Complete elliptic integral of the first kind
        k = conductorWidth / (conductorWidth + 2 * gapWidth)
        k2 = math.sqrt(1 - math.pow(k,2))

        # Total conductor geometric inductance
        return (freeSpacePermeability / 4) * ellipk(k2) / ellipk(k)

    def kineticInductance(self, conductorWidth, gapWidth, conductorThickness,
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

        # Kinetic Inductance
        return freeSpacePermeability * (math.pow(londonPenetrationDepthT, 2) / (conductorWidth * conductorThickness)) \
               * geometricFactor


    def londonPenetrationDepthT(self, temperature, criticalTemperature, londonPenetrationDepthZero):
        return londonPenetrationDepthZero / math.sqrt(1 - math.pow((temperature / criticalTemperature), 4) )

    def characteristicImpedance(self, inductance, capacitance, resistance=0, conductance=0, frequency=1):
        return cmath.sqrt(
            (resistance + 1j*2*math.pi*frequency*inductance ) /
            (conductance + 1j*2*math.pi*frequency*capacitance))

    def inputImpedance(self, length, characteristicImpedance, loadBoundaryCondition, inductance, capacitance, frequency,
                       resistance=0, conductance=0):
        gamma = cmath.sqrt(
            (resistance + 1j*2*math.pi*frequency*inductance ) *
            (conductance + 1j*2*math.pi*frequency*capacitance))

        if loadBoundaryCondition == 'Short':
            return characteristicImpedance * cmath.tanh(gamma*length)
        elif loadBoundaryCondition == 'Open':
            return characteristicImpedance / cmath.tanh(gamma * length)
        else:
            print('Error: Load boundary condition no valid!')
            return -1

    def resonantFrequency(self, inductance, capacitance, length, resonatorType):

        if resonatorType == 'half':
            resonatorTypeFactor = 2
        elif resonatorType == 'quarter':
            resonatorTypeFactor = 4
        else:
            print('Error: Incorrect resonator type provided!')
            return -1

        return  1 / (math.sqrt(inductance*capacitance) * resonatorTypeFactor * length)

    def coupledResonantFrequency(self, inductance, capacitance, length, resonatorType, couplingCapacitance, loadImpedane, mode=1):

        if resonatorType == 'half':
            resonatorTypeFactor = 2
        elif resonatorType == 'quarter':
            resonatorTypeFactor = 4
        else:
            print('Error: Incorrect resonator type provided!')
            return -1

        # Pre-coupled
        resonantFrequency = 1 / (math.sqrt(inductance * capacitance) * resonatorTypeFactor * length)

        # Post-coupled
        inductance_n = resonatorTypeFactor * length * inductance / (math.pow(math.pi, 2) * math.pow(mode, 2))
        capacitance_n = resonatorTypeFactor * capacitance * length / 4

        effectiveCouplingCapacitance = self.nortonCapacitance(couplingCapacitance, resonantFrequency, loadImpedane)

        return 1 / (math.sqrt(inductance_n * (capacitance_n + effectiveCouplingCapacitance)) * 2*math.pi)

    def nortonCapacitance(self, couplingCapacitance, resonantFrequency, loadImpedane):
        return couplingCapacitance / \
               (1 + math.pow(resonantFrequency, 2) * math.pow(couplingCapacitance, 2) * math.pow(loadImpedane, 2))

    def conductorProperties(self, material):
        return{
            'Niobium': Conductor(material=material, criticalTemperature=9.2, londonPenetrationDepthZero=33.3E-9,
                resistance=0, conductance=0),
            'Niobium Nitride': Conductor(material=material, criticalTemperature=16.2, londonPenetrationDepthZero=40E-9,
                resistance=0, conductance=0),
            'Niobium Titanium Nitride': Conductor(material=material, criticalTemperature=16.2, londonPenetrationDepthZero=40E-9,
                resistance=0, conductance=0),
            'Aluminium': Conductor(material=material, criticalTemperature=1.2, londonPenetrationDepthZero=15.4E-9,
                                   resistance=0, conductance=0),
        }[material]

    def substrateProperties(self, material):
        return{
            'Silicon': Substrate(material=material, relativePermittivity=11.9),
            'Silicon Oxide': Substrate(material=material, relativePermittivity=4.9),
            'Sapphire': Substrate(material=material, relativePermittivity=10.2),
        }[material]


class Conductor:
    'CPW conductor material class'

    def __init__(self, material, criticalTemperature, londonPenetrationDepthZero, resistance, conductance):
        self.material = material
        self.criticalTemperature = criticalTemperature
        self.londonPenetrationDepthZero = londonPenetrationDepthZero
        self.resistance = resistance
        self.conductance = conductance


class Substrate:
    'Substrate material class'

    def __init__(self, material, relativePermittivity):
        self.material = material
        self.relativePermittivity = relativePermittivity


if __name__ == '__main__':
    mCPW = CPWResonator(length=7524E-6, conductorWidth=9.06E-6, gapWidth=5E-6, conductorThickness=160E-9,
                        resonatorType='quarter', conductorMaterial='Niobium', substrateMaterial='Silicon',
                        temperature=4, couplingCapacitance=50E-15, loadBoundaryCondition='Short')
    print(tabulate(
        [['Effective permittivity', mCPW.effectivePermittivity, ''],
         ['Substrate capacitance', mCPW.substrateCapacitance * math.pow(10, 12), 'pF/m'],
         ['Geometric Inductance', mCPW.geometricInducatance, 'H/m'],
         ['Kinetic Inductance', mCPW.kineticInductance, 'H/m'],
         ['Characteristic Impedance', mCPW.characteristicImpedance, 'Ohms'],
         ['Input Impedance', mCPW.inputImpedance, 'Ohms'],
         ['Resonant frequency (Uncoupled)', mCPW.resonantFrequency/math.pow(10,9), 'Ghz'],
         ['Resonant frequency (Coupled)', mCPW.coupledResonantFrequency / math.pow(10, 9), 'Ghz']],
        headers=['Property', 'Value', 'Units']))
