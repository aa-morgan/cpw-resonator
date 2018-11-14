from scipy.special import ellipk
from scipy.constants import epsilon_0, mu_0, pi
import numpy as np
import pandas as pd
from .materials import *

class CPWResonator:
    """
    Coplanar wave-guide resonator class
    """

    def __init__(self, length, conductorWidth, gapWidth, conductorThickness, resonatorType, 
                 conductorMaterial, substrateMaterial,
                 temperature=4, couplingCapacitance=0, loadImpedance=50, loadBoundaryCondition='Short', mode=1):
        # Supplied parameters
        self.length = np.array(length)
        self.conductorWidth = np.array(conductorWidth)
        self.gapWidth = np.array(gapWidth)
        self.conductorThickness = np.array(conductorThickness)
        self.resonatorType = np.array(resonatorType)
        self.conductor = Conductor(conductorMaterial)
        self.substrate = Substrate(substrateMaterial)
        self.temperature = np.array(temperature)
        self.couplingCapacitance = np.array(couplingCapacitance)
        self.loadImpedance = np.array(loadImpedance)
        self.loadBoundaryCondition = np.array(loadBoundaryCondition)
        self.mode = np.array(mode)

    def effectivePermittivity(self):
        return (1 + self.substrate.relativePermittivity)/2

    def capacitancePerUnitLength(self):
        # Complete elliptic integral of the first kind
        k = self.conductorWidth / (self.conductorWidth + 2 * self.gapWidth)
        k2 = np.sqrt(1 - k**2)

        # Total CPW capacitance p.u.l.
        return 4 * epsilon_0 * (self.effectivePermittivity() + 0) * (ellipk(k) / ellipk(k2))
    
    def totalInductancePerUnitLength(self):
        if self.conductor.superconductor:
            return self.geometricInductancePerUnitLength() + self.kineticInductancePerUnitLength()
        else:
            return self.geometricInductancePerUnitLength()

    def geometricInductancePerUnitLength(self):
        # Complete elliptic integral of the first kind
        k = self.conductorWidth / (self.conductorWidth + 2 * self.gapWidth)
        k2 = np.sqrt(1 - k**2)

        # Total conductor geometric inductance p.u.l.
        return (mu_0 / 4) * (ellipk(k2) / ellipk(k))

    def kineticInductancePerUnitLength(self):
        # Complete elliptic integral of the first kind
        k = self.conductorWidth / (self.conductorWidth + 2 * self.gapWidth)
        K = ellipk(k)

        # Geometrical factor
        s, W, T = self.gapWidth, self.conductorWidth, self.conductorThickness
        geometricFactor = (1 / (2 * k**2 * K**2)) * (- np.log(T / (4 * W)) + ((2 * (W + s))
            / (W + 2 * s)) * np.log(s / (W + s)) - (W / (W + 2 * s)) * np.log(T / (4 * (W + 2 * s))))

        # Kinetic Inductance p.u.l.
        return mu_0 * (self.londonPenetrationDepthT()**2 / (W * T)) * geometricFactor


    def londonPenetrationDepthT(self):
        return self.conductor.londonPenetrationDepthZero /                np.sqrt(1 - (self.temperature / self.conductor.criticalTemperature)**4)

    def characteristicImpedance(self, resistance=0, conductance=0, frequency=1):
        return np.sqrt(
            (resistance  + 1j*2 * pi * frequency * self.totalInductancePerUnitLength() ) /
            (conductance + 1j*2 * pi * frequency * self.capacitancePerUnitLength()))

    def inputImpedance(self):
        gamma = np.sqrt(
            (self.conductor.resistancePerUnitLength + 
             1j*2*pi*self.coupledResonantFrequency()*self.totalInductancePerUnitLength() ) *
            (self.conductor.conductancePerUnitLength + 
             1j*2*pi*self.coupledResonantFrequency()*self.capacitancePerUnitLength()))

        if self.loadBoundaryCondition == 'Short':
            return self.characteristicImpedance() * np.tanh(gamma * self.length)
        elif self.loadBoundaryCondition == 'Open':
            return self.characteristicImpedance() / np.tanh(gamma * self.length)
        else:
            print('Error: Load boundary condition no valid!')
            return -1

    def uncoupledResonantFrequency(self):
        m = self.getModeFactor()
        return 1 / (np.sqrt(self.totalInductancePerUnitLength()*self.capacitancePerUnitLength()) * m * self.length)

    def coupledResonantFrequency(self):
        m = self.getModeFactor()        
        return 1 / (np.sqrt((self.totalInductancePerUnitLength() * self.length) * (
            (self.capacitancePerUnitLength() * self.length) + self.effectiveCouplingCapacitance())) * m)

    def effectiveCouplingCapacitance(self):
        return self.couplingCapacitance / (1 + 
               (self.uncoupledResonantFrequency() * self.couplingCapacitance * self.loadImpedance * pi)**2)

    def internalQualityFactor(self):
        m = self.getModeFactor()
        return (1/m) * (pi/(self.conductor.amplitudeAttenuationPerUnitLength * self.length))
    
    def externalQualityFactor(self, method=0):
        if method == 0:
            return self.externalQualityFactorMain()
        elif method == 1:
            return self.externalQualityFactorApprox()
        elif method == 2:
            return self.externalQualityFactorQWref()
    
    def externalQualityFactorMain(self, loadResistance=50):
        omega_n = 2 * pi * self.uncoupledResonantFrequency()
        r_star = (1+(omega_n*self.couplingCapacitance*loadResistance)**2) /                  ((omega_n*self.couplingCapacitance)**2 * loadResistance)
        C = (self.capacitancePerUnitLength() * self.length)/2
        return omega_n * r_star * C

    def externalQualityFactorApprox(self):
        m = self.getModeFactor()
        q_in = 2 * pi * self.uncoupledResonantFrequency() * self.couplingCapacitance * self.characteristicImpedance()
        return (1/m) * (pi/(q_in**2))

    def externalQualityFactorQWref(self, inputPortImpedance = 50):
        m = self.getModeFactor()
        omega_0 = 2 * pi * self.uncoupledResonantFrequency()
        mBody = 1/(omega_0**2 * self.couplingCapacitance**2 * self.characteristicImpedance() * inputPortImpedance)
        return (pi/m) * mBody
    
    def loadedQualityFactor(self):
        return 1 / ( (1/self.internalQualityFactor()) + (1/self.externalQualityFactor()) )

    def getModeFactor(self):
        if self.resonatorType == 'half':
            m = 4.0 / (2.0 * self.mode)
        elif self.resonatorType == 'quarter':
            m = 4.0 / ((2.0 * self.mode) - 1)
        else:
            print('Error: Incorrect resonator type provided!')
            return -1
        return m

    def insertionLoss(self):
        g = self.internalQualityFactor()/self.externalQualityFactor()
        return -20 * np.log10(g/(g+1))
    
    def beta(self):
        omega_n = 2 * pi * self.uncoupledResonantFrequency()
        return omega_n * np.sqrt(self.totalInductancePerUnitLength()*self.capacitancePerUnitLength())
    
    def info(self):
        supplied_parameters = {
            'length': self.length,
            'conductorWidth': self.conductorWidth,
            'gapWidth': self.gapWidth,
            'conductorThickness': self.conductorThickness,
            'resonatorType': [self.resonatorType],
            'conductor': [self.conductor.material],
            'substrate': [self.substrate.material],
            'temperature': self.temperature,
            'couplingCapacitance': self.couplingCapacitance,
            'loadImpedance': self.loadImpedance,
            'loadBoundaryCondition': [self.loadBoundaryCondition],
            'mode': self.mode
        }
        calculated_parameters = {
            'effectivePermittivity': self.effectivePermittivity(),
            'capacitancePerUnitLength': self.capacitancePerUnitLength(),
            'totalInductancePerUnitLength': self.totalInductancePerUnitLength(),
            'geometricInductancePerUnitLength': self.geometricInductancePerUnitLength(),
            'kineticInductancePerUnitLength': self.kineticInductancePerUnitLength(),
            'londonPenetrationDepthT': self.londonPenetrationDepthT(),
            'characteristicImpedance': self.characteristicImpedance(),
            'inputImpedance': self.inputImpedance(),
            'uncoupledResonantFrequency': self.uncoupledResonantFrequency(),
            'coupledResonantFrequency': self.coupledResonantFrequency(),
            'effectiveCouplingCapacitance': self.effectiveCouplingCapacitance(),
            'internalQualityFactor': self.internalQualityFactor(),
            'externalQualityFactor': self.externalQualityFactor(),
            'loadedQualityFactor': self.loadedQualityFactor(),
            'insertionLoss': self.insertionLoss(),
            'beta': self.beta()
        }
        return [pd.DataFrame.transpose(pd.DataFrame(supplied_parameters, index=['Supplied parameters'])),
                pd.DataFrame.transpose(pd.DataFrame(calculated_parameters, index=['Calculated parameters']))]
