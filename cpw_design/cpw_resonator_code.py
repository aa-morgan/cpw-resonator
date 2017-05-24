
# coding: utf-8

# ### Imports

# In[1]:

from scipy.special import ellipk
import math
from tabulate import tabulate
import numpy as np
import cmath

import matplotlib.pyplot as plt
get_ipython().magic('matplotlib inline')

from IPython.display import Image
from IPython.display import HTML, display


# ### Classes

# In[2]:

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


# In[3]:

class Substrate:
    """Substrate material class"""

    def __init__(self, material, relativePermittivity):
        self.material = material
        self.relativePermittivity = relativePermittivity


# In[4]:

class CPWResonator:
    """Coplanar wave-guide resonator class"""

    def __init__(self, length, conductorWidth, gapWidth, conductorThickness, resonatorType,
                 conductorMaterial='Niobium', substrateMaterial='Silicon',
                 temperature=4, couplingCapacitance=0E0, loadImpedance=50, loadBoundaryCondition='Short', modes=1):
        # Supplied parameters
        self.length = np.array(length)
        self.conductorWidth = np.array(conductorWidth)
        self.gapWidth = np.array(gapWidth)
        self.conductorThickness = np.array(conductorThickness)
        self.resonatorType = np.array(resonatorType)
        self.conductor = self.conductorProperties(conductorMaterial)
        self.substrate = self.substrateProperties(substrateMaterial)
        self.temperature = np.array(temperature)
        self.couplingCapacitance = np.array(couplingCapacitance)
        self.loadImpedance = np.array(loadImpedance)
        self.loadBoundaryCondition = np.array(loadBoundaryCondition)
        self.modes = np.array(modes)

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
        self.externalQualityFactorApprox = self.externalQualityFactorApproxFn(self.resonatorType, self.modes, self.uncoupledResonantFrequency,
            self.couplingCapacitance, self.characteristicImpedance)
        self.externalQualityFactor = self.externalQualityFactorFn(self.uncoupledResonantFrequency,
            self.capacitancePerUnitLength, self.length, self.couplingCapacitance)
        self.externalQualityFactorQWref = self.externalQualityFactorQWrefFn(self.resonatorType, self.modes, self.uncoupledResonantFrequency, self.capacitancePerUnitLength, self.length, 
            self.couplingCapacitance, self.characteristicImpedance)
        self.loadedQualityFactor = self.loadedQualityFactor(self.internalQualityFactor, self.externalQualityFactor)
        self.insertionLoss = self.insertionLoss(self.internalQualityFactor, self.externalQualityFactor)
        self.beta = self.beta(self.uncoupledResonantFrequency, self.totalInductancePerUnitLength, self.capacitancePerUnitLength)

    def effectivePermittivity(self, relativePermittivity):
        return (1 + relativePermittivity)/2

    def capacitancePerUnitLength(self, conductorWidth, gapWidth, effectivePermittivity):
        # Permittivity of free space
        freeSpacePermittivity = 8.85418782E-12

        # Complete elliptic integral of the first kind
        k = conductorWidth / (conductorWidth + 2 * gapWidth)
        k2 = np.sqrt(1 - np.power(k,2))

        # Total CPW capacitance p.u.l.
        return 4 * freeSpacePermittivity * (effectivePermittivity + 0) * (ellipk(k) / ellipk(k2))

    def geometricInducatancePerUnitLength(self, conductorWidth, gapWidth):
        # Permeability of freespace
        freeSpacePermeability = 1.25663706E-6

        # Complete elliptic integral of the first kind
        k = conductorWidth / (conductorWidth + 2 * gapWidth)
        k2 = np.sqrt(1 - np.power(k,2))

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
        geometricFactor = (1 / (2 * np.power(k, 2) * np.power(K, 2))) * (
            - np.log(conductorThickness / (4 * conductorWidth)) + ((2 * (conductorWidth + gapWidth))
            / (conductorWidth + 2 * gapWidth)) * np.log(gapWidth / (conductorWidth + gapWidth)) - (
            conductorWidth / (conductorWidth + 2 * gapWidth)) * np.log(conductorThickness
            / (4 * (conductorWidth + 2 * gapWidth))))

        # Kinetic Inductance p.u.l.
        return freeSpacePermeability * (np.power(londonPenetrationDepthT, 2) / (conductorWidth * conductorThickness)) * geometricFactor


    def londonPenetrationDepthT(self, temperature, criticalTemperature, londonPenetrationDepthZero):
        return londonPenetrationDepthZero / np.sqrt(1 - np.power((temperature / criticalTemperature), 4))

    def characteristicImpedance(self, inductance, capacitance, resistance=0, conductance=0, frequency=1):
        return np.sqrt(
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
        return 1 / (np.sqrt(totalInductancePerUnitLength*capacitancePerUnitLength) * m * length)

    def coupledResonantFrequency(self, totalInductancePerUnitLength, capacitancePerUnitLength, length, resonatorType, couplingCapacitance, loadImpedane, modes):
        m = self.getModeFactor(resonatorType, modes)

        # Pre-coupled
        uncoupledResonantFrequency = 1 / (np.sqrt(totalInductancePerUnitLength*capacitancePerUnitLength) * m * length)

        # Post-coupled
        effectiveCouplingCapacitance = self.effectiveCouplingCapacitance(couplingCapacitance, uncoupledResonantFrequency, loadImpedane)

        return 1 / (np.sqrt((totalInductancePerUnitLength*length) * ((capacitancePerUnitLength*length) + effectiveCouplingCapacitance)) * m)

    def effectiveCouplingCapacitance(self, couplingCapacitance, frequency, loadImpedane):
        return couplingCapacitance /                (1 + np.power(frequency * couplingCapacitance * loadImpedane * math.pi, 2))

    def internalQualityFactor(self, length, resonatorType, modes, conductor):
        m = self.getModeFactor(resonatorType, modes)
        return (1/m) * (math.pi/(conductor.amplitudeAttenuationPerUnitLength*length))
    
    def externalQualityFactorFn(self, uncoupledResonantFrequency, capacitancePerUnitLength, length, couplingCapacitance, loadResistance=50):
        omega_n = 2 * math.pi * uncoupledResonantFrequency
        r_star = (1+(omega_n*couplingCapacitance*loadResistance)**2) / ((omega_n*couplingCapacitance)**2 * loadResistance)
        C = (capacitancePerUnitLength * length)/2
        return omega_n * r_star * C

    def externalQualityFactorApproxFn(self, resonatorType, modes, uncoupledResonantFrequency, couplingCapacitance, characteristicImpedance):
        m = self.getModeFactor(resonatorType, modes)
        q_in = 2 * math.pi * uncoupledResonantFrequency * couplingCapacitance * characteristicImpedance
        return (1/m) * (math.pi/(q_in**2))

    def externalQualityFactorQWrefFn(self, resonatorType, modes, uncoupledResonantFrequency, capacitancePerUnitLength, length, couplingCapacitance, characteristicImpedance):
        m = self.getModeFactor(resonatorType, modes)
        omega_0 = 2 * math.pi * uncoupledResonantFrequency
        inputPortImpedance = 50
        mBody = 1/(omega_0**2 * couplingCapacitance**2 * characteristicImpedance * inputPortImpedance)
        return (math.pi/m) * mBody
    
    def loadedQualityFactor(self, internalQualityFactor, externalQualityFactor):
        return 1/((1/internalQualityFactor) + (1/externalQualityFactor))

    def getModeFactor(self, resonatorType, modes):
        if resonatorType == 'half':
            m = 4.0 / (2.0 * modes)
        elif resonatorType == 'quarter':
            m = 4.0 / ((2.0 * modes) - 1)
        else:
            print('Error: Incorrect resonator type provided!')
            return -1
        return m

    def insertionLoss(self, internalQualityFactor, externalQualityFactor):
        g = internalQualityFactor/externalQualityFactor
        return -20 * np.log10(g/(g+1))
    
    def beta(self, uncoupledResonantFrequency, inductancePerUnitLength, capacitancePerUnitLength):
        omega_n = 2 * math.pi * uncoupledResonantFrequency
        return omega_n * np.sqrt(inductancePerUnitLength*capacitancePerUnitLength)

    def conductorProperties(self, material):
        return{
            'Niobium': Conductor(material=material, criticalTemperature=9.2, londonPenetrationDepthZero=33.3E-9,
                resistancePerUnitLength=0, conductancePerUnitLength=0, amplitudeAttenuationPerUnitLength=2.4E-4),
            'Niobium Nitride': Conductor(material=material, criticalTemperature=16.2, londonPenetrationDepthZero=40E-9,
                resistancePerUnitLength=0, conductancePerUnitLength=0, amplitudeAttenuationPerUnitLength=2.4E-4),
        }[material]

    def substrateProperties(self, material):
        return{
            'Silicon': Substrate(material=material, relativePermittivity=11.9),
            'Sapphire': Substrate(material=material, relativePermittivity=10.2),
        }[material]


# ## ABCD Transmission matrix

# The ABCD matrix can be used to fit the transmission/reflection resonance peaks of a resonator. The matrix is defined as,
# 
# \begin{equation}
# \begin{bmatrix}
#     A & B \\
#     C & D
# \end{bmatrix}
# = 
# \begin{bmatrix}
#     1 & Z_{in} \\
#     0 & 0
# \end{bmatrix}
# \begin{bmatrix}
#     t_{11} & t_{12} \\
#     t_{21} & t_{22}
# \end{bmatrix}
# \begin{bmatrix}
#     1 & Z_{out} \\
#     0 & 0
# \end{bmatrix}
# \end{equation}
# 
# where, 
# 
# \begin{equation}
#     t_{11} = t_{22} = \cosh{\gamma l}
# \end{equation}
# 
# \begin{equation}
#     t_{12} = Z_0 \sinh(\gamma l)
# \end{equation}
# 
# \begin{equation}
#     t_{21} = \frac{1}{Z_0} \sinh(\gamma l)
# \end{equation}
# 
# where $Z_0$ is the characteristic impedance, and $\gamma = \alpha + i \beta$ where $\alpha$ is the amplitude attenuation coefficient and $\beta$ is given by,
# 
# \begin{equation}
#     \beta = \frac{\omega_n}{v_{ph}} = \omega_n \sqrt{L_l C_l} = \omega_n \frac{\sqrt{\epsilon_0}}{c}
# \end{equation}
# 
# The transmission parameter, $S_{21}$ is then given by,
# 
# \begin{equation}
#     S_{21} = \frac{2}{A + \frac{B}{R_L} + C R_L + D}
# \end{equation}
# 
# where $R_L$ is the load resistance.

# In[5]:

class ABCD:
    'ABCD matrix method'

    def __init__(self, freq, length, alpha, beta, couplingCapacitance, charImpedance, loadResistance=50):
        self.freq = freq
        self.length = length
        self.gamma = alpha + 1j*beta
        self.couplingCapacitance = couplingCapacitance
        self.charImpedance = charImpedance
        self.loadResistance = loadResistance
        
        self.s21 = self.s21(self.abcd(
            self.input(self.freq, self.couplingCapacitance),
            self.transmission(self.freq, self.length, self.gamma, self.charImpedance),
            self.output(self.freq, self.couplingCapacitance)),
            self.loadResistance)

    def input(self, freq, couplingCapacitance):
        n=np.size(freq)
        Zin = 1/(1j*freq*couplingCapacitance)
        return np.append(np.ones(n),[np.zeros(n),Zin,np.ones(n)]).reshape(n,2,2,order='F')
        
    def output(self, freq, couplingCapacitance):
        n=np.size(freq)
        Zin = 1/(1j*freq*couplingCapacitance)
        return np.append(np.ones(n),[np.zeros(n),Zin,np.ones(n)]).reshape(n,2,2,order='F')
        
    def transmission(self, freq, length, gamma, charImpedance):
        n=np.size(freq)
        t11 = np.cosh(gamma * length)
        t12 = charImpedance * np.sinh(gamma * length)
        t21 = (1/charImpedance) * np.sinh(gamma * length)
        t22 = np.cosh(gamma * length)
        return np.append(t11*np.ones(n),[t21*np.ones(n),t12*np.ones(n),t22*np.ones(n)])             .reshape(n,2,2,order='F')
    
    def abcd(self, input, transmission, output):
        return input*transmission*output
    
    def s21(self, pABCD, loadResistance):
        A = pABCD[:,0,0]
        B = pABCD[:,0,1]
        C = pABCD[:,1,0]
        D = pABCD[:,1,1]
        RL = loadResistance
        return 2/( A + (B/RL)+ (C*RL) + D )


# # 2 - Calculating Rydberg atom transition frequencies

# The wavelength of the transition between the $n_1$th and $n_2$th levels is given by,
# 
# \begin{equation}
#     \frac{1}{\lambda} = R_{M} \left( \frac{1}{(n_1-\delta_1)^2} - \frac{1}{(n_2-\delta_2)^2} \right)
# \end{equation}
# 
# where $\delta_x$ are the quantum defects, and $R_{M}$ is the reduced mass,
# 
# \begin{equation}
#     R_{M} = \frac{R_{\infty}}{1+\frac{m_e}{M}}
# \end{equation}
# 
# where $R_{\infty}$ is the Rydberg constant with an infinite mass nucleus, $m_e$ is the electron mass, and $M$ is the mass of the nucleus. $R_{\infty}$ is given by,
# 
# \begin{equation}
#     R_{\infty} = \frac{m_e e^4}{8 \epsilon_0^2 h^3 c} = 1.0973731568508 \times 10^7 m^{-1}
# \end{equation} 
# 
# The frequency of the transition is then,
# 
# \begin{equation}
#     f = \frac{c}{\lambda}
# \end{equation}
# 
# where $c$ is the speed of light.

# In[42]:

class RydbergHelium:
    def __init__(self):
        self.z = 2
        self.r_inf = 1.0973731568508 * 10**7
        self.c = 2.99792458 * 10**8
        self.h = 6.62607004 * 10**-34
        self.e = 1.60217662 * 10**-19
        self.m_e = 9.10938356 * 10**-31
        self.m_p = 1.6726219 * 10**-27
    
    def energy_level(self, n, l):
        r_m = self.r_inf / (1 + (self.m_e/(2*self.z*self.m_p)))
        defect = self.quantum_defect(n, l)
        wavelength = 1 / ( r_m * ( 1/float(n-defect)**2 ) )
        return self.h * self.c / wavelength
    
    def energy_transition(self, n_from, n_to, l_from=6, l_to=6):
        return np.abs(self.energy_level(n_from, l_from) - self.energy_level(n_to, l_to))
    
    def quantum_defect(self, n, l):
        # Routine to calculate the quantum defects of the triplet Rydberg states of helium
        # From Martin, Phys. Rev. A, vol. 36, pp. 3575-3589 (1987)

        #    s            p            d           f           g           h           +
        a = [0.29665486,  0.06835886,  0.00289043, 0.00043924, 0.00012568, 0.00004756, 0]
        b = [0.03824614, -0.01870111, -0.0064691, -0.0017850, -0.0008992, -0.000552  , 0]
        c = [0.0082574,  -0.0117730,   0.001362,   0.000465,   0.0007,     0.00112   , 0]
        d = [0.000359,   -0.008540,   -0.00325,    0,          0,          0         , 0]

        if l <= 5:
            idx = l;
        else:
            idx = 6
            
        m = n - a[idx];
        return a[idx] + b[idx]*m**(-2) + c[idx]*m**(-4) + d[idx]*m**(-6);
    
    def energy_ionisation(self):
        # E/hc = 1/lambda (cm^-1)
        return (self.h * self.c) * (198310.6663720 * 100)
    
    def energy_1s3p(self):
        # E/hc = 1/lambda (cm^-1)
        return (self.h * self.c) * (185564.561920 * 100) # J = 2
    
    def energy_1s2s(self):
        # E/hc = 1/lambda (cm^-1)
        return (self.h * self.c) * (159855.9743297 * 100)
    
    def energy_1s3p_nl(self, n, l):
        return (self.energy_ionisation() - self.energy_1s3p()) - self.energy_level(n, l)
    
    def frequency(self, E):
        return E / self.h
    
    def wavelength(self, E):
        return self.h * self.c / E


# In[ ]:



