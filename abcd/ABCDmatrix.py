import numpy as np

class ABCD:
    'ABCD matrix method'

    def __init__(self, length, alpha, beta, couplingCapacitance, freqRange, freqStep):
        self.length = length
        self.alpha = alpha
        self.beta = beta
        self.couplingCapacitance = couplingCapacitance
        self.freqRange = freqRange
        self.freqStep = freqStep

    def input(self, couplingCapacitance, freqRange, freqStep):
        freq = np.array([freqRange(0):])
        Zin = 1/1j*
