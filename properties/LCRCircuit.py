import cmath
import math

class LCRCircuit:
    'LCR oscillator circuit class'

    def __init__(self, inductance, capacitance, resistance, arrangement, couplingCapacitance=0):
        self.inductance = inductance
        self.capacitance = capacitance
        self.resistance = resistance
        self.arrangement = arrangement
        self.couplingCapacitance = couplingCapacitance
        self.frequency = self.resonantFrequency(inductance, capacitance)

        self.resonantFrequency = self.resonantFrequency(inductance, capacitance)
        self.impedance = self.impedance(self.inductance, self.capacitance, self.resistance, self.frequency, self.arrangement)
        self.qualityFactor = self.qualityFactor(self.resonantFrequency, self.resistance, self.capacitance, self.arrangement)

    def resonantFrequency(self, inductance, capacitance):
        return 1/(math.sqrt(inductance*capacitance))

    def impedance(self, inductance, capacitance, resistance, frequency, arrangement):
        if arrangement == 'Parallel':
            return 1 / resistance + 1 / (1j * frequency * inductance) + (1j * frequency * capacitance)
        elif arrangement == 'Series':
            return resistance + (1j * frequency * inductance) * (1 - 1/(1j * math.pow(frequency, 2) * capacitance * inductance))
        else:
            print('Error: Invalid arrangement type!')
            return -1

    def qualityFactor(self, resonantFrequency, resistance, capacitance, arrangement):
        if arrangement == 'Parallel':
            return resonantFrequency * resistance * capacitance
        elif arrangement == 'Series':
            return 1/(resonantFrequency * resistance * capacitance)
        else:
            print('Error: Invalid arrangement type!')
            return -1

    def nortonCapacitance(self, couplingCapacitance, loadResistance, resonantFrequency):

    def nortonResistance(self, ):

    def setFrequency(self, frequency):
        self.frequency = frequency


if __name__ == '__main__':
    mLCR = LCRCircuit(inductance=1, capacitance=1, resistance=1, arrangement='Parallel')

    print(mLCR)