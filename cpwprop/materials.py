import pandas as pd

class Conductor:
    """
    CPW conductor material class
    """
    
    materials = {'Niobium': {
                     'superconductor':                    True,
                     'criticalTemperature':               9.2,
                     'londonPenetrationDepthZero':        33.3E-9,
                     'resistancePerUnitLength':           0,
                     'conductancePerUnitLength':          0,
                     'amplitudeAttenuationPerUnitLength': 2.4E-4},
                 'Niobium Nitride': {
                     'superconductor':                    True,
                     'criticalTemperature':               16.2,
                     'londonPenetrationDepthZero':        40.0E-9,
                     'resistancePerUnitLength':           0,
                     'conductancePerUnitLength':          0,
                     'amplitudeAttenuationPerUnitLength': 2.4E-4},
                'Copper': {
                     'superconductor':                    False}}

    def __init__(self, material):
        self.material                          = material
        self.superconductor                    = self.materials[material]['superconductor']
        if self.superconductor:
            self.criticalTemperature               = self.materials[material]['criticalTemperature']
            self.londonPenetrationDepthZero        = self.materials[material]['londonPenetrationDepthZero']
            self.resistancePerUnitLength           = self.materials[material]['resistancePerUnitLength']
            self.conductancePerUnitLength          = self.materials[material]['conductancePerUnitLength']
            self.amplitudeAttenuationPerUnitLength = self.materials[material]['amplitudeAttenuationPerUnitLength']
    
    def info():
        return pd.DataFrame(Conductor.materials)

class Substrate:
    """
    Substrate material class
    """
    
    materials = {'Silicon': {
                     'relativePermittivity': 11.9},
                 'Sapphire': {
                     'relativePermittivity': 10.2},
                'ArlonAD1000': {
                     'relativePermittivity': 10.2}}

    def __init__(self, material):
        self.material = material
        self.relativePermittivity = self.materials[material]['relativePermittivity']
        
    def info():
        return pd.DataFrame(Substrate.materials)
