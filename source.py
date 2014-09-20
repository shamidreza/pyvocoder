"""
Source class for SourceFilterVocoder
"""
import numpy as np

from vocoder import Source

def getf0(filename, fs=16000, frame_len=0.010, f0_min=50, f0_max=400):
    from subprocess import call
    import os
    tclsh = 'tclsh8.5'
    path = os.path.dirname(__file__)    
    call([tclsh, os.path.join(path, 'getf0.tcl'), filename,
          "-L", str(f0_min), "-H", str(f0_max),
          "-o", filename+'.f0',
          "-r", str(fs),
          "-s", str(frame_len)])    
    f = open(filename+'.f0', 'r')
    f0 = []
    while True:
        line = f.readline()
        if line == '':
            break
        cur_f0 = float(line.split(' ')[0])
        f0.append(cur_f0)
    f.close()
    return np.array(f0)
class PulseNoiseSource(Source):    
    def __init__(self):
        pass
    
    def encode(self, wav):
        self.f0 = getf0(fname)
    
    def decode(self):
        pass
    
class MixedExcitationSource(Source):    
    def __init__(self):
        pass
    
    def encode(self, wav):
        pass 
    
    def decode(self):
        pass
    
class ResidualExcitationSource(Source):    
    def __init__(self):
        pass
    
    def encode(self, wav):
        pass 
    
    def decode(self):
        pass