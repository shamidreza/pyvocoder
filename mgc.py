"""
MGC class
"""
from vocoder import SourceFilterVocoder, Filter

def mgc_python(x, order):
    pass
def mlsa_python(param, src_signal):
    pass

mgc_func = mgc_python
mlsa_func = mlsa_python

class MGCVocoder(SourceFilterVocoder):       
    def __init__(self, order, fs):
        self.src = source.PulseNoiseSource(fs)
        self.filt = MGCFilter(order, fs)
        
class MGCFilter(Filter):
    def __init__(self, order, fs):
        Filter.__init__(self, order, fs)
    
    def _encode_frame(self, signal):
        m = mgc_func(signal, self.order)  
        return m  
    def _decode_frame(self, param, src_signal):
        return mlsa_func(param, src_signal)
    def filter2spectrum():
        pass
    def spectrum2filter():
        pass    
    
    