"""
LPC class
"""
from vocoder import SourceFilterVocoder, Filter

def levinson_durbin():
    pass

class LPCVocoder(SourceFilterVocoder):
    def __init__(self):
        self.src = None
        self.filt = LPCFilter()
        
class LPCFilter(Filter):
    def __init__(self):
        pass
    
    def encode(self, wav):
        pass    
    
    def decode(self, src_signal):
        pass
    def _filter_frame(self, src_signal):
        pass
    def filter2spectrum():
        pass
    def spectrum2filter():
        pass    
    
    