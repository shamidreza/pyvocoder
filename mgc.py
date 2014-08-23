"""
MGC class
"""
from vocoder import SourceFilterVocoder, Filter

def mgc_iteration():
    pass
def mlsa():
    pass

class MGCVocoder(SourceFilterVocoder):       
    def __init__(self):
        self.src = None
        self.filt = MGCFilter()
        
class MGCFilter(Filter):
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
    
    