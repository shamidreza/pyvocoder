"""
basic classes
"""
import abc

class Vocoder():
    __metaclass__ = abc.ABCMeta
    @abc.abstractmethod    
    def __init__(self):
        pass    
    @abc.abstractmethod
    def encode(self, wav):
        pass
    @abc.abstractmethod
    def decode(self, parameters):
        pass
    
class SourceFilterVocoder(Vocoder):
    __metaclass__ = abc.ABCMeta
    def __init__(self):
        self.src = None
        self.filt = None
    @abc.abstractmethod    
    def encode_source(self, wav):
        pass
    
    @abc.abstractmethod    
    def encode_filter(self, wav):
        pass
    
    @abc.abstractmethod    
    def decode_source(self, parameters):
        pass
    
    @abc.abstractmethod    
    def decode_filter(self, parameters):
        pass    
    
    def filter(self, src, filt):
        # No abstract
        # frame src
        # _filter_frame(every frame)
        pass
    
    @abc.abstractmethod            
    def _filter_frame(self, src, filt):
        pass
    
    def encode(self, wav):
        self. src = self.encode_source(wav)
        self. trg = self.encode_filter(wav)
        
    def decode(self):                 
        src = self.decode_source(self.src)
        filt = self.decode_filter(self.filt)
        wav = self.filter(src, filt)
        return wav
    
