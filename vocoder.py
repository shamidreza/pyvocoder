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
    def decode(self):
        pass
    

class Source():
    __metaclass__ = abc.ABCMeta
    @abc.abstractmethod      
    def __init__(self):
        pass
    
    @abc.abstractmethod          
    def encode(self, wav):
        pass 
    
    @abc.abstractmethod      
    def decode(self):
        pass
    
class Filter():
    __metaclass__ = abc.ABCMeta
    
    @abc.abstractmethod          
    def __init__(self):
        pass
    
    @abc.abstractmethod          
    def encode(self, wav):
        pass    
    
    @abc.abstractmethod      
    def decode(self, src_signal):
        pass
    @abc.abstractmethod              
    def _filter_frame(self, src_signal):
        pass
    @abc.abstractmethod              
    def filter2spectrum():
        pass
    @abc.abstractmethod              
    def spectrum2filter():
        pass    
    
    
class SourceFilterVocoder(Vocoder):
    __metaclass__ = abc.ABCMeta
       
    @abc.abstractmethod     
    def __init__(self):
        self.src = Source()
        self.filt = Filter()
    
    def encode(self, wav):
        self. src = self.encode_source(wav)
        self. trg = self.encode_filter(wav)
        
    def decode(self):                 
        src =  self.src.decode()
        wav = self.filt.decode(src)
        return wav
    
