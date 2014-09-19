"""
basic classes
"""
import abc

import frame

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
    
    def encode(self, wav, order, frame_len, frame_step):
        self.frames = frame.framesig(wav, frame_len, frame_step, winfunc=lambda x:numpy.ones((1,x)))
        self.params = np.zeros((self.frames.shape[0], order))
        for i in range(self.frames.shape[0]):
            self.params[i, :] = self._encode_frame(self.frames[i])
            
    def decode(self, src_signal):
        src_frames = frame.framesig(src_signal, frame_len, frame_step, winfunc=lambda x:numpy.ones((1,x)))
        for i in range(self.frames.shape[0]):
            self.frames[i, :] = self._decode_frame(self.src_frames[i])
        wav = frame.deframesig(self.frames, frame_len, frame_step, winfunc=lambda x:numpy.ones((1,x)))
        return wav
    
    @abc.abstractmethod              
    def _encode_frame(self, src_signal):
        pass    
    @abc.abstractmethod              
    def _decode_frame(self, src_signal):
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
        self.src = self.encode_source(wav)
        self.trg = self.encode_filter(wav)
        
    def decode(self):                 
        src =  self.src.decode()
        wav = self.filt.decode(src)
        return wav
    
