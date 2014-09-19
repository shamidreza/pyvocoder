"""
basic classes
"""
import abc

import numpy as np
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
    
    def __init__(self, order, fs):
        self.order = order
        self.fs = fs
    def encode(self, wav, frame_len=0.025, frame_step=0.010):
        self.frame_len = int(frame_len * self.fs)
        self.frame_step = int(frame_step * self.fs)
        
        self.frames = frame.framesig(wav, self.frame_len, self.frame_step, winfunc=lambda x:np.ones((1,x)))
        param_size = len(self._encode_frame(self.frames[0]))
        self.params = np.zeros((self.frames.shape[0], param_size))
        for i in range(self.frames.shape[0]):
            self.params[i, :] = self._encode_frame(self.frames[i])
            
    def decode(self, src_signal):
        src_frames = frame.framesig(src_signal, self.frame_len, self.frame_step, winfunc=lambda x:np.ones((1,x)))
        for i in range(self.frames.shape[0]):
            self.frames[i, :] = self._decode_frame(self.params[i, :], self.src_frames[i])
        wav = frame.deframesig(self.frames, self.frame_len, self.frame_step, winfunc=lambda x:np.ones((1,x)))
        return wav
    
    @abc.abstractmethod              
    def _encode_frame(self, signal):
        pass    
    @abc.abstractmethod              
    def _decode_frame(self, params, src_signal):
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
        #self.src = Source()
        #self.filt = Filter()
        pass
    
    def encode(self, wav):
        self.src = self.src.encode(wav)
        self.trg = self.filt.encode(wav)
        
    def decode(self):                 
        src =  self.src.decode()
        wav = self.filt.decode(src)
        return wav
    
