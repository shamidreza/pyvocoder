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
    @abc.abstractmethod
    def spectrogram(self):
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
        winfunc=lambda x:np.hamming(x).reshape((1,x))
        self.frame_len = int(frame_len * self.fs)
        self.frame_step = int(frame_step * self.fs)        
        wav = frame.preemphasis(wav,coeff=0.97)##       
        self.energy = frame.get_energy(wav, self.frame_len, self.frame_step, winfunc=winfunc)        
        
        self.frames = frame.framesig(wav, self.frame_len, self.frame_step, winfunc=winfunc)
        param_size = len(self._encode_frame(self.frames[0]))
        self.params = np.zeros((self.frames.shape[0], param_size))
        for i in range(self.frames.shape[0]):
            self.params[i, :] = self._encode_frame(self.frames[i])
    def decode(self, src_signal):
        winfunc=lambda x:np.hamming(x).reshape((1,x))
        src_energy = frame.get_energy(src_signal, self.frame_len, self.frame_step)
        gain = self.energy / src_energy
        gain_interp = np.interp(np.linspace(0,1,src_signal.shape[0]), np.linspace(0, 1, gain.shape[0]), gain)
        src_signal *= gain_interp
        src_frames = frame.framesig(src_signal, self.frame_len, self.frame_step)
        
        for i in range(self.frames.shape[0]):
            self.frames[i, :] = self._decode_frame(self.params[i, :], np.r_[src_frames[max(0,i-1)],src_frames[i]])[self.frame_len:]
        wav = frame.deframesig(self.frames, src_signal.shape[0], self.frame_len, self.frame_step, winfunc=winfunc)
        wav = frame.deemphasis(wav,coeff=0.97)##
        #wav *= gain_interp
        return wav
    
    @abc.abstractmethod              
    def _encode_frame(self, signal):
        pass    
    @abc.abstractmethod              
    def _decode_frame(self, params, src_signal):
        pass
    @abc.abstractmethod              
    def filter2spectrum(self, params):
        pass
    @abc.abstractmethod              
    def spectrum2filter(self):
        pass    
    
    
class SourceFilterVocoder(Vocoder):
    __metaclass__ = abc.ABCMeta
        
    @abc.abstractmethod     
    def __init__(self, order, fs):
        #self.src = Source()
        #self.filt = Filter()
        pass
    
    def encode(self, wav):
        self.src.encode(wav)
        self.filt.encode(wav)
        
    def decode(self):                 
        src =  self.src.decode()
        wav = self.filt.decode(src)
        wav /= wav.max()
        wav *= 30000.0
        wav = wav.astype(np.int16)
        return wav
    
    def spectrogram(self):
        spec = np.zeros((self.filt.params.shape[0], 512))
        for i in range(self.filt.params.shape[0]):
            spec[i, :] = self.filt.filter2spectrum(self.filt.params[i,:])    
        return spec