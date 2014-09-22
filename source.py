"""
Source class for SourceFilterVocoder
"""
import numpy as np
import scipy.signal
import copy
from vocoder import Source

def getf0_python(wav, fs, frame_len, frame_step):
    def lowpassfilter(x, cutoff, gain):
        h = scipy.signal.firwin(31, cutoff)
        h /= h.sum()
        return gain * scipy.signal.filtfilt(h, np.array([1.0]), x)    
    # convert dtype ro float64
    wav = wav.astype(np.float64)
    wav /= 30000.0    
    # convert dtype ro float64
    wav = wav.astype(np.float64)    
    # lowpass at 900Hz
    wav = lowpassfilter(wav, 2.0*900.0/fs, 1.0)
    #frame the signal
    ENERGY_THRESHOLD = 0.1
    wsize = frame_len
    wrate = frame_step
    st = 0
    en = wsize
    correlogram = np.zeros((0, wsize//2))
    pit = np.zeros(0)
    while True:
        if en > wav.shape[0]:
            break
        frame = wav[st:en].copy()
        # classify as speech/non-speech using enerfy
        speech = True
        rms = np.sqrt(np.mean(frame**2))
        if rms < ENERGY_THRESHOLD:
            speech = False
        # center clipping
        C = frame.max() * 0.4
        frame[np.abs(frame) < C] = 0.0
        frame[frame > C] = 1.0
        frame[frame < -C] = -1.0
        
        # autocorrelation
        r = np.correlate(frame, frame, 'same')[wsize//2:]
        if r[0] > 0:
            r /= r[0]
        correlogram = np.r_[correlogram, r.reshape((1,r.shape[0]))]
        
        # find peak
        r_limit = r[50:150]
        r_limit -= r_limit.min()
        r_limit /= r_limit.max()
        peak = np.argmax(r_limit)+50
        if r_limit.max() < 0.2:
            cur_pit = 0
        else:
            cur_pit = fs / peak
        pit = np.r_[pit, cur_pit]
        st += wrate
        en += wrate
    return pit


def getf0_tcl(filename, fs=16000, frame_len=0.010, f0_min=50, f0_max=400):
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
    def __init__(self, fs):
        self.fs = fs
        
    def encode(self, wav, frame_len=0.025, frame_step=0.010):
        self.frame_len = int(frame_len * self.fs)
        self.frame_step = int(frame_step * self.fs) 
        self.duration = wav.shape[0]
        self.f0 = getf0_python(wav, self.fs, frame_len=self.frame_len, frame_step=self.frame_step)
    
    def decode(self):
        ret = np.zeros(self.duration)
        cur = 0
        cur_inx = 0
        while True:
            cur += 1
            if self.f0[cur] > 0.0:
                ret[self.frame_step]=1
            
    
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