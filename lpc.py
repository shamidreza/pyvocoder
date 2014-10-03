"""
LPC class
"""
import numpy as np
from matplotlib import pyplot as pp

from vocoder import SourceFilterVocoder, Filter
import source

from scikits.talkbox import lpc
from scipy.signal import lfilter, hamming, freqz, deconvolve, convolve
        
def levinson_durbin():
    pass
def poly2lsf(a):
    a = a / a[0]        
    A = np.r_[a, 0.0]
    B = A[::-1]
    P = A - B  
    Q = A + B  
    
    P = deconvolve(P, np.array([1.0, -1.0]))[0]
    Q = deconvolve(Q, np.array([1.0, 1.0]))[0]
    
    roots_P = np.roots(P)
    roots_Q = np.roots(Q)
    
    angles_P = np.angle(roots_P[::2])
    angles_Q = np.angle(roots_Q[::2])
    angles_P[angles_P < 0.0] += np.pi
    angles_Q[angles_Q < 0.0] += np.pi
    lsf = np.sort(np.r_[angles_P, angles_Q])
    return lsf

def lsf2poly(L):
    order = len(L)
    Q = L[::2]
    P = L[1::2]
    poles_P = np.r_[np.exp(1j*P),np.exp(-1j*P)]
    poles_Q = np.r_[np.exp(1j*Q),np.exp(-1j*Q)]
    
    P = np.poly(poles_P)
    Q = np.poly(poles_Q)
    
    P = convolve(P, np.array([1.0, -1.0]))
    Q = convolve(Q, np.array([1.0, 1.0]))
    
    a = 0.5*(P+Q)
    return a[:-1]
class LPCVocoder(SourceFilterVocoder):
    def __init__(self, order, fs):
        self.src = source.PulseNoiseSource(fs)
        self.filt = LPCFilter(order, fs)

class LSFVocoder(SourceFilterVocoder):
    def __init__(self, order, fs):
        self.src = source.PulseNoiseSource(fs)
        self.filt = LSFFilter(order, fs)
    
class LPCFilter(Filter):
    def __init__(self, order, fs):
        Filter.__init__(self, order, fs)
   
    def _encode_frame(self, signal):
        A, e, k = lpc(signal, self.order) 
        return A
    def _decode_frame(self, param, src_signal):
        return lfilter([1], param, src_signal)
    def filter2spectrum(self, param):
        return np.abs(freqz([1], param)[1])
    def spectrum2filter():
        pass   
   
class LSFFilter(LPCFilter):
    def __init__(self, order, fs):
        LPCFilter.__init__(self, order, fs)
   
    def _encode_frame(self, signal):
        A = LPCFilter._encode_frame(self, signal)
        return poly2lsf(A)
    def _decode_frame(self, param, src_signal):
        A = lsf2poly(param)
        return LPCFilter._decode_frame(self, A, src_signal)
    def filter2spectrum(self, param):
        A = lsf2poly(param)        
        return LPCFilter.filter2spectrum(self, A)
    def spectrum2filter():
        pass    
 
if __name__ == '__main__':
    fname = '/Users/hamid/Code/gitlab/voice-conversion2/src/test/wav/ga_8_Kevin.wav'
    order = 18
    import wave
    spf = wave.open(fname, 'r') # http://www.linguistics.ucla.edu/people/hayes/103/Charts/VChart/ae.wav
    fs = spf.getframerate()
    # Get file as numpy array.
    x = spf.readframes(-1)
    x = np.fromstring(x, 'Int16')
    import source
    #f0 = source.getf0_python(x, fs, frame_len=0.025, frame_step=0.010 )
    #voc = LPCVocoder(order, fs)
    voc = LSFVocoder(order, fs)
    voc.encode(x)
    spec = voc.spectrogram()
    wav = voc.decode()
    import scipy
    scipy.io.wavfile.write('testlsf.wav', fs, wav)
    pass