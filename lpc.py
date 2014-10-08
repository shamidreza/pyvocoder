"""
LPC class
"""
# python libraries
import logging
logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING) 

# third-party librareis
import numpy as np
from matplotlib import pyplot as pp
from scipy.signal import lfilter, hamming, freqz, deconvolve, convolve
from numpy.fft import rfft, irfft       

# pyvocoder
from vocoder import SourceFilterVocoder, Filter
import source

def lpc_python(x, order):
    def levinson(R,order):
        """ input: autocorrelation and order, output: LPC coefficients """
        a   = np.zeros(order+2)
        a[0] = -1
        k = np.zeros(order+1)
        # step 1: initialize prediction error "e" to R[0]
        e = R[0]
        # step 2: iterate over [1:order]
        for i in range(1,order+1):
            # step 2-1: calculate PARCOR coefficients
            k[i]= (R[i] - np.sum(a[1:i] * R[i-1:0:-1])) / e
            # step 2-2: update LPCs
            a[i] = np.copy(k[i])
            a_old = np.copy(a) 
            for j in range(1,i):
                a[j] -=  k[i]* a_old[i-j] 
            # step 2-3: update prediction error "e" 
            e = e * (1.0 - k[i]**2)
        return -1*a[0:order+1], e, -1*k[1:]
    # N: compute next power of 2
    n = x.shape[0]
    N = int(np.power(2, np.ceil(np.log2(2 * n - 1))))
    # calculate autocorrelation using FFT
    X = rfft(x, N)
    r = irfft(abs(X) ** 2)[:order+1]
    return levinson(r, order)

try:
    from scikits.talkbox import lpc as lpc
    lpc_func = lpc
except:
    logging.warning('scikits.talkbox could not be imported. Using python version of levinson recursion.')
    lpc_func = lpc_python
    
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
        A, e, k = lpc_func(signal, self.order)  
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