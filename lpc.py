"""
LPC class
"""
import numpy as np
from matplotlib import pyplot as pp

from vocoder import SourceFilterVocoder, Filter
import source

from scikits.talkbox import lpc
from scipy.signal import lfilter, hamming, freqz
        
def levinson_durbin():
    pass

class LPCVocoder(SourceFilterVocoder):
    def __init__(self, order, fs):
        self.src = source.PulseNoiseSource(fs)
        self.filt = LPCFilter(order, fs)
    
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
    voc = LPCVocoder(order, fs)
    voc.encode(x)
    spec = voc.spectrogram()
    wav = voc.decode()
    import scipy
    scipy.io.wavfile.write('test.wav', fs, wav)
    pass