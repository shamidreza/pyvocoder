"""
LPC class
"""
import numpy 
from vocoder import SourceFilterVocoder, Filter

def levinson_durbin():
    pass

class LPCVocoder(SourceFilterVocoder):
    def __init__(self):
        self.src = None
        self.filt = LPCFilter()
        
class LPCFilter(Filter):
    
    def __init__(self, order):
        Filter.__init__(order)
        from scikits.talkbox import lpc
        from scipy.signal import lfilter, hamming, freqz
   
    def _encode_frame(self, signal):
        A, e, k = lpc(signal, self.order) 
        return A
    def _decode_frame(self, param, src_signal):
        return lfilter([1], param, src_signal)
    def filter2spectrum(param):
        freqz([1], param)
    def spectrum2filter():
        pass    
 
if __name__ == '__main__':
    fname = '/Users/hamid/Code/gitlab/voice-conversion/src/test/wav/ga_8_Kevin.wav'
    import wave
    spf = wave.open(fname, 'r') # http://www.linguistics.ucla.edu/people/hayes/103/Charts/VChart/ae.wav
    # Get file as numpy array.
    x = spf.readframes(-1)
    x = numpy.fromstring(x, 'Int16')
     