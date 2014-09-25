# This file includes routines for basic signal processing including framing and computing power spectra.
# Author: James Lyons 2012
# Website: https://github.com/jameslyons/python_speech_features

import numpy
import math

def framesig(sig,frame_len,frame_step,winfunc=lambda x:numpy.ones((1,x))):
    """Frame a signal into overlapping frames.

    :param sig: the audio signal to frame.
    :param frame_len: length of each frame measured in samples.
    :param frame_step: number of samples after the start of the previous frame that the next frame should begin.
    :param winfunc: the analysis window to apply to each frame. By default no window is applied.    
    :returns: an array of frames. Size is NUMFRAMES by frame_len.
    """
    slen = len(sig)
    frame_len = int(round(frame_len))
    frame_step = int(round(frame_step))
    if slen <= frame_len: 
        numframes = 1
    else:
        numframes = 1 + int(math.ceil((1.0*slen - frame_len)/frame_step))
        
    padlen = int((numframes-1)*frame_step + frame_len)
    
    zeros = numpy.zeros((padlen - slen,))
    padsignal = numpy.concatenate((sig,zeros))
    
    indices = numpy.tile(numpy.arange(0,frame_len),(numframes,1)) + numpy.tile(numpy.arange(0,numframes*frame_step,frame_step),(frame_len,1)).T
    indices = numpy.array(indices,dtype=numpy.int32)
    frames = padsignal[indices]
    win = numpy.tile(winfunc(frame_len),(numframes,1))
    return frames*win
    
    
def deframesig(frames,siglen,frame_len,frame_step,winfunc=lambda x:numpy.ones((1,x))):
    """Does overlap-add procedure to undo the action of framesig. 

    :param frames: the array of frames.
    :param siglen: the length of the desired signal, use 0 if unknown. Output will be truncated to siglen samples.    
    :param frame_len: length of each frame measured in samples.
    :param frame_step: number of samples after the start of the previous frame that the next frame should begin.
    :param winfunc: the analysis window to apply to each frame. By default no window is applied.    
    :returns: a 1-D signal.
    """
    frame_len = round(frame_len)
    frame_step = round(frame_step)
    numframes = numpy.shape(frames)[0]
    assert numpy.shape(frames)[1] == frame_len, '"frames" matrix is wrong size, 2nd dim is not equal to frame_len'
 
    indices = numpy.tile(numpy.arange(0,frame_len),(numframes,1)) + numpy.tile(numpy.arange(0,numframes*frame_step,frame_step),(frame_len,1)).T
    indices = numpy.array(indices,dtype=numpy.int32)
    padlen = (numframes-1)*frame_step + frame_len   
    
    if siglen <= 0: siglen = padlen
    
    rec_signal = numpy.zeros((1,padlen))
    window_correction = numpy.zeros((1,padlen))
    win = winfunc(frame_len)
    
    for i in range(0,numframes):
        window_correction[0,indices[i,:]] = window_correction[0,indices[i,:]] + win[0,:] + 1e-15 #add a little bit so it is never zero
        rec_signal[0, indices[i,:]] = rec_signal[0, indices[i,:]] + frames[i,:]
        
    rec_signal = rec_signal/window_correction
    return rec_signal[0:siglen]
    
def preemphasis(signal,coeff=0.95):
    """perform preemphasis on the input signal.
    
    :param signal: The signal to filter.
    :param coeff: The preemphasis coefficient. 0 is no filter, default is 0.95.
    :returns: the filtered signal.
    """    
    return numpy.append(signal[0],signal[1:]-coeff*signal[:-1])


