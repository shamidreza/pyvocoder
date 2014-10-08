"""
Mel Generalized Cepstral Analysis/Synthesis class
Based on SPTK-3.7/bin/mgcep/_mgcep.c
"""
# python libraries
import logging
logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING) 

# third-party librareis
import numpy as np
from matplotlib import pyplot as pp

# pyvocoder
from vocoder import SourceFilterVocoder, Filter
import source

def _gain(er, c, m, g):
    """
    /*  gain(epsilon) calculation  */
    static double gain(double *er, double *c, int m, double g)
    {
       int i;
       double t;
    
       if (g != 0.0) {
          for (t = 0.0, i = 1; i <= m; i++)
             t += er[i] * c[i];
          return (er[0] + g * t);
       } else
          return (er[0]);
    }
    """    
    if g != 0.0:
        t = 0.0
        for i in range(1,m+1):
            t += er[i] * c[i]
        return er[0] + g * t
    else:
        return er[0]
def _b2c(b, m1, c, m2, a):
    """
    static void b2c(double *b, int m1, double *c, int m2, double a)
    {
       int i, j;
       static double *d = NULL, *g;
       static int size;
       double k;
    
       if (d == NULL) {
          size = m2;
          d = dgetmem(size + size + 2);
          g = d + size + 1;
       }
       if (m2 > size) {
          free(d);
          size = m2;
          d = dgetmem(size + size + 2);
          g = d + size + 1;
       }
    
       k = 1 - a * a;
    
       fillz(g, sizeof(*g), m2 + 1);
    
       for (i = -m1; i <= 0; i++) {
          d[0] = g[0];
          g[0] = b[-i];
    
          if (1 <= m2)
             g[1] = k * d[0] + a * (d[1] = g[1]);
    
          for (j = 2; j <= m2; j++)
             g[j] = d[j - 1] + a * ((d[j] = g[j]) - g[j - 1]);
       }
       movem(g, c, sizeof(*g), m2 + 1);
    
       return;
    }
    
    """
    i = j = 0
    d = None
    g = None
    size = 0
    k = 0.0
    
    if d is None:
        size = m2
        d = np.zeros(size + size + 2)
        g = d + size + 1
    if m2 > size:
        d = None
        size = m2
        d = np.zeros(size + size + 2)
        g = d + size + 1   
    
    k = 1 - a * a
    
    for i in range(-m1, 0+1):
        d[0] = g[0]
        g[0] = b[-i]
        if 1 <= m2:
            d[1] = g[1] # not sure about the order of this line and next line
            g[1] = k * d[0] + a * d[1]
        for j in range(2, m2+1):
            d[j] = g[j] # not sure about the order of this line and next line
            g[j] = d[j - 1] + a * (d[j] - g[j - 1])
    return g # compare with movem
     
     
           
    
def _ptrans(p, m, a):
    """
    /*  recursion for p(m)  */
    static void ptrans(double *p, int m, double a)
    {
       double d, o;
    
       d = p[m];
       for (m--; m > 0; m--) {
          o = p[m] + a * d;
          d = p[m];
          p[m] = o;
       }
       o = a * d;
       p[m] = (1. - a * a) * p[m] + o + o;
    
       return;
    }
    """
    d = o = 0.0
    d = p[m]
    for mp in range(m-1, 0, -1):
        o = p[mp] + a * d
        d= p[mp]
        p[mp] = o
    o = a * d
    # not sure about the following line
    p[0] =  (1.0 - a * a) * p[0] + o + o
    
def _qtrans(q, m, a):
    """
    /*  recursion for q(m)  */
    static void qtrans(double *q, int m, double a)
    {
       int i;
       double d, o;
    
       m += m;
       i = 1;
       d = q[i];
       for (i++; i <= m; i++) {
          o = q[i] + a * d;
          d = q[i];
          q[i] = o;
       }
    
       return;
    }
    
    """
    i = 0
    d = o = 0.0
    m += m
    i = 1
    d = q[i]
    for ii in range(2, m+1):
        o = q[ii] + a * d
        d = q[ii]
        q[ii] = o
        
def mgc_python(x, order):
    pass
def mlsa_python(param, src_signal):
    pass

mgc_func = mgc_python
mlsa_func = mlsa_python

class MGCVocoder(SourceFilterVocoder):       
    def __init__(self, order, fs):
        self.src = source.PulseNoiseSource(fs)
        self.filt = MGCFilter(order, fs)
        
class MGCFilter(Filter):
    def __init__(self, order, fs):
        Filter.__init__(self, order, fs)
    
    def _encode_frame(self, signal):
        m = mgc_func(signal, self.order)  
        return m  
    def _decode_frame(self, param, src_signal):
        return mlsa_func(param, src_signal)
    def filter2spectrum():
        pass
    def spectrum2filter():
        pass    
    
    