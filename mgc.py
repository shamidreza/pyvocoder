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
        
def _newton(x, flng, c, m, a, g, n, j, f):
    """
    double newton(double *x, const int flng, double *c, const int m, const double a,
                  const double g, const int n, const int j, const double f)
    {
       int i, m2;
       double t = 0, s, tr, ti, trr, tii;
       static double *cr = NULL, *ci, *pr, *qr, *qi, *rr, *ri, *b;
       static int size_cr, size_b;
    
       if (cr == NULL) {
          cr = dgetmem(7 * flng);
          ci = cr + flng;
          pr = ci + flng;
          qr = pr + flng;
          qi = qr + flng;
          rr = qi + flng;
          ri = rr + flng;
          size_cr = flng;
    
          b = dgetmem(m + 1);
          size_b = m;
       }
       if (flng > size_cr) {
          free(cr);
          cr = dgetmem(7 * flng);
          ci = cr + flng;
          pr = ci + flng;
          qr = pr + flng;
          qi = qr + flng;
          rr = qi + flng;
          ri = rr + flng;
          size_cr = flng;
       }
       if (m > size_b) {
          free(b);
          b = dgetmem(m + 1);
          size_b = m;
       }
    
       m2 = m + m;
    
       fillz(cr, sizeof(*cr), flng);
       movem(&c[1], &cr[1], sizeof(*c), m);
    
       if (a != 0.0)
          b2c(cr, m, cr, n, -a);
    
       fftr(cr, ci, flng);          /* cr +j ci : FFT[c]  */
    
       if (g == -1.0)
          movem(x, pr, sizeof(*x), flng);
       else if (g == 0.0)
          for (i = 0; i < flng; i++)
             pr[i] = x[i] / exp(cr[i] + cr[i]);
       else
          for (i = 0; i < flng; i++) {
             tr = 1 + g * cr[i];
             ti = g * ci[i];
             s = (trr = tr * tr) + (tii = ti * ti);
             t = x[i] * pow(s, -1.0 / g);
             pr[i] = (t /= s);
             rr[i] = tr * t;
             ri[i] = ti * t;
             t /= s;
             qr[i] = (trr - tii) * t;
             s = tr * ti * t;
             qi[i] = s + s;
          }
    
       ifftr(pr, ci, flng);
    
       if (a != 0.0)
          b2c(pr, n, pr, m2, a);
    
       if (g == 0.0 || g == -1.0) {
          movem(pr, qr, sizeof(*pr), m2 + 1);
          movem(pr, rr, sizeof(*pr), m + 1);
       } else {
          ifft(qr, qi, flng);
          ifft(rr, ri, flng);
    
          if (a != 0.0) {
             b2c(qr, n, qr, n, a);
             b2c(rr, n, rr, m, a);
          }
       }
    
       if (a != 0.0) {
          ptrans(pr, m, a);
          qtrans(qr, m, a);
       }
    
       /*  c[0] : gain, t : epsilon  */
       if (g != -1.0)
          c[0] = sqrt(t = gain(rr, c, m, g));
    
       if (g == -1.0)
          fillz(qr, sizeof(*qr), m2 + 1);
       else if (g != 0.0)
          for (i = 2; i <= m2; i++)
             qr[i] *= 1.0 + g;
    
       if (theq(pr, &qr[2], &b[1], &rr[1], m, f)) {
          fprintf(stderr, "mgcep : Error in theq() at %dth iteration!\n", j);
          exit(1);
       }
    
       for (i = 1; i <= m; i++)
          c[i] += b[i];
    
       /*  c[0] : gain, t : epsilon  */
       if (g == -1.0)
          c[0] = sqrt(t = gain(rr, c, m, g));
    
       return (log(t));
    }
    
    """
    i, m2 = 0
    t = s = tr = ti = trr = tii = 0.0
    cr = ci = pr = qr = qi = rr = ri = b = None
    size_cr = size_b = 0
    if cr is None:
        cr = np.zeros(flng)
        ci = np.zeros(flng)
        pr = np.zeros(flng)
        qr = np.zeros(flng)
        qi = np.zeros(flng)
        rr = np.zeros(flng)
        ri = np.zeros(flng)
        size_cr = flng   
   
        b = np.zeros(m+1)
        size_b = m
    
    if flng > size_cr: # what?!
        cr = np.zeros(flng)
        ci = np.zeros(flng)
        pr = np.zeros(flng)
        qr = np.zeros(flng)
        qi = np.zeros(flng)
        rr = np.zeros(flng)
        ri = np.zeros(flng)
        size_cr = flng  
   
    if m > size_b:
        b = np.zeros(m+1)
        size_b = m
   
    m2 = m + m
    cr[:] = 0.0
    
    cr[:] = c[:]
    if a != 0.0:
        _b2c(cr, m, cr, n, -a)
    # ci = np.fft.rfft(cr)
    # or the following 3 lines
    cr = np.fft.rfft(c)
    ci = cr.imag
    cr = cr.real
    
    if g == -1.0:
        pr[:] = x[:]
    elif g == 0.0:
        for i in range(0, flng):
            pr[i] = x[i] / np.exp(cr[i] + cr[i])
    else:
        for i in range(0, flng):
            tr = 1 + g * cr[i]
            ti = g * ci[i]
            trr = tr * tr
            tii = ti * ti
            s = (trr) + (tii)
            t = x[i] * np.power(s, -1.0 / g);
            t /= s # order
            pr[i] = (t);
            rr[i] = tr * t
            ri[i] = ti * t
            t /= s
            qr[i] = (trr - tii) * t
            s = tr * ti * t
            qi[i] = s + s          

      

    pr_tmp = np.fft.irfft(pr)
    pr = pr_tmp.real
    ci = pr_tmp.imag
    
    if a != 0.0:
        _b2c(pr, n, pr, m2, a)
    
    if g == 0.0 or g == -1.0:
        qr[:] = pr[:] #?
        rr[:] = pr[:]
    else:
        qr_tmp = np.fft.irfft(qr)
        qr = qr_tmp.real
        qi = qr_tmp.imag  
        rr_tmp = np.fft.irfft(rr)
        rr = rr_tmp.real
        ri = rr_tmp.imag            
        
        if a != 0.0:
            _b2c(qr, n, qr, n, a)
            _b2c(rr, n, rr, m, a)
    
    if a != 0.0:
        _ptrans(pr, m, a)
        _qtrans(qr, m, a)
          
    if g != -1.0:
        t = _gain(rr, c, m, g)
        c[0] = np.sqrt(t)
    if g == -1.0:
        qr[:] = 0.0
    elif g != 0.0:
        for i in range(2, m2+1):
            qr[i] *= 1.0 + g
            
    if _theq(pr, qr[2:], b[1:], rr[1:], m, f) != 0:
        raise Exception("mgcep : Error in _theq() at %dth iteration!\n") 

    for i in range(1, m+1):
        c[i] += b[i]
    if g == -1.0:
        t = _gain(rr, c, m, g)
        c[0] = np.sqrt(t)
    
    return np.log(t)
  
def _mgcep(xw, flng, b, m, a, g, n, itr1, itr2, dd, etype, e, f, itype):
    i = j = flag = 0
    x = y = d = None 
    size_x = size_c = 0
    ep = epo = eps = 0.0
    _min = _max = 0.0
       
    if etype == 1 and e < 0.0:
        raise Exception("mgcep : value of e must be e>=0!") 
    
    if etype == 2 and e >= 0.0:
        raise Exception("mgcep : value of E must be E<0!") 
    
    if etype == 1:
        eps = e
        
    if x == None:
        x = np.zeros(flng)
        y = np.zeros(flng)
        size_x = flng
        d = np.zeros(m+1)
        size_c = m
    if flng > size_x:
        x = np.zeros(flng)
        y = np.zeros(flng)
        size_x = flng  

    if m > size_c:
        d = np.zeros(m+1)
        size_c = m
    
    x[:] = xw[:]

    if itype == 0: # time-domain signal
        xt = np.fft.rfft(x)
        x = xt.real
        y = xt.imag
        for i in range(0, flng):
            x[i] = x[i] * x[i] + y[i] * y[i] + eps
    elif itype == 1: # dB
        for i in range(0, 1+flng/2):
            x[i] = np.exp((x[i]/20.0)*np.log(10.0))
            x[i] = x[i] * x[i] + eps
    elif itype == 2: # log
        for i in range(0, 1+flng/2):
            x[i] = np.exp(x[i])
            x[i] = x[i] * x[i] + eps    
    elif itype == 3: # amplitude
        for i in range(0, 1+flng/2):
            x[i] = x[i] * x[i] + eps            
    elif itype == 4: # periodogram
        for i in range(0, 1+flng/2):
            x[i] = x[i] + eps     
    else:
        raise Exception("mgcep : Input type %d is not supported!") 
      
    if itype > 0:
        for i in range(1, flng/2):
            x[flng-i] = x[i]

    if etype == 2 and e < 0.0:
        _max = x[0]
        for i in range(1, flng):
            if _max < x[i]:
                _max = x[i]
        _max = np.sqrt(_max)
        _min = _max * np.power(10.0, e/20.0)
        _min = _min * _min
        for i in range(0, flng):
            if x[i] < _min:
                x[i] = _min
        
    b[:] = 0.0
    ep = _newton(x, flng, b, m, a, -1.0, n, 0, f)
    if g != -1.0:
        if a != 0.0:
            _ignorm(b, b, m, -1.0)
            _b2mc(b, b, m, a)
            _gnorm(b, d, m, -1.0)
        else:
            d[:] = b[:]
             
        _gc2gc(d, m, -1.0, b, m, g)
        
        if 1 != 0.0:
            _ignorm(b, b, m, g)
            _mc2b(b, b, m, a)
            _gnorm(b, b, m, g) 
    
    # Newton-Raphson method
    if g != -1.0:
        for j in range(1, itr2+1):
            epo = ep
            ep = _newton(x, flng, b, m, a, g, n, j, f)
            
            if j >= itr1:
                if np.abs((epo - ep) / ep) < dd:
                    flag = 1
            if flag == 1:
                break
    
    if flag:
        return 0
    else:
        return -1
            
def _ignorm(c1, c2, m, g):
    """
    void ignorm(double *c1, double *c2, int m, const double g)
    {
       double k;
    
       k = pow(c1[0], g);
       if (g != 0.0) {
          for (; m >= 1; m--)
             c2[m] = k * c1[m];
          c2[0] = (k - 1.0) / g;
       } else {
          movem(&c1[1], &c2[1], sizeof(*c1), m);
          c2[0] = log(c1[0]);
       }
    
       return;
    }
    
    """
    k = 0.0
    k = np.power(c1[0], g)
    if g != 0.0:
        for mm in range(m, 0, -1):
            c2[mm] = k * c1[mm]
        c2[0] = (k-1.0)/g
    else:
        c2[1:] = c1[1:]
        c2[0] = np.log(c1[0])
        
def _b2mc(b, mc, m, a):
    """
    void b2mc(double *b, double *mc, int m, const double a)
    {
       double d, o;
    
       d = mc[m] = b[m];
       for (m--; m >= 0; m--) {
          o = b[m] + a * d;
          d = b[m];
          mc[m] = o;
       }
    
       return;
    }    
    """
    d = o = 0.0
    d = mc[m] = b[m]
    for mm in range(m-1, -1, -1): # sure?
        o = b[m] + a * d
        d = b[m]
        mc[m] = o
        
def _gnorm(c1, c2, m, g):
    """
    void gnorm(double *c1, double *c2, int m, const double g)
    {
       double k;
    
       if (g != 0.0) {
          k = 1.0 + g * c1[0];
          for (; m >= 1; m--)
             c2[m] = c1[m] / k;
          c2[0] = pow(k, 1.0 / g);
       } else {
          movem(&c1[1], &c2[1], sizeof(*c1), m);
          c2[0] = exp(c1[0]);
       }
    
       return;
    }
    
    """
    k = 0.0
    if g != 0.0:
        k = 1.0 + g * c1[0]
        for mm in range(m, 0, -1):
            c2[mm] = c1[mm] / k
        c2[0] = np.power(k, 1.0/g)
    else:
        c2[1:] = c1[1:]
        c2[0] = np.exp(c1[0])
    
       
def _gc2gc(c1, m1, g1, c2, m2, g2):
    """
    void gc2gc(double *c1, const int m1, const double g1, double *c2, const int m2,
               const double g2)
    {
       int i, min, k, mk;
       double ss1, ss2, cc;
       static double *ca = NULL;
       static int size;
    
       if (ca == NULL) {
          ca = dgetmem(m1 + 1);
          size = m1;
       }
       if (m1 > size) {
          free(ca);
          ca = dgetmem(m1 + 1);
          size = m1;
       }
    
       movem(c1, ca, sizeof(*c1), m1 + 1);
    
       c2[0] = ca[0];
       for (i = 1; i <= m2; i++) {
          ss1 = ss2 = 0.0;
          min = (m1 < i) ? m1 : i - 1;
          for (k = 1; k <= min; k++) {
             mk = i - k;
             cc = ca[k] * c2[mk];
             ss2 += k * cc;
             ss1 += mk * cc;
          }
    
          if (i <= m1)
             c2[i] = ca[i] + (g2 * ss2 - g1 * ss1) / i;
          else
             c2[i] = (g2 * ss2 - g1 * ss1) / i;
       }
    
       return;
    }    
    """
    i = k = mk = 0
    ss1 = ss2 = cc = 0.0
    ca = None
    _size = 0
    
    if ca is None:
        ca = np.zeros(m1+1)
        _size = m1
    if m1 > _size:
        ca = np.zeros(m1+1)
        _size = m1   
   
    ca[:] = c1[:]
    for i in range(1, m2+1):
        ss1 = ss2 = 0.0
        if m1 < i:
            _min = m1
        else:
            _min = i - 1
        for k in range(1, _min+1):
            mk = i - k
            cc = ca[k] * c2[mk]
            ss2 += k * cc
            ss1 += mk * cc
        if i <= m1:
            c2[i] = ca[i] + (g2 * ss2 - g1 * ss1) / i
        else:
            c2[i] = (g2 * ss2 - g1 * ss1) / i
            
def _theq():
    pass

def _mv_mul(t, x, y):
    t[0] = x[0] * y[0] + x[1] * y[1]
    t[1] = x[2] * y[0] + x[3] * y[1]
    
def _mm_ml(t, x, y):
    t[0] = x[0] * y[0] + x[1] * y[2]
    t[1] = x[0] * y[1] + x[1] * y[3]
    t[2] = x[2] * y[0] + x[3] * y[2]
    t[3] = x[2] * y[1] + x[3] * y[3]    
 
def _inverse(x, y, eps):
    det = 0.0
    det = y[0] * y[3] - y[1] * y[2]
    if np.abs(det) or np.isnan(det):
        raise Exception('theq() : determinant of the normal matrix is too small!')
    x[0] = y[3] / det
    x[1] = -y[1] / det
    x[2] = -y[2] / det
    x[3] = y[0] / det   

def _crstrns(x, y):
    x[0] = y[3]
    x[1] = y[2]
    x[2] = y[1]
    x[3] = y[0]
    
def _cal_p0(p, r, b, n, eps):
    _inverse(t, r[0], eps)
    s[0] = b[0]
    s[1] = b[n - 1]
    _mv_mul(p[0], t, s)
    
"""




"""
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
    
    