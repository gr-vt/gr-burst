#!/usr/bin/env python

from gnuradio import gr
import pmt, numpy, time, pprint, math
# import scipy.io as sio      # just for debugging, remove when tested

class synchronizer_v3(gr.sync_block):
    def __init__(self, sps = 2, Fs = 100000.0):
        gr.sync_block.__init__(self,"synchronizer_v3",[],[])
        self.message_port_register_in(pmt.intern("cpdus"))
        self.message_port_register_out(pmt.intern("cpdus"))
        
        self.set_msg_handler(pmt.intern("cpdus"), self.handler)
        
        self.sqrtTwoDivTwo = math.sqrt(2)/2

        self.sps = float(sps)    # samples per symbol
        assert(self.sps == 2)
        self.Fs = Fs
        
        self.preSyms = self.genPreamble()
        self.preSyms_x2 = self.genPreamble_x2()
        
        self.burstidx = 0;
        self.cfo = 0;
#         self.debugMode = False
        self.task_times = {};
    
#     def enableDebugMode(self, debugFilename):
#         self.debugMode = True
#         self.debugFilename = debugFilename
    
    def qpskBurstCFOCorrect(self, x, Fs):
        # the expected input x should be a numpy.array
        # this function quadruples the signal, takes the FFT, and determines
        # the baud-rate lines to correct the CFO that may exist

        x_pow = numpy.power(x, 4)
        X4 = numpy.fft.fftshift(numpy.absolute(numpy.fft.fft(x_pow)))
        f = numpy.linspace(-Fs/2,Fs/2,num=len(x_pow))
        maxIdx = numpy.argmax(X4)
        
        cfoEstimate = f[maxIdx]/4.0
        self.cfo = cfoEstimate;
        #print "cfo_estimate: %f"%(cfoEstimate);

        # perform the frequency correction
        rot = 1.0+0j;
        rdiff = numpy.exp(-1j*2*numpy.pi*cfoEstimate/float(Fs))
        y = numpy.multiply(x, numpy.cumprod( numpy.ones(x.size)*rdiff ))
         
        return y
    
    def determineOptimalFilter(self, x):
        # performs the direct solution (Wiener-Hopf solution)
        # to determine optimal filter weights for equalization (optimal w/ respect to MMSE)
        # this function expects data @ 1sps, not 2 sps
        # if the data is available @ a higher samples/sym than 1,
        # it is OK to decimate the data and pass it into this block because
        # the equalization should take care of some of the timing.  Of course,
        # if the sps is 4 or greater, something smarter could probably be done
        # than just throwing data away, so think about that
        numTrainingSyms = len(self.preSyms) 
        x_n = x[0:numTrainingSyms]       # slice the appropriate preamble data
        
        # generate the input correlation matrix
        m = numTrainingSyms       
        X = numpy.fft.fft(x_n, 128)  # the FFT Size is 128 - this will change if the preamble length changes, so beware! 
        X_magSq = numpy.square(numpy.absolute(X))
        rxx = numpy.fft.ifft(X_magSq)
        rxx = rxx/m
        
        toeplitzMatCol = rxx[0:m]
        R = self.toeplitz(toeplitzMatCol)
        
        # generate the P vector
        xc = numpy.correlate(self.preSyms, x_n, 'full')
        P = xc[0:m]
                
        w = numpy.linalg.solve(R,P)
        w /= numpy.amax(numpy.absolute(w))      # scale the filter coefficients
        
        return w
    
    def toeplitz(self, c, r=None):
        """
        *** STOLEN FROM SCIPY SOURCE CODE *** b/c this module cannot have SCIPY dependencies,
        only NUMPY dependencies are allowed
        
        Construct a Toeplitz matrix.
    
        The Toeplitz matrix has constant diagonals, with c as its first column
        and r as its first row.  If r is not given, ``r == conjugate(c)`` is
        assumed.
    
        Parameters
        ----------
        c : array_like
            First column of the matrix.  Whatever the actual shape of `c`, it
            will be converted to a 1-D array.
        r : array_like
            First row of the matrix. If None, ``r = conjugate(c)`` is assumed;
            in this case, if c[0] is real, the result is a Hermitian matrix.
            r[0] is ignored; the first row of the returned matrix is
            ``[c[0], r[1:]]``.  Whatever the actual shape of `r`, it will be
            converted to a 1-D array.
    
        Returns
        -------
        A : (len(c), len(r)) ndarray
            The Toeplitz matrix. Dtype is the same as ``(c[0] + r[0]).dtype``.
    
        See also
        --------
        circulant : circulant matrix
        hankel : Hankel matrix
    
        Notes
        -----
        The behavior when `c` or `r` is a scalar, or when `c` is complex and
        `r` is None, was changed in version 0.8.0.  The behavior in previous
        versions was undocumented and is no longer supported.
    
        Examples
        --------
        >>> from scipy.linalg import toeplitz
        >>> toeplitz([1,2,3], [1,4,5,6])
        array([[1, 4, 5, 6],
               [2, 1, 4, 5],
               [3, 2, 1, 4]])
        >>> toeplitz([1.0, 2+3j, 4-1j])
        array([[ 1.+0.j,  2.-3.j,  4.+1.j],
               [ 2.+3.j,  1.+0.j,  2.-3.j],
               [ 4.-1.j,  2.+3.j,  1.+0.j]])
    
        """
        c = numpy.asarray(c).ravel()
        if r is None:
            r = c.conjugate()
        else:
            r = numpy.asarray(r).ravel()
        # Form a 1D array of values to be used in the matrix, containing a reversed
        # copy of r[1:], followed by c.
        vals = numpy.concatenate((r[-1:0:-1], c))
        a, b = numpy.ogrid[0:len(c), len(r) - 1:-1:-1]
        indx = a + b
        # `indx` is a 2D array of indices into the 1D array `vals`, arranged so
        # that `vals[indx]` is the Toeplitz matrix.
        return vals[indx]
    
    def genPreamble(self):
        return numpy.array([-0.70710677+0.70710677j,  0.70710677+0.70710677j,
                             0.70710677-0.70710677j, -0.70710677-0.70710677j,
                            -0.70710677-0.70710677j, -0.70710677+0.70710677j,
                            -0.70710677+0.70710677j,  0.70710677-0.70710677j,
                            -0.70710677+0.70710677j, -0.70710677-0.70710677j,
                            -0.70710677-0.70710677j,  0.70710677-0.70710677j,
                             0.70710677+0.70710677j,  0.70710677+0.70710677j,
                             0.70710677-0.70710677j, -0.70710677+0.70710677j,
                             0.70710677-0.70710677j, -0.70710677+0.70710677j,
                             0.70710677-0.70710677j, -0.70710677+0.70710677j,
                             0.70710677-0.70710677j, -0.70710677-0.70710677j,
                             0.70710677+0.70710677j,  0.70710677-0.70710677j,
                             0.70710677+0.70710677j,  0.70710677+0.70710677j,
                             0.70710677-0.70710677j, -0.70710677+0.70710677j,
                            -0.70710677+0.70710677j, -0.70710677-0.70710677j,
                            -0.70710677-0.70710677j, -0.70710677+0.70710677j,
                             0.70710677+0.70710677j, -0.70710677-0.70710677j,
                            -0.70710677-0.70710677j,  0.70710677+0.70710677j,
                            -0.70710677-0.70710677j,  0.70710677+0.70710677j,
                             0.70710677+0.70710677j,  0.70710677+0.70710677j,
                            -0.70710677-0.70710677j, -0.70710677-0.70710677j,
                            -0.70710677-0.70710677j, -0.70710677+0.70710677j,
                            -0.70710677+0.70710677j, -0.70710677+0.70710677j,
                             0.70710677-0.70710677j,  0.70710677+0.70710677j], dtype=numpy.complex64)
    
    def genPreamble_x2(self):
        # generated using scipy's interpolate function, but we are writing this module without any
        # scipy dependencies, so hence the hardcoding
        return numpy.array([-0.70710682+0.70710678j, -0.38806983+0.57187419j,
                             0.70710683+0.70710681j,  1.29566346+0.25711846j,
                             0.70710667-0.70710677j, -0.27108426-1.10157359j,
                            -0.70710670-0.70710674j, -0.69809613-0.48355645j,
                            -0.70710682-0.70710674j, -0.73197603-0.41725175j,
                            -0.70710677+0.70710678j, -0.77846512+1.41305459j,
                            -0.70710671+0.70710679j, -0.04705258-0.51324481j,
                             0.70710666-0.70710686j,  0.43727133+0.15147361j,
                            -0.70710682+0.70710686j, -1.24164102+0.18950648j,
                            -0.70710679-0.70710679j, -0.34947090-0.99419741j,
                            -0.70710672-0.70710672j, -0.53795634-0.57596541j,
                             0.70710666-0.70710682j,  1.46415169-0.35750557j,
                             0.70710680+0.70710668j, -0.03159594+1.39756935j,
                             0.70710681+0.70710671j,  1.56957037-0.55650835j,
                             0.70710673-0.70710672j, -0.80945058+0.28088036j,
                            -0.70710678+0.70710675j,  0.56852669-0.10137521j,
                             0.70710691-0.70710671j, -0.41400600-0.04977918j,
                            -0.70710678+0.70710681j,  0.28946681+0.20185674j,
                             0.70710675-0.70710677j, -0.17450830-0.38476532j,
                            -0.70710670+0.70710669j,  0.05508756+0.66894427j,
                             0.70710671-0.70710668j,  0.08827468-1.53255144j,
                            -0.70710683-0.70710676j, -0.30967342+0.567434j  ,
                             0.70710684+0.7071067j ,  1.03661071-0.15496387j,
                             0.70710668-0.7071067j ,  0.56178943-0.24885255j,
                             0.70710684+0.7071067j ,  0.73506473+1.18435336j,
                             0.70710666+0.70710669j,  0.81873886-0.25460533j,
                             0.70710667-0.70710673j,  0.04223688-0.18974672j,
                            -0.70710671+0.70710682j, -0.90893525+1.08450901j,
                            -0.70710677+0.70710676j, -0.62667697-0.03305677j,
                            -0.70710682-0.7071068j , -0.70443785-1.00534061j,
                            -0.70710671-0.70710656j, -0.83524833+0.04208908j,
                            -0.70710677+0.70710682j,  0.00463670+0.91634622j,
                             0.70710676+0.70710672j,  0.43081511+0.13785199j,
                            -0.70710677-0.70710678j, -1.37390293-1.22124058j,
                            -0.70710668-0.70710678j,  0.46924280+0.38897232j,
                             0.70710666+0.70710677j, -0.11851995-0.08858704j,
                            -0.70710671-0.70710674j, -0.20128489-0.18980072j,
                             0.70710682+0.70710676j,  0.97398927+0.92401198j,
                             0.70710671+0.70710667j,  0.64430445+0.73423678j,
                             0.70710684+0.7071067j ,  0.21010288+0.073215j  ,
                            -0.70710683-0.70710674j, -1.08900163-0.88655834j,
                            -0.70710669-0.70710666j, -0.42076342-0.74770238j,
                            -0.70710679-0.70710672j, -0.95984520-0.06683378j,
                            -0.70710672+0.70710683j, -0.45783546+0.87726958j,
                            -0.70710682+0.70710677j, -1.00477580+0.77428944j,
                            -0.70710675+0.70710679j,  0.04994498-0.02322075j,
                             0.70710680-0.70710688j,  0.96787699-0.31263848j,
                             0.70710664+0.70710667j, -0.05751960+1.06877927j], dtype=numpy.complex64)
    
    def qpskFirstOrderPLL(self, x, alpha):
        phiHat = 0;
        #phiHat = 1+0j;
        y = numpy.zeros(len(x),dtype=numpy.complex64)
        for ii in range(0,len(x)):
            # correct w/ estimated phase
            y[ii] = x[ii]*numpy.exp(-1j*phiHat)

            # demodulating circuit
            xHat = (int(math.copysign(1, y[ii].real))<<1 - 1)  + (int(math.copysign(1, y[ii].imag))<<1 - 1)*1j;
            
            # loop filter to update phase estimate
            phiHatT = numpy.multiply(numpy.conjugate(xHat),y[ii])
            phiHatT = phiHatT.imag / phiHatT.real
            if(math.isnan(phiHatT)): phiHatT = 0;
            phiHat = alpha*phiHatT + phiHat

        return y

#     def qpskFirstOrderPLL(self, x, alpha):
#         phiHat = 0;
#         y = numpy.zeros(len(x),dtype=numpy.complex64)
#         for ii in range(0,len(x)):
#             # correct w/ estimated phase
#             y[ii] = x[ii]*numpy.exp(-1j*phiHat)
# 
#             # demodulating circuit
#             xHat = (int(numpy.sign(y[ii].real))<<1 - 1)  + (int(numpy.sign(y[ii].imag))<<1 - 1)*1j
# 
#             # loop filter to update phase estimate
#             phiHatT = numpy.multiply(numpy.conjugate(xHat),y[ii])
#             phiHatT = phiHatT.imag / phiHatT.real
#             phiHat = alpha*phiHatT + phiHat
#         return y

    def start_timer(self):
        self.time = time.time();
 
    def diff_timer(self, task):
        now = time.time()
        diff = now - self.time;
        self.time=now;
        #print "TASK TIME: %s (%f ms)"%(task, diff*1000.0);
        
        if(not self.task_times.has_key(task)):
            self.task_times[task] = (0,0,0);

        nmeas = self.task_times[task][0];
        ttime = self.task_times[task][1];

        self.task_times[task] = (nmeas+1, (ttime+diff), 1000.0*(ttime+diff)/(nmeas+1));
  
    def print_times(self):
        print " ********* AVG TASK TIMES (ms/item) ***********"
        avgtime = {};
        for t in self.task_times.keys():
            avgtime[t] = self.task_times[t][2];
        pprint.pprint(avgtime);

    def handler(self, msg):
        #print "synch starting new sync"
        self.start_timer();

        # get input
        meta = pmt.car(msg);
        samples = pmt.cdr(msg);
        x = numpy.array(pmt.c32vector_elements(samples), dtype=numpy.complex64)
        self.diff_timer("get_input");

        # correct for CFO 
        eqBurst = self.qpskBurstCFOCorrect(x, self.Fs)
        self.diff_timer("cfo")

        # determine preamble start
        preCrossCorr = numpy.absolute(numpy.correlate(self.preSyms_x2, eqBurst, 'full'))
        maxIdx = numpy.argmax(preCrossCorr)
        preambleIdxStart = len(eqBurst) - maxIdx - 1
        eqIn = eqBurst[preambleIdxStart:]
        # decimate the signal to 1sps (right now the constructor enforces teh 2 sps)
        eqIn_decimated = eqIn[0:len(eqIn):2]
        self.diff_timer("timing")

        wOpt = self.determineOptimalFilter(eqIn_decimated)
        self.diff_timer("est_filter")
        
        # filter the input signal w/ the optimal filter (in MMSE sense)
        # todo, do we need to remove the extra syms at the end produced by conv?
        whFilt = numpy.convolve(wOpt, eqIn_decimated)
        whFiltMean = numpy.mean(numpy.absolute(whFilt))
        whFilt /= whFiltMean
        self.diff_timer("apply filt")
        
        # correct any residual phase offset w/ the 1st order PLL
        alpha = 0.002
        phRecoveredSyms = self.qpskFirstOrderPLL(whFilt, alpha)
        self.diff_timer("residual phase")

#         if(self.debugMode):
#             dictToSave = {}
#             dictToSave['gr_eqBurst'] = eqBurst
#             dictToSave['gr_preCrossCorr'] = preCrossCorr
#             dictToSave['gr_eqIn'] = eqIn            
#             dictToSave['gr_wOpt'] = wOpt
#             dictToSave['gr_whFilt'] = whFilt
#             dictToSave['gr_phRecoveredSyms'] = phRecoveredSyms
#             sio.savemat(self.debugFilename, dictToSave)

        # publish equalized symbols
        c = map(lambda i: complex(i), phRecoveredSyms);
        xcm = pmt.init_c32vector(len(c), c)
        #xcm = pmt.init_c32vector(phRecoveredSyms.size, map(lambda i: complex(i), phRecoveredSyms))
        xcm_cpdu = pmt.cons(meta,xcm)
        self.message_port_pub(pmt.intern("cpdus"), xcm_cpdu);
        self.diff_timer("publish result")

        self.burstidx = self.burstidx + 1;
        #self.print_times();

    def work(self, input_items, output_items):
        pass

