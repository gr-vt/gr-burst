#!/usr/bin/env python

from gnuradio import gr
import pmt,numpy
import math
from scipy import signal, linalg
import scipy.io as sio

class synchronizer_v2(gr.sync_block):
    def __init__(self, sps = 2, Fs = 100000.0):
        gr.sync_block.__init__(self,"synchronizer_v2",[],[])
        self.message_port_register_in(pmt.intern("cpdus"))
        self.message_port_register_out(pmt.intern("cpdus"))
        
        self.set_msg_handler(pmt.intern("cpdus"), self.handler)
        
        self.sqrtTwoDivTwo = math.sqrt(2)/2

        self.sps = float(sps)    # samples per symbol
        assert(self.sps == 2)
        self.Fs = Fs
        
        self.preSyms = self.genPreamble()
        self.preSyms_x2 = signal.resample(self.preSyms, len(self.preSyms)*2)      # interpolate by 2x
        
        self.debugMode = False
    
    def enableDebugMode(self, debugFilename):
        self.debugMode = True
        self.debugFilename = debugFilename
    
    def qpskBurstCFOCorrect(self, x, Fs):
        # the expected input x should be a numpy.array
        # this function quadruples the signal, takes the FFT, and determines
        # the baud-rate lines to correct the CFO that may exist
        x_pow = numpy.power(x, 4)
        X4 = numpy.fft.fftshift(numpy.absolute(numpy.fft.fft(x_pow)))
        f = numpy.linspace(-Fs/2,Fs/2,num=len(x_pow))
        maxIdx = numpy.argmax(X4)
        
        cfoEstimate = f[maxIdx]/4.0
        
        # perform the frequency correction
        t0 = numpy.arange(0,len(x),1)/float(Fs)
        freqCorrVector = numpy.exp(-1j*2*numpy.pi*cfoEstimate*t0)
        y = numpy.multiply(x, freqCorrVector)
        
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
        R = linalg.toeplitz(toeplitzMatCol)
        
        # generate the P vector
        xc = numpy.correlate(self.preSyms, x_n, 'full')
        P = xc[0:m]
                
        w = numpy.linalg.solve(R,P)
        w /= numpy.amax(numpy.absolute(w))      # scale the filter coefficients
        
        return w
    
    def genPreamble(self):
        # TODO: should probably pass in the preamble as a parameter huh?
        preSeqBits = [0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, \
                      0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, \
                      0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0]
        y = self.qpskModulate(preSeqBits)
        return y
    
    def qpskModulate(self, x):
        # modulates a bitstream to QPSK symbols according to the 
        # mapping that was chosen for thsi project
        y = numpy.empty(len(x)/2, dtype=numpy.complex64)
        for ii in range(0,len(y)):
            idxLo = ii*2
            idxHi = idxLo+1 
            bitGrp = x[idxLo:idxHi+1]       # need to add the +1 in python
            if(bitGrp[0]==0 and bitGrp[1]==0):
                # [0 0]
                y[ii] = self.sqrtTwoDivTwo + 1j*self.sqrtTwoDivTwo
            elif(bitGrp[0]==0 and bitGrp[1]==1):
                # [0 1]
                y[ii] = -1*self.sqrtTwoDivTwo + 1j*self.sqrtTwoDivTwo
            elif(bitGrp[0]==1 and bitGrp[1]==0):
                # [1 0]
                y[ii] = self.sqrtTwoDivTwo - 1j*self.sqrtTwoDivTwo
            else:
                # [1 1]
                y[ii] = -1*self.sqrtTwoDivTwo - 1j*self.sqrtTwoDivTwo
                
        return y
            
    
    def qpskFirstOrderPLL(self, x, alpha):
        phiHat = 0
        y = numpy.zeros(len(x),dtype=numpy.complex64)
        for ii in range(0,len(x)):
            # correct w/ estimated phase
            y[ii] = x[ii]*numpy.exp(-1j*phiHat)
            
            # demodulating circuit
            if(y[ii].real>=0 and y[ii].imag>=0):
                xHat = numpy.exp(1j*numpy.pi/2)
            elif(y[ii].real>=0 and y[ii].imag<0):
                xHat = numpy.exp(1j*3*numpy.pi/2)
            elif(y[ii].real<0 and y[ii].imag<0):
                xHat = numpy.exp(1j*5*numpy.pi/2)
            else:
                xHat = numpy.exp(1j*7*numpy.pi/2)
            
            # loop filter to update phase estimate
            phiHatT = numpy.angle(numpy.multiply(numpy.conjugate(xHat),y[ii]))
            phiHat = alpha*phiHatT + phiHat
            
        return y
        
    def handler(self, msg):
        # get input
        meta = pmt.car(msg);
        samples = pmt.cdr(msg);
        x = numpy.array(pmt.c32vector_elements(samples), dtype=numpy.complex64)
 
        # correct for CFO 
        eqBurst = self.qpskBurstCFOCorrect(x, self.Fs)
        
        # determine preamble start
        preCrossCorr = numpy.absolute(numpy.correlate(self.preSyms_x2, eqBurst, 'full'))
        maxIdx = numpy.argmax(preCrossCorr)
        preambleIdxStart = len(eqBurst) - maxIdx - 1
        eqIn = eqBurst[preambleIdxStart:]
        # decimate the signal to 1sps (right now the constructor enforces teh 2 sps)
        eqIn_decimated = eqIn[0:len(eqIn):2]

        wOpt = self.determineOptimalFilter(eqIn_decimated)
        
        # filter the input signal w/ the optimal filter (in MMSE sense)
        whFilt = signal.lfilter(wOpt, [1], eqIn_decimated)
        whFiltMean = numpy.mean(numpy.absolute(whFilt))
        whFilt /= whFiltMean
        
        # correct any residual phase offset w/ the 1st order PLL
        alpha = 0.002
        phRecoveredSyms = self.qpskFirstOrderPLL(whFilt, alpha)

        if(self.debugMode):
            dictToSave = {}
            dictToSave['gr_eqBurst'] = eqBurst
            dictToSave['gr_preCrossCorr'] = preCrossCorr
            dictToSave['gr_eqIn'] = eqIn            
            dictToSave['gr_wOpt'] = wOpt
            dictToSave['gr_whFilt'] = whFilt
            dictToSave['gr_phRecoveredSyms'] = phRecoveredSyms
            sio.savemat(self.debugFilename, dictToSave)

        # publish equalized symbols
        xcm = pmt.init_c32vector(phRecoveredSyms.size, map(lambda i: complex(i), phRecoveredSyms))
        xcm_cpdu = pmt.cons(meta,xcm)
        self.message_port_pub(pmt.intern("cpdus"), xcm_cpdu);


    def work(self, input_items, output_items):
        pass

