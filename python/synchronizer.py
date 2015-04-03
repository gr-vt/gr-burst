#!/usr/bin/env python

from gnuradio import gr;
import pmt,numpy;
from scipy import signal

# timing algorithm from 00585803.pdf
class synchronizer(gr.sync_block):
    def __init__(self, sps = 2):
        gr.sync_block.__init__(self,"synchronizer",[],[])
        self.message_port_register_in(pmt.intern("cpdus"));
        self.message_port_register_out(pmt.intern("cpdus"));
        
        self.message_port_register_out(pmt.intern("timing_metric"));
        self.message_port_register_out(pmt.intern("sym_timed"));

        self.message_port_register_out(pmt.intern("phase_ramp"));

        self.set_msg_handler(pmt.intern("cpdus"), self.handler);

        self.T = float(sps);    # samples per symbol
        assert(self.T == 2);
        #self.N = 8;            # interpolation factor
        #self.N = 32;            # interpolation factor
        self.N = 16;            # interpolation factor
        #self.N = 4;            # interpolation factor
        self.O = 4;             # power for freq corr

        self.nburst = 0;
        self.nburst_ok = 0;

    def handler(self, msg):
        # get input
        meta = pmt.car(msg);
        samples = pmt.cdr(msg);
        x = numpy.array(pmt.c32vector_elements(samples), dtype=numpy.complex64)
 
        # upsample and normalize power
        xi = signal.resample(x, len(x)* (self.N / self.T));
        # compute the symbol timing
        xt = numpy.real(xi*xi.conjugate()) * numpy.exp( (-1j*2.0*numpy.pi/self.N) * numpy.linspace(1,len(xi),len(xi)) );
        s = numpy.sum(x);    
        tau = (-self.T/(2*numpy.pi)) * numpy.arctan2(numpy.imag(s), numpy.real(s));
        
        # publish timing metric for debugging
        tm = pmt.init_c32vector(xt.size, map(lambda i: complex(i), xt))
        tm_cpdu = pmt.cons(meta,tm)
        self.message_port_pub(pmt.intern("timing_metric"), tm_cpdu);

        # extract symbols
        offset = round(self.N*tau/self.T);
        fo = (offset + self.N)%self.N;
        sym = xi[fo:-1:self.N];
        
        # normalize power to 1
        sym = sym / numpy.mean(numpy.real(sym * sym.conjugate()));

        # publish timing correct symbols (with frequency offset)
        sm = pmt.init_c32vector(sym.size, map(lambda i: complex(i), sym))
        sm_cpdu = pmt.cons(meta,sm)
        self.message_port_pub(pmt.intern("sym_timed"), sm_cpdu);

        # compute symbol frequency offset (linear phase offset within block)
        x_n = numpy.power(sym[200:1000], self.O);
        phase_ramp = numpy.unwrap(numpy.angle( x_n ));
        f_off_O = numpy.mean(numpy.diff(phase_ramp));
        goodstat = numpy.std(numpy.diff(phase_ramp));
        f_off = f_off_O / self.O;

        # check percentages
        self.nburst = self.nburst + 1;
        if(goodstat < 1.0):
            self.nburst_ok = self.nburst_ok + 1;
        else:
            print "WARNING: feedforward synchronizer discarding burst, goodness metric %f < 1.0 (likely poor timing recovery occurred, the CFO phase ramp looks like garbage)"%(goodstat)
            return
        print "sync: "+str((goodstat, self.nburst, self.nburst_ok, self.nburst_ok*100.0 / self.nburst));       

        # export phase ramp
        pr = pmt.init_f32vector(phase_ramp.size, map(lambda i: float(i), phase_ramp))
        pr_fpdu = pmt.cons(meta,pr)
        self.message_port_pub(pmt.intern("phase_ramp"), pr_fpdu);

        # apply frequency offset correction
        xc = numpy.multiply(sym, numpy.exp(-1j * f_off * numpy.linspace(1,sym.size,sym.size)));
        
        # compute and correct static symbol phase offset
        xcp = numpy.power(xc[400:1000], self.O);

# linear mean
        theta = numpy.mean( numpy.angle( xcp ) ) / self.O + numpy.pi/4;
# weighted mean
#        theta = numpy.sum(numpy.angle(xcp) * numpy.abs(xcp)) / numpy.sum(numpy.abs(xcp));
#        theta = theta / self.O + numpy.pi/4;

        xc = xc * numpy.exp(-1j*theta);

        # show time, frequency and phase estimates
        #print "tau = %f, f_off = %f, theta = %f"%(tau, f_off, theta);

        # add our estimates to the metadata dictionary
        meta = pmt.dict_add(meta, pmt.intern("tau"), pmt.from_double(tau));
        meta = pmt.dict_add(meta, pmt.intern("f_off"), pmt.from_double(f_off));
        meta = pmt.dict_add(meta, pmt.intern("theta"), pmt.from_double(theta));

        # publish freq corrected symbols
        xcm = pmt.init_c32vector(xc.size, map(lambda i: complex(i), xc))
        xcm_cpdu = pmt.cons(meta,xcm)
        self.message_port_pub(pmt.intern("cpdus"), xcm_cpdu);


    def work(self, input_items, output_items):
        pass




