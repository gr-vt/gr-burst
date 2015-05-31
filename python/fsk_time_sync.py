#!/usr/bin/env python

from gnuradio import gr;
import pmt,numpy;
from scipy import signal

class fsk_time_sync(gr.sync_block):
    def __init__(self, sample_rate=32e3, sps=None, upsample_rate=100.0, n_offsets = 10.0):
        gr.sync_block.__init__(self,"fsk_synchronizer",[],[])
        self.message_port_register_in(pmt.intern("pdus"));
        self.message_port_register_out(pmt.intern("pdus"));
        self.message_port_register_out(pmt.intern("autocorr"));
        self.message_port_register_out(pmt.intern("timing"));
        self.set_msg_handler(pmt.intern("pdus"), self.handler);
        self.sample_rate = sample_rate
        self.sps = sps

        self.upsample_rate = upsample_rate
        self.sps = sps
        self.nburst = 0;
        self.nburst_ok = 0;
        self.sps_samps = []
        self.n_offsets = n_offsets

    def handler(self, msg):
        # get input
        meta = pmt.car(msg);
        x = pmt.to_python(pmt.cdr(msg))

        if( self.sps == None and len(self.sps_samps) < 10):
            # compute the cross correlation metric first peak (to find baud rate)
            (clen, ncut) = (500,500)
            e = numpy.zeros(clen*2-1)
            for i in range(1,ncut):
                c = x[i:i+clen]
                d = numpy.correlate(c,c, mode='full')
                e += d

            # upsample to xcorr to interpolate fractional sym rate
            e = e[clen-1:]
            e_upsamp = signal.resample(e, self.upsample_rate*len(e))
            e_upsamp = e_upsamp[:len(e_upsamp)/2]
            #e_upsamp = e_upsamp[0:len(e_upsamp)/2]
            self.message_port_pub(pmt.intern("autocorr"), pmt.cons(meta, pmt.to_pmt(e_upsamp)))
    
            # locate first minimum and next peak (need to see how generalizable this is ... )
            firstmin = numpy.argmin(e_upsamp)
            firstmax = firstmin + numpy.argmax(e_upsamp[firstmin:])

            # determine samples per symbol
            sps = firstmax/self.upsample_rate
            self.sps_samps.append(sps)
            self.sps = numpy.mean(self.sps_samps)
        else:
            sps = self.sps

        meta = pmt.dict_add(meta, pmt.intern("meta"), pmt.from_double(sps));
        print "sps = %f"%(sps)
        
        ovf = []
        ovals = {}
        n_offsets = float(self.n_offsets)
        nsyms = (len(x)/sps)-1
        best = 0
        best_syms = None
        for o in numpy.arange(0,sps,sps/n_offsets):
            syms = signal.resample(x[o:o+nsyms*sps], len(x[o:])*10/(sps))
            syms = syms[0:len(syms)-len(syms)%10]
            syms = syms.reshape( [len(syms)/10, 10] )
            #syms = syms[:,4]
            syms = syms[:,3:7].mean(1)
            dist = numpy.mean(numpy.abs(syms))
            if dist > best:
                best = dist
                best_syms = syms
            ovals[o] = dist
            ovf.append(dist)

        # output timing metric (should look sinusoidal-ish)
        self.message_port_pub(pmt.intern("timing"), pmt.cons(meta, pmt.to_pmt(ovf)))
        best_offset = ovals.keys()[ numpy.argmax(ovals.values()) ]
        meta = pmt.dict_add(meta, pmt.intern("tau"), pmt.from_double(best_offset));

        # publish our recovered symbols
        self.message_port_pub(pmt.intern("pdus"), pmt.cons(meta, pmt.to_pmt(best_syms)))

    def work(self, input_items, output_items):
        pass




