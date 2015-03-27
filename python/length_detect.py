#!/usr/bin/env python

from gnuradio import gr;
import pmt,numpy;

class length_detect(gr.sync_block):
    def __init__(self):
        gr.sync_block.__init__(self,"length_detect",[],[])
        self.message_port_register_in(pmt.intern("cpdus"));
        self.message_port_register_out(pmt.intern("cpdus"));
        self.set_msg_handler(pmt.intern("cpdus"), self.handler);

    def handler(self, msg):
        meta = pmt.car(msg);
        samples = pmt.cdr(msg);
        
        x = numpy.array(pmt.c32vector_elements(samples), dtype=numpy.complex64)
        x_2 = numpy.real(x * x.conjugate())
        
        # smoothing filter
        x_2f = numpy.convolve(50*[1], x_2);

        # find max power to compute power thresh
        maxidx = numpy.argmax(x_2f);
        maxpow = x_2f[maxidx];
        thr = maxpow / 16; # 6db down

        # find where we are below thresh
        start_offset = 1000;
        idx = numpy.argmax(x_2f[start_offset:] < thr) + start_offset + 1000;
#        print "below = (%d, %f)"%(idx, x_2f[idx])

        # discard bursts where we dont find an end
        if idx == start_offset:
            print "WARNING: length detect: discarding burst"
            return

        # tack on some metadata
        meta = pmt.dict_add(meta, pmt.intern("est_len"), pmt.from_double(int(idx)));
        meta = pmt.dict_add(meta, pmt.intern("orig_len"), pmt.from_double(len(x)));
        
        # extract the useful signal
        x = x[0:idx];

        # send it on its way
        samples_out = pmt.init_c32vector(len(x), map(lambda i: complex(i), x))
        cpdu = pmt.cons(meta,samples_out)
        self.message_port_pub(pmt.intern("cpdus"), cpdu);

    def work(self, input_items, output_items):
        pass




