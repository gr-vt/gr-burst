#!/usr/bin/env python

from gnuradio import gr;
import pmt;
import sys;
import bitarray, array, numpy

class packer(gr.sync_block):
    def __init__(self):
        gr.sync_block.__init__(self,"packer",[],[])
        self.message_port_register_in(pmt.intern("pdus"));
        self.message_port_register_out(pmt.intern("pdus"));
        self.set_msg_handler(pmt.intern("pdus"), self.handler);

    def handler(self, msg):
        meta = pmt.car(msg);
        x = pmt.to_python(pmt.cdr(msg))
        
        # put bits into bitarray ..
        ba = bitarray.bitarray();          
        for i in x:
            ba.append(i)

        # send packed
        vec2 = numpy.frombuffer( ba.tobytes(), dtype='uint8' )
        self.message_port_pub(pmt.intern("pdus"), 
                pmt.cons(meta, pmt.to_pmt(vec2))
                    )

    def work(self, input_items, output_items):
        pass



