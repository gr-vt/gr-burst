#!/usr/bin/env python

from gnuradio import gr;
import numpy, pmt, sys, pprint, bitarray, array, struct, binascii

class tofpdu(gr.sync_block):
    def __init__(self, blocksize=1024):
        gr.sync_block.__init__(self,"binary_tofpdu",[],[])
        self.message_port_register_in(pmt.intern("pdus"));
        self.message_port_register_out(pmt.intern("fpdus"));
        self.set_msg_handler(pmt.intern("pdus"), self.handler);

    def handler(self, msg):
        ba = bitarray.bitarray();
        meta = pmt.car(msg);
        data = pmt.cdr(msg);
        
        f_vec = pmt.to_pmt(  numpy.array(pmt.to_python(data), dtype="float32") * 2 - 1 )
        pdu = pmt.cons(meta, f_vec);

        # send it on its way
        self.message_port_pub(pmt.intern("fpdus"), pdu);

    def work(self, input_items, output_items):
        pass



