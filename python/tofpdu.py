#!/usr/bin/env python

from gnuradio import gr;
import pmt, sys, pprint, bitarray, array, struct, binascii

class tofpdu(gr.sync_block):
    def __init__(self, blocksize=1024):
        gr.sync_block.__init__(self,"binary_slicer_async",[],[])
        self.message_port_register_in(pmt.intern("pdus"));
        self.message_port_register_out(pmt.intern("fpdus"));
        self.set_msg_handler(pmt.intern("pdus"), self.handler);

    def handler(self, msg):
        ba = bitarray.bitarray();
        meta = pmt.car(msg);
        data = pmt.cdr(msg);
        
        # convert pmt -> int list (of bits)
        data = array.array('B', pmt.u8vector_elements(data))
        bindata = map(lambda x: float(x)*2.0-1.0, data);
        bindatap = pmt.init_f32vector(len(bindata), bindata);

        # make the new pdu
        pdu = pmt.cons(meta, bindatap);

        # send it on its way
        self.message_port_pub(pmt.intern("fpdus"), pdu);

    def work(self, input_items, output_items):
        pass



