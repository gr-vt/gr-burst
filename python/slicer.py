#!/usr/bin/env python

from gnuradio import gr;
import pmt, sys, pprint, bitarray, array, struct, binascii

class slicer(gr.sync_block):
    def __init__(self, blocksize=1024):
        gr.sync_block.__init__(self,"binary_slicer_async",[],[])
        self.message_port_register_in(pmt.intern("fpdus"));
        self.message_port_register_out(pmt.intern("pdus"));
        self.set_msg_handler(pmt.intern("fpdus"), self.handler);

    def handler(self, msg):
        ba = bitarray.bitarray();
        meta = pmt.car(msg);
        data = pmt.cdr(msg);
        
        # convert pmt -> int list (of bits)
        data = array.array('f', pmt.f32vector_elements(data))
        bindata = map(lambda x: x > 0, data);
        bindatap = pmt.init_u8vector(len(bindata), bindata);

        # make the new pdu
        pdu = pmt.cons(meta, bindatap);

        # send it on its way
        self.message_port_pub(pmt.intern("pdus"), pdu);

    def work(self, input_items, output_items):
        pass



