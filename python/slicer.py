#!/usr/bin/env python

from gnuradio import gr;
import numpy, pmt, sys, pprint, bitarray, array, struct, binascii

class slicer(gr.sync_block):
    def __init__(self, blocksize=1024):
        gr.sync_block.__init__(self,"binary_slicer_async",[],[])
        self.message_port_register_in(pmt.intern("fpdus"));
        self.message_port_register_out(pmt.intern("pdus"));
        self.set_msg_handler(pmt.intern("fpdus"), self.handler);

    def handler(self, pdu):
        # grab float vector from pdu
        meta = pmt.car(pdu);
        x = pmt.to_python(pmt.cdr(pdu))

        # convert to uint8 vector of mapped bits
        bindata = numpy.array(map(lambda x: x > 0, x), dtype='uint8')

        # send it on its way
        self.message_port_pub(pmt.intern("pdus"), pmt.cons(meta, pmt.to_pmt(bindata)))

    def work(self, input_items, output_items):
        pass



