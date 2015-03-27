#!/usr/bin/env python

from gnuradio import gr;
import pmt;
import sys;
import bitarray, array

class unpacker(gr.sync_block):
    def __init__(self):
        gr.sync_block.__init__(self,"unpacker",[],[])
        self.message_port_register_in(pmt.intern("packed_pdus"));
        self.message_port_register_out(pmt.intern("unpacked_pdus"));
        self.set_msg_handler(pmt.intern("packed_pdus"), self.handler);

    def handler(self, msg):
        meta = pmt.car(msg);
        packed_data = pmt.cdr(msg);
        
        # convert pmt -> int list -> string -> bit array
        ba = bitarray.bitarray();          
        ba.frombytes(array.array('B', pmt.u8vector_elements(packed_data)).tostring())

        # convert the unpacked bits to a list and back into a pmt u8vector
        unpacked_data = pmt.init_u8vector(len(ba), ba.tolist());

        # make the new pdu
        pdu = pmt.cons(meta, unpacked_data);

        # send it on its way
        self.message_port_pub(pmt.intern("unpacked_pdus"), pdu);

    def work(self, input_items, output_items):
        pass



