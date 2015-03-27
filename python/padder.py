#!/usr/bin/env python

from gnuradio import gr;
import pmt, sys, pprint, bitarray, array, struct, binascii, math

class padder(gr.sync_block):
    def __init__(self, blocksize=144):
        gr.sync_block.__init__(self,"padder",[],[])
        self.message_port_register_in(pmt.intern("pdus"));
        self.message_port_register_out(pmt.intern("pdus"));
        self.set_msg_handler(pmt.intern("pdus"), self.handler);
        self.blocksize = blocksize;

    def handler(self, msg):
        meta = pmt.car(msg);
        packed_data = pmt.cdr(msg);
        data = array.array('B', pmt.u8vector_elements(packed_data))

        # break up into blocks 
        nbits = len(data); # 32 for the payload crc      
        final_len = int(math.ceil( 1.0*nbits / self.blocksize) * self.blocksize);
        padding = final_len - nbits;

        # add some padding
        for i in range(0,padding):
            data.append(0);

        padded_bits = pmt.init_u8vector(len(data), data);
        pdu = pmt.cons(meta, padded_bits);

        # send it on its way
        self.message_port_pub(pmt.intern("pdus"), pdu);

    def work(self, input_items, output_items):
        pass



