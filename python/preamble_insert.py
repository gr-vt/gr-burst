#!/usr/bin/env python

from gnuradio import gr;
import pmt, array

class preamble_insert(gr.sync_block):
    def __init__(self, preamble=[0,0,0,0,0,0,1,1,0,1,1,0,1,1,0,0,1,1,1,1,0,0,0,0]):
        gr.sync_block.__init__(self,"unpacker",[],[])
        self.message_port_register_in(pmt.intern("pdus"));
        self.message_port_register_out(pmt.intern("pdus"));
        self.set_msg_handler(pmt.intern("pdus"), self.handler);
        self.preamble = array.array("B", preamble);

    def handler(self, msg):
        meta = pmt.car(msg);
        data_in = pmt.cdr(msg);
        data = array.array('B', pmt.u8vector_elements(data_in))
        #data = array.array('B', pmt.u8vector_elements(data_in))
        pre_data = self.preamble + data;
#        print pre_data[:100];
        #burst_bits = pmt.to_pmt(pre_data);
        burst_bits = pmt.init_u8vector(len(pre_data), pre_data.tolist());
        pdu = pmt.cons(meta, burst_bits);
        self.message_port_pub(pmt.intern("pdus"), pdu);
    
    def work(self, input_items, output_items):
        pass



