#!/usr/bin/env python
from gnuradio import gr
import pmt
import pdulib
import random
class burst_scheduler2(gr.sync_block):
    def __init__(self, fs):
        gr.sync_block.__init__(
            self,
            name="burst_scheduler2",
            in_sig = None,
            out_sig = None)
        self.fs = fs;
        self.nproduced_val = 0;
        self.message_port_register_in(pmt.intern("sched_pdu"));
        self.message_port_register_in(pmt.intern("nproduced"));
        self.message_port_register_out(pmt.intern("sched_pdu"));
        self.set_msg_handler(pmt.intern("sched_pdu"), self.sched_pdu);
        self.set_msg_handler(pmt.intern("nproduced"), self.nproduced);

    def sched_pdu(self, pdu):
        sched_time = (self.nproduced_val + 10000); # pick a time in the future
        sched_time = sched_time - sched_time%5000; # round to nearest slot
        pdu = pdulib.pdu_arg_add(pdu, pmt.intern("event_time"), pmt.from_uint64(sched_time));
        pdu = pdulib.pdu_arg_add(pdu, pmt.intern("interp"), pmt.from_long(8));

        hop_offset = 10000;
        offset = (random.random()*self.fs-self.fs/2);
        offset = round(offset/hop_offset)*hop_offset; # quantize to nearest 1k offset
        pdu = pdulib.pdu_arg_add(pdu, pmt.intern("freq_offset"), pmt.from_double(offset));
        self.message_port_pub(pmt.intern("sched_pdu"), pdu);

    def nproduced(self, produced_pmt):
        nproduced = pmt.to_uint64(produced_pmt);
        self.nproduced_val = nproduced;



