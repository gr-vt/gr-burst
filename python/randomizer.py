#!/usr/bin/env python

from gnuradio import gr;
import time, pmt, sys, pprint, bitarray, array, struct, binascii


class randomizer(gr.sync_block):
    def __init__(self, poly=[0,14,15], init=[1,0,0,1,0,1,0,1,0,0,0,0,0,0,0], gensteps=1000000):
        gr.sync_block.__init__(self,"unpacker",[],[])
        self.message_port_register_in(pmt.intern("pdus"));
        self.message_port_register_out(pmt.intern("pdus"));
        self.set_msg_handler(pmt.intern("pdus"), self.handler);

        # make sure the poly and init match happily
        for p in poly:
            assert(len(init) >= p)

        print "pregenerating randomizer"
        a = time.time();
        self.seq = [];
        for i in range(0,gensteps):
            out = 0;
            for v in poly:              # compute output
                if(v==0):
                    continue;
                out = out ^ init[v-1];
            init = [out] + init[:-1];   # shift       
            self.seq.append(out);
        b = time.time();
        print "pregenerating randomizer done (%f seconds)"%(b-a);
        #print self.seq;

    def handler(self, msg):
        ba = bitarray.bitarray();
        meta = pmt.car(msg);
        bits_in = pmt.cdr(msg);
       
        # convert pmt -> int list (of bits)
        data = array.array('B', pmt.u8vector_elements(bits_in))

        # randomize!
        rndata = map(lambda x,y: (x&0x01)^y, data, self.seq[0:len(data)]);
#        print "randomizer in %d, out %d"%(len(pmt.u8vector_elements(bits_in)), len(rndata))
        
        # convert the unpacked bits to a list and back into a pmt u8vector
        burst_bits = pmt.init_u8vector(len(rndata), rndata);

        # make the new pdu
        pdu = pmt.cons(meta, burst_bits);

        # send it on its way
        self.message_port_pub(pmt.intern("pdus"), pdu);

    def work(self, input_items, output_items):
        pass



