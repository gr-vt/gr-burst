#!/usr/bin/env python

from gnuradio import gr;
import Queue, time, threading, random, pmt, sys, pprint, bitarray, array, struct, binascii, math

class random_drop(gr.sync_block):
    def __init__(self, p_drop=0.1, delay=1.0):
        gr.sync_block.__init__(self,"random_drop",[],[])
        self.message_port_register_in(pmt.intern("pdus"));
        self.message_port_register_out(pmt.intern("pdus"));
        self.set_msg_handler(pmt.intern("pdus"), self.handler);
        self.p_drop = p_drop;
        self.delay = delay;
        self.q = Queue.Queue();

    def run(self):
        time.sleep(self.delay);
        msg = self.q.get();
        print msg
        self.message_port_pub(pmt.intern("pdus"), msg);

    def handler(self, msg):
        if(random.uniform(0,1) >= self.p_drop):
            #print "RANDOM_DROP: PASSING"
            t = threading.Thread(target=self.run);
            self.q.put(pmt.cons(pmt.car(msg), pmt.cdr(msg)));
            t.start();
        #    , args=(msg)).start();
            #self.message_port_pub(pmt.intern("pdus"), msg);
        else:
            pass
            #print "RANDOM_DROP: DROPPING"

    def work(self, input_items, output_items):
        pass



