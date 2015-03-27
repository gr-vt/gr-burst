#!/usr/bin/env python

from gnuradio import gr;
import time, pmt, sys, pprint, bitarray, array, struct, binascii


class cpdu_matlab_writer(gr.sync_block):
    def __init__(self):
        gr.sync_block.__init__(self,"cpdu_matlab_writer",[],[])
        self.message_port_register_in(pmt.intern("cpdus"));
        self.set_msg_handler(pmt.intern("cpdus"), self.handler);
        self.idx = 0;
        self.maxidx = 1000;

    def handler(self, msg):
        if(self.idx > self.maxidx):
            return;

        ba = bitarray.bitarray();
        meta = pmt.car(msg);
        cdata = pmt.cdr(msg);
 
        pydict = None;
        try:   
            pydict = pmt.to_python(meta);       
        except:
            pass;

        x = pmt.c32vector_elements(cdata)

 
        #fn = "/tmp/cpdu_burst_%f.txt"%(time.time());       
        fn = "/tmp/cpdu_burst_%f.txt"%(self.idx);       
        print "writing %s"%(fn);
        f = open(fn,"w");
        x = ", ".join( map( lambda x: "%f+%fj"%(x.real, x.imag), x) );
        f.write("x = [%s]\n"%(str(x))); 
#        f.write("meta = %s"%(str(pydict)));
        f.close();
        f = None;

        self.idx = self.idx + 1;

    def work(self, input_items, output_items):
        pass



