#!/usr/bin/env python

from gnuradio import gr;
import pmt, sys, pprint, bitarray, array, struct, binascii, numpy, time

# this block does a few things
#
#  1. add a thing header that looks like
#          [ 0x1337:16 | data_len_bits:16 | header_crc:32 ]
#  2. it then inserts the payload bits
#  3. a crc 32 is placed after the payload
#  4. the burst bits are then padded out to a multiple of block size
#
class deframer(gr.sync_block):
    def __init__(self):
        gr.sync_block.__init__(self,"burst_deframer",[],[])
        self.message_port_register_in(pmt.intern("pdus"));
        self.message_port_register_out(pmt.intern("pdus"));
        self.set_msg_handler(pmt.intern("pdus"), self.handler);
        self.npkt = 0;
        self.npkt_ok = 0;
        self.npkt_hok = 0;

    def handler(self, msg):
        meta = pmt.car(msg);
        bits = pmt.cdr(msg);
 
        self.npkt = self.npkt + 1;

        # convert pmt -> int list (of bits)
        data = pmt.u8vector_elements(bits);
        ba = bitarray.bitarray(data);
        datab = ba.tobytes();

#        print map(lambda x: hex(ord(x)), datab[0:4]);
#        print 'Received Len Bits= ' + str(len(datab)*8)

#        print "Rx Packet:" + str(data[0:10]);
#        print "Rx Packet: "+":".join("{:02x}".format(ord(c)) for c in datab[0:8])
        try:
            (prefix, pktlen) = struct.unpack('<hh', datab[0:4]);
            pprint.pprint({"prefix":hex(prefix), "pktlen":pktlen, "pktlen_bytes":pktlen/8});
            if(not (prefix == 0x1337)):
                print "Deframer: BAD PREFIX!"
                return;
    
            # check header crc
            c2 = binascii.crc32(datab[0:4]);
            hcrc = struct.unpack('i', datab[4:8])[0]
    #        print "CRC: " + str((c2,hcrc))
            if not (c2 == hcrc):
                print "Deframer: bad header crc"
                return;
            self.npkt_hok = self.npkt_hok + 1;
    
            # make sure we got enough bits for the given header len
            if(len(data) < (pktlen/8 + 8 + 4)):
                print "Deframer: not enough bits received for full payload!"
                print "Deframer: pktlen field = %d, received = %d\n"%(pktlen, len(data))
                return;
    
            if(not (pktlen % 8 == 0)):
                print "Deframer: payload should be a multiple of 8"
                return
    
            # extract header bytes
            c1 = binascii.crc32(datab[0:8+pktlen/8]);
            payload = datab[8:8+pktlen/8];
            #print "RX Payload len = %d"%(len(payload))
            #print  ":".join("{:02x}".format(ord(c)) for c in datab)
    #        print payload;
            ex_crc2 = datab[(8+pktlen/8):(8+pktlen/8+4)];
        except:
            print 'Not enough data to read! dropping'
            return
        try:
            c1h = struct.unpack('i', ex_crc2)[0]
            #print "rx payload CRC = %d (%s)"%(c1h, ":".join("{:02x}".format(ord(c)) for c in ex_crc2))
        except:
            print "shortened packet length dropping"
            return
#        print "CRC2:" + str((c1, c1h));          
        if(not c1 == c1h):
            print "Failed payload CRC"
            return

#        print "BURST OK!"
        self.npkt_ok = self.npkt_ok + 1;
        pct_ok = 100.0*self.npkt_ok / self.npkt;
        pct_hok = 100.0*self.npkt_hok / self.npkt;
        print "Deframer: Percent ok = %f (%f header)%%"%(pct_ok, pct_hok);

        # send it on its way
        payload = numpy.fromstring(payload, dtype=numpy.uint8);
        v = pmt.to_pmt(payload);
        meta = pmt.dict_add(meta, pmt.intern("timestamp"), pmt.from_double(time.time()));
        meta = pmt.dict_add(meta, pmt.intern("npkt"), pmt.from_long(self.npkt));
        meta = pmt.dict_add(meta, pmt.intern("npkt_hok"), pmt.from_long(self.npkt_hok));
        meta = pmt.dict_add(meta, pmt.intern("npkt_ok"), pmt.from_long(self.npkt_ok));
        meta = pmt.dict_add(meta, pmt.intern("header_pass_rate"), pmt.from_double(pct_hok));
        meta = pmt.dict_add(meta, pmt.intern("payload_pass_rate"), pmt.from_double(pct_ok));
        pdu = pmt.cons( meta, v )
        self.message_port_pub(pmt.intern("pdus"), pdu);

    def work(self, input_items, output_items):
        pass



