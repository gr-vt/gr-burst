#!/usr/bin/env python

from gnuradio import gr;
import pmt, sys, pprint, bitarray, array, struct, binascii, math


# this block does a few things
#
#  1. add a thing header that looks like
#          [ 0x1337:16 | data_len_bits:16 | header_crc:32 ]
#  2. it then inserts the payload bits
#  3. a crc 32 is placed after the payload
#  4. the burst bits are then padded out to a multiple of block size
#
class framer(gr.sync_block):
    def __init__(self, blocksize=1024):
        gr.sync_block.__init__(self,"unpacker",[],[])
        self.message_port_register_in(pmt.intern("packed_pdus"));
        self.message_port_register_out(pmt.intern("unpacked_pdus"));
        self.set_msg_handler(pmt.intern("packed_pdus"), self.handler);
        self.blocksize = blocksize;

    def handler(self, msg):
        ba = bitarray.bitarray();
        meta = pmt.car(msg);
        packed_data = pmt.cdr(msg);
        
        # convert pmt -> int list (of packed bytes)
        data = array.array('B', pmt.u8vector_elements(packed_data))

        # add header on front
        header = struct.pack('hh', 0x1337, 8*len(data));
        ba.frombytes(header);

        # compute header crc
        c2 = binascii.crc32(ba.tobytes());
        hcrc = struct.pack('i', c2);
        ba.frombytes(hcrc);

        # add the data payload
        ba.frombytes(data.tostring())
        
        # compute payload crc
        c2 = binascii.crc32(ba.tobytes());
        pcrc = struct.pack('i', c2);
        ba.frombytes(pcrc);
        
        # convert the unpacked bits to a list and back into a pmt u8vector
        burst_bits = pmt.init_u8vector(len(ba), ba.tolist());
        print "Tx Packet: " + ":".join("{:02x}".format(ord(c)) for c in ba.tobytes()[0:8])
        
        # make the new pdu
        pdu = pmt.cons(meta, burst_bits);

        # send it on its way
        self.message_port_pub(pmt.intern("unpacked_pdus"), pdu);

    def work(self, input_items, output_items):
        pass



