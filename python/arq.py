#!/usr/bin/env python

from gnuradio import gr;
import time, threading, random, pmt, sys, pprint, bitarray, array, struct, binascii, math

class arq(gr.sync_block):
    def __init__(self, timeout=2.5, max_attempts=3):
        gr.sync_block.__init__(self,"arq",[],[])
        self.message_port_register_in(pmt.intern("rx_in"));
        self.message_port_register_in(pmt.intern("tx_in"));
        self.message_port_register_out(pmt.intern("rx_out"));
        self.message_port_register_out(pmt.intern("tx_out"));
        self.set_msg_handler(pmt.intern("tx_in"), self.tx_handler);
        self.set_msg_handler(pmt.intern("rx_in"), self.rx_handler);
        self.timeout = timeout
        self.max_attempts = max_attempts
        self.check_interval = 0.25;

        self.tx_cnt = 0;
        self.rx_cnt = 0;
        self.rtx_cnt = 0;
        self.ack_tx_cnt = 0;
        self.ack_rx_cnt = 0;
        self.max_at = 0;
       
        # pick a starting packet index
        #self.pkt_idx = random.randint(0,100000);
        self.pkt_idx = 0;

        # op codes
        self.operations = {
                "data":[0,0,0,0,0,0,0,0],
                "ack":[0,0,0,0,0,0,0,1],
                };

        self.retrans = {};
        self.rx_record = {};

    def run(self):
        while True:
            print "Rx: %d  Tx: %d  ReTx: %d  AckTx: %d AckRx: %d MaxAttempts: %d Watching: %d"% (
                     self.rx_cnt, self.tx_cnt, self.rtx_cnt, 
                     self.ack_tx_cnt, self.ack_rx_cnt, self.max_at,
                     len(self.rx_record.keys()));

            now = time.time();
            # send retransmissions after non-ack for a while
            for seq in self.retrans.keys():
                if(now - self.retrans[seq]["time_last"]  > self.timeout):
                    print "timeout on seq %d"%(seq)
                    if(self.retrans[seq]["attempts"] > self.max_attempts):
                        print "******** WARNING *********"
                        print " Tx Pkt Seq = %d expired after max attempts!"%(seq);
                        del self.retrans[seq];
                        self.max_at = self.max_at + 1;
                    else:
                        self.rtx_cnt = self.rtx_cnt + 1;
                        self.send_data(seq, self.retrans[seq]["pdu"]);
                        self.retrans[seq]["attempts"] = self.retrans[seq]["attempts"] + 1;
                        self.retrans[seq]["last_time"] = time.time();

            # remove receive records after a while 
            for seq in self.rx_record.keys():
                if( now - self.rx_record[seq]  > self.timeout * self.max_attempts):
                    del self.rx_record[seq];

            time.sleep(self.check_interval);

        

    def start(self):
        self.finished=False;
        self.thread = threading.Thread(target=self.run);
        self.thread.start();

    def stop(self):
        self.finished=True;
        self.thread.join();

    def unpack_bytes(self, b):
        ba = bitarray.bitarray();
        ba.frombytes(b);
        return map(lambda x: int(x), ba.tolist());

    def pack_bits(self, b):
        b = bitarray.bitarray(b);
        return b.tobytes();

    def tx_handler(self, msg):
        self.tx_cnt = self.tx_cnt + 1;
        if(len(self.retrans) > 1000):
            print "WARNING!!!!!!!!! DISCARDING PACKET BECAUSE RE-TX QUEUE IS TOO LONG!"
            return
        meta = pmt.car(msg);
        data_in = array.array('B', pmt.u8vector_elements(pmt.cdr(msg)))
        now = time.time();
        self.retrans[self.pkt_idx] = {"time_orig":now, 
                                      "time_last":now,
                                      "attempts":0, 
                                      "pdu":(meta,data_in)};
        self.send_data(self.pkt_idx,(meta,data_in));
        self.pkt_idx = self.pkt_idx + 1;

    def rx_handler(self, msg):
        meta = pmt.car(msg);
        data_in = array.array('B', pmt.u8vector_elements(pmt.cdr(msg)))
        data_list = data_in.tolist();
        if data_list[0:8] == self.operations["data"]:
            self.rx_cnt = self.rx_cnt + 1;
            seq = struct.unpack("<i", self.pack_bits(data_list[8:8+32]))[0];
            #print "rx sequence: %d"%(seq);
            # send ACK bak
            self.send_ack(seq);

            # pass along sequence number for fun
            if(pmt.is_null(meta)):
                meta = pmt.make_dict();
            meta = pmt.dict_add(meta, pmt.intern("arq_seq"), pmt.from_long(seq));

            if(self.rx_record.has_key(seq)):
                print "duplicate recieve data pkt seq! %d"%(seq)
                return;
            self.rx_record[seq] = time.time();

            # send payload to next layer
            data_upper = data_list[8+32:];
            self.message_port_pub(pmt.intern("rx_out"), pmt.cons(meta,pmt.init_u8vector(len(data_upper),data_upper)));
            return;

        if data_list[0:8] == self.operations["ack"]:
            self.ack_rx_cnt = self.ack_rx_cnt + 1;
            # set pkt ack'd locally
            seq = struct.unpack("<i", self.pack_bits(data_list[8:8+32]))[0];
            print "got pkt ack (%d)"%(seq)
            self.ack(seq);
            return;           

        # fall through fail
        print "got invalid ARQ header, discarding!"

    def send_data(self, pkt_idx, (meta,data)):
        print "SENDING DATA %d"%(pkt_idx);
        data = self.operations["data"] + self.unpack_bytes(struct.pack("<i", pkt_idx)) + data.tolist();
        self.message_port_pub(pmt.intern("tx_out"), pmt.cons(meta,pmt.init_u8vector(len(data),data)));

    # send selective ack message
    def send_ack(self, seq):
        print "SENDING ACK: %d"%(seq)
        self.ack_tx_cnt = self.ack_tx_cnt + 1;
        pkt = self.operations["ack"] + self.unpack_bytes(struct.pack("<i", seq));
        self.message_port_pub(pmt.intern("tx_out"), pmt.cons(pmt.PMT_NIL,pmt.init_u8vector(len(pkt),pkt)));

    def ack(self, seq):
        print "setting %d to received"%(seq);
        if(self.retrans.has_key(seq)):
            del self.retrans[seq];
            #print "removed ok"
        else:
            print "ack'ed pkt not found :: duplicate ack?"

    def work(self, input_items, output_items):
        pass



