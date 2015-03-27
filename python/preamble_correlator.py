#!/usr/bin/env python

from gnuradio import gr, filter, blocks
from gnuradio.filter import firdes
import mapper

class preamble_correlator(gr.hier_block2):   
    def __init__(self, sps=2.0, rolloff=0.35, preamble=[0,0,0,0,0,0,1,1,0,1,1,0,1,1,0,0,1,1,1,1,0,0,0,0],
                    modtype=mapper.QPSK, greymap=[0,1,3,2] ):
        gr.hier_block2.__init__(self, "preamble_correlator",
                                    gr.io_signature(1,1,gr.sizeof_gr_complex),
                                    gr.io_signature(1,1,gr.sizeof_float))
                                    #gr.io_signature(0,0,0))

        # vet preamble bits
        for b in preamble:
            assert(b >= 0 and b<=1);

        tb = gr.top_block();
        vs = blocks.vector_source_b( preamble );
        mp = mapper.mapper(modtype, greymap);  
        it = filter.interp_fir_filter_ccf(2, firdes.root_raised_cosine(1, 1.0, 1.0/sps, rolloff, 21))
        vk = blocks.vector_sink_c();
        tb.connect(vs,mp,it,vk);
        tb.run();
        self.taps = list(vk.data());
        self.taps.reverse();
        self.taps = map(lambda x: x.conjugate(), self.taps);

        self.flt = filter.fft_filter_ccc(1, self.taps);
        self.mag = blocks.complex_to_mag_squared();
        self.connect(self, self.flt, self.mag);

        # connect output
        self.connect(self.mag, self);



a = preamble_correlator();

