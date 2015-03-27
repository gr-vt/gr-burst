#!/usr/bin/env python
import pmt;
def pdu_arg_add(pdu, k, v):
    meta = pmt.car(pdu);
    data = pmt.cdr(pdu);
    if(pmt.is_null(meta)):
        meta = pmt.make_dict();
    assert(pmt.is_dict(meta));
    meta = pmt.dict_add(meta, k, v);
    return pmt.cons(meta,data);

