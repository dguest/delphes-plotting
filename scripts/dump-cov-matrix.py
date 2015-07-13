#!/usr/bin/env python2
from __future__ import print_function

import ROOT
import sys
import re, json

def _cov_as_list(cov_matrix):
    n_rows = cov_matrix.GetNrows()
    n_cols = cov_matrix.GetNcols()
    assert (n_rows, n_cols) == (5,5)
    out = []
    for rown in xrange(n_rows):
        row = []
        for coln in xrange(n_cols):
            row.append(cov_matrix(rown, coln))
        out.append(row)
    return out

if __name__ == '__main__':
    regtail = "_ptbin([0-9]{2})_etabin([0-9]{2})"
    covre = re.compile('^covmat' + regtail)
    meanre = re.compile('^meanvec' + regtail)
    rfile = ROOT.TFile(sys.argv[1], "read")
    cov_dict = {}
    for key in ROOT.TIter(rfile.GetListOfKeys()):
        name = key.GetName()
        cov_match = covre.match(name)
        mean_match = meanre.match(name)
        if cov_match:
            ptbin, etabin = cov_match.group(1,2)
            obj = key.ReadObj()
            covlist = _cov_as_list(obj)
            pt_dict = cov_dict.setdefault(int(ptbin), {})
            # print(ptbin, len(pt_dict), len(cov_dict))
            pt_dict[int(etabin)] = covlist

    # dump the info
    print(json.dumps(cov_dict, indent=2, sort_keys=True))
