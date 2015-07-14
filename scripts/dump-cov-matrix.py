#!/usr/bin/env python2
from __future__ import print_function

"""outer dict is pt bin, inner is eta"""


import ROOT
import sys
import re, json, argparse

def _get_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input')
    parser.add_argument('-b', '--bins', help='only dump numbers',
                        action='store_true')
    return parser.parse_args()

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

regtail = "_ptbin([0-9]{2})_etabin([0-9]{2})"
covre = re.compile('^covmat' + regtail)
meanre = re.compile('^meanvec' + regtail)

def _dump_bins(root_file):
    pt_bins = {}
    for key in ROOT.TIter(rfile.GetListOfKeys()):
        cov_match = covre.match(key.GetName())
        if cov_match:
            pt_bin, eta_bin = cov_match.group(1,2)
            eta_bins = pt_bins.setdefault(pt_bin, set())
            eta_bins.add(eta_bin)

    for ptbin in sorted(pt_bins):
        print('pt bin {}, eta: {}'.format(
            ptbin, ', '.join(sorted(pt_bins[ptbin]))))

if __name__ == '__main__':
    args = _get_args()
    rfile = ROOT.TFile(args.input, "read")
    if args.bins:
        _dump_bins(rfile)
        sys.exit()
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
