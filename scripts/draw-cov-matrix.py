#!/usr/bin/env python3

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.cm import get_cmap

import numpy as np

import argparse
import os
import json
from enum import Enum

_parameter = Enum('parameter', 'd0 a0 phi theta qoverp')

def _get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file')
    parser.add_argument('output_dir', nargs='?', default='test')
    args = parser.parse_args()
    return args

def run():
    args = _get_args()
    with open(args.input_file, 'r') as jfile:
        cov_pars = json.load(jfile)
    d0d0_cov = _get_par_eta_vs_pt(cov_pars, _parameter.d0)
    _print_par_eta_vs_pt(d0d0_cov, output_dir=args.output_dir)

_pt_bins = dict(enumerate([10, 20, 50, 100, 200, 250, 500, 750]))
_eta_bins = dict(enumerate([0.0, 0.4, 0.8, 1.05, 1.5, 1.7, 2.0, 2.25, 2.7]))

def _print_par_eta_vs_pt(cov_array, output_dir):
    fig = Figure(figsize=(5.0, 5.0*3/4))
    canvas = FigureCanvas(fig)
    ax = fig.add_subplot(1,1,1)
    zeroed_pt = np.r_[[0], list(_pt_bins.values())]
    # zeroed_eta = np.r_[[0], list(_eta_bins.values())]
    zeroed_eta = np.array(list(_eta_bins.values()) + [3.0])
    valid_div = cov_array > 0.0
    stdiv_array = np.ones(cov_array.shape) * -1
    stdiv_array[valid_div] = cov_array[valid_div]**0.5
    stdiv_array[~valid_div] = -1
    cmap = get_cmap('jet')
    cmap.set_under('k')
    img = ax.pcolormesh(zeroed_pt, zeroed_eta, stdiv_array.T, vmin=0,
                        cmap=cmap)
    cb = fig.colorbar(img, label='$d_0$ width [mm]', )
    # img.set_under('k')
    ax.set_xscale('log')
    ax.set_xlim(5, zeroed_pt[-1])
    ax.set_xlabel(r'$p_{\rm T}$ [GeV]', x=0.98, ha='right')
    ax.set_ylabel(r'$|\eta|$', y=0.98, ha='right')
    outname = 'eta_vs_pt.pdf'
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    canvas.print_figure(os.path.join(output_dir,outname), bbox_inches='tight')


def _get_par_eta_vs_pt(cov_pars, parameter):
    pt_eta_array = np.ones((len(_pt_bins), len(_eta_bins))) * -1
    par = parameter.value
    for pt_bin, pt_upper in _pt_bins.items():
        for eta_bin, eta_upper in _eta_bins.items():
            spt, seta = str(pt_bin), str(eta_bin)
            pt_eta_list = cov_pars[spt].get(seta)
            if not pt_eta_list: continue
            pt_eta_array[pt_bin, eta_bin] = pt_eta_list[par][par]
    return pt_eta_array


if __name__ == '__main__':
    run()

