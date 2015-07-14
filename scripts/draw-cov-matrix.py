#!/usr/bin/env python3

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.cm import get_cmap
from matplotlib.colors import SymLogNorm, LogNorm
from matplotlib.ticker import FuncFormatter, LogFormatterMathtext


import numpy as np

import math
import argparse
import os
import json
import warnings
from enum import Enum

class _trkpar(Enum):
    d0 = (0, r'$d_0$', 'mm')
    z0 = (1, r'$z_0$', 'mm')
    phi= (2, r'$\phi$', '')
    theta = (3, r'$\theta$', '')
    qoverp = (4, r'$q/p$', r'MeV$^{-1}$')
    def __init__(self, index, string, units):
        self.index = index
        self.string = string
        self.units = units
    # def __len__(self):
    #     return len(__class__.__members__)


def _get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file')
    parser.add_argument('output_dir', nargs='?', default='test')
    args = parser.parse_args()
    return args

def run():
    args = _get_args()
    with open(args.input_file, 'r') as jfile:
        cov_pars = _int_keys(json.load(jfile))
    _print_cor_matrix(cov_pars, 0, 0, output_dir=args.output_dir)
    for par in [_trkpar.d0, _trkpar.z0, _trkpar.phi, _trkpar.theta, _trkpar.qoverp]:
        cov = _get_par_eta_vs_pt(cov_pars, par)
        _print_par_eta_vs_pt(cov, output_dir=args.output_dir,
                             parameter=par)

_pt_bins = dict(enumerate([10, 20, 50, 100, 200, 250, 500, 750]))
_eta_bins = dict(enumerate([0.0, 0.4, 0.8, 1.05, 1.5, 1.7, 2.0, 2.25, 2.7]))

def _int_keys(json_dict):
    """recursive translation to int keyed dict"""
    try:
        odict = {int(k): _int_keys(v) for k,v in json_dict.items()}
    except AttributeError:
        return json_dict
    return odict

def _label_ticks(axis, labels=_trkpar):
    tickpos = np.linspace(0, len(labels) - 1, len(labels))
    axis.set_ticks(tickpos)
    axis.set_ticklabels([x.string for x in labels])

def _print_cor_matrix(matlist, pt_bin, eta_bin, output_dir):
    ar = np.asarray(matlist[pt_bin][eta_bin])
    invstdev = np.matrix(np.diag(np.diag(ar)**(-0.5)))
    cov = invstdev * ar * invstdev
    fig = Figure(figsize=(5.0, 5.0*3/4))
    canvas = FigureCanvas(fig)
    ax = fig.add_subplot(1,1,1)
    cmap = get_cmap('RdBu')
    lin_order = 3               # make log flat after this
    norm = SymLogNorm(linthresh=10**(-lin_order), vmin=-1.001, vmax=1-0.001)
    img = ax.matshow(np.asarray(cov), norm=norm, cmap=cmap)
    # cmap.set_over('k')
    _label_ticks(ax.xaxis)
    _label_ticks(ax.yaxis)
    ticks = 0.1**np.r_[lin_order:0:-1]
    full_ticks = (np.r_[-10:11].T * ticks[:,None]).flatten()
    # formatter = LogFormatterMathtext()
    formatter = FuncFormatter(_log_formatting)
    cb = fig.colorbar(img, label='Correlation Coefficient', ticks=full_ticks,
                      format=formatter)
    outname = 'cor_mat_ptbin{}_etabin{}.pdf'.format(pt_bin, eta_bin)
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    canvas.print_figure(os.path.join(output_dir,outname), bbox_inches='tight')


def _print_par_eta_vs_pt(cov_array, output_dir, parameter):
    fig = Figure(figsize=(5.0, 5.0*3/4))
    canvas = FigureCanvas(fig)
    ax = fig.add_subplot(1,1,1)
    # WARNING: right now I don't know if these bins are correct...
    warnings.warn('we need need confirmation from Shih-Chieh that these bins '
                  'are correct', stacklevel=2)
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
    lab = parameter.string + ' width'
    units = parameter.units
    if units:
        lab += ' [{}]'.format(units)
    cb = fig.colorbar(img, label=lab)
    # img.set_under('k')
    ax.set_xscale('log')
    ax.set_xlim(5, zeroed_pt[-1])
    ax.set_xlabel(r'$p_{\rm T}$ [GeV]', x=0.98, ha='right')
    ax.set_ylabel(r'$|\eta|$', y=0.98, ha='right')
    outname = '{}_eta_vs_pt.pdf'.format(parameter.name)
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    canvas.print_figure(os.path.join(output_dir,outname), bbox_inches='tight')


def _get_par_eta_vs_pt(cov_pars, parameter):
    pt_eta_array = np.ones((len(_pt_bins), len(_eta_bins))) * -1
    par = parameter.index
    for pt_bin, pt_upper in _pt_bins.items():
        for eta_bin, eta_upper in _eta_bins.items():
            # spt, seta = str(pt_bin), str(eta_bin)
            pt_eta_list = cov_pars[pt_bin].get(eta_bin)
            if not pt_eta_list: continue
            pt_eta_array[pt_bin, eta_bin] = pt_eta_list[par][par]
    return pt_eta_array

# ________________________________________________________________________

def _log_formatting(value, pos):
    if value == 0:
        return r'0'
    roundlog = round(math.log(abs(value), 10))
    roundval = 10**roundlog
    if not np.isclose(roundval, abs(value)):
        return ''
    if np.isclose(roundval, 0.1):
        return r'0.1'
    if roundval == 1:
        base = math.copysign(1, value)
        exp = ''
    elif roundval == 10:
        base = math.copysign(10, value)
        exp = ''
    else:
        base = math.copysign(10, value)
        exp = round(math.log(abs(value),10))
    return r'{:.0f}$^{{\mathdefault{{ {} }} }}$'.format(base, exp)

if __name__ == '__main__':
    run()

