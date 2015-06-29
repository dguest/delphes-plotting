#!/usr/bin/env python3

import os, sys
import numpy as np
import h5py

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
# from matplotlib.ticker import FuncFormatter

def run():
    h5f = h5py.File(sys.argv[1])
    d0_plots = [('raw',h5f['particle_d0']), ('smeared',h5f['track_d0'])]
    _make_plot(d0_plots, 'd0_comp.pdf', 'd0')
    pt_plots = [(_pt, h5f['track_pt'])]
    _make_plot(pt_plots, 'track_pt.pdf', _pt, log=True)
    z0_plots = [('raw',h5f['particle_z0']), ('smeared', h5f['track_z0'])]
    _make_plot(z0_plots, 'z0_comp.pdf', 'z0')


_pt = r'$p_{\rm T}$'
_ax_size = 12
_text_size = 12

def _make_plot(names, outname, xname, log=False):
    """return the canvas with everything drawn on it"""
    fig = Figure(figsize=(5.0, 5.0*3/4))
    canvas = FigureCanvas(fig)
    ax = fig.add_subplot(1,1,1)
    color_itr = iter(['black', 'red'])
    units = set()
    for leg, dataset in names:
        extent = [dataset.attrs[x][0] for x in ['min', 'max']]
        yvals = np.array(dataset)[1:]
        xvals = np.linspace(*extent, num=(yvals.shape[0]))
        color = next(color_itr)
        # print(xvals.shape, yvals.shape)
        ax.plot(xvals, yvals, drawstyle='steps-post', label=leg, color=color)
        units.add(dataset.attrs['units'][0])
    assert len(units) == 1, "incompatible units: {}".format(units)

    xlab = '{} [{}]'.format(xname, units.pop())
    ax.set_xlabel(xlab, x=0.98, ha='right', size=_ax_size)
    ax.set_ylabel('tracks', y=0.98, ha='right', size=_ax_size)
    ax.set_xlim(*extent)
    if log:
        ax.set_ylim(1, ax.get_ylim()[1])
        ax.set_yscale('log')
    ax.legend(framealpha=0, prop=dict(size=_text_size))
    canvas.print_figure(outname, bbox_inches='tight')

if __name__ == '__main__':
    run()
