"""
Wrappers for basic matplotlib figures.
"""

import os

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.colors import LogNorm

import numpy as np

class Canvas:
    default_name = 'test.pdf'
    def __init__(self, out_path=None, figsize=(5.0,5.0*3/4), ext=None):
        self.fig = Figure(figsize)
        self.canvas = FigureCanvas(self.fig)
        self.ax = self.fig.add_subplot(1,1,1)
        self.out_path = out_path
        self.ext = ext

    def save(self, out_path=None, ext=None):
        output = out_path or self.out_path
        assert output, "an output file name is required"
        out_dir, out_file = os.path.split(output)
        if ext:
            out_file = '{}.{}'.format(out_file, ext.lstrip('.'))
        if out_dir and not os.path.isdir(out_dir):
            os.makedirs(out_dir)
        self.canvas.print_figure(output, bbox_inches='tight')

    def __enter__(self):
        if not self.out_path:
            self.out_path = self.default_name
        return self
    def __exit__(self, extype, exval, extb):
        if extype:
            return None
        self.save(self.out_path, ext=self.ext)
        return True


# _________________________________________________________________________
# specific draw routines
_ax_size = 12
_text_size = 12

def draw2d(can, hist, log=False, **kwargs):
    """
    Simple draw routine for 2d hist. Assumes an ndhist, and a
    canvas with attributes `ax` (the axes) and `fig` (the figure).
    """
    ax = can.ax
    fig = can.fig

    axes = hist.axes
    xlims, ylims = axes[0].lims, axes[1].lims
    imextent = list(xlims) + list(ylims)

    args = dict(aspect='auto', origin='lower', extent=imextent)
    args.update(**kwargs)
    if log:
        args['norm'] = LogNorm()

    image = hist.hist[1:-1,1:-1].T
    if image.sum() > 0:
        im = ax.imshow(image, interpolation='nearest', **args)
        cb = fig.colorbar(im)

    ax.set_xlim(*xlims)
    ax.set_ylim(*ylims)
    ax.set_xlabel(_ax_name(axes[0]), x=0.98, ha='right', size=_ax_size)
    ax.set_ylabel(_ax_name(axes[1]), y=0.98, ha='right', size=_ax_size)

def draw1d(can, hists, ylabel='entries', log=False):
    """
    Draw a list of 1d histograms on canvas `can`.

    The Histograms should come from `ndhist`, and can contain can
    include a few additional attributes:
     - color: the color to plot them as
     - label: legend label
     - norm: multiply by this
    """
    draw_opts = dict(drawstyle='steps-post')
    xax = hists[0].axes[0]
    for hist in hists:
        assert len(hist.axes) == 1, "only works for 1d hists (for now)"
        assert np.all(np.isclose(xax.lims, hist.axes[0].lims))
        opts = draw_opts.copy()

        def add(attribute):
            """Pull draw options from the hist objects"""
            if hasattr(hist, attribute):
                opts[attribute] = getattr(hist, attribute)
        add('color')
        add('label')
        x, y = getxy(hist)
        normy = getattr(hist, 'norm', 1.0) * y
        can.ax.plot(x, normy, **opts)
    can.ax.set_xlim(*xax.lims)
    can.ax.set_xlabel(_ax_name(xax), x=0.98, ha='right', size=_ax_size)
    can.ax.set_ylabel(ylabel, y=0.98, ha='right', size=_ax_size)
    handles, labels = can.ax.get_legend_handles_labels()
    if handles:
        can.ax.legend(framealpha=0)
    if log:
        can.ax.set_ylim(1, ax.get_ylim()[1]*2)
        can.ax.set_yscale('log')


# ________________________________________________________________________
# lower level draw utilities

def getxy(hist, rebin_target=None):
    """
    Return values for `plot` by clipping overflow bins and duplicating
    the final y-value.
    Assume we'll use `drawstyle='steps-post'`
    """
    x_ax = hist.axes[0]
    y_vals = hist.hist[1:-1]

    # TODO: make rebinning work
    assert rebin_target is None, "not implemented"
    # while len(y_vals) > rebin_target:
    #     y_vals = y_vals.reshape([-1,10]).sum(axis=1)
    x_vals = np.linspace(*x_ax.lims, num=(y_vals.shape[0] + 1))
    return x_vals, np.r_[y_vals, y_vals[-1]]

def _ax_name(ax):
    nm, un = ax.name, ax.units
    return '{} [{}]'.format(nm, un) if un else nm

