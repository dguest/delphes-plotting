from matplotlib.ticker import Locator, FuncFormatter
import numpy as np
import math

# __________________________________________________________________________
# axes class from ndhist
def get_axes(ds):
    """returns a list of axes from a Dataset produced via ndhist"""
    axes_ar = ds.attrs['axes']
    ax_props = axes_ar.dtype.names
    axes = []
    for ax in axes_ar:
        the_ax = Axis(ax_props, ax)
        axes.append(the_ax)
    return axes

class Axis:
    def __init__(self, prop_list, array):
        self.name = array[prop_list.index('name')]
        self.lims = [array[prop_list.index(x)] for x in ['min', 'max']]
        self.units = array[prop_list.index('units')]
    def __str__(self):
        prints = [self.name] + list(self.lims) + [self.units]
        return 'name: {}, range: {}-{}, units {}'.format(
            *(str(x) for x in prints))

# _________________________________________________________________________
# stuff to deal with log1p axes
def fauxify(axis, units=None, minval=None):
    axis.set_minor_locator(FauxLogLocator(subs=np.r_[1:5], numdecs=2))
    axis.set_major_locator(FauxLogLocator(minval=minval))
    formatter = FuncFormatter(_faux_log_formatting)
    axis.set_major_formatter(formatter)
    if units:
        less = [x for x in units.split() if x != 'log1p']
        faux_units = ' '.join(less)
        return faux_units


def _exp1m(val):
    return math.exp(val) - 1

def _faux_log_formatting(value, pos):
    value = _exp1m(value)
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
    if not exp:
        return r'{:.0f}'.format(base)
    return r'{:.0f}$^{{\mathdefault{{ {} }} }}$'.format(base, exp)

# see more here: https://github.com/matplotlib/matplotlib/blob/master/lib/matplotlib/ticker.py#L930

class FauxLogLocator(Locator):
    """
    Determine the tick locations for log axes
    """

    def __init__(self, base=10.0, subs=[1.0], numdecs=4, numticks=15,
                 minval=None):
        """
        place ticks on the location= base**i*subs[j]
        """
        self.base(base)
        self.subs(subs)
        self.numticks = numticks
        self.numdecs = numdecs
        self.minval = minval

    def set_params(self, base=None, subs=None, numdecs=None, numticks=None):
        """Set parameters within this locator."""
        if base is not None:
            self.base = base
        if subs is not None:
            self.subs = subs
        if numdecs is not None:
            self.numdecs = numdecs
        if numticks is not None:
            self.numticks = numticks

    def base(self, base):
        """
        set the base of the log scaling (major tick every base**i, i integer)
        """
        self._base = base + 0.0

    def subs(self, subs):
        """
        set the minor ticks the log scaling every base**i*subs[j]
        """
        if subs is None:
            self._subs = None  # autosub
        else:
            self._subs = np.asarray(subs) + 0.0

    def __call__(self):
        'Return the locations of the ticks'
        vmin, vmax = self.axis.get_view_interval()
        return self.tick_values(vmin, vmax)

    def tick_values(self, vmin, vmax):
        vmin, vmax = (_exp1m(x) for x in (vmin, vmax))
        b = self._base
        # dummy axis has no axes attribute
        if hasattr(self.axis, 'axes') and self.axis.axes.name == 'polar':
            vmax = math.ceil(math.log(vmax) / math.log(b))
            decades = np.arange(vmax - self.numdecs, vmax)
            ticklocs = b ** decades

            return ticklocs

        if vmin <= 0.0:
            if self.axis is not None:
                vmin = self.axis.get_minpos()

            if vmin <= 0.0 or not np.isfinite(vmin):
                raise ValueError(
                    "Data has no positive values, and therefore can not be "
                    "log-scaled.")

        vmin = math.log(vmin) / math.log(b)
        vmax = math.log(vmax) / math.log(b)

        if vmax < vmin:
            vmin, vmax = vmax, vmin

        numdec = math.floor(vmax) - math.ceil(vmin)

        if self._subs is None:  # autosub
            if numdec > 10:
                subs = np.array([1.0])
            elif numdec > 6:
                subs = np.arange(2.0, b, 2.0)
            else:
                subs = np.arange(2.0, b)
        else:
            subs = self._subs

        stride = 1
        while numdec / stride + 1 > self.numticks:
            stride += 1

        decades = np.arange(math.floor(vmin) - stride,
                            math.ceil(vmax) + 2 * stride, stride)
        if hasattr(self, '_transform'):
            ticklocs = self._transform.inverted().transform(decades)
            if len(subs) > 1 or (len(subs == 1) and subs[0] != 1.0):
                ticklocs = np.ravel(np.outer(subs, ticklocs))
        else:
            if len(subs) > 1 or (len(subs == 1) and subs[0] != 1.0):
                ticklocs = []
                for decadeStart in b ** decades:
                    ticklocs.extend(subs * decadeStart)
            else:
                ticklocs = b ** decades

        if self.minval is not None:
            ticklocs = ticklocs[ticklocs > math.log1p(self.minval)]

        return self.raise_if_exceeds(np.log1p(np.asarray(ticklocs)))

    def view_limits(self, vmin, vmax):
        'Try to choose the view limits intelligently'
        vmin, vmax = (_exp1m(x) for x in (vmin, vmax))
        b = self._base

        if vmax < vmin:
            vmin, vmax = vmax, vmin

        if self.axis.axes.name == 'polar':
            vmax = math.ceil(math.log(vmax) / math.log(b))
            vmin = b ** (vmax - self.numdecs)
            return vmin, vmax

        minpos = self.axis.get_minpos()

        if minpos <= 0 or not np.isfinite(minpos):
            raise ValueError(
                "Data has no positive values, and therefore can not be "
                "log-scaled.")

        if vmin <= minpos:
            vmin = minpos

        if not is_decade(vmin, self._base):
            vmin = decade_down(vmin, self._base)
        if not is_decade(vmax, self._base):
            vmax = decade_up(vmax, self._base)

        if vmin == vmax:
            vmin = decade_down(vmin, self._base)
            vmax = decade_up(vmax, self._base)
        result = mtransforms.nonsingular(vmin, vmax)
        return np.log1p(result)
