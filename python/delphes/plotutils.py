

# note: want to override this class to give locactions on log1p axes
# see more here: https://github.com/matplotlib/matplotlib/blob/master/lib/matplotlib/ticker.py#L930

class LogLocator(Locator):
    """
    Determine the tick locations for log axes
    """

    def __init__(self, base=10.0, subs=[1.0], numdecs=4, numticks=15):
        """
        place ticks on the location= base**i*subs[j]
        """
        self.base(base)
        self.subs(subs)
        self.numticks = numticks
        self.numdecs = numdecs

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

        return self.raise_if_exceeds(np.asarray(ticklocs))

    def view_limits(self, vmin, vmax):
        'Try to choose the view limits intelligently'
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
        return result
