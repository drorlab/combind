import numpy as np

class DensityEstimate:
    """
    Computes and stores density estimates f(x) for input values x.
    
    DensityEstimates are stored in the form of a probability density
    to get a value in the form of counts at a position, just multiply
    by n_samples.

    DensityEstimates can be averaged.
    """
    def __init__(self, points=100, domain=None, sd=1.0,
                 reflect=True, out_of_bounds=None):
        """
        points (int): number of values at which to compute density.
        domain ((float, float)): range of values at which to compute density
            if left as None, use min and max of input data.
        sd (float): standard deviation of gaussian kernel.
        reflect (bool): If True compute density for domain + left and right
            flanks, then reflect flanks and add to center. (This provides
            better behaviour at boundaries.)
        out_of_bounds (None or float): If asked to return a density for a
            value outside the domain, if None, return closest density, else
            return this value.
        """
        self.points = points
        self.out_of_bounds = out_of_bounds
        self.sd = sd
        self.reflect = reflect
        self.domain = domain
        self.n_samples = 0

    # I/O
    # File format: first line -> n_samples, sd, reflect, remaining -> x, fx
    def __str__(self):
        return '\n'.join([','.join([str(self.n_samples), str(self.sd),
                                    str(self.reflect)])]
                         + ['{},{}'.format(_x, _fx)
                            for _x, _fx in zip(self.x, self.fx)])

    def write(self, fname):
        with open(fname, 'w') as fp:
            fp.write(str(self))

    @classmethod
    def read(cls, fname):
        x, fx = [], []
        with open(fname) as fp:
            n_samples, sd, reflect = fp.readline().strip().split(',')
            for line in fp:
                _x, _fx = line.strip().split(',')
                x += [float(_x)]
                fx += [float(_fx)]
        de = DensityEstimate(points = len(x),
                             sd = float(sd),
                             reflect = bool(reflect))
        de.n_samples=float(n_samples)
        de.x, de.fx = np.array(x), np.array(fx)
        return de

    # Core methods
    def __call__(self, x):
        """
        Returns f(x) for the given value of x by linear interpolation.
        If x is out of functions domain, return the closest response
        and print a warning if self.out_of_bounds is None or else
        self.out_of_bounds.
        """
        return np.interp(x, self.x, self.fx)

    def fit(self, X, weights=1):
        """
        Given an array of values X and weights weights,
        compute a density estimate with standard deviation self.sd.
        If reflect, compute densities for each flank and add
        computed densities back to the center.
        If hist, normalize so that area under the curve is equal to 
        """
        if self.domain is None:
            self.x = np.linspace(X.min(), X.max(), self.points)
        else:
            self.x = np.linspace(self.domain[0], self.domain[1], self.points)

        if not X.shape[0]:
            return self._uniform()

        if self.reflect:
            if X.max() > self.x[-1] or X.min() < self.x[0]:
                print('Warning: Data out of domain of density estimate'
                      ' with reflected boundary conditions. Squishing'
                      ' data to be on specified domain.')
                X[X > self.x[-1]] = self.x[-1]
                X[X < self.x[0]] = self.x[0]
            r = self.x[-1] - self.x[0]
            self.x = np.hstack([self.x-r, self.x, self.x+r])

        self._kde(X, weights)
        
        if self.reflect:
            # left, center, right
            self.fx = (  self.fx[self.points:0:-1]
                       + self.fx[self.points:2*self.points]
                       + self.fx[-1:2*self.points-1:-1])
            self.x = self.x[self.points:2*self.points]

        self.fx *= (self.x.shape[0] / (self.x[-1]-self.x[0])) / self.fx.sum()
        self.n_samples = (weights*np.ones(X.shape)).sum()
        return self

    def data_loglikelihood(self, X, weights=1):
        return np.sum(np.log(self(X))*weights)

    def _gauss(self, mean, x):
        """
        Return PDF of N(mean, sd**2) at x.
        """
        return (np.exp(-0.5*((x - mean)/self.sd)**2)
                / (self.sd*np.sqrt(2*np.pi)))

    def _kde(self, X, weights):
        """
        Returns density estimate at each point in self.x for the input data X
        weighted by weights.
        """
        self.fx = []
        for mean in self.x:
            kernel = self._gauss(mean, X)
            self.fx += [(weights*kernel).sum()]
        self.fx = np.array(self.fx)

    def _uniform(self):
        self.n_samples = 0
        self.fx = np.ones(self.x.shape)
        self.fx /= self.x[-1]-self.x[0]
        return self

    def _average(self, other):
        """
        Returns a new function representing the average of the self and other
        functions. The domain of the new function covers the domain of both
        input functions with the same number of points as self.
        """
        assert self.reflect == other.reflect, "Either reflect or don't."
        if not other.n_samples:
            return self
        if not self.n_samples:
            return other

        de = DensityEstimate(points = self.points,
                             out_of_bounds=self.out_of_bounds,
                             reflect = self.reflect)
        de.x = np.linspace(min(self.x[0], other.x[0]),
                           max(self.x[-1], other.x[-1]), self.points)
        de.n_samples = self.n_samples + other.n_samples
        de.fx = (self(de.x)*self.n_samples + other(de.x)*other.n_samples) / de.n_samples
        return de

    @classmethod
    def merge(cls, des):
        out = des[0]
        for de in des[1:]:
            out = out._average(de)
        return out
