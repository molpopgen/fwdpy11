import sys

import numpy as np
import pandas as pd

import fwdpy11

N = int(sys.argv[1])
rho = float(sys.argv[2])

moving_optimum_deme_0 = fwdpy11.GSSmo(
    [
        (0, 0, 1),  # Generation 0, optimum = 0, VS = 1
        (10 * N, 1, 1),  # Generation 10N, optimum = 1, VS = 1
    ]
)
moving_optimum_deme_1 = fwdpy11.GSSmo([(0, 0, 1), (10 * N, -1, 1)])

# The covariance matrix for effect sizes.
# The marginal Gaussians will have mean zero and sd = 0.1
# The deviates will be highly positively correlated.
sd = 0.1
covariance_matrix = np.matrix([0.999] * 4).reshape((2, 2))
np.fill_diagonal(covariance_matrix, 1)
covariance_matrix *= np.power(sd, 2)

pdict = {
    "nregions": [],
    "sregions": [
        # Multivariate Gaussian distribution of effect sizes
        fwdpy11.mvDES(
            fwdpy11.MultivariateGaussianEffects(
                beg=0, end=1, weight=1, h=1, matrix=covariance_matrix
            ),
            # Means of zero for each marginal Gaussian
            np.zeros(2),
        )
    ],
    "recregions": [fwdpy11.PoissonInterval(0, 1, rho / (4 * N))],
    "rates": (0, 1e-3, None),
    "demography": fwdpy11.DiscreteDemography(
        migmatrix=np.array([0.99, 0.01, 0.01, 0.99]).reshape((2, 2))
    ),
    "simlen": 10 * N + 200,  # 200 gens past optimum shift
    # Specify one gvalue object per deme
    "gvalue": [
        fwdpy11.Additive(
            ndemes=2,  # Number of demes
            scaling=2,  # 0, 0+sh, 0+2s for AA, Aa, and aa, respectively
            gv2w=moving_optimum_deme_0,  # Mapping of trait (genetic) value to fitness
        ),
        fwdpy11.Additive(ndemes=2, scaling=2, gv2w=moving_optimum_deme_1),
    ],
}

params = fwdpy11.ModelParams(**pdict)
pop = fwdpy11.DiploidPopulation([N, N], 1.0)
rng = fwdpy11.GSLrng(1010)
fwdpy11.evolvets(rng, pop, params, 100)
md = np.array(pop.diploid_metadata, copy=False)
df = pd.DataFrame.from_records(md[["deme", "g", "w"]])
g = df.groupby(["deme"])
print(g.mean())
