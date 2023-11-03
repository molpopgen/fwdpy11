import sys

import numpy as np
import pandas as pd

import fwdpy11

N = int(sys.argv[1])
rho = float(sys.argv[2])

yaml = f"""
time_units: generations
demes:
  - name: alpha
    epochs:
      - start_size: {N}
  - name: beta
    epochs:
      - start_size: {N}
migrations:
  - demes: [alpha, beta]
    rate: 1e-2
"""

demography = fwdpy11.ForwardDemesGraph.from_demes(
    yaml, burnin=10 * N + 200, burnin_is_exact=True
)

moving_optimum_deme_0 = fwdpy11.GaussianStabilizingSelection.single_trait(
    optima=[
        fwdpy11.Optimum(when=0, optimum=0.0, VS=1.0),
        fwdpy11.Optimum(when=10 * N, optimum=1.0, VS=1.0),
    ],
)
moving_optimum_deme_1 = fwdpy11.GaussianStabilizingSelection.single_trait(
    optima=[
        fwdpy11.Optimum(when=0, optimum=0.0, VS=1.0),
        fwdpy11.Optimum(when=10 * N, optimum=-1.0, VS=1.0),
    ],
)

# The covariance matrix for effect sizes.
# The marginal Gaussians will have mean zero and sd = 0.1
# The deviates will be highly positively correlated.
sd = 0.1
covariance_matrix = np.array([0.999] * 4).reshape((2, 2))
np.fill_diagonal(covariance_matrix, 1)
covariance_matrix *= np.power(sd, 2)

pdict = {
    "nregions": [],
    "sregions": [
        # Multivariate Gaussian distribution of effect sizes
        fwdpy11.mvDES(
            fwdpy11.MultivariateGaussianEffects(
                beg=0, end=1, weight=1, h=1, cov_matrix=covariance_matrix
            ),
            # Means of zero for each marginal Gaussian
            np.zeros(2),
        )
    ],
    "recregions": [fwdpy11.PoissonInterval(0, 1, rho / (4 * N))],
    "rates": (0, 1e-3, None),
    "demography": demography,
    "simlen": 10 * N + 200,  # 200 gens past optimum shift
    # Specify one gvalue object per deme
    "gvalue": [
        fwdpy11.Additive(
            ndemes=2,  # Number of demes
            scaling=2,  # 0, 0+sh, 0+2s for AA, Aa, and aa, respectively
            # Mapping of trait (genetic) value to fitness
            gvalue_to_fitness=moving_optimum_deme_0,
        ),
        fwdpy11.Additive(ndemes=2, scaling=2, gvalue_to_fitness=moving_optimum_deme_1),
    ],
}

params = fwdpy11.ModelParams(**pdict)
pop = fwdpy11.DiploidPopulation([N, N], 1.0)
rng = fwdpy11.GSLrng(1010)
fwdpy11.evolvets(rng, pop, params, 100)
assert pop.generation == 10 * N + 200
md = np.array(pop.diploid_metadata, copy=False)
df = pd.DataFrame.from_records(md[["deme", "g", "w"]])
g = df.groupby(["deme"])
print(g.mean())
