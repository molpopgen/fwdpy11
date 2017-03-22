import fwdpy11
x=fwdpy11.ConstantS(0,1,1,0.1)
y=x.callback()
print(y)

x=fwdpy11.GammaS(0,1,1,0.1,1)
y=x.callback()
rng=fwdpy11.GSLrng(101)
print(y(rng))
print(y)

nregions = [fwdpy11.Region(0,1,0.5),fwdpy11.Region(0.5,0.6,0.5)]
sregions = [fwdpy11.ExpS(0,1,1,1)]

x = fwdpy11.makeMutationRegions(nregions,sregions)

print(x)
print(x.nbeg)
print(x.shmodels)
