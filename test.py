import fwdpy11 as fp11
import fwdpy11.evolve
from collections import Counter

class RecordSFS:
    def __init__(self):
        self.data=[]
    def __call__(self,pop,generation):
        c=Counter()
        for m in pop.mcounts:
            if m > 0:
                c[m]+=1
        self.data.append((generation,c))

pop = fp11.Spop(1000)
rec=RecordSFS()
rng=fp11.GSLrng(101)
#This function is written in C++ and uses fwdpp,
#and we are SENDING A PYTHON CLASS IN TO PROCESS
#THE POP EACH GENERATION OMG.
fwdpy11.evolve.evolve(pop,rng,1000,10000,0.001,0.001,rec)

print (rec.data)
