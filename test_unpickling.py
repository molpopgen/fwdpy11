import pickle
import fwdpy11
with open("test.pickle","rb") as f:
    for i in range(100):
        p=pickle.load(f)
        print(type(p[0]),type(p[1]),len(p[1].mutations))
