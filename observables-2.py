import MDAnalysis as mda
import numpy as np
import pytim
import matplotlib.pyplot as plt
from   pytim.datafiles   import *
from pytim.observables import Profile,Number

u       = mda.Universe(WATER_GRO,WATER_XTC)
oxygens = u.select_atoms("name OW")

obs     = Number()

interface = pytim.ITIM(u, alpha=2.0, max_layers=4,cluster_cut=3.5,centered=True,molecular=False)

profile = Profile(group=oxygens,observable=obs)

for ts in u.trajectory[::]:
    profile.sample()

low, up, avg = profile.get_values(binwidth=1.0)
plt.plot((low+up)/2., avg)
plt.show()