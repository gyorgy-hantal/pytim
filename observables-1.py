import MDAnalysis as mda
import numpy as np
import pytim
import matplotlib.pyplot as plt
from   pytim.datafiles import *
u = mda.Universe(WATER_GRO,WATER_XTC)
L = np.min(u.dimensions[:3])
oxygens = u.select_atoms("name OW")
radii=pytim_data.vdwradii(G43A1_TOP)
interface = pytim.ITIM(u,alpha=2.,itim_group=oxygens,                               max_layers=4,radii_dict=radii,                               cluster_cut=3.5)
rdf=pytim.observables.RDF2D(u,nbins=120)
for ts in u.trajectory[::50] :
    layer=interface.layers[0,1]
    rdf.sample(layer,layer)
rdf.rdf[0]=0.0
plt.plot(rdf.bins, rdf.rdf)
plt.show()