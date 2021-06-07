import numpy as np
from glob import glob

# Read in all the files and find the time region in which they overlap
all_links=[]
latest_start=0
earliest_end=2**64

for f in glob("/data/lar/dunedaq/rodrigues/hit-dumps/neutron-source-runs-2020-07-09/felix5*_off.txt"):
    tmp=np.loadtxt(f, dtype=int,
                        converters={0: lambda x : int(x, base=16)}, max_rows=500000)
    tstart=tmp[:,2]
    latest_start=max(latest_start, np.min(tstart))
    earliest_end=min(earliest_end, np.max(tstart))
    all_links.append(tmp)


all_links_overlap=[]

for i,link in enumerate(all_links):
    # Reduce to the overlapping time region
    tstart=link[:,2]
    overlap_mask=np.logical_and(tstart>latest_start, tstart<earliest_end)
    all_links_overlap.append(link[overlap_mask])


a=np.vstack(all_links_overlap)
tstart=a[:,2]
# Sort by start time
inds=np.argsort(tstart)
a=a[inds]

np.savetxt(f"full-apa.txt", a[:,1:3], fmt="% 10d")
