import sys; sys.path.insert(0, "/Users/leighpierce/Git/246/analysis/nl_lowband30/cand_pipeline")
from sim import *
np.set_printoptions(linewidth=200, precision=4)
key=tuple(float(v) for v in sys.argv[1:4]); r=align(key)
x,fs=sf.read(file_path(key),always_2d=True); nf=len(x)
cache=nl.ir(key)
h,d=pipeline(x); print("pipeline(real) vs cache max abs diff", np.abs(h[:len(cache)]-cache[:len(h)]).max(), "cache max", np.abs(cache).max(), h.shape, cache.shape)
nn=int(nf/RATIO)+40
m=model_cache(key,nn,r)
raw=forward_warp(m,r["ppm"]); dc=to_capture_rate(raw); xs=chain(dc,nf)
print("xs max", np.abs(xs).max(0), np.abs(xs).argmax(0)); print("x real max", np.abs(x).max(0), np.abs(x).argmax(0))
hs,ds=pipeline(xs); print("pipeline(sim) max", np.abs(hs).max(0), "cache", np.abs(cache).max(0))
A=[(x[:,c]@xs[:,c])/(xs[:,c]@xs[:,c]) for c in (0,1)]; print("A",A)
for c in (0,1):
    e=x[:,c]-A[c]*xs[:,c]; print("x null ch",c,10*np.log10((e@e)/(x[:,c]@x[:,c])))
