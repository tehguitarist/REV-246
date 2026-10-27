import sys; sys.path.insert(0, "/Users/leighpierce/Git/246/analysis/nl_lowband30/cand_pipeline")
from sim import *
from scipy.signal import butter, sosfiltfilt
np.set_printoptions(linewidth=220, precision=3, suppress=True)
key=tuple(float(v) for v in sys.argv[1:4]); r=align(key)
x,fs=sf.read(file_path(key),always_2d=True); nf=len(x)
nn=int(nf/RATIO)+40
m=model_cache(key,nn,r); raw=forward_warp(m,r["ppm"]); dc=to_capture_rate(raw); xs=chain(dc,nf)*0.2586
sos=butter(4,[300,8000],btype="band",fs=FS96,output="sos")
xb,sb=sosfiltfilt(sos,x,axis=0),sosfiltfilt(sos,xs,axis=0)
B=2400
print("file len s",nf/FS96)
print("block rms ratio real/sim (ch0), blocks of 25 ms from end backward:")
nb=nf//B
rat=[]
for i in range(nb):
    a=nf-(i+1)*B; b=a+B
    rat.append((np.sqrt((xb[a:b,0]**2).mean()), np.sqrt((sb[a:b,0]**2).mean())))
rat=np.array(rat)
for i in range(min(nb,48)): print(f"t_end-{(i+1)*B/FS96:6.3f}s real {rat[i,0]:.3e} sim {rat[i,1]:.3e} ratio {rat[i,0]/rat[i,1]:.3f}")
