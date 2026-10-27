import soundfile as sf, numpy as np, sys
L="/Users/leighpierce/Git/246/irs/EMT 246 captures/96kHz/Non-Lin/"
np.set_printoptions(linewidth=200, precision=3)
for nm in ["12_0_05_10_00","12_0_03_10_04","12_0_04_100_20"]:
    x,fs=sf.read(L+nm+"_dc.wav",always_2d=True)
    print(nm, x.shape, fs, "len s", len(x)/fs, "mean", x.mean(0), "peak", np.abs(x).max(0), sf.info(L+nm+"_dc.wav").subtype)
    a=np.abs(x).max(1)
    print(" first nonzero", np.nonzero(a)[0][:3], "last nonzero", np.nonzero(a)[0][-1])
    print(" x[:8]", x[:8,0]); print(" first 300 rms by 50:", [float(np.sqrt((x[i:i+50,0]**2).mean())) for i in range(0,400,50)])
    n=len(x); seg=[float(np.sqrt((x[n-int(t*fs)-4800:n-int(t*fs),0]**2).mean())) for t in (0,0.2,0.5,0.8,1.0,1.5)]
    print(" rms 50ms ending at end-t for t=0,.2,.5,.8,1,1.5:", seg)
    # mean of last 1s
    print(" mean last 0.5s", x[-int(.5*fs):].mean(0), "mean first 2ms of data", x[:200].mean(0))
