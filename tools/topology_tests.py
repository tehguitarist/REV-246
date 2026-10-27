"""Topology tests T1 (L/R per-loop circular shift, PHAT) and T2 (per-harmonic mode-magnitude spread)
from docs/algorithm_research.md section 7.2, on the CLA IR plus two synthetic controls:
  A = common all-pass diffuser -> 4 loops -> one read tap per channel
  B = direct injection -> 4 loops -> 16 exponentially weighted random read taps per channel

    .venv/bin/python tools/topology_tests.py   # needs analysis/shortRT/h_dev.npy from tools/analyze_ir.py
"""
import numpy as np, scipy.sparse as sp
from scipy.sparse.linalg import lsqr
from scipy.signal import butter, sosfiltfilt, lfilter
fsd=19.2e6/728; rt=0.98; P=[1321,1389,1461,1587]
def decompose(h,k0,flp=6000):
    sos=butter(8,flp,fs=fsd,output='sos'); comps={}
    for c in range(2):
        y=sosfiltfilt(sos,h[:,c])*10**(3*np.arange(len(h))/fsd/rt)
        y=y[k0+700:k0+int(0.75*fsd)]; n=np.arange(len(y))+700
        rows=[];cols=[];off=0
        for p in P: rows.append(np.arange(len(n))); cols.append(off+n%p); off+=p
        A=sp.csr_matrix((np.ones(len(n)*4),(np.concatenate(rows),np.concatenate(cols))),shape=(len(n),off))
        sol=lsqr(A,y,iter_lim=3000,atol=1e-10,btol=1e-10)[0]; o=0
        for p in P: comps[(c,p)]=sol[o:o+p]-sol[o:o+p].mean(); o+=p
    return comps
def t1(comps):
    out=[]
    for p in P:
        A=np.fft.rfft(comps[(0,p)]); B=np.fft.rfft(comps[(1,p)]); f=np.arange(len(A))*fsd/p
        X=A*np.conj(B); X=np.where((f>150)&(f<5000), X/(np.abs(X)+1e-30), 0)
        x=np.fft.irfft(X,p); x/=x.max() if abs(x.max())>abs(x.min()) else -x.min()
        i=int(np.argmax(np.abs(x))); srt=np.sort(np.abs(x))[::-1]
        out.append((p, int(i if i<p/2 else i-p), round(float(srt[1]/srt[0]),2)))
    return out
def t2(comps):
    res=[]
    for c in range(2):
        row=[]
        for p in P:
            W=np.abs(np.fft.rfft(comps[(c,p)])); f=np.arange(len(W))*fsd/p
            others=np.concatenate([np.arange(1,700)*fsd/q for q in P if q!=p])
            sel=(f>150)&(f<5000)&np.array([np.min(np.abs(others-fr))>1.0 for fr in f])
            L=20*np.log10(W[sel]+1e-12); d=(L-np.convolve(L,np.ones(15)/15,'same'))[7:-7]
            row.append(round(float(np.std(d)),2))
        res.append(row)
    return res
h=np.load(__import__('os').path.join(__import__('os').path.dirname(__file__), '..', 'analysis', 'shortRT', 'h_dev.npy'))
c=decompose(h,263)
print('EMT  T1 (loop, best lag, 2nd/1st PHAT peak):',t1(c)); print('EMT  T2 std dB:',t2(c))
# controls: synthetic 4-comb reverbs, 2 s at fsd, same loops and RT, one-pole damping
rng=np.random.default_rng(3); n=int(2*fsd)
def comb(x,p):
    g=10**(-3*p/fsd/rt); a=np.zeros(p+1); a[0]=1; a[p]=-g; b=np.zeros(p+1); b[p]=1
    return lfilter(b,a,x)
def allpass_chain(x):
    for d,g in [(142,0.7),(107,0.7),(379,0.7),(277,0.7),(613,0.6),(997,0.6)]:
        b=np.zeros(d+1); b[0]=-g; b[d]=1; a=np.zeros(d+1); a[0]=1; a[d]=-g; x=lfilter(b,a,x)
    return x
imp=np.zeros(n); imp[263]=1
diff=allpass_chain(imp)
def build(kind):
    out=np.zeros((n,2))
    for i,p in enumerate(P):
        if kind=='A':   # common allpass diffuser -> loop -> one read tap per channel
            y=comb(diff,p); out[:,0]+=np.roll(y,37*i); out[:,1]+=np.roll(y,211+53*i)*(-1)**i
        else:           # B: direct injection, 16 random read taps per channel, exponential energy
            y=comb(imp,p)
            for ch in range(2):
                taps=rng.integers(0,p,16); g=rng.standard_normal(16)*10**(-1.5*taps/fsd/rt)
                for t,gg in zip(taps,g): out[:,ch]+=gg*np.roll(y,t)
    return out+rng.standard_normal((n,2))*1e-5
for kind in 'AB':
    s=build(kind); cc=decompose(s,263)
    print(f'ctl {kind} T1:',t1(cc)); print(f'ctl {kind} T2 std dB:',t2(cc))
