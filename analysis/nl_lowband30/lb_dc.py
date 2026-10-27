from lb_common import *
print("READ_SIGNS (sum per loop => DC gain of the late reads is sum*1/(1-g)):")
for ch in (0, 1):
    print(" ch", ch, {P: (list(s), sum(s)) for P, s in nlt.READ_SIGNS[ch].items()})
# model DC / LF gain of the branch (late only) vs LATE, relative to main-tap unit: |B(f)| at 5, 10, 20, 50 Hz for SPAN 0.5 EARLY 10
N = 1 << 18
f = np.fft.rfftfreq(N, 1 / FS)
for late in (0.4, 1.0, 2.0):
    for ch in (0, 1):
        F, Bb = M.spectra((0.5, 10.0, late), ch, N)
        Fn, Bn = M.spectra((0.5, 10.0, 0.0), ch, N)
        late_only = Bb - Bn
        i = [np.argmin(abs(f - x)) for x in (1, 5, 10, 20, 50)]
        print(f"LATE {late} ch{ch}: |late branch| at 1/5/10/20/50 Hz (dB re 1): ", np.round(20 * np.log10(abs(late_only[i]) + 1e-12), 1), " |front| ", np.round(20 * np.log10(abs(F[i]) + 1e-12), 1))
