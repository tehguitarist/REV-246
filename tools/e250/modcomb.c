// Feedback comb with a modulated (fractional) delay and a one-pole loop filter, for the E-program modulation-depth fit.
//   y[n] = x[n] + g * lp[n],  lp[n] = (1-p) * yd[n] + p * lp[n-1],  yd[n] = y(n - D(n)) by 4-point Lagrange interpolation,
//   D(n) = L + sum_i A_i sin(2 pi f_i n / rate + ph_i)   (A_i, L in samples)
#include <math.h>
#include <stdlib.h>
void modcomb(const double *x, double *y, long n, double L, double g, double p, int nl, const double *f, const double *A,
             const double *ph, double rate) {
    const double tp = 2.0 * M_PI / rate;
    double lp = 0.0;
    for (long i = 0; i < n; i++) {
        double D = L;
        for (int k = 0; k < nl; k++) D += A[k] * sin(tp * f[k] * (double)i + ph[k]);
        double pos = (double)i - D;
        long i0 = (long)floor(pos); double d = pos - (double)i0;
        double yd = 0.0;
        if (i0 - 1 >= 0 && i0 + 2 < i) {
            double ym1 = y[i0 - 1], y0 = y[i0], y1 = y[i0 + 1], y2 = y[i0 + 2];
            double c0 = -d * (d - 1) * (d - 2) / 6.0, c1 = (d + 1) * (d - 1) * (d - 2) / 2.0;
            double c2 = -(d + 1) * d * (d - 2) / 2.0, c3 = (d + 1) * d * (d - 1) / 6.0;
            yd = c0 * ym1 + c1 * y0 + c2 * y1 + c3 * y2;
        } else if (i0 >= 0 && i0 + 1 < i) {
            yd = y[i0] * (1 - d) + y[i0 + 1] * d;
        }
        lp = (1.0 - p) * yd + p * lp;
        y[i] = x[i] + g * lp;
    }
}
