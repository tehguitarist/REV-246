/* S's four recirculating loops in float or 16-bit fixed point (model_s_results 36, Q10).
   Per loop i: r = D_i[n - P_i]; s = c0 r + c1 r[n-1] + a s[n-1] (shelf with the loop gain); y = (1 - p) s + p y[n-1] (HF one-pole);
   D_i[n] = b x[n] + y; out[n] = sum_i sign_i r.  Fixed point: Q15 coefficients, 16-bit data, products summed exactly in a 32-bit
   accumulator (int64 here, no overflow at these levels), each stored value (s, y, the delay write, the output) = acc >> 15 with
   truncation (floor, the TMS32010's SACH) or rounding (+2^14 first), saturated to 16 bits; the input is quantised to 16 bits too.
   usage: s_loops_fx <params.txt> <mode 0 float | 1 trunc | 2 round> <in.f64> <out.f64>
   params.txt: n_samples, then 4 lines "P c0 c1 a p b sign".  in.f64: n_samples doubles, full scale 1.0 (impulse or noise). */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

static int mode;
static long nclip = 0;
static int64_t q15(double c) { return (int64_t) llround(c * 32768.0); }
static int32_t store(int64_t acc) {            /* acc in Q30 -> Q15 data */
    int64_t v = mode == 2 ? (acc + (1 << 14)) >> 15 : acc >> 15;   /* >> on negative int64 is arithmetic (floor) with clang/gcc */
    if (v > 32767) { v = 32767; nclip++; }
    if (v < -32768) { v = -32768; nclip++; }
    return (int32_t) v;
}

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage\n"); return 1; }
    FILE *fp = fopen(argv[1], "r");
    long N; int P[4]; double c0[4], c1[4], a[4], p[4], b[4], sg[4];
    if (fscanf(fp, "%ld", &N) != 1) return 2;
    for (int i = 0; i < 4; i++) if (fscanf(fp, "%d %lf %lf %lf %lf %lf %lf", &P[i], &c0[i], &c1[i], &a[i], &p[i], &b[i], &sg[i]) != 7) return 2;
    fclose(fp);
    mode = atoi(argv[2]);
    double *x = malloc(N * sizeof(double)), *out = malloc(N * sizeof(double));
    FILE *fi = fopen(argv[3], "rb");
    if (fread(x, sizeof(double), N, fi) != (size_t) N) return 3;
    fclose(fi);
    int L = 2048;                                   /* power of two > every P */
    if (mode == 0) {
        double *D[4], s[4] = {0}, y[4] = {0}, rp[4] = {0};
        for (int i = 0; i < 4; i++) D[i] = calloc(L, sizeof(double));
        for (long n = 0; n < N; n++) {
            double o = 0;
            for (int i = 0; i < 4; i++) {
                double r = D[i][(n - P[i]) & (L - 1)];
                s[i] = c0[i] * r + c1[i] * rp[i] + a[i] * s[i];
                rp[i] = r;
                y[i] = (1 - p[i]) * s[i] + p[i] * y[i];
                D[i][n & (L - 1)] = b[i] * x[n] + y[i];
                o += sg[i] * r;
            }
            out[n] = o;
        }
    } else {
        int32_t *D[4], s[4] = {0}, y[4] = {0}, rp[4] = {0};
        int64_t C0[4], C1[4], A[4], Pp[4], P1[4], B[4];
        for (int i = 0; i < 4; i++) {
            D[i] = calloc(L, sizeof(int32_t));
            C0[i] = q15(c0[i]); C1[i] = q15(c1[i]); A[i] = q15(a[i]); Pp[i] = q15(p[i]); P1[i] = q15(1 - p[i]); B[i] = q15(b[i]);
        }
        for (long n = 0; n < N; n++) {
            int64_t xa = (int64_t) llround(x[n] * 32768.0 * 32768.0);       /* the ADC: 16 bits */
            int32_t xq = store(xa);
            int64_t o = 0;
            for (int i = 0; i < 4; i++) {
                int32_t r = D[i][(n - P[i]) & (L - 1)];
                s[i] = store(C0[i] * r + C1[i] * rp[i] + A[i] * s[i]);
                rp[i] = r;
                y[i] = store(P1[i] * s[i] + Pp[i] * y[i]);
                D[i][n & (L - 1)] = store(B[i] * xq + ((int64_t) y[i] << 15));
                o += (int64_t) (sg[i] * 32768.0) * r;
            }
            out[n] = store(o) / 32768.0;
        }
    }
    printf("%ld\n", nclip);
    FILE *fo = fopen(argv[4], "wb");
    fwrite(out, sizeof(double), N, fo);
    fclose(fo);
    return 0;
}
