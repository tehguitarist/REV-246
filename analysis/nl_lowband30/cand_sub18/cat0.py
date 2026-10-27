from common import *
import collections
B, ks = catalog(0.0, 0.5, 2.0, -10)
c = collections.Counter((sess(B[k]), k[0], 'L0' if k[2]==0 else ('Ls' if k[2]<=0.4 else 'Ll')) for k in ks)
for kk in sorted(c): print(kk, c[kk])
print(len(ks))
import numpy as np
print(sorted(set(sess(B[k]) for k in B)))
