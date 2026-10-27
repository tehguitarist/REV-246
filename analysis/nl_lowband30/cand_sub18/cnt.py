from common import *
import collections
B, ks = catalog(0.0, 0.5, 2.0, 0.3)
print(len(ks)); c = collections.Counter((sess(B[k]), k[0]) for k in ks); print(sorted(c.items()))
B2, k2 = catalog(0.0, 0.5, 2.0, -10); print(len(k2))
