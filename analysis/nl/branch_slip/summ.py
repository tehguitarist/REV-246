import json,glob,collections
import numpy as np
S="/Users/leighpierce/Git/246/analysis/nl/branch_slip"
step=96000/26373.6
rows=collections.defaultdict(list); bands=collections.defaultdict(list); before=collections.defaultdict(list)
def cls(n): return 'n=0' if n==0 else ('n=-1..-7' if n>=-7 else 'n=-233..-320')
orig={}
for f in glob.glob('/Users/leighpierce/Git/246/analysis/plugin_check/n_fix/*.json'):
    d=json.load(open(f)); orig[f.split('/')[-1]]=d
for src,files in (('orig',None),('deslip',glob.glob(S+'/deslip/n_fix/*.json'))):
    it=[(k,d) for k,d in orig.items() if d['bench'] and round((d['bench']['o']-1)*step)==0] if src=='orig' else [(f.split('/')[-1],json.load(open(f))) for f in files]
    for name,d in it:
        n=round((d['bench']['o']-1)*step); c=cls(n)+(' trimmed' if d['trimmed'] else '')
        for t,v in d['chain3']['windows'].items():
            rows[(c,src,t)].append((v['null1'],v['null0'],v['null1']-v['floor1'],d['bench']['windows'].get(t,{}).get('null1',np.nan)))
            if src=='deslip': before[(c,t)].append(orig[name]['chain3']['windows'].get(t,{}).get('null1',np.nan))
        for b,v in d['chain3']['bands'].items(): bands[(c,src,b)].append(v[1])
print("class | src | window | N | plugin R (L) | excess R | bench R | <=-50 R | before de-slip R")
for k in sorted(rows,key=lambda k:(k[0],k[1],float(k[2].split('-')[0]))):
    a=np.array(rows[k]); bf=np.nanmedian(before[(k[0],k[2])]) if k[1]=='deslip' else np.nan
    print(" | ".join([k[0],k[1],k[2],str(len(a)),'%.1f (%.1f)'%(np.median(a[:,0]),np.median(a[:,1])),'%+.1f'%np.median(a[:,2]),'%.1f'%np.nanmedian(a[:,3]),'%d/%d'%((a[:,0]<=-50).sum(),len(a)),'%.1f'%bf]))
print()
for k in sorted(bands): print(k,len(bands[k]),'%.1f'%np.median(bands[k]))
