#!/usr/bin/env python3
"""Exact-genome gapseq inventory and response-blind balanced pilot selection."""
import csv,gzip,hashlib,json,sys
from pathlib import Path
from collections import defaultdict,Counter
R=Path(sys.argv[1]);W=R/'06_atmospheric/gapseq_traits_pilot_clean'
def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,z):
 with p.open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(z[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(z)
def seqhash(p):
 op=gzip.open if str(p).endswith('.gz') else open
 hashes=[];seq=[]
 with op(p,'rt') as f:
  for line in f:
   if line.startswith('>'):
    if seq:hashes.append(hashlib.sha256(''.join(seq).upper().encode()).hexdigest())
    seq=[]
   else:seq.append(line.strip())
 if seq:hashes.append(hashlib.sha256(''.join(seq).upper().encode()).hexdigest())
 return hashlib.sha256('\n'.join(sorted(hashes)).encode()).hexdigest()
sources={}
for line in (R/'05_metabolism/hq_mag_list.tsv').read_text().splitlines():
 p,s,b,*_=line.split('\t');sources[s+'__'+b]=Path(p)
rows=[]
for r in read(W/'inputs/model_inventory.tsv'):
 if '__single_' not in r['rep'] or float(r['completeness'])<90 or float(r['contamination'])>5 or r['flagged']=='1' or r['subset']=='1' or not r['order'] or r['group']=='other':continue
 rep=r['rep'];tag=rep.replace('__single_','__');d=R/'05_metabolism/gapseq_out'/tag;orig=sources.get(tag)
 filled=d/(tag+'.xml');draft=d/(tag+'-draft.xml');current=R/'05_refine/_drep/out/dereplicated_genomes'/(rep+'.fa')
 z={k:v for k,v in r.items() if k not in ['filled','draft','eligible']}
 z.update(source=str(orig or ''),filled_path=str(filled),draft_path=str(draft),original_hash='',current_hash='',identity='missing_source_or_models')
 if orig and orig.exists() and filled.exists() and draft.exists() and current.exists():
  z['original_hash']=seqhash(orig);z['current_hash']=seqhash(current);z['identity']='exact' if z['original_hash']==z['current_hash'] else 'mismatch'
 rows.append(z)
write(W/'results/identity_inventory.tsv',rows)
by=defaultdict(lambda:defaultdict(list))
for r in rows:
 if r['identity']=='exact':by[r['order']][r['group']].append(r)
orders=sorted((o for o,g in by.items() if all(g[k] for k in ['1l','1h','neither'])),key=lambda o:hashlib.sha256(o.encode()).hexdigest())
panel=[];used=set()
for repeat in range(3):
 for order in orders:
  if len(panel)>=60:break
  g=by[order];cases=sorted((r for r in g['1l'] if r['rep'] not in used),key=lambda r:(-float(r['completeness']),float(r['contamination']),r['rep']))
  if not cases:continue
  case=cases[0];chosen=[]
  for group in ['1h','neither']:
   avail=[r for r in g[group] if r['rep'] not in used]
   if not avail:break
   chosen.append(min(avail,key=lambda r:(r['genus']!=case['genus'],r['family']!=case['family'],abs(float(r['completeness'])-float(case['completeness'])),float(r['contamination']),r['rep'])))
  if len(chosen)!=2:continue
  block=len(panel)//3+1
  for r in [case]+chosen:panel.append(dict(block=block,**r));used.add(r['rep'])
assert panel
write(W/'results/pilot_panel.tsv',panel)
s=dict(candidates=len(rows),identity=dict(Counter(r['identity'] for r in rows)),exact_groups=dict(Counter(r['group'] for r in rows if r['identity']=='exact')),panel=len(panel),orders=dict(Counter(r['order'] for r in panel)))
(W/'results/inventory_summary.json').write_text(json.dumps(s,indent=2)+'\n');print(json.dumps(s,indent=2))
