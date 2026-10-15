#!/usr/bin/env python3
"""Descriptive reconstruction sensitivity on exact-sequence paired models."""
import csv,json
from collections import Counter
from pathlib import Path
from statistics import median
P=Path(__file__).resolve().parents[1]/'results/metabolic_pilot'
def read(n):
 with (P/n).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write(n,z):
 with (P/n).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(z[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(z)
ids=read('gapseq_identity.tsv');valid={r['rep']:r for r in ids if r['status']=='identical'}
assert len(valid)==25
assert all(r['original_contigs_sha256']==r['current_contigs_sha256'] for r in valid.values())
g=read('gapseq_scenarios.tsv');assert len(g)==750 and all(r['rep'] in valid for r in g)
assert len({(r['rep'],r['variant'],r['scenario']) for r in g})==750
assert all(r['status']=='optimal' for r in g)
t={(r['rep'],r['variant']):r for r in read('model_traits.tsv')}
pairs=[]
for rep,info in valid.items():
 for variant in ['draft','filled']:
  z={r['scenario']:float(r['value']) for r in g if r['rep']==rep and r['variant']==variant}
  assert abs(z['closed'])<1e-6 and abs(z['inorganic'])<1e-6
  n=sum(v>1e-6 for k,v in z.items() if k not in {'closed','inorganic','H2'})
  q=t[rep,variant]
  pairs.append(dict(rep=rep,group=info['group'],variant=variant,gapseq_carbon_ATP_breadth=n,
                    gapsmith_carbon_ATP_breadth=q['carbon_ATP_breadth'],
                    gapseq_hydrogen_ATP_benefit=z['H2']-z['inorganic'],
                    gapsmith_hydrogen_ATP_benefit=q['hydrogen_ATP_benefit']))
write('reconstruction_pairs.tsv',pairs)
summary=[]
for variant in ['draft','filled']:
 for group in ['1l','1h','neither']:
  z=[r for r in pairs if r['variant']==variant and r['group']==group]
  summary.append(dict(variant=variant,group=group,n=len(z),
   gapseq_carbon_ATP_breadth_median=median(r['gapseq_carbon_ATP_breadth'] for r in z),
   gapsmith_carbon_ATP_breadth_median=median(int(r['gapsmith_carbon_ATP_breadth']) for r in z),
   gapseq_hydrogen_positive=sum(r['gapseq_hydrogen_ATP_benefit']>1e-6 for r in z),
   gapsmith_hydrogen_positive=sum(float(r['gapsmith_hydrogen_ATP_benefit'])>1e-6 for r in z)))
write('reconstruction_summary.tsv',summary)
m=read('gapseq_marker_mapping.tsv')
s=dict(exact_sequence_models=len(valid),scenario_LPs=len(g),statuses=dict(Counter(r['status'] for r in g)),
       gpr_translation_statuses=dict(Counter(r['status'] for r in m)),
       filled_hydrogen_positive_models=[r['rep'] for r in pairs if r['variant']=='filled' and r['gapseq_hydrogen_ATP_benefit']>1e-6])
(P/'reconstruction_summary.json').write_text(json.dumps(s,indent=2)+'\n')
print(json.dumps(s,indent=2))
