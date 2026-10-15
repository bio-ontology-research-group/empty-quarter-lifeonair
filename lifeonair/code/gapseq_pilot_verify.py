#!/usr/bin/env python3
"""Verify compact pilot outputs without loading source models or refitting."""
import csv,json,math
from collections import defaultdict
from pathlib import Path
P=Path(__file__).resolve().parents[1]/'results/gapseq_traits_pilot'
def read(n):
 with (P/n).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def mean(z):return sum(z)/len(z)
panel=read('pilot_panel.tsv');reps={r['rep'] for r in panel};n=len(reps);assert n==len(panel)
assert all(r['identity']=='exact' and r['original_hash']==r['current_hash'] for r in panel)
for block in {r['block'] for r in panel}:assert sorted(r['group'] for r in panel if r['block']==block)==['1h','1l','neither']
s=read('scenarios.tsv');a=read('audit.tsv');t=read('traits.tsv');assert len(s)==92*n and len(a)==len(t)==2*n
assert len({(r['rep'],r['variant'],r['kind'],r['name'],r['carbon'],r['oxygen'],r['objective']) for r in s})==len(s)
for r in s:assert r['rep'] in reps
monotonic=0
for rep in reps:
 for variant in ['filled','draft']:
  z={(float(r['carbon']),float(r['oxygen']),r['objective']):r for r in s if r['rep']==rep and r['variant']==variant and r['kind']=='mixture'}
  for objective in ['ATPM','biomass']:
   for c in [.06,.6,6]:
    for lo,hi in [(.1,1),(1,10)]:
     x,y=z[c,lo,objective],z[c,hi,objective]
     if x['status']==y['status']=='optimal':assert float(y['value'])>=float(x['value'])-1e-6;monotonic+=1
   for o in [.1,1,10]:
    for lo,hi in [(.06,.6),(.6,6)]:
     x,y=z[lo,o,objective],z[hi,o,objective]
     if x['status']==y['status']=='optimal':assert float(y['value'])>=float(x['value'])-1e-6;monotonic+=1
pairs=read('paired_abundance.tsv');assert all(1<=int(r['site'])<=60 for r in pairs)
actual_visits={(r['site'],r['trip']) for r in pairs}
with (P.parent/'h123_v2/h1_visit_log2fc_root_vs_bulk.tsv').open() as f:
 reference_visits={(r['site'],r['trip']) for r in csv.DictReader(f,delimiter='\t')}
assert actual_visits==reference_visits and len(actual_visits)==106
pred=read('predictions.tsv');metrics=read('prediction_metrics.tsv');checked=0
for m in metrics:
 z=[r for r in pred if r['variant']==m['variant'] and r['epsilon']==m['epsilon'] and r['kind']==m['kind'] and (m['fold']=='all' or r['fold']==m['fold'])]
 byrep=defaultdict(list)
 for r in z:byrep[r['rep']].append(r)
 assert len(z)==int(m['n_predictions']) and len(byrep)==int(m['n_models'])
 for col in ['baseline','augmented']:
  mse=mean([mean([(float(r['observed'])-float(r[col]))**2 for r in q]) for q in byrep.values()])
  mae=mean([mean([abs(float(r['observed'])-float(r[col])) for r in q]) for q in byrep.values()])
  assert math.isclose(mse,float(m[col+'_mse']),rel_tol=1e-10,abs_tol=1e-10)
  assert math.isclose(mae,float(m[col+'_mae']),rel_tol=1e-10,abs_tol=1e-10)
 assert math.isclose(float(m['mse_improvement']),1-float(m['augmented_mse'])/float(m['baseline_mse']),abs_tol=1e-10)
 checked+=1
summary=json.loads((P/'summary.json').read_text());assert summary['scale_up']==all(summary['checks'].values())
assert summary['scenario_LPs']==len(s) and summary['panel']==n
null=read('permutation_null.tsv');assert len(null)==199
print(f'PASS {n} exact-genome models, {len(s)} unique scenario LPs, {monotonic} optional-resource monotonicity checks, {checked} independently recomputed prediction metrics, 199 null shuffles; scale_up={summary["scale_up"]}')
