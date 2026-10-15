#!/usr/bin/env python3
"""Summarise the complete pilot with closure and LP monotonicity checks."""
import argparse,csv,hashlib,json,math
from pathlib import Path
from collections import defaultdict,Counter
from statistics import median

ROOT=Path(__file__).resolve().parents[1]/'results/metabolic_pilot'
def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows):
 with p.open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def main():
 parser=argparse.ArgumentParser();parser.add_argument('--from-frozen',action='store_true');args=parser.parse_args()
 panel=read(ROOT/'pilot_panel.tsv');assert len(panel)==72 and len({r['rep'] for r in panel})==72
 assert all(set(r['group'] for r in panel if r['block']==b)=={'1l','1h','neither'} for b in {r['block'] for r in panel})
 script=hashlib.sha256((ROOT.parents[1]/'code/metabolic_simulate.py').read_bytes()).hexdigest()
 tables={}
 for kind in ['audit','scenarios','hydrogen_reactions','marker_mapping','input_hashes','substrate_transport']:
  if args.from_frozen:
   tables[kind]=read(ROOT/(kind+'.tsv'));continue
  rows=[]
  for i,p in enumerate(panel):
   fn=ROOT/'simulations'/f'{kind}_{i}.tsv'
   if fn.exists():
    rr=read(fn);assert all(r['rep']==p['rep'] for r in rr);rows+=rr
   else:
    assert kind in {'marker_mapping','hydrogen_reactions'},str(fn)
    if kind=='hydrogen_reactions':
     assert all(int(r['h2_reactions'])==0 for r in tables['audit'] if r['rep']==p['rep']),str(fn)
  tables[kind]=rows
  if rows:write(ROOT/(kind+'.tsv'),rows)
 runtimes=json.loads((ROOT/'runtime_manifest.json').read_text()) if args.from_frozen else [dict(index=i,rep=panel[i]['rep'],**json.loads((ROOT/'simulations'/f'runtime_{i}.json').read_text())) for i in range(72)]
 assert len(runtimes)==72
 for i,r in enumerate(runtimes):
  assert r['index']==i and r['rep']==panel[i]['rep']
  assert r['script_sha256']==script,(i,'stale script')
 if not args.from_frozen:(ROOT/'runtime_manifest.json').write_text(json.dumps(runtimes,indent=2)+'\n')
 audits=tables['audit'];scenarios=tables['scenarios'];assert len(audits)==216 and len(scenarios)==10584
 assert len({(r['rep'],r['variant']) for r in audits})==216
 ix=defaultdict(list)
 for r in scenarios:ix[(r['rep'],r['variant'])].append(r)
 byrep={r['rep']:r for r in panel};traits=[];monotonic_checks=0
 for a in audits:
  z=ix[(a['rep'],a['variant'])];assert len(z)==49
  grad={(float(r['carbon_input']),float(r['h2_bound'])):r for r in z if r['kind']=='gradient'}
  for c in [0,.006,.06,.6,6]:
   vals=[grad[(c,h)] for h in [0,.001,.01,.1]]
   for x,y in zip(vals,vals[1:]):
    if x['status']==y['status']=='optimal':
     assert float(y['value'])>=float(x['value'])-1e-6,(a['rep'],'H2 monotonicity');monotonic_checks+=1
  def val(c,h):
   r=grad[(c,h)];return float(r['value']) if r['status']=='optimal' else None
  benefit=val(0,.1)-val(0,0) if val(0,.1) is not None and val(0,0) is not None else None
  breadth={obj:sum(r['status']=='optimal' and float(r['value'])>1e-6 for r in z if r['kind']=='substrate' and r['objective']==obj) for obj in ['ATPM','biomass']}
  info=byrep[a['rep']]
  traits.append(dict(rep=a['rep'],block=info['block'],order=info['order'],group=a['group'],variant=a['variant'],energy_qc_pass=a['energy_qc_pass'],
                     rich_growth=a['rich_growth'],native_typed_h2_reactions=a['native_typed_h2_reactions'],
                     carbon_ATP_breadth=breadth['ATPM'],carbon_growth_breadth=breadth['biomass'],hydrogen_ATP_benefit=benefit,
                     max_hydrogen_ATP_benefit=max(val(c,.1)-val(c,0) for c in [0,.006,.06,.6,6]) if all(val(c,h) is not None for c in [0,.006,.06,.6,6] for h in [0,.1]) else None))
 write(ROOT/'model_traits.tsv',traits)
 summaries=[]
 for v in ['draft','filled','filled_grounded_h2']:
  for g in ['1l','1h','neither']:
   z=[r for r in traits if r['variant']==v and r['group']==g]
   good=[r for r in z if r['energy_qc_pass']=='1']
   summaries.append(dict(variant=v,group=g,n=len(z),energy_qc_pass=len(good),rich_growth_positive=sum(float(r['rich_growth'] or 0)>1e-6 for r in z),
                         native_typed_h2_positive=sum(int(r['native_typed_h2_reactions'])>0 for r in z),
                         carbon_ATP_breadth_median=median(r['carbon_ATP_breadth'] for r in good) if good else None,
                         carbon_growth_breadth_median=median(r['carbon_growth_breadth'] for r in good) if good else None,
                         hydrogen_ATP_benefit_positive=sum(r['hydrogen_ATP_benefit'] is not None and r['hydrogen_ATP_benefit']>1e-6 for r in good)))
 write(ROOT/'group_summary.tsv',summaries)
 # Descriptive matched differences; no biological inference from QC failures.
 contrasts=[]
 for v in ['draft','filled','filled_grounded_h2']:
  z=[r for r in traits if r['variant']==v]
  for b in sorted({r['block'] for r in z},key=int):
   q={r['group']:r for r in z if r['block']==b}
   if not all(r['energy_qc_pass']=='1' for r in q.values()):continue
   for control in ['1h','neither']:
    for metric in ['carbon_ATP_breadth','carbon_growth_breadth','hydrogen_ATP_benefit']:
     if q['1l'][metric] is not None and q[control][metric] is not None:
      contrasts.append(dict(variant=v,block=b,order=q['1l']['order'],control=control,metric=metric,difference_1l_minus_control=q['1l'][metric]-q[control][metric]))
 write(ROOT/'matched_differences.tsv',contrasts)
 summary=dict(models=72,orders=len({r['order'] for r in panel}),variants=3,scenario_LPs=len(scenarios),
              statuses=dict(Counter(r['status'] for r in scenarios)),monotonicity_checks=monotonic_checks,
              missing_model_gene_counts=sorted({int(r['model_genes_absent_from_input']) for r in audits}),
              mapped_carrier_models=sum(int(r['exact_mapped_genes'])>0 for r in audits if r['variant']=='filled' and r['group']!='neither'),
              native_linked_carrier_models=sum(int(r['native_typed_h2_reactions'])>0 for r in audits if r['variant']=='filled' and r['group']!='neither'),
              script_sha256=script)
 (ROOT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
