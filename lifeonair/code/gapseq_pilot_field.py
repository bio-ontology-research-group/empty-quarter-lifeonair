#!/usr/bin/env python3
"""Prespecified gapseq trait prediction with grouped out-of-fold evaluation."""
import csv,gzip,hashlib,json,math,sys
from pathlib import Path
from collections import defaultdict,Counter
import numpy as np
R=Path(sys.argv[1]);W=R/'06_atmospheric/gapseq_traits_pilot_clean';P=W/'results'
def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,z):
 if not z:return
 with p.open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(z[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(z)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def numeric(x):
 try:return float(x)
 except (TypeError,ValueError):return np.nan
panel=read(P/'pilot_panel.tsv');byrep={r['rep']:r for r in panel};aud=[];sc=[];hashes=[];runtimes=[]
for i,r in enumerate(panel):
 a=read(P/f'audit_{i}.tsv');s=read(P/f'scenarios_{i}.tsv')
 assert len(a)==2 and len(s)==92 and all(x['rep']==r['rep'] for x in a+s)
 run=json.loads((P/f'runtime_{i}.json').read_text());assert run['script_sha256']==sha(W/'scripts/gapseq_pilot_simulate.py')
 assert run['helper_sha256']==sha(W/'scripts/metabolic_simulate.py')
 assert run['rep']==r['rep'];runtimes.append(run);aud+=a;sc+=s;hashes+=read(P/f'hashes_{i}.tsv')
write(P/'audit.tsv',aud);write(P/'scenarios.tsv',sc);write(P/'model_hashes.tsv',hashes)
(P/'runtime_manifest.json').write_text(json.dumps(runtimes,indent=2)+'\n')
traits=[]
for a in aud:
 z=[s for s in sc if s['rep']==a['rep'] and s['variant']==a['variant']]
 sub=[s for s in z if s['kind']=='substrate' and s['objective']=='ATPM'];assert len(sub)==12
 mix=next(s for s in z if s['kind']=='mixture' and float(s['carbon'])==6 and float(s['oxygen'])==10 and s['objective']=='ATPM')
 good=a['energy_qc']=='1' and all(s['status']=='optimal' for s in sub+[mix])
 breadth=sum(float(s['value'])>1e-6 for s in sub if s['status']=='optimal')
 traits.append(dict(rep=a['rep'],variant=a['variant'],group=a['group'],energy_qc=a['energy_qc'],trait_qc=int(good),
  atp_breadth=breadth,mixture_atp=numeric(mix['value']),growth_breadth=sum(s['status']=='optimal' and float(s['value'])>1e-6 for s in z if s['kind']=='substrate' and s['objective']=='biomass'),
  rich_growth=numeric(a['rich_growth'])))
write(P/'traits.tsv',traits)
meta=read(W/'inputs/samples.tsv')
qc={r['sample'] for r in read(W/'inputs/guild_genome_equivalents_v2_r3noflagged.tsv') if int(r['id_floor'])==60 and int(r['scg_genes_detected'])>=10}
meta=[r for r in meta if r['sample_id'] in qc and r['compartment'] in ['Deep','PlantRoot'] and r['trip'] and r['site'].isdigit() and 1<=int(r['site'])<=60]
assert len({r['sample_id'] for r in meta})==len(meta)
samples={r['sample_id'] for r in meta};ab=defaultdict(dict);measured=set()
with gzip.open(W/'inputs/coverm_rep_abundance_long.tsv.gz','rt') as f:
 for sample,rep,value in csv.reader(f,delimiter='\t'):
  measured.add(sample)
  if rep in byrep and sample in samples:
   assert sample not in ab[rep];ab[rep][sample]=float(value)
missing_coverm=sorted(samples-measured)
meta=[r for r in meta if r['sample_id'] in measured]
samples={r['sample_id'] for r in meta}
cells=defaultdict(list)
for m in meta:cells[m['site'],m['trip'],m['compartment']].append(m['sample_id'])
visits=sorted({(s,t) for s,t,c in cells if (s,t,'Deep') in cells and (s,t,'PlantRoot') in cells})
pairs=[]
for rep in byrep:
 for site,trip in visits:
  root=float(np.median([ab[rep].get(s,0) for s in cells[site,trip,'PlantRoot']]))
  bulk=float(np.median([ab[rep].get(s,0) for s in cells[site,trip,'Deep']]))
  pairs.append(dict(rep=rep,site=site,trip=trip,root_pct=root,bulk_pct=bulk,detected=int(max(root,bulk)>=.01)))
write(P/'paired_abundance.tsv',pairs)
coverage=[];eligible=[]
for rep in byrep:
 z=[r for r in pairs if r['rep']==rep and r['detected']]
 ok=len(z)>=5 and len({r['site'] for r in z})>=3
 coverage.append(dict(rep=rep,detected_visits=len(z),detected_sites=len({r['site'] for r in z}),field_eligible=int(ok)))
 if ok:eligible.append(rep)
write(P/'field_coverage.tsv',coverage)
growth={r['rep']:numeric(r['d_hours']) for r in read(W/'inputs/rep_growth_v2.tsv')}
def genus(rep):
 r=byrep[rep];return r['order']+';'+(r['genus'] or 'unassigned_family:'+r['family'])
def folds(values):
 return {g:i%5 for i,g in enumerate(sorted(set(values),key=lambda x:hashlib.sha256(str(x).encode()).hexdigest()))}
gfold=folds(genus(r) for r in eligible);sfold=folds(r['site'] for r in pairs)
write(P/'fold_assignments.tsv',[dict(kind='genus',group=g,fold=f) for g,f in gfold.items()]+[dict(kind='site',group=g,fold=f) for g,f in sfold.items()])
# Fixed covariates: fit imputation/scaling/levels ONLY on training MAGs.
def design(train,test,feature):
 def continuous(reps):
  z=[]
  for rep in reps:
   r=byrep[rep];d=growth.get(rep,np.nan);ld=math.log(d) if np.isfinite(d) and d>0 else np.nan
   v=[ld,int(not np.isfinite(ld)),math.log(float(r['length_bp'])),float(r['completeness']),float(r['contamination'])]
   if feature is not None:v+=list(feature[rep])
   z.append(v)
  return np.asarray(z,dtype=float)
 a,b=continuous(train),continuous(test)
 for j in range(a.shape[1]):
  finite=np.isfinite(a[:,j]);med=np.median(a[finite,j]) if finite.any() else 0
  a[~finite,j]=med;b[~np.isfinite(b[:,j]),j]=med
 mean=a.mean(axis=0);std=a.std(axis=0);std[std<1e-12]=1
 a=(a-mean)/std;b=(b-mean)/std
 for col in ['order','group']:
  levels=sorted({byrep[r][col] for r in train})
  a=np.column_stack([a]+[np.array([byrep[r][col]==level for r in train],float) for level in levels])
  b=np.column_stack([b]+[np.array([byrep[r][col]==level for r in test],float) for level in levels])
 return np.column_stack([np.ones(len(a)),a]),np.column_stack([np.ones(len(b)),b])
def predict(train,test,y,features):
 a,b=design(train,test,features);pen=np.eye(a.shape[1])*10;pen[0,0]=0
 return b@np.linalg.solve(a.T@a+pen,a.T@np.asarray(y))
def outcomes(eps):
 bysite=defaultdict(list)
 for r in pairs:
  if max(r['root_pct'],r['bulk_pct'])>0:
   bysite[r['rep'],r['site']].append(math.log2((r['root_pct']+eps)/(r['bulk_pct']+eps)))
 return {k:float(np.mean(v)) for k,v in bysite.items()}
def evaluate(reps,features,ys):
 allmean={r:float(np.mean([v for (p,s),v in ys.items() if p==r])) for r in reps}
 rows=[]
 for kind in ['genus','site']:
  for f in range(5):
   if kind=='genus':
    train=[r for r in reps if gfold[genus(r)]!=f];test=[r for r in reps if gfold[genus(r)]==f]
    train_y=[allmean[r] for r in train];test_y=[allmean[r] for r in test]
   else:
    tr={r:[v for (p,s),v in ys.items() if p==r and sfold[s]!=f] for r in reps}
    te={r:[v for (p,s),v in ys.items() if p==r and sfold[s]==f] for r in reps}
    train=[r for r in reps if len(tr[r])>=3];test=[r for r in train if len(te[r])>=3]
    train_y=[float(np.mean(tr[r])) for r in train];test_y=[float(np.mean(te[r])) for r in test]
   if len(train)<10 or not test:continue
   p0=predict(train,test,train_y,None);p1=predict(train,test,train_y,features)
   for r,y,b,a in zip(test,test_y,p0,p1):rows.append(dict(kind=kind,fold=f,rep=r,observed=y,baseline=float(b),augmented=float(a)))
 return rows
def metrics(rows):
 result=[]
 for kind in ['genus','site']:
  for f in ['all']+list(range(5)):
   z=[r for r in rows if r['kind']==kind and (f=='all' or r['fold']==f)]
   if not z:continue
   # Repeated site-fold predictions get equal total weight per MAG.
   unit=defaultdict(list)
   for r in z:unit[r['rep']].append(r)
   def error(pred,absolute=False):return float(np.mean([np.mean([abs(r['observed']-r[pred]) if absolute else (r['observed']-r[pred])**2 for r in a]) for a in unit.values()]))
   b,a=error('baseline'),error('augmented')
   result.append(dict(kind=kind,fold=f,n_models=len(unit),n_predictions=len(z),baseline_mse=b,augmented_mse=a,mse_improvement=1-a/b if b>0 else 0,baseline_mae=error('baseline',True),augmented_mae=error('augmented',True)))
 return result
all_metrics=[];all_preds=[];nulls=[];analysis_counts=[]
for variant in ['filled','draft']:
 t={r['rep']:r for r in traits if r['variant']==variant and r['trait_qc']}
 reps=sorted(set(eligible)&set(t));features={r:np.array([t[r]['atp_breadth']/12,math.log1p(max(0,t[r]['mixture_atp']))]) for r in reps}
 analysis_counts.append(dict(variant=variant,models=len(reps),genera=len({genus(r) for r in reps}),orders=len({byrep[r]['order'] for r in reps})))
 if len(reps)<10:continue
 for eps in [.001,.0001,.01]:
  ys=outcomes(eps);pred=evaluate(reps,features,ys);met=metrics(pred)
  all_preds += [dict(variant=variant,epsilon=eps,**r) for r in pred]
  all_metrics += [dict(variant=variant,epsilon=eps,**r) for r in met]
  if variant=='filled' and eps==.001:
   rng=np.random.default_rng(20261012);orders=defaultdict(list)
   for r in reps:orders[byrep[r]['order']].append(r)
   for perm in range(199):
    shuffled={}
    for z in orders.values():
     for dst,src in zip(z,rng.permutation(z)):shuffled[dst]=features[src]
    pm=metrics(evaluate(reps,shuffled,ys));v={r['kind']:r['mse_improvement'] for r in pm if r['fold']=='all'}
    nulls.append(dict(permutation=perm,genus=v.get('genus',float('nan')),site=v.get('site',float('nan')),combined=np.mean(list(v.values()))))
write(P/'prediction_metrics.tsv',all_metrics);write(P/'predictions.tsv',all_preds);write(P/'permutation_null.tsv',nulls);write(P/'analysis_counts.tsv',analysis_counts)
main=[r for r in all_metrics if r['variant']=='filled' and r['epsilon']==.001];mainall={r['kind']:r for r in main if r['fold']=='all'}
count=next(r for r in analysis_counts if r['variant']=='filled')
qc_fraction=sum(r['trait_qc'] and r['atp_breadth']>0 for r in traits if r['variant']=='filled')/len(panel)
checks=dict(model_qc=qc_fraction>=.8,field_coverage=count['models']>=30 and count['genera']>=15 and count['orders']>=6,
 cv_coverage=all(k in mainall and mainall[k]['n_models']>=30 and sum(r['kind']==k and r['fold']!='all' for r in main)==5 for k in ['genus','site']),
 predictive_gain=all(k in mainall and mainall[k]['mse_improvement']>=.05 for k in ['genus','site']),
 fold_consistency=all(sum(r['kind']==k and r['fold']!='all' and r['mse_improvement']>0 for r in main)>=3 for k in ['genus','site']))
combined=float(np.mean([r['mse_improvement'] for r in mainall.values()])) if len(mainall)==2 else None
q90=float(np.quantile([r['combined'] for r in nulls],.9)) if nulls else None
checks['permutation']=combined is not None and q90 is not None and combined>q90
sensitivity=[r for r in all_metrics if r['variant']=='filled' and r['epsilon']!=.001 and r['fold']=='all']
checks['pseudocount']=len(sensitivity)==4 and all(r['mse_improvement']>0 for r in sensitivity)
s=dict(panel=len(panel),paired_visits=len(visits),samples=len(samples),missing_coverm_samples=missing_coverm,sites=len({s for s,t in visits}),scenario_LPs=len(sc),statuses=dict(Counter(r['status'] for r in sc)),qc_fraction=qc_fraction,
 analysis_counts=analysis_counts,checks=checks,scale_up=all(checks.values()),combined_improvement=combined,null_90th_percentile=q90,
 permutation_tail_fraction=(1+sum(r['combined']>=combined for r in nulls))/(1+len(nulls)) if nulls and combined is not None else None)
(P/'summary.json').write_text(json.dumps(s,indent=2)+'\n');print(json.dumps(s,indent=2))
write(P/'provenance.tsv',[dict(path=str(p),sha256=sha(p)) for p in sorted((W/'inputs').iterdir()) if p.is_file()]+[dict(path=str(p),sha256=sha(p)) for p in sorted((W/'scripts').glob('*.py'))])
