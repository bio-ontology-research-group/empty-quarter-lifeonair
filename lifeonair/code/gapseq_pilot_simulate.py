#!/usr/bin/env python3
"""Gapseq-only energy and carbon-use pilot. Run inside a Slurm CPU allocation."""
import sys,json,hashlib
from pathlib import Path
from cobra.io import read_sbml_model
from cobra.util.solver import linear_reaction_coefficients
from metabolic_simulate import read,write,digest,base,exchanges,set_medium,atpm,solve,MINERALS,GROWTH_MINERALS,H2,SUBSTRATES
R=Path(sys.argv[1]);W=R/'06_atmospheric/gapseq_traits_pilot_clean';index=int(sys.argv[2])
panel=read(W/'results/pilot_panel.tsv')
if index>=len(panel):sys.exit(0)
r=panel[index];out=[];audit=[];hashes=[];rx=[]
for variant in ['draft','filled']:
 path=Path(r[variant+'_path']);hashes.append(dict(rep=r['rep'],variant=variant,path=str(path),sha256=digest(path)))
 m=read_sbml_model(str(path));m.solver='glpk';m.solver.configuration.timeout=30
 original_status,original_growth=solve(m)
 obj=linear_reaction_coefficients(m);assert len(obj)==1;bio=next(iter(obj));ex=exchanges(m);eid={x.id for x in ex}
 biomass=[x for x in m.reactions if 'biomass' in x.name.lower() or x.id.startswith('bio') or x==bio]
 bounds={x.id:(max(0,x.lower_bound),max(0,x.upper_bound)) for x in biomass}
 for x in m.reactions:
  if x.boundary and x.id not in eid:x.bounds=(0,0)
 def growth(on):
  for x in biomass:x.bounds=bounds[x.id] if on else (0,0)
 hyd=[x for x in m.reactions if x.id not in eid and any(base(z)==H2 for z in x.metabolites) and len({base(z) for z in x.metabolites})>1]
 for x in hyd:rx.append(dict(rep=r['rep'],variant=variant,reaction=x.id,gpr=x.gene_reaction_rule,stoichiometry=x.reaction))
 maintenance=atpm(m)
 def rec(kind,name,carbon,oxygen,objective):
  s,v=solve(m);out.append(dict(rep=r['rep'],group=r['group'],variant=variant,kind=kind,name=name,carbon=carbon,oxygen=oxygen,objective=objective,status=s,value='' if v is None else v));return s,v
 growth(False);m.objective=maintenance
 for x in ex:x.bounds=(0,0)
 closed_s,closed=rec('control','closed',0,0,'ATPM')
 set_medium(m,ex,MINERALS);min_s,min_v=rec('control','inorganic',0,10,'ATPM')
 set_medium(m,ex,{**MINERALS,H2:.1});rec('control','hydrogen',0,10,'ATPM')
 # Permissive growth is a feasibility check, never a biological phenotype.
 set_medium(m,ex,{base(next(iter(x.metabolites))):10 for x in ex});growth(True);m.objective=bio
 rich_s,rich=rec('control','permissive',-1,-1,'biomass');growth(False)
 for name,cpd,nc in SUBSTRATES:
  set_medium(m,ex,{**MINERALS,cpd:6/nc});m.objective=maintenance
  rec('substrate',name,6,10,'ATPM')
  set_medium(m,ex,{**GROWTH_MINERALS,cpd:6/nc});growth(True);m.objective=bio
  rec('substrate',name,6,10,'biomass');growth(False)
 for carbon in [.06,.6,6]:
  for oxygen in [.1,1,10]:
   medium={**MINERALS,'cpd00007':oxygen,**{cpd:carbon/(len(SUBSTRATES)*nc) for _,cpd,nc in SUBSTRATES}}
   set_medium(m,ex,medium);m.objective=maintenance;rec('mixture','equal_C',carbon,oxygen,'ATPM')
   set_medium(m,ex,{**GROWTH_MINERALS,**medium});growth(True);m.objective=bio;rec('mixture','equal_C',carbon,oxygen,'biomass');growth(False)
 audit.append(dict(rep=r['rep'],group=r['group'],variant=variant,reactions=len(m.reactions),genes=len(m.genes),
  formulas=sum(bool(x.formula) for x in m.metabolites),metabolites=len(m.metabolites),original_status=original_status,original_growth=original_growth,
  closed_status=closed_s,closed_atp=closed,inorganic_status=min_s,inorganic_atp=min_v,rich_status=rich_s,rich_growth=rich,
  energy_qc=int(closed_s==min_s=='optimal' and abs(closed)<1e-6 and abs(min_v)<1e-6)))
for name,z in [('scenarios',out),('audit',audit),('hashes',hashes),('hydrogen_reactions',rx)]:write(W/'results'/f'{name}_{index}.tsv',z)
(W/'results'/f'runtime_{index}.json').write_text(json.dumps(dict(rep=r['rep'],script_sha256=digest(Path(__file__)),helper_sha256=digest(Path(__file__).parent/'metabolic_simulate.py'),python=sys.version,solver='glpk'),indent=2)+'\n')
print(index,r['rep'],len(out),flush=True)
