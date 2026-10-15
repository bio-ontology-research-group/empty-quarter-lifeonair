#!/usr/bin/env python3
"""Secondary check of existing gapseq models, gated by exact contig sequences."""
import sys,json,hashlib
from pathlib import Path
from metabolic_simulate import read,write,fasta,digest,base,exchanges,set_medium,atpm,solve,MINERALS,H2,SUBSTRATES
from cobra.io import read_sbml_model
from cobra.util.solver import linear_reaction_coefficients
R=Path(sys.argv[1]);W=R/'06_atmospheric/metabolic_pilot_clean'
panel=read(W/'results/pilot_panel.tsv')
sources={}
for line in (R/'05_metabolism/hq_mag_list.tsv').read_text().splitlines():
 p,s,b,*_=line.split('\t');sources[s+'__'+b]=Path(p)
def seqhash(p):
 hashes=sorted(hashlib.sha256(s.encode()).hexdigest() for _,s in fasta(p))
 return hashlib.sha256('\n'.join(hashes).encode()).hexdigest()
audits=[];rxs=[];out=[]
for row in panel:
 rep=row['rep']
 if '__single_' not in rep:continue
 tag=rep.replace('__single_','__');d=R/'05_metabolism/gapseq_out'/tag
 xml=d/(tag+'.xml');orig=sources.get(tag);current=R/'05_refine/_drep/out/dereplicated_genomes'/(rep+'.fa')
 if not xml.exists() or orig is None:
  audits.append(dict(rep=rep,group=row['group'],status='missing_source_or_model',source='',original_contigs_sha256='',current_contigs_sha256='',model_sha256=''));continue
 h1,h2=seqhash(orig),seqhash(current)
 audits.append(dict(rep=rep,group=row['group'],status='identical' if h1==h2 else 'sequence_mismatch',source=str(orig),original_contigs_sha256=h1,current_contigs_sha256=h2,model_sha256=digest(xml)))
 if h1!=h2:continue
 for variant,suffix in [('draft','-draft.xml'),('filled','.xml')]:
  m=read_sbml_model(str(d/(tag+suffix)));m.solver='glpk';m.solver.configuration.timeout=30
  bio=next(iter(linear_reaction_coefficients(m)));ex=exchanges(m);exids={r.id for r in ex}
  hyd=[r for r in m.reactions if r.id not in exids and any(base(x)==H2 for x in r.metabolites) and len({base(x) for x in r.metabolites})>1]
  for r in hyd:rxs.append(dict(rep=rep,group=row['group'],variant=variant,reaction=r.id,name=r.name,stoichiometry=r.reaction,gpr=r.gene_reaction_rule))
  for r in m.reactions:
   if (r.boundary and r.id not in exids) or 'biomass' in r.name.lower() or r.id.startswith('bio') or r==bio:r.bounds=(0,0)
  m.objective=atpm(m)
  def rec(scenario):
   s,v=solve(m);out.append(dict(rep=rep,group=row['group'],variant=variant,scenario=scenario,status=s,value='' if v is None else v))
  for r in ex:r.bounds=(0,0)
  rec('closed')
  set_medium(m,ex,MINERALS);rec('inorganic')
  set_medium(m,ex,{**MINERALS,H2:.1});rec('H2')
  for name,cpd,nc in SUBSTRATES:
   set_medium(m,ex,{**MINERALS,cpd:6/nc});rec(name)
 print(rep,flush=True)
for name,rows in [('gapseq_identity',audits),('gapseq_hydrogen_reactions',rxs),('gapseq_scenarios',out)]:write(W/'results'/(name+'.tsv'),rows)
