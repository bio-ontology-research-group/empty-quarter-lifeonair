#!/usr/bin/env python3
"""Matched-model pilot: provenance, conservation controls and substrate/H2 scenarios.

Run in a CPU allocation. Source reconstructions are never modified on disk.
Uptake bounds are scenarios in mmol/gDW/h, not atmospheric kinetic estimates.
"""
import argparse,csv,gzip,hashlib,json,logging,math,re,sys,time
from pathlib import Path
from collections import defaultdict
import cobra
from cobra.io import read_sbml_model
from cobra.util.solver import linear_reaction_coefficients

logging.getLogger('cobra').setLevel(logging.ERROR)
SUBSTRATES=[('glucose','cpd00027',6),('fructose','cpd00082',6),('sucrose','cpd00076',12),
 ('acetate','cpd00029',2),('pyruvate','cpd00020',3),('succinate','cpd00036',4),
 ('malate','cpd00130',4),('citrate','cpd00137',6),('alanine','cpd00035',3),
 ('glutamate','cpd00023',5),('aspartate','cpd00041',4),('glycine','cpd00033',2)]
MINERALS={'cpd00001':100,'cpd00067':100,'cpd00009':10,'cpd00013':10,'cpd00048':10,
 'cpd00099':10,'cpd00205':10,'cpd00971':10,'cpd00063':10,'cpd00254':10,'cpd00011':1,'cpd00007':10}
H2='cpd11640'; TOL=1e-7
GROWTH_MINERALS={**MINERALS,'cpd00030':.01,'cpd00034':.01,'cpd00058':.01,
 'cpd00149':.01,'cpd00244':.01,'cpd10515':.01,'cpd10516':.01}

def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def fasta(path):
 op=gzip.open if str(path).endswith('.gz') else open
 with op(path,'rt') as f:
  name=None;seq=[]
  for line in f:
   if line.startswith('>'):
    if name is not None:yield name,''.join(seq).rstrip('*').upper()
    name=line[1:].split()[0];seq=[]
   else:seq.append(line.strip())
  if name is not None:yield name,''.join(seq).rstrip('*').upper()
def base(m):
 z=re.search(r'cpd\d+',m.id);return z.group() if z else m.id
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,rows):
 if not rows:return
 with p.open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def solve(m):
 try:
  v=m.slim_optimize(error_value=float('nan'));status=m.solver.status
  return status,float(v) if status=='optimal' and math.isfinite(v) else None
 except Exception as e:return 'error:'+type(e).__name__,None
def exchanges(m):
 # Only single-metabolite reactions in extracellular compartments are exchanges.
 return [r for r in m.reactions if len(r.metabolites)==1 and next(iter(r.metabolites)).compartment.startswith('e')]
def set_medium(m,ex,medium):
 for r in ex:
  met,coef=next(iter(r.metabolites.items()));u=medium.get(base(met),0)
  r.bounds=(-u/abs(coef),1000) if coef<0 else (-1000,u/abs(coef))
def atpm(m):
 # Add only an energy-dissipating diagnostic objective; never a source reaction.
 sto={'cpd00001':-1,'cpd00002':-1,'cpd00008':1,'cpd00009':1,'cpd00067':1}
 d={base(x):x for x in m.metabolites if x.compartment=='c0' and base(x) in sto}
 if len(d)!=5:raise ValueError('Missing ATP hydrolysis metabolites')
 r=cobra.Reaction('P3_ATPM');r.name='Diagnostic ATP maintenance';r.bounds=(0,1000)
 r.add_metabolites({d[k]:v for k,v in sto.items()});m.add_reactions([r]);return r

def run(a):
 panel=read(a.work/'results/pilot_panel.tsv')
 if a.index is not None:panel=[panel[a.index]]
 ids={r['rep'] for r in panel}
 # Exact amino-acid identity links curated Bakta markers to reconstruction Prodigal genes.
 typed=defaultdict(dict)
 for r in read(a.work/'inputs/atmo_classification_v2.tsv'):
  if r['rep_id'] in ids and r['marker']=='hydB' and r['assignment'] in ['[NiFe]_Group_1l','[NiFe]_Group_1h']:
   typed[r['rep_id']][r['protein_id']]=r['assignment']
 seqs=defaultdict(lambda:defaultdict(set))
 for name,seq in fasta(a.marker_fasta):
  z=name.split('|')
  if len(z)>=4 and z[0]=='EQ' and z[1]=='hydB' and z[2] in typed and z[3] in typed[z[2]]:seqs[z[2]][seq].add(z[3])
 results=[];audits=[];rxrows=[];mapping=[];hashes=[];transports=[]
 for row in panel:
  rep=row['rep'];print(rep,flush=True);d=a.models/rep;start=time.time()
  faa=a.proteomes/(rep+'.faa.gz');genes=set();protein_ids=set();nproteins=0
  for gene,seq in fasta(faa):
   nproteins+=1;protein_ids.add(gene)
   if seq in seqs[rep]:
    genes.add(gene)
    for pid in sorted(seqs[rep][seq]):mapping.append(dict(rep=rep,model_gene=gene,typed_protein=pid,assignment=typed[rep][pid],match='exact_aa'))
  ledger=read(d/(rep+'-filled-added.tsv'));added={r['rxn_id']:r['step'] for r in ledger}
  for p in [faa,d/(rep+'-draft.xml'),d/(rep+'-filled.xml'),d/(rep+'-filled-added.tsv')]:hashes.append(dict(rep=rep,path=str(p),sha256=digest(p)))
  for variant in ['draft','filled','filled_grounded_h2']:
   src='draft' if variant=='draft' else 'filled'
   m=read_sbml_model(str(d/(rep+'-'+src+'.xml')));m.solver='glpk';m.solver.configuration.timeout=30
   biomass=[r for r in m.reactions if 'biomass' in r.name.lower() or r.id.startswith('bio')]
   objective=linear_reaction_coefficients(m)
   if len(objective)!=1:raise ValueError('Expected one biomass objective')
   bio=next(iter(objective));original_status,original_growth=solve(m)
   # gapsmith maximises a drain of biomass pseudometabolite cpd11416, fed by
   # bio1. Both the synthesis reaction and its drain must reopen for growth.
   if bio not in biomass:biomass.append(bio)
   bio_bounds={r.id:(max(0,r.lower_bound),max(0,r.upper_bound)) for r in biomass}
   def growth_mode(on):
    for r in biomass:r.bounds=bio_bounds[r.id] if on else (0,0)
   ex=exchanges(m);exids={r.id for r in ex}
   for name,cpd,nc in SUBSTRATES:
    tx=[r for r in m.reactions if r.id not in exids and len({x.compartment for x in r.metabolites if base(x)==cpd})>1]
    transports.append(dict(rep=rep,variant=variant,substrate=name,compound=cpd,carbon_atoms=nc,
                           exchange_present=int(any(base(next(iter(r.metabolites)))==cpd for r in ex)),
                           transport_reactions=len(tx),transport_with_gpr=sum(bool(r.genes) for r in tx),
                           transport_not_gapfilled=sum(r.id not in added for r in tx),reaction_ids=';'.join(r.id for r in tx)))
   hyd=[r for r in m.reactions if r.id not in exids and any(base(x)==H2 for x in r.metabolites) and len({base(x) for x in r.metabolites})>1]
   native=[r for r in hyd if {g.id for g in r.genes}&genes and r.id not in added]
   for r in hyd:
    rxrows.append(dict(rep=rep,variant=variant,reaction=r.id,name=r.name,stoichiometry=r.reaction,lower=r.lower_bound,upper=r.upper_bound,
                      gpr=r.gene_reaction_rule,gapfill_step=added.get(r.id,''),exact_marker_gpr=int(bool({g.id for g in r.genes}&genes))))
   if variant=='filled_grounded_h2':
    keep={r.id for r in native}
    for r in hyd:
     if r.id not in keep:r.bounds=(0,0)
   # Internal sinks/demands are closed; biomass cannot run in reverse.
   for r in m.reactions:
    if r.boundary and r.id not in exids:r.bounds=(0,0)
   for r in biomass:r.bounds=(0,0)
   maintenance=atpm(m);m.objective=maintenance
   for r in ex:r.bounds=(0,0)
   closed_status,closed_atp=solve(m)
   set_medium(m,ex,MINERALS);inorg_status,inorg_atp=solve(m)
   def record(kind,substrate,carbon,h2,obj,status,value):
    results.append(dict(rep=rep,block=row['block'],group=row['group'],variant=variant,kind=kind,substrate=substrate,carbon_input=carbon,h2_bound=h2,objective=obj,status=status,value='' if value is None else value))
   record('control','closed',0,0,'ATPM',closed_status,closed_atp)
   record('control','minerals',0,0,'ATPM',inorg_status,inorg_atp)
   # Permissive nutrient positive control checks structural biomass functionality.
   rich={base(next(iter(r.metabolites))):10 for r in ex};rich.pop(H2,None)
   set_medium(m,ex,rich);growth_mode(True);m.objective=bio
   rich_status,rich_growth=solve(m);record('control','all_exchanges',0,0,'biomass',rich_status,rich_growth)
   growth_mode(False);m.objective=maintenance
   audit=dict(rep=rep,group=row['group'],variant=variant,n_reactions=len(m.reactions)-1,n_genes=len(m.genes),n_proteins=nproteins,
              model_genes_absent_from_input=len({g.id for g in m.genes}-protein_ids),metabolites_with_formula=sum(bool(x.formula) for x in m.metabolites),
              typed_markers=len(typed[rep]),markers_with_sequence=sum(map(len,seqs[rep].values())),exact_mapped_genes=len(genes),h2_reactions=len(hyd),native_typed_h2_reactions=len(native),
              h2_exchange=int(any(base(next(iter(r.metabolites)))==H2 for r in ex)),gapfill_reactions=len(added),
              original_status=original_status,original_growth=original_growth,closed_status=closed_status,closed_atp=closed_atp,
              inorganic_status=inorg_status,inorganic_atp=inorg_atp,rich_status=rich_status,rich_growth=rich_growth,
              energy_qc_pass=int(closed_status=='optimal' and inorg_status=='optimal' and closed_atp is not None and inorg_atp is not None and abs(closed_atp)<TOL and abs(inorg_atp)<TOL))
   audits.append(audit)
   # Carbon panel at 6 mmol C/gDW/h, no hydrogen; both energy and growth outputs.
   for name,cpd,nc in SUBSTRATES:
    med={**MINERALS,cpd:6/nc};set_medium(m,ex,med)
    m.objective=maintenance;status,val=solve(m);record('substrate',name,6,0,'ATPM',status,val)
    set_medium(m,ex,{**GROWTH_MINERALS,cpd:6/nc})
    m.objective=bio;growth_mode(True);status,val=solve(m);record('substrate',name,6,0,'biomass',status,val);growth_mode(False)
   # Fixed glucose carbon gradient with optional hydrogen bounds.
   m.objective=maintenance
   for carbon in [0,.006,.06,.6,6]:
    for h in [0,.001,.01,.1]:
     med={**MINERALS,'cpd00027':carbon/6,H2:h};set_medium(m,ex,med)
     status,val=solve(m);record('gradient','glucose',carbon,h,'ATPM',status,val)
   # Oxygen-dependence and all-H2-reaction knockout controls at zero organic carbon.
   med={**MINERALS,H2:.1,'cpd00007':0};set_medium(m,ex,med)
   status,val=solve(m);record('control','H2_anoxic',0,.1,'ATPM',status,val)
   set_medium(m,ex,{**MINERALS,H2:.1})
   with m:
    for r in hyd:r.bounds=(0,0)
    status,val=solve(m);record('control','H2_reaction_knockout',0,.1,'ATPM',status,val)
  print('seconds',round(time.time()-start,2),flush=True)
 out=a.work/'results/simulations';out.mkdir(exist_ok=True)
 tag='all' if a.index is None else str(a.index)
 for name,rs in [('scenarios',results),('audit',audits),('hydrogen_reactions',rxrows),('marker_mapping',mapping),('input_hashes',hashes),('substrate_transport',transports)]:write(out/(name+'_'+tag+'.tsv'),rs)
 (out/('runtime_'+tag+'.json')).write_text(json.dumps(dict(cobra=cobra.__version__,python=sys.version,solver='glpk',panel=len(panel),script_sha256=digest(Path(__file__))),indent=2)+'\n')

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--work',type=Path,required=True);p.add_argument('--models',type=Path,required=True)
 p.add_argument('--proteomes',type=Path,required=True);p.add_argument('--marker-fasta',type=Path,required=True);p.add_argument('--index',type=int)
 run(p.parse_args())
