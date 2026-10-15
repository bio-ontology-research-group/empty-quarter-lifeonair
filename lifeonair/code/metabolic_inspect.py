#!/usr/bin/env python3
"""Read a few matched models; report exact reaction conventions before simulations."""
import csv,json,sys,hashlib
from pathlib import Path
import cobra
from cobra.io import read_sbml_model
root=Path(sys.argv[1]); models=Path(sys.argv[2]); out=[]
with (root/'results/pilot_panel.tsv').open() as f:panel=list(csv.DictReader(f,delimiter='\t'))
for row in panel[:3]:
 rep=row['rep'];d=models/rep
 m=read_sbml_model(str(d/(rep+'-filled.xml')))
 z=dict(rep=rep,group=row['group'],n_reactions=len(m.reactions),n_genes=len(m.genes),objective=str(m.objective.expression),compartments=m.compartments,
        notes=m.notes,hydrogen_metabolites=[],hydrogen_reactions=[],atp_reactions=[],boundaries=[],genes=[g.id for g in list(m.genes)[:5]])
 for met in m.metabolites:
  if 'cpd11640' in met.id or met.formula=='H2':
   z['hydrogen_metabolites'].append(dict(id=met.id,name=met.name,formula=met.formula))
   for r in met.reactions:z['hydrogen_reactions'].append(dict(id=r.id,name=r.name,reaction=r.reaction,bounds=r.bounds,gpr=r.gene_reaction_rule,notes=r.notes))
 for r in m.reactions:
  if 'rxn00062' in r.id or 'maintenance' in r.name.lower():z['atp_reactions'].append(dict(id=r.id,name=r.name,reaction=r.reaction,bounds=r.bounds,gpr=r.gene_reaction_rule))
  if r.boundary:z['boundaries'].append(dict(id=r.id,name=r.name,reaction=r.reaction,bounds=r.bounds))
 out.append(z)
(root/'results/model_inspection.json').write_text(json.dumps(dict(cobra=cobra.__version__,models=out),indent=2)+'\n')
print('Inspected',len(out),'models')
