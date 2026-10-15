#!/usr/bin/env python3
"""Standalone pilot QC figure; no changes to manuscript figures."""
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root=Path(__file__).resolve().parents[1]/'results/metabolic_pilot'
with (root/'model_traits.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
with (root/'audit.tsv').open() as f:audit=list(csv.DictReader(f,delimiter='\t'))
groups=['1l','1h','neither'];colors=['#aa3e39','#297a9c','#68736d']
fig,axs=plt.subplots(1,3,figsize=(11,3.6),layout='constrained')
rng=np.random.default_rng(20261012)
for i,g in enumerate(groups):
 z=[r for r in rows if r['group']==g and r['variant']=='filled' and r['energy_qc_pass']=='1']
 axs[0].scatter(i+rng.uniform(-.16,.16,len(z)),[float(r['carbon_ATP_breadth']) for r in z],s=25,color=colors[i],alpha=.7)
 med=np.median([float(r['carbon_ATP_breadth']) for r in z]) if z else np.nan
 axs[0].plot([i-.22,i+.22],[med,med],color='black',lw=2)
 a=[r for r in audit if r['group']==g and r['variant']=='filled']
 mapped=sum(int(r['exact_mapped_genes'])>0 for r in a)
 linked=sum(int(r['native_typed_h2_reactions'])>0 for r in a)
 axs[1].bar(i-.17,mapped,width=.32,color=colors[i],label='Marker matched to model input' if i==0 else None)
 axs[1].bar(i+.17,linked,width=.32,color=colors[i],hatch='///',alpha=.45,label='Linked to native H₂ reaction' if i==0 else None)
 axs[1].text(i+.17,linked+.5,str(linked),ha='center',fontsize=9)
 for j,variant in enumerate(['draft','filled','filled_grounded_h2']):
  q=[r for r in rows if r['group']==g and r['variant']==variant and r['energy_qc_pass']=='1']
  n=sum(float(r['hydrogen_ATP_benefit'] or 0)>1e-6 for r in q)
  x=i+(j-1)*.23
  axs[2].bar(x,n,width=.21,color=colors[i],alpha=[.4,1,.65][j],hatch=['','', '//'][j])
  axs[2].text(x,n+.3,str(n),ha='center',fontsize=8)
for ax in axs:
 ax.set_xticks(range(3),['Group 1l','Group 1h','Neither'])
 ax.spines[['top','right']].set_visible(False)
axs[0].set(title='A  Carbon-supported ATP maintenance',ylabel='Substrates supporting ATP / 12',ylim=(-.5,12.8))
axs[1].set(title='B  Hydrogenase representation',ylabel='Models / 24',ylim=(0,30))
axs[1].legend(frameon=False,fontsize=7,loc='upper left')
axs[2].set(title='C  Hydrogen benefit in simulations',ylabel='Models with positive ΔATP maintenance',ylim=(0,27))
axs[2].text(.5,-.2,'Within each group: draft, filled, restricted H₂',transform=axs[2].transAxes,ha='center',fontsize=8)
fig.savefig(root/'pilot_summary.pdf');fig.savefig(root/'pilot_summary.png',dpi=180)
