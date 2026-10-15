#!/usr/bin/env python3
"""Plot the frozen small gapseq predictive pilot, without fitting models."""
import csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parents[1]/'results/gapseq_traits_pilot'
def read(n):
 with (P/n).open() as f:return list(csv.DictReader(f,delimiter='\t'))
t=read('traits.tsv');m=read('prediction_metrics.tsv');null=read('permutation_null.tsv');s=json.loads((P/'summary.json').read_text())
fig,ax=plt.subplots(1,3,figsize=(11,3.6),layout='constrained');rng=np.random.default_rng(42)
for i,g in enumerate(['1l','1h','neither']):
 v=[int(r['atp_breadth']) for r in t if r['group']==g and r['variant']=='filled' and r['trait_qc']=='1']
 ax[0].scatter(i+rng.uniform(-.15,.15,len(v)),v,s=24,color=['#b25442','#337d9c','#66786a'][i],alpha=.7)
 if v:ax[0].plot([i-.2,i+.2],[np.median(v)]*2,color='black')
ax[0].set(xticks=[0,1,2],xticklabels=['Group 1l','Group 1h','Neither'],ylabel='Substrates supporting ATP / 12',ylim=(-.5,12.5),title='A  Gapseq carbon-supported ATP')
for j,e in enumerate([.0001,.001,.01]):
 vals=[]
 for kind in ['genus','site']:
  r=next(r for r in m if r['variant']=='filled' and float(r['epsilon'])==e and r['kind']==kind and r['fold']=='all');vals.append(float(r['mse_improvement'])*100)
 ax[1].bar(np.arange(2)+(j-1)*.25,vals,.23,label=f'ε = {e:g}',color=['#92abb5','#28617a','#bb775c'][j])
ax[1].axhline(0,color='black',lw=.8);ax[1].axhline(5,color='black',ls=':',lw=.8)
ax[1].set(xticks=[0,1],xticklabels=['Held-out genera','Held-out sites'],ylabel='MSE reduction versus baseline (%)',title='B  Predictive gain from two traits')
ax[1].legend(frameon=False,fontsize=8)
a=[float(r['combined'])*100 for r in null];ax[2].hist(a,bins=20,color='#a5b8ad',edgecolor='white')
ax[2].axvline(s['combined_improvement']*100,color='#a83830',lw=2,label='Observed')
ax[2].axvline(s['null_90th_percentile']*100,color='black',ls=':',label='Shuffle 90th percentile')
ax[2].set(xlabel='Mean predictive improvement (%)',ylabel='Within-order trait shuffles',title='C  Screening reference (199 shuffles)')
ax[2].legend(frameon=False,fontsize=8)
for a in ax:
 a.spines[['top','right']].set_visible(False)
 a.title.set_fontsize(10)
fig.savefig(P/'pilot_prediction.pdf');fig.savefig(P/'pilot_prediction.png',dpi=180)
