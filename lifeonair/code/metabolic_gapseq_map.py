#!/usr/bin/env python3
"""Map coordinate-based gapseq H2 GPRs to curated proteins by translated sequence."""
import itertools,re,sys
from pathlib import Path
from collections import defaultdict
from metabolic_simulate import read,write,fasta
R=Path(sys.argv[1]);W=R/'06_atmospheric/metabolic_pilot_clean'
typed={(r['rep_id'],r['protein_id']):r['assignment'] for r in read(W/'inputs/atmo_classification_v2.tsv') if r['marker']=='hydB'}
seqs=defaultdict(lambda:defaultdict(list))
for name,seq in fasta(R/'06_atmospheric/lifeonair_reads/db/markers.faa'):
 z=name.split('|')
 if len(z)>=4 and (z[2],z[3]) in typed:seqs[z[2]][seq].append((z[3],typed[(z[2],z[3])]))
table=dict(zip(map(''.join,itertools.product('TCAG',repeat=3)),'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))
rows=read(W/'results/gapseq_hydrogen_reactions.tsv');out=[];cache={}
for r in rows:
 rep=r['rep']
 if rep not in cache:cache[rep]=dict(fasta(R/'05_refine/_drep/out/dereplicated_genomes'/(rep+'.fa')))
 genes=re.findall(r'[^\s()]+',r['gpr']);genes=[g for g in genes if g not in ['and','or']]
 for gene in genes:
  try:
   contig,start,end=gene.rsplit('_',2);start,end=int(start),int(end)
   nt=cache[rep][contig][min(start,end)-1:max(start,end)]
   if start>end:nt=nt.translate(str.maketrans('ACGT','TGCA'))[::-1]
   seq=''.join(table.get(nt[i:i+3],'X') for i in range(0,len(nt)-2,3)).rstrip('*')
   if nt[:3] in ('ATG','GTG','TTG','CTG','ATT','ATC','ATA') and seq:seq='M'+seq[1:]
   hits=seqs[rep].get(seq,[]);status='exact_translated_match' if hits else 'no_curated_hydrogenase_match'
  except (ValueError,KeyError):hits=[];status='unresolved_coordinates'
  out.append(dict(rep=rep,group=r['group'],variant=r['variant'],reaction=r['reaction'],model_gene=gene,status=status,
                  typed_proteins=';'.join(x[0] for x in hits),assignments=';'.join(x[1] for x in hits)))
write(W/'results/gapseq_marker_mapping.tsv',out)
