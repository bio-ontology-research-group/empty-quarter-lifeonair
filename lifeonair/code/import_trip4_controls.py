#!/usr/bin/env python3
"""Import the September control replay from the companion repository.

Run from the repository root. This optional provenance import needs the companion
repository; normal P3 reanalysis uses the committed candidate files and union.
"""
from pathlib import Path
import csv,hashlib,json,re
src=Path('../empty-quarter-data-paper/evidence/control-replay-20260909/screen/primary_contaminant_calls.tsv')
rows=[r for r in csv.DictReader(src.open(),delimiter='\t') if r['screen']=='Trip4']
out=Path('lifeonair/results/review_r3');out.mkdir(exist_ok=True)
with (out/'trip4_candidate_calls.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=rows[0],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
ids={r['feature_id'] for r in rows};fasta=Path('../empty-quarter-data-paper/metadata/taxonomy/ASV_seqs-trips1-5.fasta')
assert len(ids)==7
found=set();keep=False
with fasta.open() as f, (out/'trip4_candidate_sequences.fasta').open('w') as o:
 for line in f:
  if line.startswith('>'):
   id=line[1:].split()[0];keep=id in ids
   if keep:found.add(id)
  if keep:o.write(line)
assert found==ids, 'Missing candidate sequences'
genera={r['taxonomy'].split(';')[5] for r in rows}-{ 'NA',''}
old=set(Path('lifeonair/results/reads/excl_v2_noflagged.txt').read_text().split())
matched=[]
for fn in ['r232_bac120.summary.tsv','r232_ar53.summary.tsv']:
 for r in csv.DictReader(open('lifeonair/results/inputs/'+fn),delimiter='\t'):
  g=next((x[3:] for x in r['classification'].split(';') if x.startswith('g__')),'');g=re.sub(r'_[A-Z]+$','',g)
  if g in genera:matched.append((r['user_genome'],g,r['user_genome'] in old))
with (out/'trip4_genus_mag_matches.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['rep','genus','already_in_exclusion_union']);w.writerows(matched)
union=old|{r[0] for r in matched}
(out/'mag_exclusion_union.txt').write_text(''.join(x+'\n' for x in sorted(union)))
info={'source':str(src.resolve()),'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'source_screen':'September 9 replay, Trip4 only','candidate_ids':len(ids),'sequences_found':len(found),'missing_sequence_ids':sorted(ids-found),'named_genera':sorted(genera),'matched_MAGs':len(matched),'new_MAGs':len(union-old),'union_n':len(union)}
(out/'trip4_provenance.json').write_text(json.dumps(info,indent=2)+'\n');print(info)
