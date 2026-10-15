#!/usr/bin/env python3
"""Pin secondary model inputs after the read-only diagnostic run."""
import csv,hashlib,sys
from pathlib import Path
R=Path(sys.argv[1]);W=R/'06_atmospheric/metabolic_pilot_clean'
with (W/'results/gapseq_identity.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
paths={R/'05_metabolism/hq_mag_list.tsv',R/'06_atmospheric/lifeonair_reads/db/markers.faa'}
for r in rows:
 if r['status']!='identical':continue
 rep=r['rep'];tag=rep.replace('__single_','__');d=R/'05_metabolism/gapseq_out'/tag
 paths.update([Path(r['source']),R/'05_refine/_drep/out/dereplicated_genomes'/(rep+'.fa'),d/(tag+'.xml'),d/(tag+'-draft.xml')])
with (W/'results/gapseq_input_hashes.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['path','sha256'])
 for p in sorted(paths):w.writerow([str(p),hashlib.sha256(p.read_bytes()).hexdigest()])
