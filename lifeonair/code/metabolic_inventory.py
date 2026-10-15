#!/usr/bin/env python3
"""Inventory exact-name reconstructions and select a phenotype-blind matched pilot."""
import argparse, csv, json
from pathlib import Path
from collections import Counter, defaultdict

def read(p):
    with p.open() as f: return list(csv.DictReader(f, delimiter='\t'))

def write(p, rows):
    with p.open('w') as f:
        w=csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(rows)

def main():
    p=argparse.ArgumentParser(); p.add_argument('--inputs',type=Path,required=True)
    p.add_argument('--models',type=Path,required=True); p.add_argument('--out',type=Path,required=True)
    a=p.parse_args(); a.out.mkdir(parents=True,exist_ok=True)
    flags=set((a.inputs/'mag_exclusion_union.txt').read_text().split())
    cls=defaultdict(set)
    for r in read(a.inputs/'atmo_classification_v2.tsv'):
        if r['marker']=='hydB': cls[r['rep_id']].add(r['assignment'].split('_')[-1])
    co=set((a.inputs/'guild_coxL_reps_v2.txt').read_text().split())
    tax={r['user_genome']:r['classification'] for fn in ['r232_bac120.summary.tsv','r232_ar53.summary.tsv'] for r in read(a.inputs/fn)}
    rows=[]
    for r in read(a.inputs/'rep_genome_size.tsv'):
        rep=r['genome_id']; gs=cls[rep]; d=a.models/rep
        g='1l' if '1l' in gs and not gs&{'1h','2a'} else ('1h' if '1h' in gs and not gs&{'1l','2a'} else ('neither' if not gs&{'1l','1h','2a'} and rep not in co else 'other'))
        ranks=dict(x.split('__',1) for x in tax[rep].split(';'))
        model=(d/(rep+'-filled.xml')).is_file();draft=(d/(rep+'-draft.xml')).is_file()
        quality=float(r['completeness'])>=90 and float(r['contamination'])<=5
        eligible=model and draft and quality and not rep.endswith('_sub') and rep not in flags and ranks.get('o','') and g!='other'
        rows.append(dict(rep=rep,group=g,completeness=r['completeness'],contamination=r['contamination'],length_bp=r['length_bp'],
                         phylum=ranks.get('p',''),order=ranks.get('o',''),family=ranks.get('f',''),genus=ranks.get('g',''),
                         coxI=int(rep in co),filled=int(model),draft=int(draft),subset=int(rep.endswith('_sub')),flagged=int(rep in flags),eligible=int(bool(eligible))))
    write(a.out/'model_inventory.tsv',rows)
    # At most 24 non-reused triplets, <=3 per order; closest family/genus first,
    # then genome-quality differences. No field abundance or model phenotype used.
    by=defaultdict(lambda:defaultdict(list))
    for r in rows:
        if r['eligible']:by[r['order']][r['group']].append(r)
    panel=[]; used=set(); counts=Counter()
    for order,z in sorted(by.items()):
        if not all(z[g] for g in ['1l','1h','neither']):continue
        cases=sorted(z['1l'],key=lambda r:(-float(r['completeness']),float(r['contamination']),r['rep']))
        for case in cases:
            if counts[order]>=3 or len(panel)>=72:break
            choices=[]
            for g in ['1h','neither']:
                avail=[r for r in z[g] if r['rep'] not in used]
                if not avail:break
                def distance(r):return (int(r['genus']!=case['genus']),int(r['family']!=case['family']),abs(float(r['completeness'])-float(case['completeness'])),float(r['contamination']),r['rep'])
                choices.append(min(avail,key=distance))
            if len(choices)!=2:continue
            block=len(panel)//3+1
            for r in [case]+choices:
                panel.append(dict(block=block,**r));used.add(r['rep'])
            counts[order]+=1
    if not panel:raise RuntimeError('No complete matched triplets')
    write(a.out/'pilot_panel.tsv',panel)
    summary=dict(n=len(rows),filled=sum(r['filled'] for r in rows),draft=sum(r['draft'] for r in rows),eligible=sum(r['eligible'] for r in rows),
                 eligible_groups=dict(Counter(r['group'] for r in rows if r['eligible'])),matched_models=len(panel),matched_orders=dict(counts))
    (a.out/'inventory_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
