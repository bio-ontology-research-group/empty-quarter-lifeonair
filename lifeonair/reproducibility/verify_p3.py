#!/usr/bin/env python3
"""Check frozen P3 evidence and independently recompute primary paired contrasts.

TABLE checks compare reported, rounded claims with saved model outputs; they do
not refit those models. COMPUTE checks derive cohorts and paired contrasts from
sample metadata and frozen GE, independently of the R analysis. No network.
"""
import csv
import hashlib
import json
import math
from pathlib import Path
from statistics import median
from collections import defaultdict

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
RES = ROOT / 'lifeonair/results'
MAN = (ROOT / 'lifeonair/manuscript/main.tex').read_text()

def read(path):
    with Path(path).open() as f:
        return list(csv.DictReader(f, delimiter='\t'))

failures = []
nchecks = 0

def check(ok, label):
    global nchecks
    nchecks += 1
    if not ok:
        failures.append(label)
        print('FAIL:', label)

for c in read(HERE / 'claims_p3.tsv'):
    rows = read(RES / c['source'])
    select = json.loads(c['selector'])
    chosen = [r for r in rows if all(r[k] == v for k, v in select.items())]
    check(len(chosen) == 1, c['id'] + ' unique source row')
    if len(chosen) != 1:
        continue
    value = float(chosen[0][c['field']]) * float(c['scale'])
    check(abs(value - float(c['reported'])) <= float(c['tolerance']), c['id'] + ' value')
    if c['text_anchor'] != '-':
        check(c['text_anchor'] in MAN, c['id'] + ' manuscript anchor')

metadata = {r['sample_id']: r for r in read(RES / 'inputs/samples.tsv')
            if r['site'].isdigit() and 1 <= int(r['site']) <= 60 and r['trip'] and r['compartment']}
ge = [r for r in read(RES / 'reads/guild_genome_equivalents_v2.tsv')
      if r['id_floor'] == '60' and int(r['scg_genes_detected']) >= 10 and r['sample'] in metadata]
check(len(ge) == 546 and len({r['sample'] for r in ge}) == 546, 'COMPUTE primary sample closure')
check(len({metadata[r['sample']]['site'] for r in ge}) == 59, 'COMPUTE primary site closure')
markers = ['hydB_1l','hydB_1h','hydB_2a','coxL_formI','hydB_other_nife','coxL_other']
groups = defaultdict(list)
for r in ge:
    m = metadata[r['sample']]
    groups[(m['site'], m['trip'], m['compartment'])].append(r)
expected = read(RES / 'h123_v2/contrasts_by_position.tsv')
for name, numerator, denominator in [('root/bulk','PlantRoot','Deep'),('surface/bulk','Surface','Deep'),('root/surface','PlantRoot','Surface')]:
    values = defaultdict(list)
    sites = set()
    for site, trip, pos in groups:
        if pos != numerator or (site, trip, denominator) not in groups:
            continue
        a = groups[(site,trip,numerator)]; b = groups[(site,trip,denominator)]
        ratios = {v: math.log2((median(float(r[v]) for r in a)+0.01)/(median(float(r[v]) for r in b)+0.01)) for v in markers}
        for v,x in ratios.items():values[v].append(x)
        for v in markers[:4]:values['excess_'+v].append(ratios[v]-(ratios[markers[4]]+ratios[markers[5]])/2)
        sites.add(site)
    for r in expected:
        if r['contrast'] != name:continue
        vals = values[r['quantity']]
        check(len(vals) == int(r['n_visits']) and len(sites) == int(r['n_sites']), 'COMPUTE '+name+' cohort '+r['quantity'])
        check(abs(median(vals)-float(r['median_log2'])) < 1e-10, 'COMPUTE '+name+' '+r['quantity'])
        check(abs(sum(v<0 for v in vals)/len(vals)-float(r['frac_below0'])) < 1e-10, 'COMPUTE '+name+' direction '+r['quantity'])

# Independently regenerate the additional round-3 medians and benchmark audit.
import importlib.util
spec = importlib.util.spec_from_file_location('r3_evidence', ROOT/'lifeonair/code/r3_evidence_checks.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
r3_rows, audit = module.compute()
saved = read(RES/'review_r3/read_sensitivity_contrasts.tsv')
check(len(saved) == len(r3_rows), 'COMPUTE R3 contrast row closure')
for a,b in zip(r3_rows,saved):
    check(all(str(a[k]) == b[k] for k in a if k != 'median_log2') and
          abs(a['median_log2']-float(b['median_log2'])) < 1e-10, 'COMPUTE R3 contrast '+str((a['variant'],a['contrast'],a['quantity'])))
check(audit == json.loads((RES/'review_r3/benchmark_residual_homology.json').read_text()), 'COMPUTE benchmark residual homology')
candidate_ids = {r['feature_id'] for r in read(RES/'review_r3/trip4_candidate_calls.tsv')}
sequence_ids = {s[1:].split()[0] for s in (RES/'review_r3/trip4_candidate_sequences.fasta').read_text().splitlines() if s.startswith('>')}
check(len(candidate_ids) == 7 and sequence_ids == candidate_ids, 'COMPUTE Trip4 candidate sequence closure')
old_flags = set((RES/'reads/excl_v2_noflagged.txt').read_text().split())
new_matches = {r['rep'] for r in read(RES/'review_r3/trip4_genus_mag_matches.tsv')}
union = set((RES/'review_r3/mag_exclusion_union.txt').read_text().split())
check(union == old_flags | new_matches and len(union) == 387 and len(union-old_flags) == 4, 'COMPUTE expanded exclusion union')

for path, keys in [('inputs/env_site.tsv',['site']),('inputs/env_site_trip.tsv',['site','trip'])]:
    rows=read(RES/path)
    check(len(rows)==len({tuple(r[k] for k in keys) for r in rows}), path+' unique keys')
labels=read(RES/'inputs/atmo_classification_v2.tsv')
for group,assignment,n in [('1l','[NiFe]_Group_1l',597),('1h','[NiFe]_Group_1h',410),('2a','[NiFe]_Group_2a',134)]:
    check(len({r['rep_id'] for r in labels if r['assignment']==assignment})==n,'COMPUTE '+group+' carriers')
coxl=[r for r in labels if r['marker']=='coxL' and r['assignment'].startswith('form_I|')]
check(len(coxl)==284 and len({r['rep_id'] for r in coxl})==245,'COMPUTE retained form-I proteins/carriers')
check(not any(r['assignment'].startswith('form_I|AYRGAGR') for r in labels),'form-I motif labels')
for line in (HERE/'inputs.sha256').read_text().splitlines():
    digest,name=line.split('  ',1)
    path=ROOT/name
    check(path.is_file() and hashlib.sha256(path.read_bytes()).hexdigest()==digest, 'INPUT '+name)
# LaTeX source tables: catch stale column counts and undefined citation keys.
for name,ncol in [('benchmark',7),('mock',8),('sensitivity',6),('sensitivity_secondary',6),('lineages',5)]:
    rows=(ROOT/f'lifeonair/manuscript/tab_{name}.tex').read_text().splitlines()
    check(all(len(r.split('&'))==ncol for r in rows if r.strip()), 'table columns '+name)
import re
bib=(ROOT/'lifeonair/manuscript/references.bib').read_text()
keys=set(re.findall(r'@\w+\s*\{\s*([^,]+),',bib))
cited={k.strip() for match in re.findall(r'\\cite\w*\{([^}]+)\}',MAN) for k in match.split(',')}
check(cited <= keys, 'citation keys: '+str(sorted(cited-keys)))
print(f'P3 verification: {nchecks} checks, {len(failures)} failures; model-output checks are TABLE, not refits.')
raise SystemExit(bool(failures))
