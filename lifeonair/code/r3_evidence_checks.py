#!/usr/bin/env python3
"""Recompute R3 paired medians and audit residual benchmark homology; no model fits."""
import csv
import json
import math
from collections import defaultdict
from pathlib import Path
from statistics import median

RES = Path(__file__).resolve().parents[1] / 'results'
OUT = RES / 'review_r3'

def read(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter='\t'))

def compute():
    meta = {r['sample_id']: r for r in read(RES / 'inputs/samples.tsv')
            if r['site'].isdigit() and 1 <= int(r['site']) <= 60 and r['trip'] and r['compartment']}
    markers = ['hydB_1l', 'hydB_1h', 'hydB_2a', 'coxL_formI', 'hydB_other_nife', 'coxL_other']
    result = []
    for variant, suffix in [('baseline', ''), ('excluded', '_r3noflagged')]:
        groups = defaultdict(list)
        for r in read(RES / f'reads/guild_genome_equivalents_v2{suffix}.tsv'):
            if r['id_floor'] != '60' or int(r['scg_genes_detected']) < 10 or r['sample'] not in meta:
                continue
            m = meta[r['sample']]
            groups[(m['site'], m['trip'], m['compartment'])].append(r)
        for contrast, numerator, denominator in [('root/bulk', 'PlantRoot', 'Deep'), ('surface/bulk', 'Surface', 'Deep'), ('root/surface', 'PlantRoot', 'Surface')]:
            values = defaultdict(list)
            sites = set()
            for site, trip, pos in groups:
                if pos != numerator or (site, trip, denominator) not in groups:
                    continue
                a, b = groups[(site, trip, numerator)], groups[(site, trip, denominator)]
                ratios = {v: math.log2((median(float(r[v]) for r in a)+.01)/(median(float(r[v]) for r in b)+.01)) for v in markers}
                for v, x in ratios.items():
                    values[v].append(x)
                for v in markers[:4]:
                    values['excess_'+v].append(ratios[v]-(ratios[markers[4]]+ratios[markers[5]])/2)
                sites.add(site)
            for quantity, vals in values.items():
                result.append(dict(variant=variant, contrast=contrast, quantity=quantity,
                                   n_visits=len(vals), n_sites=len(sites), median_log2=median(vals)))
    original = {(r['contrast'], r['quantity']): r for r in read(RES / 'h123_v2/contrasts_by_position.tsv')}
    for r in result:
        if r['variant'] == 'baseline':
            old = original[(r['contrast'], r['quantity'])]
            assert abs(r['median_log2']-float(old['median_log2'])) < 1e-10
            assert r['n_visits'] == int(old['n_visits'])
    rows = read(RES / 'benchmark_r2/disjoint_verification.tsv')
    audit = {}
    for scheme, threshold in [('comp50', 50), ('comp70', 70)]:
        z = [r for r in rows if r['scheme'] == scheme]
        hits = [r for r in z if r['max_identity_cov50'] not in ('NA', '', 'nan')]
        above = [r for r in hits if float(r['max_identity_cov50']) > threshold]
        audit[scheme] = dict(total_queries=len(z), queries_with_cov50_hits=len(hits),
                             above_threshold_queries=len(above), max_cov50_identity=max(float(r['max_identity_cov50']) for r in hits),
                             affected_groups=sorted({r['group'] for r in above}))
    return result, audit

if __name__ == '__main__':
    result, audit = compute()
    with (OUT / 'read_sensitivity_contrasts.tsv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=result[0], delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(result)
    (OUT / 'benchmark_residual_homology.json').write_text(json.dumps(audit, indent=2)+'\n')
    print('R3 primary contrasts reproduced; exclusion contrasts and benchmark audit saved.')
