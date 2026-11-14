"""Paired oracle/blind decisions and cross-model orbital diagnostics; no new cuts."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from xds.moving_shadow import period_limits


def write(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--tag',default='pilot');args=parser.parse_args()
    root=Path(__file__).resolve().parents[1]
    with (root/f'docs/co_blind_{args.tag}_decisions.csv').open() as f:data=list(csv.DictReader(f))
    comparisons=[];agreement=[];coverage=[];lo,hi=period_limits()
    for stream in ['validation','injection']:
        labels=sorted({r['family'] for r in data if r['stream']==stream})
        for label in labels:
            groups={n:[r for r in data if r['stream']==stream and r['family']==label and int(r['n_controls'])==n] for n in [3,5]}
            for n,rows in groups.items():
                if label!='uniform_phase':
                    fixed=np.array([float(r['oracle_p_robust'])<=.05 for r in rows]);blind=np.array([r['flag_robust']=='True' for r in rows])
                    comparisons.append(dict(stream=stream,family=label,n_controls=n,n=len(rows),
                        fixed_orbit_flags=int(fixed.sum()),blind_flags=int(blind.sum()),
                        fixed_only_flags=int(np.sum(fixed&~blind)),blind_only_flags=int(np.sum(blind&~fixed))))
            a={r['draw']:r for r in groups[3]};b={r['draw']:r for r in groups[5]}
            for draw,row in a.items():
                peer=b[draw]
                agreement.append(dict(stream=stream,family=label,draw=int(draw),
                    consensus_flag=row['flag_consensus']=='True',
                    absolute_log_period_difference=float(abs(np.log(float(row['period'])/float(peer['period'])))),
                    circular_phase_difference=float(abs((float(row['phase'])-float(peer['phase'])+.5)%1-.5))))
            if stream=='injection':
                rows=groups[3]
                coverage.append(dict(family=label,n=len(rows),
                    nonzero_physical_signal=sum(float(r['physical_signal_norm'])>1e-12 for r in rows),
                    zero_physical_signal=sum(float(r['physical_signal_norm'])<=1e-12 for r in rows),
                    min_signal_norm=min(float(r['physical_signal_norm']) for r in rows),
                    median_signal_norm=float(np.median([float(r['physical_signal_norm']) for r in rows])),
                    max_signal_norm=max(float(r['physical_signal_norm']) for r in rows)))
    write(root/f'docs/co_blind_{args.tag}_oracle_comparison.csv',comparisons)
    write(root/f'docs/co_blind_{args.tag}_model_agreement.csv',agreement)
    write(root/f'docs/co_blind_{args.tag}_phase_coverage.csv',coverage)
    fits=json.loads((root/f'docs/co_blind_{args.tag}_candidate_fits.json').read_text())
    runs=[r for case in fits for candidate in case['nonlinear'] for r in candidate['runs']]
    diagnostics=dict(n_selected_nonlinear_starts=len(runs),failed_selected_nonlinear_starts=sum(not r['success'] for r in runs),
        note='oracle comparison only for the predeclared fixed showcase orbit; uniform-phase group excluded',
        numerical_zero_visibility_threshold=1e-12,visibility_is_diagnostic_not_a_selection_cut=True)
    (root/f'docs/co_blind_{args.tag}_extra_checks.json').write_text(json.dumps(diagnostics,indent=2)+'\n')
    print(json.dumps(dict(coverage=coverage,nonlinear_checks=diagnostics),indent=2))


if __name__=='__main__':main()
