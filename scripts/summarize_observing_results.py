"""Validate the complete ensemble and export a reproducible referee report."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import scipy
from observing_contamination_pilot import FAMILIES
from xds.orbit_search import binomial_interval


def read(path):
    with path.open() as handle:return list(csv.DictReader(handle))


def main(tag):
    root=Path(__file__).resolve().parents[1]
    config=json.loads((root/f'docs/co_observing_{tag}_config.json').read_text())
    expected=config.pop('fingerprint')
    files=['scripts/observing_contamination_pilot.py','src/xds/observing_contamination.py',
        'src/xds/orbit_search.py','src/xds/nonlinear_atmosphere.py','src/xds/moving_shadow.py',
        'data/observing/skycalc_paranal.csv','data/observing/phoenix_5700_4p5.csv']
    fingerprint=hashlib.sha256(json.dumps(config,sort_keys=True).encode())
    for name in files:fingerprint.update((root/name).read_bytes())
    if fingerprint.hexdigest()!=expected:raise RuntimeError('source differs from calibrated procedure')
    decisions=read(root/f'docs/co_observing_{tag}_decisions.csv')
    calibration=read(root/f'docs/co_observing_{tag}_calibration.csv')
    rates=read(root/f'docs/co_observing_{tag}_rates.csv')
    all_rows=calibration+decisions
    keys=[tuple(r[k] for k in ['level','method','stream','family','draw']) for r in all_rows]
    assert len(set(keys))==len(keys)
    for level in config['levels']:
        for method in ['unprotected','protected']:
            for stream,families,n in [('calibration',FAMILIES,config['n_calibration']),
                ('validation',FAMILIES,config['n_validation']),
                ('injection',['clean_fixed','combined_fixed','combined_random'],config['n_injections'])]:
                for family in families:
                    group=[r for r in all_rows if r['level']==level and r['method']==method and r['stream']==stream and r['family']==family]
                    assert len(group)==n,(level,method,stream,family,len(group),n)
    library_files=['outputs/co_nonlinear_emulator.npz','outputs/co_blind_shadow.npz','docs/co_blind_library.json','scripts/blind_co_recovery.py']
    manifest=dict(procedure_fingerprint=expected,scenario_data_sets=len(all_rows)//2,
        underlying_paired_realizations=len(all_rows)//4,searches=len(all_rows),
        versions=dict(numpy=np.__version__,scipy=scipy.__version__),
        sha256={name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in files+library_files})
    checks=dict(seed_refinements=3*len(all_rows),
        seed_refinements_without_reported_convergence=sum(int(r['refinement_failures']) for r in all_rows),
        maximum_null_start_spread=max(float(r['null_start_spread']) for r in all_rows),
        null_start_spreads_above_0p001=sum(float(r['null_start_spread'])>1e-3 for r in all_rows))
    (root/f'docs/co_observing_{tag}_checks.json').write_text(json.dumps(checks,indent=2)+'\n')
    (root/f'docs/co_observing_{tag}_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    audit=json.loads((root/'docs/co_observing_metric_audit.json').read_text())
    budget=read(root/'docs/co_observing_budgets.csv')
    information=read(root/f'docs/co_observing_{tag}_information.csv')
    start_audit=json.loads((root/f'docs/co_observing_{tag}_null_start_audit.json').read_text())
    start_checks=start_audit['checks']
    assert len(start_checks)==5 and all(r['all_audit_starts_success'] for r in start_checks)
    score_change=max(abs(r['audited_score']-r['original_score'])/max(abs(r['original_score']),1e-30) for r in start_checks)
    text=['# Stellar/telluric contamination and instrumental count-budget pilot','',
        '## Scope and main interpretation','',
        'This is a reproducible, conditional observation-space sensitivity experiment, not an end-to-end instrument forecast or a claim of detection. The physical benchmark is unchanged: a 0.7 Mj/1.1 Rj planet, 640 K reference atmosphere, 12 km/s spin, and a 0.09 Rp moon on a 21.2383 h circular aligned orbit; phase and period are searched, not supplied to the recovery. The 30-day planet orbits a roughly solar host, represented by a 5700 K/logg 4.5/solar PHOENIX flux spectrum. K=6 is an assumed apparent brightness, not an identified target.','',
        'The previous idealized 11.39465ppm control is not supported by the new count budget. Even with the orbit known, the narrow-band physical moon has kinematic information far below unity at the count-based noise. False-alarm calibration is necessary, but cannot replace missing photons. No population search was repeated.','',
        '## Physical inputs and observing budget','',
        'Stellar line positions and shapes come from PHOENIX-ACES; terrestrial transmission and sky radiance come from ESO SkyCalc, Paranal, airmass 1.3 and PWV 2.5 mm. The original model requests, raw-file hashes and cropped data are in `data/observing/`. Both templates are convolved and binned to R=100000. The fixed mask retains 3455/3468 bins (T>=0.3 and normalized stellar flux>=0.2).','',
        'The source counts per resolution element are approximated by','',
        r'$$N_{\star,k}=\frac{F_{\nu,k}}{hR}\,A\eta\,t_{\rm live}\,T_k,$$','',
        'with a monochromatic 666.7 Jy Ks zero point at 2.159 um, a 5700 K continuum slope and the stellar line spectrum. The fractional noise is','',
        r'$$\sigma_k^2=\frac{N_{\star,k}+2(N_{{\rm sky},k}+N_{{\rm th},k}+N_{{\rm dark},k}+n_{\rm pix}n_{\rm read}\sigma_{\rm read}^2)}{N_{\star,k}^2}.$$','',
        'The factor 2 represents paired background subtraction. Assumptions are efficiency 8%, extraction 0.08 arcsec2, 8 pixels per resolution bin, 12 e read noise, 0.003 e/s/pixel dark, and a 280 K/emissivity 0.1 thermal surrogate. No freely multiplied number of lines or transits is used.','',
        '| Scenario | Live time / 5min | Single-exposure noise q16--median--q84 [ppm] | Peak e/read/pixel | Linearity screen |',
        '| --- | ---: | ---: | ---: | --- |']
    for r in budget:
        text.append(f"| {r['label']} | {float(r['live_seconds']):.3f}s | {float(r['q16_sigma_ppm']):.1f}--{float(r['median_sigma_ppm']):.1f}--{float(r['q84_sigma_ppm']):.1f} | {float(r['peak_source_e_per_read_pixel']):.0f} | {'Pass' if r['valid_linearity']=='True' else 'FAIL: photon upper bound only'} |")
    text.extend(['',
        'The 8 m case uses 8x30 s integrations; 39 m uses 84x1.427 s. A peak/mean pixel factor 2 and 12900 e linearity screen are adopted. The 39 m/K3 row FAILS this screen and is not an executable observation or recovery forecast. Attenuation, lower efficiency, shorter integrations, spatial spreading and real overheads would need an instrument-specific design. The 39 m area is 978 m2, but its 8% efficiency and detector assumptions are NOT certified ELT spectrograph performance. A Ks passband-integrated flux calibration, real order coverage and an ETC remain necessary.','',
        'A one-hour out-of-transit reference comprises 12 equivalent exposures. In the spectrally whitened coordinates, its exact time covariance is C=I+11^T/12. Thus the diagonal variance is 13/12 and off-diagonal correlation 1/13. The same C is simulated and inverted in every nonlinear atmospheric fit and every score. It is not hidden inside an independent "floor".','',
        '## Contamination and the reduction comparison','',
        'Every generator includes the previous smooth bounded atmospheric histories plus random wind pulses. Additional stellar residuals retain 3--15% of a disk-minus-local-spectrum RM surrogate (vsini=2 km/s), with a phenomenological 20%(1-mu) line-depth change. PHOENIX surface flux is only a local-profile surrogate: this is NOT a mu-resolved CLV calculation, full spot model or physical stellar surface integration.','',
        'Telluric residuals use optical-depth errors of 0.1--0.5% and wavelength errors of 3--15 m/s, with 0.7 h correlated Gaussian-process histories. Planetary orbital acceleration and an assumed 0.1 km/s/h terrestrial drift shift the templates into the planet frame. Systemic velocity is fixed at zero; these results do not marginalize over observing epochs, barycentric velocity or stellar velocity. Residual amplitudes are explicit sensitivity assumptions, NOT measured correction performance.','',
        'The protected reduction removes six fixed spectral directions with unconstrained coefficients per exposure: tau, its derivative, tau squared, stellar contrast, and its first and second derivatives. Shapes are externally fixed, but actual shifted nonlinear residuals need not lie in their span. No injected lunar orbit, atmospheric history or contaminant coefficients are supplied to the fit. Protected and unprotected noise realizations are paired.','',
        '## Calibration design and numerical checks','',
        f"There are {manifest['scenario_data_sets']} noise-level data sets ({manifest['underlying_paired_realizations']} underlying paired realizations) and {manifest['searches']} paired-reduction searches. Each noise level uses 64 calibration nulls and 32 independent held-out nulls for EACH of four generators, plus 24 injections in EACH of three groups. Noise levels are paired; their ensembles are not independent repeats. Each draw refits the nonlinear atmosphere with three temporal controls. The full 73x1024 period/phase grid and the previous top-three local score refinements are retained. Period and separation remain Kepler-coupled.",'',
        'The two noise levels are the 39 m/K6 budget and a laboratory rescaling with median 11.39465 ppm but the SAME wavelength-dependent noise shape. The control also includes the new mask and finite shared reference, so its rates cannot be subtracted directly from the previous flat-noise/perfect-reference pilot as a pure contaminant effect.','',
        'Every level/reduction receives its own search-maximum calibration. The atmosphere-only p is compared to the conservative maximum p over all four null generators. p=(1+count[Tnull>=Tobs])/65; the diagnostic threshold is 0.05 and the minimum reportable p is 1/65=0.01538, not zero. The finite ensemble cannot support five-sigma claims or guarantee control for arbitrary contaminants. No score p-value is assigned to an uncalibrated nonlinear likelihood ratio.','',
        f"The direct-map audit gives maximum absolute point-score errors below {max(c['maximum_point_score_absolute_error'] for c in audit['checks']):.2g}, grid-score errors below {max(c['maximum_grid_score_absolute_error'] for c in audit['checks']):.2g}, derivative relative errors below {max(c['maximum_derivative_relative_error'] for c in audit['checks']):.2g}, and shared-reference whitening error {audit['shared_reference_whitening_error']:.2g}. See `co_observing_metric_audit.json`.",'',
        f"There are {checks['seed_refinements_without_reported_convergence']} local seed refinements without reported convergence out of {checks['seed_refinements']}. None is discarded: the grid fallback and best computed score are retained by the same algorithm in every calibration/validation draw. Atmospheric starts differ by more than 1e-3 in chi2 in {checks['null_start_spreads_above_0p001']} searches, with maximum spread {checks['maximum_null_start_spread']:.3f}. The better null fit is always retained; convergence of each start is not a proof of a global optimum.",'',
        f"A separate eight-start numerical audit checks the worst start spread in each level/reduction and one clean control (five cases). All 40 audit starts converge. The maximum further improvement of the selected null chi2 is {max(r['improvement'] for r in start_checks):.2g}; the maximum relative search-score change is {100*score_change:.5f}%. The two-start BEST solutions are therefore stable in these stress cases despite the large spread of a worse start. This finite check does not prove global convergence. Audited trials do NOT replace the frozen calibration or receive its p-values; see `co_observing_{tag}_null_start_audit.json`.",'',
        '## Lunar information before contaminant mismatches','',
        '| Count scenario | Unprotected known-orbit kinematic norm | Protected norm |',
        '| --- | ---: | ---: |'])
    for scenario in ['8m_K6','39m_K6']:
        vals=[next(r['known_orbit_kinematic_norm'] for r in audit['information'] if r['scenario']==scenario and r['protected']==p) for p in [False,True]]
        text.append(f'| {scenario} | {vals[0]:.5f} | {vals[1]:.5f} |')
    text.extend(['','These norms are noise-weighted physical model separations after a local atmospheric tangent removal, not discovery significances. They are measured at a reference atmosphere and favorable phase. The protected/unprotected norm ratio in 39 m/K6 is'+f" {audit['information'][3]['known_orbit_kinematic_norm']/audit['information'][2]['known_orbit_kinematic_norm']:.4f}"+'; this particular projection removes about 1% of the signal norm, but does not guarantee atmospheric or stellar uniqueness.','',
        '| Laboratory / count level | Reduction | Full signal norm | Known-orbit kinematic norm |',
        '| --- | --- | ---: | ---: |'])
    for r in information:text.append(f"| {r['level']} | {r['method']} | {float(r['full_signal_norm']):.4f} | {float(r['known_orbit_kinematic_norm']):.4f} |")
    for level in config['levels']:
        text.extend(['',f'## Held-out results: {level}','','False candidates without a moon; counts are out of 32.','',
            '| Generator | Unprotected, atmospheric calibration | Unprotected, all-null calibration | Protected, atmospheric calibration | Protected, all-null calibration |',
            '| --- | ---: | ---: | ---: | ---: |'])
        for family in FAMILIES:
            counts=[next(int(r['flags']) for r in rates if r['level']==level and r['stream']=='validation' and r['family']==family and r['method']==m and r['calibration']==c) for m,c in [('unprotected','atmosphere'),('unprotected','robust'),('protected','atmosphere'),('protected','robust')]]
            text.append('| '+family+' | '+' | '.join(f'{k}/32' for k in counts)+' |')
        text.extend(['','Injection candidate flags with all-null calibration; these are NOT correct-orbit recoveries. Counts are out of 24.','',
            '| Injection | Unprotected flags (95% CI) | Protected flags (95% CI) |',
            '| --- | ---: | ---: |'])
        for family in ['clean_fixed','combined_fixed','combined_random']:
            entries=[]
            for method in ['unprotected','protected']:
                r=next(r for r in rates if r['level']==level and r['stream']=='injection' and r['family']==family and r['method']==method and r['calibration']=='robust')
                entries.append(f"{r['flags']}/24 ({100*float(r['ci95_low']):.1f}--{100*float(r['ci95_high']):.1f}%)")
            text.append('| '+family+' | '+' | '.join(entries)+' |')
    upper=100*binomial_interval(0,32)[1]
    text.extend(['',f'With 32 held-out nulls per family, zero false candidates still has a 95% upper bound of {upper:.1f}%. These small samples do not establish a precisely known 5% false-alarm rate. All rate confidence intervals are in the numerical tables and validation figure.','',
        '## What can and cannot be claimed','',
        'Calibration can reject many spurious candidates by raising the threshold; that is NOT successful contaminant subtraction and can also erase sensitivity. The protected reduction and contaminated-null calibration solve different problems and must be evaluated together with injection rates. Candidate flags in the photon-budget scenario must be compared to its own no-moon false-candidate rates, not interpreted as recovery because a moon happened to be injected.','',
        'The preserved planet model is already instrument-binned. We add a small-signal normalized-residual observation surrogate; we do not multiply high-resolution stellar, terrestrial and planetary spectra before the LSF. The signal is averaged over the full cadence, approximating evenly distributed short detector integrations. Noise uses mean observing conditions and a Gaussian high-count approximation. Grey transit attenuation, detailed DIT overheads/ramp nonlinearity, order gaps, detector systematics, line-spread variability and real reference-condition changes are not simulated. The fixed planetary spectral subspace also discards external information that could constrain contaminants.','',
        'The unknown orbital inclination, eccentricity, moon radius, host/spin uncertainties, actual stellar specific intensities/spots and empirical pipeline residual distributions remain outside this conditional experiment. Shared-reference photon/detector covariance is exact under the adopted budget; unknown systematic errors in the reference are not comprehensively modelled. The moon\'s own photometric and stellar-RM occultation signals are absent: the detectability conclusion concerns its atmospheric EXOVEL signature, not every exomoon detection channel.','',
        'The conclusion is not that EXOVEL is impossible, but that this single-transit Earth-scale moon benchmark in this CO window is not supported by this count budget. Prioritize concrete bright targets, physically larger signals or broader molecular coverage before another population sweep; do not restore detectability by silently amplifying the moon or dividing the noise. Each reduction has a separate calibration; selecting the smaller p after inspecting both reductions would require calibrating that additional selection.','',
        '## Referee-response wording','',
        'We have added an observation-space sensitivity experiment using physical PHOENIX stellar spectra and ESO SkyCalc transmission and radiance, with an explicit photon/background/detector budget and finite shared-reference covariance. We repeated the phase- and period-blind search on independent contaminated no-moon controls and injected systems, both with and without fixed spectral nuisance removal. The false-alarm calibration is performed separately for each reduction and observing-noise scenario. We distinguish candidate-flag rates from correct orbital recovery and discovery significance. The count-based experiment shows that the adopted Earth-scale moon benchmark is below the information required for a single-transit detection in the modelled CO band. Accordingly, the noise-reduced experiments are retained only as methodological controls, not as observational forecasts. Residual correction amplitudes, the stellar local-profile surrogate and instrumental efficiencies are explicitly stated assumptions; a full stellar-intensity calculation and an instrument-specific end-to-end reduction remain beyond this pilot.','',
        '## Reproduction and figures','',
        'See `observing_contamination_protocol.md`, the README command sequence, the frozen config/manifest and tracked numerical CSV tables. No original population selections were modified. Normal regeneration uses the tracked cropped spectra, not a live ESO query.','',
        f'- `outputs/co_observing_{tag}_budget.png/pdf`: physical input spectra, count budget, representative residuals, shared-reference covariance.',
        f'- `outputs/co_observing_{tag}_null_tails.png/pdf`: global search tails, separately calibrated by reduction and noise.',
        f'- `outputs/co_observing_{tag}_validation.png/pdf`: held-out false alarms and injection flags with exact95% binomial intervals.',
        f'- `outputs/co_observing_{tag}_orbital_errors.png/pdf`: continuous orbital errors; filled points indicate candidate flags only.','',
        f'- `outputs/co_observing_{tag}_line_maps.png/pdf`: noiseless physical moon and representative contaminants near a model-selected strong CO line, with INDEPENDENT color scales. Nuisance removal is shown before planetary-subspace compression. This is not a detected residual map.','',
        'Sources: [ESO CRIRES detector performance](https://www.eso.org/sci/facilities/paranal/instruments/crires/inst.html), [ELT collecting area](https://elt.eso.org/about/facts/), [2MASS zero point](https://www.ipac.caltech.edu/2mass/releases/allsky/doc/sec6_4a.html), [ESO SkyCalc](https://www.eso.org/observing/etc/doc/skycalc/helpskycalc.html), [PHOENIX library](https://doi.org/10.1051/0004-6361/201219058).',''])
    (root/'docs/co_observing_contamination_report.md').write_text('\n'.join(text))
    print('Verified ensemble:',manifest['searches'],'searches;',manifest['scenario_data_sets'],'noise-level data sets')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--tag',default='pilot');args=parser.parse_args();main(args.tag)
