# Observation-template provenance

These are cropped physical-model spectra, not observations of a target.

- `phoenix_5700_4p5.csv`: PHOENIX-ACES surface flux at Teff5700K, logg4.5,
  solar metallicity, vacuum2.26--2.37um. Preserve attribution to
  [Husser et al. 2013](https://doi.org/10.1051/0004-6361/201219058) and the
  [Gottingen archive](https://phoenix.astro.physik.uni-goettingen.de/data/v2.0/HiResFITS/).
  This is a flux spectrum, NOT a mu-dependent specific-intensity grid.
- `skycalc_paranal.csv`: ESO Cerro Paranal Advanced Sky Model, airmass1.3,
  PWV2.5mm, vacuum2260--2370nm, sampling lambda/dlambda1e6, no LSF. Columns
  are wavelength(nm), transmission, sky radiance(ph/s/m2/um/arcsec2).
  Full request parameters and raw response hash are in `skycalc_provenance.json`.
  Cite [Noll et al. 2012](https://doi.org/10.1051/0004-6361/201219040),
  [Jones et al. 2013](https://doi.org/10.1051/0004-6361/201322433) and
  [ESO SkyCalc](https://www.eso.org/observing/etc/doc/skycalc/helpskycalc.html).

The fetch script records source URLs and SHA256 hashes of the downloaded FITS
responses. No observational residual amplitudes come from these downloads:
those remain explicitly specified sensitivity hypotheses in the pilot.
These third-party model data retain their source attribution; the project's
software licence does not claim ownership of these libraries.

To regenerate the crops (network access, ESO service availability and optional
FITS/download dependencies are required):

```bash
python -m pip install '.[observing-data]'
python scripts/fetch_observing_templates.py
python scripts/fetch_observing_templates.py --stellar-only
```

The normal pilot reads the tracked CSV crops and does not need a live service.
Regenerating templates changes their content fingerprint and requires a new
Monte Carlo run tag; it must not silently reuse the existing calibration.
