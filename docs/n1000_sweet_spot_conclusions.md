# XDS Sweet-Spot Search: n=220 vs n=1000

This note updates the conclusions from the first exploratory sweet-spot search.

## What Changed

The first run used:

```text
sample_candidates(n=220, seed=7)
```

and found 57 dynamically valid systems. The best candidate had:

```text
SNR = 3.78
```

No candidate exceeded SNR 5.

The new run used:

```text
sample_candidates(n=1000, seed=7)
```

and found 239 dynamically valid systems. The best candidate had:

```text
SNR = 12.41
```

Four candidates exceeded SNR 5. This shows that the favourable region of
parameter space is relatively narrow: a sparse random sample can miss it.

## Best n=1000 Candidate

```text
SNR = 12.41
peak XDS = 265.8 ppm
line depth = 5512.9 ppm
scale height = 822.1 km
Mp = 0.141 Mj
Rp = 1.642 Rj
Rstar = 0.601 Rsun
Pplanet = 116.8 d
Teq = 295 K
Prot = 23.3 h
Pbreakup = 16.6 h
v_eq = 8.78 km/s
am = 3.17 Rp
Pm = 93.8 h
Rm = 0.326 Rp
Tdur = 5.75 h
```

The high SNR is driven mainly by very low surface gravity and a large planetary
radius, which produce a very large atmospheric scale height and line depth.

## Moon-Separation Optimisation

For the best sampled candidate, scanning the moon separation gives:

```text
best am = 1.39 Rp
Pm = 27.3 h
Pm / Tdur = 4.75
SNR = 20.38
contact duration = 5.75 h
exposure limit = 15.0 min
peak XDS = 280.0 ppm
```

The scan confirms that the satellite separation is a central optimisation
parameter. In this synthetic case, the best geometry is a close moon just above
the hard inner limit set by the planet plus moon radii. The Roche limit is lower
because the planet is extremely low-density.

## Revised Scientific Conclusion

The earlier conclusion was too conservative because the n=220 run undersampled
the favourable region. With n=1000, the model finds detectable configurations
under the assumed observing setup:

```text
R = 100000
continuum SNR per resolution element = 300
N_lines = 200
N_transits = 1
systematic floor = 10 ppm
```

However, the most favourable systems are astrophysically special:

- very inflated;
- sub-Saturn to Saturn-mass;
- low-density;
- around relatively small stars;
- with close, large moons;
- rotating fast enough for Doppler leverage but not near breakup.

The manuscript should therefore avoid claiming broad detectability. A stronger
and more defensible claim is:

```text
XDS can be detectable in a narrow but physically interpretable region of
parameter space: low-gravity inflated giant planets with close massive moons,
long enough transits, and enough high-resolution spectral lines.
```

## Caution

The best candidate has an unusually low density and an enormous scale height.
Before using it as a benchmark, we should check:

- whether such a planet can retain an H/He atmosphere at the adopted age and
  irradiation;
- whether the assumed moon size and close orbit are dynamically plausible;
- whether the moon survives tides and migration;
- whether 200 useful spectral lines at the assumed SNR are realistic;
- whether systematics remain below 10 ppm in the velocity-resolved map.
