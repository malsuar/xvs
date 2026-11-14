# ExoMol CO Li2015 input data

Obtain inputs with `python scripts/fetch_co_template_data.py` from the repository
root. SHA256 hashes and source URLs are pinned in manifest.json. Raw downloads
are intentionally excluded from Git. An upstream change fails verification;
updating the scientific inputs requires a new documented manifest.

Source: https://exomol.com/data/molecules/CO/12C-16O/Li2015/
Reference: Li et al. (2015), ApJS 216,15,
https://doi.org/10.1088/0067-0049/216/1/15.
H2 broadening source: https://exomol.com/data/data-types/broadening_coefficients/CO/
Line-intensity expression: Yurchenko et al. (2018), A&A 614,A131,
https://doi.org/10.1051/0004-6361/201732531.

ExoMol data license: CC BY-SA 4.0,
https://exomol.com/data/licence/ and https://creativecommons.org/licenses/by-sa/4.0/.
This is separate from the software's MIT license. Attribution and ShareAlike
apply to the source data and derived molecular spectral data products.

Derived profiles in outputs/co_molecular_injection.npz are modified products:
LTE populations, line filtering, pressure/Doppler broadening, hydrostatic slant
columns, opacity exponentiation and instrument/geometry integration have been
applied by the XVS experiment. They are not original ExoMol spectra and are
provided under CC BY-SA 4.0 with the above attribution. No ExoMol endorsement
is implied. The experiment is not a chemical-equilibrium calculation.
