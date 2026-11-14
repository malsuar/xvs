"""Fetch the ESO SkyCalc physical sky; retain request and raw-file provenance."""
import argparse
import hashlib
import json
import sys
from importlib.metadata import version
from pathlib import Path


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--client-wheel')
    parser.add_argument('--stellar-only',action='store_true')
    args=parser.parse_args()
    if args.client_wheel:sys.path.insert(0,args.client_wheel)
    from skycalc_cli.skycalc import SkyModel
    root=Path(__file__).resolve().parents[1]
    folder=root/'data/observing';folder.mkdir(exist_ok=True)
    if args.stellar_only:
        import io
        import requests
        import numpy as np
        from astropy.io import fits
        base='https://phoenix.astro.physik.uni-goettingen.de/data/v2.0/HiResFITS/'
        names=['WAVE_PHOENIX-ACES-AGSS-COND-2011.fits',
            'PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte05700-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits']
        arrays=[];sources=[]
        for name in names:
            response=requests.get(base+name,timeout=120);response.raise_for_status()
            arrays.append(fits.getdata(io.BytesIO(response.content)))
            sources.append(dict(url=base+name,sha256=hashlib.sha256(response.content).hexdigest()))
        wave=arrays[0]/1e4;keep=(wave>=2.26)&(wave<=2.37)
        np.savetxt(folder/'phoenix_5700_4p5.csv',np.c_[wave[keep],arrays[1][keep]],
            delimiter=',',header='vacuum_wavelength_um,flux_erg_s_cm2_cm',comments='')
        (folder/'phoenix_provenance.json').write_text(json.dumps(dict(sources=sources,
            model='PHOENIX-ACES 5700 K, logg 4.5, solar metallicity; surface flux, not specific intensity'),indent=2)+'\n')
        print('PHOENIX cropped:',int(keep.sum()),'samples')
        return
    model=SkyModel()
    model.params.update(airmass=1.3,pwv=2.5,wmin=2260.,wmax=2370.,
        wgrid_mode='fixed_spectral_resolution',wres=1000000,lsf_type='none',
        incl_moon='N',incl_therm='N',observatory='paranal')
    model.call()
    raw=model.getdata()
    (folder/'skycalc_paranal.fits').write_bytes(raw)
    metadata=dict(source='https://www.eso.org/observing/etc/doc/skycalc/helpskycalc.html',
        service=model.base_url,client='skycalc_cli '+version('skycalc_cli'),parameters=model.params,
        sha256=hashlib.sha256(raw).hexdigest())
    (folder/'skycalc_provenance.json').write_text(json.dumps(metadata,indent=2)+'\n')
    import io
    import numpy as np
    from astropy.io import fits
    table=fits.getdata(io.BytesIO(raw),1)
    np.savetxt(folder/'skycalc_paranal.csv',np.c_[table['lam'],table['trans'],table['flux']],
        delimiter=',',header='vacuum_wavelength_nm,transmission,photons_s_m2_um_arcsec2',comments='')
    print('SkyCalc downloaded:',len(raw),'bytes',flush=True)


if __name__=='__main__':main()
