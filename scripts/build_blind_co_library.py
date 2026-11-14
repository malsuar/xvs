"""Reusable physical moving-shadow library; independent of search/noise trials."""
import hashlib
import inspect
import itertools
import json
from pathlib import Path
import numpy as np
from nonlinear_co_recovery import AXES,direct
from molecular_injection_recovery import build_molecular_geometry
from xds.moving_shadow import moment_geometry,stellar_moments
from xds.nonlinear_atmosphere import physical_parameters


def build(root):
    old=np.load(root/'outputs/co_nonlinear_emulator.npz');basis=old['basis']
    base=build_molecular_geometry(nr=52,n_scale_heights=16.,split_limbs=True)
    fingerprint=hashlib.sha256(basis.tobytes())
    for filename in ['src/xds/moving_shadow.py','src/xds/molecular.py','scripts/molecular_injection_recovery.py',
                     'data/exomol/CO_Li2015/manifest.json']:
        fingerprint.update((root/filename).read_bytes())
    fingerprint.update(inspect.getsource(physical_parameters).encode())
    fingerprint.update(json.dumps([a.tolist() for a in AXES]).encode())
    key=fingerprint.hexdigest()[:16];cache=root/'outputs/co_blind_cache'/key;cache.mkdir(parents=True,exist_ok=True)
    if (cache/'geometry.npz').exists():
        saved=np.load(cache/'geometry.npz');centers=saved['centers'];left=saved['left'];error=float(saved['error'])
        geometry=dict(base,full=np.zeros_like(saved['missing']),missing=saved['missing'],time=np.arange(saved['missing'].shape[1]))
    else:
        centers,left,geometry,error=moment_geometry(base)
        np.savez_compressed(cache/'geometry.npz',centers=centers,left=left,missing=geometry['missing'],error=error)
    count=geometry['missing'].shape[1]
    library=np.zeros((*map(len,AXES),count,basis.shape[1]))
    indices=list(itertools.product(*(range(len(a)) for a in AXES)))
    for n,index in enumerate(indices):
        path=cache/('node_'+'_'.join(map(str,index))+'.npz')
        if path.exists():coeff=np.load(path)['coefficients']
        else:
            p=np.array([[a[i] for a,i in zip(AXES,index)]])
            coeff=direct(root,geometry,p)['signal']@basis
            np.savez_compressed(path,coefficients=coeff)
        library[index]=coeff
        if n%15==0 or n==len(indices)-1:print(f'Moving-shadow node {n+1}/{len(indices)}',flush=True)
    physical=[]
    for t in [-2.,0.,2.]:
        for b,sign in enumerate([-1.,1.]):
            star=stellar_moments(np.full(len(centers),t),sign*centers,base['rp'])
            physical.append(np.einsum('nmk,nm->nk',left[b],star))
    physical=np.concatenate(physical)
    gram=np.zeros((basis.shape[1],basis.shape[1]))
    for coefficients in library.reshape(-1,count,basis.shape[1]):
        spectra=physical@coefficients;gram+=spectra.T@spectra
    eigen,v=np.linalg.eigh(gram)
    order=np.argsort(eigen)[::-1];eigen=np.maximum(eigen[order],0);v=v[:,order]
    rank=int(np.searchsorted(np.cumsum(eigen),eigen.sum()*(1-1e-10))+1)
    rotation=v[:,:rank];compressed=library@rotation
    target=root/'outputs/co_blind_shadow.npz'
    np.savez_compressed(target,centers=centers,left=left,library=compressed,rotation=rotation,
        basis=basis,rp=base['rp'],geometry_relative_error=error,cache_key=key)
    summary=dict(cache_key=key,n_positions=len(centers),n_geometric_modes=count,
        geometry_relative_frobenius_error=error,spectral_rank=rank,stellar_moments=4,
        spectral_compression_metric='actual shadows at three stellar positions; relative squared tail 1e-10',
        geometric_scope='circular aligned equatorial moon, central fully backed transit, both blocking channels union')
    (root/'docs/co_blind_library.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':build(Path(__file__).resolve().parents[1])
