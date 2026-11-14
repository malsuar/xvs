"""Joint refinement of the new molecular synthesis, not just Gaussian maps."""
import json
from pathlib import Path
import numpy as np
from molecular_injection_recovery import synthesize
from xds.molecular import depth_orthogonal_component


def main():
    root=Path(__file__).resolve().parents[1]
    old=np.load(root/'outputs/co_molecular_injection.npz')
    fine=synthesize(root,na=960,nr=64,nt=2401,dv=.1,shift_step=.125,wing_sigma=200.)
    assert np.allclose(old['wavenumber'],fine['wavenumber'])
    s=old['signal'];f=fine['signal']
    q=depth_orthogonal_component(s,old['reference']);fq=depth_orthogonal_component(f,fine['reference'])
    values=dict(full_snr_relative_change=float(abs(np.linalg.norm(s)/np.linalg.norm(f)-1)),
        depth_profiled_snr_relative_change=float(abs(np.linalg.norm(q)/np.linalg.norm(fq)-1)),
        residual_map_relative_norm_difference=float(np.linalg.norm(s-f)/np.linalg.norm(f)),
        reference_map_relative_norm_difference=float(np.linalg.norm(old['reference']-fine['reference'])/np.linalg.norm(fine['reference'])),
        fine_full_snr=float(np.linalg.norm(f)/fine['noise']),fine_depth_profiled_snr=float(np.linalg.norm(fq)/fine['noise']),
        passes_1percent=bool(max(abs(np.linalg.norm(s)/np.linalg.norm(f)-1),abs(np.linalg.norm(q)/np.linalg.norm(fq)-1))<.01))
    # Maintain approximately the same radial sampling per scale height.
    for nh,nr in [(5.,16),(12.,40)]:
        atmosphere=synthesize(root,nr=nr,n_scale_heights=nh)
        aq=depth_orthogonal_component(atmosphere['signal'],atmosphere['reference'])
        values[f'heights_{nh}_full_snr']=float(np.linalg.norm(atmosphere['signal'])/atmosphere['noise'])
        values[f'heights_{nh}_kinematic_snr']=float(np.linalg.norm(aq)/atmosphere['noise'])
        values[f'heights_{nh}_kinematic_change_from10']=float(abs(np.linalg.norm(aq)/np.linalg.norm(q)-1))
    (root/'docs/co_template_convergence.json').write_text(json.dumps(values,indent=2)+'\n')
    print(json.dumps(values,indent=2))


if __name__=='__main__':main()
