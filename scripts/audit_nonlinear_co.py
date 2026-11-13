"""Outer-boundary and legacy-map regression for the bounded CO experiment."""
import json
from pathlib import Path
import numpy as np
from molecular_injection_recovery import build_molecular_geometry,synthesize
from nonlinear_co_recovery import direct


def main():
    root=Path(__file__).resolve().parents[1]
    saved=np.load(root/'outputs/co_nonlinear_emulator.npz');noise=float(saved['noise'])/130.8285010091114
    original_path=root/'outputs/co_molecular_injection.npz';regression=None
    if original_path.exists():
        original=np.load(original_path);legacy=synthesize(root)
        regression={name:float(np.max(abs(legacy[name]-original[name]))) for name in ['reference','signal']}
        assert all(np.allclose(legacy[n],original[n],rtol=1e-12,atol=1e-18) for n in regression)
    geometry=build_molecular_geometry(nr=64,n_scale_heights=20.,split_limbs=True)
    base=direct(root,geometry,np.zeros((2,3)))
    production_geometry=build_molecular_geometry(nr=52,n_scale_heights=16.,split_limbs=True)
    rows=[]
    for label,p in [('reference',np.zeros((2,3))),('upper_T_CO',np.array([[1.,1.,0.],[1.,1.,0.]]))]:
        high=direct(root,geometry,p);low=direct(root,production_geometry,p)
        rows.append(dict(atmosphere=label,reference_difference_control_noise_norm=float(np.linalg.norm(high['reference']-low['reference'])/noise),
            moon_difference_control_noise_norm=float(np.linalg.norm(high['signal']-low['signal'])/noise)))
    result=dict(previous_pilot_max_absolute_differences=regression,outer_boundary_16H_vs20H=rows,
        control_noise_ppm=noise*1e6)
    (root/'docs/co_nonlinear_geometry_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
