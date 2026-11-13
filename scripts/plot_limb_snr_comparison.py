"""Plot the tracked ingress/egress S/N table without rerunning models."""
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    root=Path(__file__).resolve().parents[1]
    with (root/'docs/ingress_egress_comparison.csv').open() as f:rows=list(csv.DictReader(f))
    fig,axes=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
    for ax,name in zip(axes,['ingress','egress']):
        subset=[r for r in rows if r['encounter']==name]
        phases=[float(r['exposure_phase']) for r in subset]
        for key,label,style in [('point_star_snr','Point-sampled star','o--'),
                                ('joint_snr','Joint cell integration','s-'),
                                ('refined_snr','Joint refined','x:')]:
            ax.plot(phases,[float(r[key]) for r in subset],style,label=label)
        ax.set(title=name,xlabel='Exposure start / duration',ylabel='S/N in selected limb window')
        ax.legend(fontsize=8)
    for ext in ['png','pdf']:fig.savefig(root/f'outputs/limb_snr_comparison.{ext}',dpi=180)


if __name__=='__main__':main()
