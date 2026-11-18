"""Genera los cuadros del informe UIS y la figura del sobrevuelo de Saturno."""
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.constantes import RAD
from src.transferencia_leo import circle, nominal_saturn_leg

TABLES = ROOT / 'resultados' / 'tablas'
FIGURES = ROOT / 'resultados' / 'figuras'


def load(name):
    return json.loads((TABLES / f'{name}.json').read_text())


def number(value, decimals=0):
    return f'{value:,.{decimals}f}'.replace(',', r'\,').replace('.', '{,}')


def row(*columns):
    return ' & '.join(columns) + r' \\'


def main():
    encounters = load('encuentros')
    comparison = load('comparacion')
    design = load('diseno')
    hohmann = load('hohmann')
    macros = {
        'ImpulsoLEO': number(hohmann['dv_leo_kms'], 3),
        'DiasHohmann': number(hohmann['tof_days'], 2),
        'GiroJoviano': number(design['delta_deg'], 2),
        'ExcesoJoviano': number(design['vinf_in_kms'], 3),
        'EnergiaGanada': number(design['energy_gain_km2s2'], 2),
        'MinimoDiseno': number(design['refined_saturn_rp_km']),
        'DiasDiseno': number(design['refined_time_days'], 2),
    }
    doc = [rf'\newcommand{{\{key}}}{{\ensuremath{{{value}}}}}' for key, value in macros.items()]
    (TABLES / 'reporte-cifras.tex').write_text('\n'.join(doc)+'\n')
    rows = []
    for name in ('Efemerides', 'N', 'P', 'P0', 'K'):
        jovian = next(r for r in encounters if r['planet']=='J' and r['model']==name)
        saturn = next(r for r in encounters if r['planet']=='S' and r['model']==name)
        rows.append(row('Horizons' if name=='Efemerides' else name,
                        f"${number(jovian['rp_km'])}$", f"${number(saturn['rp_km'])}$",
                        saturn['date_tdb'][:16].replace('T',' ')))
    (TABLES / 'reporte-encuentros.tex').write_text('\n'.join([*rows, r'\bottomrule'])+'\n')
    rows = []
    for factor in (.1,.25,.5,1.):
        values = {r['model']:r['drp_S_km'] for r in comparison if r['f']==factor}
        rows.append(row(f'${number(factor,2)}$', *(f'${number(values[m],1)}$' for m in ('P','P0','K'))))
    (TABLES / 'reporte-frontera.tex').write_text('\n'.join([*rows, r'\bottomrule'])+'\n')

    leg = nominal_saturn_leg()
    xy = leg.relative_path/1e6
    mid = len(xy)//2
    fig, ax = plt.subplots(figsize=(5.5,3.1), layout='constrained')
    ax.fill(*circle(RAD['S']/1e6,np.linspace(0,2*np.pi,180)).T,color='#a67e53',label='Saturno')
    ax.plot(*xy[:mid+1].T,color='#15768b',label='Aproximación')
    ax.plot(*xy[mid:].T,color='#b84b29',label='Salida')
    ax.scatter(*xy[mid],color='black',s=16,label='Mínimo')
    ax.set(xlabel='x saturnocéntrica (millones de km)',ylabel='y saturnocéntrica (millones de km)',
           aspect='equal')
    ax.grid(alpha=.2)
    ax.legend(loc='center left',bbox_to_anchor=(1.04,.5),frameon=False)
    for suffix in ('pdf','png'):
        fig.savefig(FIGURES / f'10_sobrevuelo_saturno.{suffix}',dpi=220)
    plt.close(fig)


if __name__ == '__main__':
    main()
