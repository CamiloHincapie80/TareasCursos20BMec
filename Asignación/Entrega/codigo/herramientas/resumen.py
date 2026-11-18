"""Genera tablas de lectura y materiales de reporte desde JSON, sin copiar cifras."""
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[1]
D=ROOT/'resultados'/'tablas'
load=lambda name:json.loads((D/(name+'.json')).read_text())
e=load('encuentros');c=load('comparacion');d=load('diseno')
lines=['# Resultados calculados','',
       'Las cifras provienen de los JSON. Condiciones y límites en documentos/VALIDACION.md. Las efemérides diarias en Saturno son insuficientes para certificar el periapsis.','',
       '## Encuentros con f=1','',
       '| Planeta | Modelo | Fecha TDB | Distancia central [km] | Rapidez relativa [km/s] |',
       '|---|---|---|---:|---:|']
for r in e:lines.append(f"| {r['planet']} | {r['model']} | {r['date_tdb']} | {r['rp_km']:.1f} | {r['vp_kms']:.5f} |")
lines+=['','## Sensibilidad a frontera','',
        'Estado al 1981-10-01; distancia mínima S en enero–octubre de 1981. El estado inicial es el de N al entrar a cada frontera joviana.','',
        '| f | Modelo | Error posición [km] | Error velocidad [km/s] | Error mínimo S [km] | Equivalente (<1000 km) |',
        '|---:|---|---:|---:|---:|---|']
for r in c:lines.append(f"| {r['f']} | {r['model']} | {r['dr_km']:.2f} | {r['dv_kms']:.6g} | {r['drp_S_km']:+.2f} | {r['equivalent_1000km']} |")
lines+=['','## Diseño','',f"Hohmann + asistencia asintótica ideal: rp_J={d['rp_jupiter_km']:.0f} km; δ={d['delta_deg']:.4f}°; Δε={d['energy_gain_km2s2']:.4f} km²/s².",
        f"Con Sol+Saturno circular prescrito: fase ajustada {d['refined_phase_deg']:.6f}°, mínimo {d['refined_saturn_rp_km']:.4f} km, tiempo {d['refined_time_days']:.4f} días. No es diseño continuo con efemérides reales."]
(D/'resumen.md').write_text('\n'.join(lines)+'\n')
# Tabla exigida: errores a época fija, no diferencia entre dos posiciones en épocas diferentes.
sat=[r for r in e if r['planet']=='S']
lines=['# Comparación en Saturno (f=1)','','Referencia de Voyager: interpolante diario insuficiente; no precisión histórica certificada.','',
       '| Magnitud | Efemérides | N | P | K |','|---|---|---|---|---|']
for title,key,fmt in [('Fecha TDB','date_tdb',str),('Distancia mínima [km]','rp_km',lambda x:f'{x:.1f}'),('Rapidez relativa [km/s]','vp_kms',lambda x:f'{x:.5f}'),('Error de posición a periapsis de referencia [km]','position_error_at_reference_periapsis_km',lambda x:f'{x:.1f}'),('Error temporal [s]','periapsis_time_error_s',lambda x:f'{x:+.1f}')]:
    vals=[]
    for model in ['Efemerides','N','P','K']:
        r=next(a for a in sat if a['model']==model);vals.append(fmt(r.get(key,float('nan'))))
    lines.append('| '+title+' | '+' | '.join(vals)+' |')
(D/'saturno.md').write_text('\n'.join(lines)+'\n')
