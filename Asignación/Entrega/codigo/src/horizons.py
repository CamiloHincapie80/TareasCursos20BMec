"""Vectores geométricos Horizons: CSV por etiquetas o bloques X=/VX=.
Unidades internas km, s; tiempos JD TDB. Nunca extrapola.
"""
from dataclasses import dataclass
from pathlib import Path
import csv
import hashlib
import re
import numpy as np
from scipy.interpolate import CubicHermiteSpline

NUM = r'[-+]?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?'

@dataclass
class Ephemeris:
    path: Path
    jd: np.ndarray
    states: np.ndarray
    meta: dict

    def spline(self, epoch=None, stride=1):
        epoch = self.jd[0] if epoch is None else epoch
        return CubicHermiteSpline((self.jd[::stride]-epoch)*86400,
                                 self.states[::stride,:3], self.states[::stride,3:],
                                 extrapolate=False)

    def audit(self):
        dt = np.diff(self.jd)*86400
        s = self.spline(stride=2)
        q = (self.jd[1:-1:2]-self.jd[0])*86400
        truth = self.states[1:-1:2]
        dr = np.linalg.norm(s(q)-truth[:,:3],axis=1)
        dv = np.linalg.norm(s(q,1)-truth[:,3:],axis=1)
        return dict(file=self.path.name, sha256=hashlib.sha256(self.path.read_bytes()).hexdigest(),
                    metadata=self.meta, rows=len(self.jd), jd_start=float(self.jd[0]),
                    jd_end=float(self.jd[-1]), step_s=float(np.median(dt)),
                    gaps=int(np.sum(dt>1.5*np.median(dt))),
                    holdout_position_max_km=float(dr.max()),
                    holdout_velocity_max_kms=float(dv.max()))

def read(path):
    path = Path(path)
    text = path.read_text(encoding='utf-8')
    if text.count('$$SOE')!=1 or text.count('$$EOE')!=1:
        raise ValueError(f'{path}: se requiere un bloque $$SOE/$$EOE')
    header, rest = text.split('$$SOE')
    body = rest.split('$$EOE')[0]
    def field(key):
        m = re.search(r'^'+re.escape(key)+r'\s*:\s*(.*)$',header,re.M)
        if not m: raise ValueError(f'{path}: falta metadato {key}')
        return m[1].strip()
    meta = {k:field(k) for k in ['Target body name','Center body name','Output units',
                                'Reference frame','Output type']}
    scale = re.search(r'Start time\s*:.*\b(TDB|TT|UT)\s*$', header,re.M)
    if not scale: raise ValueError(f'{path}: escala temporal no reconocida')
    meta['time_scale'] = scale[1]
    if meta['Output units'] not in ('KM-S','KM-D','AU-D'):
        raise ValueError(f'{path}: unidades incompatibles {meta["Output units"]}')
    if not meta['Output type'].startswith('GEOMETRIC'):
        raise ValueError(f'{path}: se requieren estados geométricos sin aberración')
    cols = None
    for line in header.splitlines():
        if ',' in line and re.search(r'\bJDTDB\b|\bJDTT\b|\bJDUT\b',line):
            cols = [x.strip() for x in next(csv.reader([line]))]
    jd, states = [], []
    if cols:
        required = ['X','Y','Z','VX','VY','VZ']
        if any(k not in cols for k in required):
            raise ValueError(f'{path}: faltan posiciones o velocidades: {cols}')
        indices = [cols.index(k) for k in required]
        for line in body.splitlines():
            if not line.strip(): continue
            row = next(csv.reader([line]))
            try:
                jd.append(float(row[0])); states.append([float(row[i].replace('D','E')) for i in indices])
            except (ValueError,IndexError) as exc:
                raise ValueError(f'{path}: fila vectorial incompatible: {line[:100]}') from exc
    else:
        current = None
        values = {}
        for line in body.splitlines():
            m = re.match(r'\s*(\d{7}\.\d+)',line)
            if m:
                if current is not None:
                    if len(values)!=6: raise ValueError(f'{path}: bloque incompleto en JD {current}')
                    jd.append(current); states.append([values[k] for k in ['X','Y','Z','VX','VY','VZ']])
                current=float(m[1]); values={}
            for key,val in re.findall(r'\b(VX|VY|VZ|X|Y|Z)\s*=\s*('+NUM+r')',line):
                values[key]=float(val.replace('D','E'))
        if current is not None:
            if len(values)!=6: raise ValueError(f'{path}: último bloque incompleto')
            jd.append(current); states.append([values[k] for k in ['X','Y','Z','VX','VY','VZ']])
    jd, states = np.asarray(jd), np.asarray(states)
    if not np.all(np.isfinite(jd)): raise ValueError(f'{path}: tiempos no finitos')
    if len(jd)<3 or states.shape!=(len(jd),6) or not np.all(np.isfinite(states)):
        raise ValueError(f'{path}: tabla vacía, insuficiente o no finita')
    if np.any(np.diff(jd)<=0): raise ValueError(f'{path}: tiempos duplicados o desordenados')
    unit=meta['Output units']
    if unit=='AU-D': states[:,:3]*=149597870.7; states[:,3:]*=149597870.7/86400
    elif unit=='KM-D': states[:,3:]/=86400
    return Ephemeris(path,jd,states,meta)

def compatible(*items):
    for key in ['Center body name','Reference frame','time_scale','Output type']:
        vals=[e.meta[key].split('{')[0].strip() for e in items]
        if len(set(vals))!=1: raise ValueError(f'Metadato incompatible {key}: {vals}')
    if max(e.jd[0] for e in items)>=min(e.jd[-1] for e in items):
        raise ValueError('Las efemérides no tienen intervalo común')
