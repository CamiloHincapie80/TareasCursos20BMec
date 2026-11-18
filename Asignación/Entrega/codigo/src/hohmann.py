"""Modelo circular idealizado original; no usa estados de Voyager."""
import numpy as np
from .constantes import MU,GM,RAD,AU,DAY

R1=1.49598e8
R2=7.78479e8
R3=9.58*AU  # radio circular elegido para Saturno, no efeméride

def calculate():
    a=(R1+R2)/2; e=(R2-R1)/(R1+R2)
    v1=np.sqrt(MU/R1); v2=np.sqrt(MU/R2)
    vt1=np.sqrt(MU*(2/R1-1/a)); vt2=R1*vt1/R2
    v0=np.sqrt(398600.4418/(6378.137+500))
    u=v2-vt2; rp=10*RAD['J']; eh=1+rp*u*u/GM['J']
    tof=np.pi*np.sqrt(a**3/MU)
    return dict(a_km=a,e=e,v_leo_kms=v0,v_earth_kms=v1,v_jupiter_kms=v2,
                vt1_kms=vt1,vt2_kms=vt2,vinf_earth_kms=vt1-v1,vinf_jupiter_kms=u,
                dv_leo_kms=np.sqrt((vt1-v1)**2+2*v0*v0)-v0,
                tof_days=tof/DAY,phase_deg=np.degrees(np.pi-np.sqrt(MU/R2**3)*tof),
                rp_jupiter_km=rp,e_hyperbola=eh,delta_deg=np.degrees(2*np.arcsin(1/eh)),
                impact_km=rp*np.sqrt(1+2*GM['J']/(rp*u*u)))
