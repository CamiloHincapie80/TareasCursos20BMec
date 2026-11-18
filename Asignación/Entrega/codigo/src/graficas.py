"""Figuras académicas reproducibles, km/s y fechas TDB."""
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/asistencia-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from .constantes import MU,GM,AU,DAY
from .dinamica import norm,invariants
from .encuentros import calendar
from .hohmann import calculate,R1,R2
from .optimizacion import design

OUT=Path(__file__).resolve().parents[1]/'resultados'/'figuras'
OUT.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.size':10,'axes.grid':True,'figure.figsize':(7,4.6)})

def save(fig,name):
    fig.tight_layout();fig.savefig(OUT/(name+'.png'),dpi=220);fig.savefig(OUT/(name+'.pdf'));plt.close(fig)

def idealized():
    h=calculate();a=h['a_km'];e=h['e'];f=np.linspace(0,np.pi,600)
    r=a*(1-e*e)/(1+e*np.cos(f));theta=np.linspace(0,2*np.pi,600)
    fig,ax=plt.subplots()
    for R,label in [(R1,'Tierra circular'),(R2,'Júpiter circular')]:ax.plot(R*np.cos(theta)/AU,R*np.sin(theta)/AU,label=label)
    ax.plot(r*np.cos(f)/AU,r*np.sin(f)/AU,label='Transferencia Hohmann');ax.scatter([0],[0],label='Sol')
    ax.set(xlabel='x heliocéntrica [UA]',ylabel='y heliocéntrica [UA]',aspect='equal',title='Modelo idealizado circular coplanar');ax.legend();save(fig,'01_hohmann')
    d,sol,ps,pbad,tm=design();ts=np.linspace(0,tm+60*DAY,700);pos=sol.sol(ts)[:3].T
    fig,ax=plt.subplots();ax.plot(pos[:,0]/AU,pos[:,1]/AU,label='Sonda diseñada')
    for fn,label in [(ps,'Saturno: fase diseñada'),(pbad,'Saturno: fase inicial que falla')]:
        p=np.array([fn(t) for t in ts]);ax.plot(p[:,0]/AU,p[:,1]/AU,label=label)
        p=fn(tm);ax.scatter(p[0]/AU,p[1]/AU)
    ax.set(xlabel='x heliocéntrica [UA]',ylabel='y heliocéntrica [UA]',aspect='equal',title='Diseño idealizado: fase de Saturno libre');ax.legend();save(fig,'09_diseno')

def historical(sys,epoch,ref,n,runs,comparison,windows,cs):
    a=windows['J'][0];b=windows['S'][1];ts=np.linspace(a,b,1500)
    fig,ax=plt.subplots();y=np.array([ref(t) for t in ts]);ax.plot(y[:,0]/AU,y[:,1]/AU,label='Voyager: Horizons aproximado')
    for k in ('J','S'):
        p=sys.planets[k](ts);ax.plot(p[:,0]/AU,p[:,1]/AU,label=f'{k}: órbita durante el intervalo')
        for t in (windows[k][0],windows[k][1]):
            p=sys.planets[k](t);ax.scatter(p[0]/AU,p[1]/AU,s=12);ax.annotate(calendar(epoch+t/DAY)[:10],(p[0]/AU,p[1]/AU),fontsize=7)
    ax.set(xlabel='x eclíptica heliocéntrica [UA]',ylabel='y eclíptica heliocéntrica [UA]',aspect='equal',title='Julio de 1979 – septiembre de 1981 (TDB)');ax.legend();save(fig,'02_historia')
    ts=np.linspace(cs[0],cs[-1],1200);rel=np.array([n(t)-sys.state('J',t) for t in ts])
    fig,ax=plt.subplots();ax.plot(rel[:,0]/1e6,rel[:,1]/1e6,label='N perturbado (proyección xy)')
    r=sys.soi('J',cs[0])/1e6;theta=np.linspace(0,2*np.pi,400);ax.plot(r*np.cos(theta),r*np.sin(theta),'--',label='SOI en entrada')
    for i,label in [(0,'Entrada finita'),(-1,'Salida finita')]:
        q=rel[i,:2]/1e6;v=rel[i,3:5];ax.annotate('',xy=q+v*.5,xytext=q,arrowprops={'arrowstyle':'->'});ax.text(*(q+v*.5),label,fontsize=8)
    i=np.argmin(np.linalg.norm(rel[:,:3],axis=1));ax.scatter(*rel[i,:2]/1e6,label='Periapsis')
    ax.set(xlabel='x jovicéntrica [10⁶ km]',ylabel='y jovicéntrica [10⁶ km]',aspect='equal');ax.legend();save(fig,'03_dispersion')
    fig,ax=plt.subplots();days=(ts-cs[0])/DAY
    for name,sol in runs.items():
        ax.plot(days,[invariants(sol(t),MU)['energy_km2s2'] for t in ts],label=name)
    ax.set(xlabel='Días TDB desde entrada SOI de N',ylabel='Energía kepleriana solar [km²/s²]');ax.legend();save(fig,'04_energia')
    fig,axs=plt.subplots(2,1,sharex=True)
    for name,sol in runs.items():
        vals=[invariants(sol(t)-sys.state('J',t),GM['J']) if name!='K' or t<=sol.segments[0][1].t[-1] else {'energy_km2s2':np.nan,'h_km2s':np.nan} for t in ts]
        axs[0].plot(days,[v['energy_km2s2'] for v in vals],label=name)
        axs[1].plot(days,[v['h_km2s'] for v in vals],label=name)
    axs[0].set_ylabel('ε J [km²/s²]');axs[1].set(ylabel='|h J| [km²/s]',xlabel='Días TDB desde entrada SOI N');axs[0].legend();save(fig,'05_invariantes')
    fig,axs=plt.subplots(3,1,sharex=True,figsize=(7,7));ts2=np.linspace(cs[0],b,600)
    for name,sol in runs.items():
        if name=='N':continue
        dy=np.array([sol(t)-n(t) for t in ts2]);
        metrics=[np.linalg.norm(dy[:,:3],axis=1),np.linalg.norm(dy[:,3:],axis=1),np.abs([invariants(sol(t),MU)['energy_km2s2']-invariants(n(t),MU)['energy_km2s2'] for t in ts2])]
        for ax,metric in zip(axs,metrics):ax.semilogy((ts2-cs[0])/DAY,np.maximum(metric,1e-12),label=name)
    for ax,label in zip(axs,['|Δr| [km]','|Δv| [km/s]','|Δε solar| [km²/s²]']):ax.set_ylabel(label)
    axs[-1].set_xlabel('Días TDB desde entrada SOI N');axs[0].legend();save(fig,'06_errores')
    fig,ax=plt.subplots()
    for name in ('K','P0','P'):
        rows=[r for r in comparison if r['model']==name];ax.semilogy([r['f'] for r in rows],[max(abs(r['drp_S_km']),1e-6) for r in rows],'o-',label=name)
    ax.axhline(1000,color='gray',ls='--',label='Tolerancia declarada: 1000 km');ax.set(xlabel='Factor f del radio instantáneo de Laplace',ylabel='Error de distancia mínima a Saturno [km]');ax.legend();save(fig,'07_radio')
    fig,ax=plt.subplots();ts=np.linspace(*windows['S'],1000)
    for name,sol in dict(Efemerides=ref,**runs).items():
        ax.semilogy((ts-windows['S'][0])/DAY,[norm(sol(t)[:3]-sys.planets['S'](t)) for t in ts],label=name)
    ax.set(xlabel='Días TDB desde 1981-08-01',ylabel='Distancia al centro de Saturno [km]');ax.legend();save(fig,'08_saturno')
