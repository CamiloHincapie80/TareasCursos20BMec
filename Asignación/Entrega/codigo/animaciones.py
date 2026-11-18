#!/usr/bin/env python3
"""Animaciones científicas interactivas VisPy + Qt, con caché reproducible."""
import argparse
import subprocess
import sys
from src.animacion_comparacion import comparison

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('modo',nargs='?',default='comparacion',choices=['comparacion','hohmann','esferas','ambas'],help='Vista superpuesta; los nombres anteriores se conservan como alias')
    parser.add_argument('--factor',type=float,default=1.,choices=[.1,.25,.5,1.])
    parser.add_argument('--preparar',action='store_true',help='Calcular cachés sin abrir ventanas ni requerir OpenGL')
    parser.add_argument('--recalcular',action='store_true',help='Ignorar caché, conservando datos originales')
    parser.add_argument('--comprobar',action='store_true',help='Ejercitar controles y render, sin capturas, y cerrar')
    parser.add_argument('--sin-pantalla',action='store_true',help='Solo con --comprobar: Qt offscreen y OpenGL EGL')
    args=parser.parse_args()
    if args.sin_pantalla:
        if not args.comprobar:parser.error('--sin-pantalla requiere --comprobar')
        import os
        os.environ['QT_QPA_PLATFORM']='offscreen'
        os.environ.setdefault('EGL_PLATFORM','surfaceless')
        os.environ.setdefault('LIBGL_ALWAYS_SOFTWARE','1')
        os.environ.setdefault('MPLCONFIGDIR','/tmp/asistencia-mpl')
    if args.preparar:
        comparison(args.factor,args.recalcular,launch=True)
        print('Cachés listas en resultados/animaciones.');return
    try:
        from src.animacion_vispy import launch
    except ImportError as exc:
        parser.exit(1,f'Falta dependencia de visualización: {exc}.\nInstale: python -m pip install -r requirements-vispy.txt\n')
    # Qt aborta si no abre el servidor gráfico. Comprobarlo en un proceso aparte.
    probe = subprocess.run(
        [sys.executable, '-c', 'from PyQt6.QtWidgets import QApplication; app = QApplication([])'],
        capture_output=True, text=True,
    )
    if probe.returncode:
        parser.exit(1, 'Qt no pudo abrir la interfaz gráfica. Ejecute el lanzador desde una terminal de la sesión gráfica.\n'
                    + probe.stderr)
    data=comparison(args.factor,args.recalcular,launch=True)
    raise SystemExit(launch(data,smoke=args.comprobar,headless=args.sin_pantalla))

if __name__=='__main__':main()
