# Asistencia gravitacional de Júpiter a Saturno — UIS

Reto III.3 de Mecánica Clásica. El proyecto diseña un sobrevuelo de Júpiter para llegar a Saturno y mide el error de las cónicas empalmadas frente a la integración continua.

## Cómo leer la entrega

1. **Reporte técnico.** Explica la física del encuentro, el diseño circular desde la Tierra y el diseño ajustado con efemérides. La comparación del error usa el mismo estado de entrada para los dos modelos en cada frontera joviana.
2. **Repositorio Git.** El archivo `asistencia.bundle` contiene las fuentes y las seis efemérides en una sola revisión reproducible. Para abrirlo desde la carpeta de entrega, ejecuta `git clone ./asistencia.bundle asistencia`.
3. **Hoja de hipótesis.** Reúne las decisiones, alternativas y efectos esperados. [Hoja de hipótesis](documentos/HIPOTESIS.md).
4. **Bitácora 4D.** Registra el uso de IA, sus correcciones y las comprobaciones hechas. [Bitácora](documentos/BITACORA_4D.md).
5. **Presentación.** Las diapositivas se compilan con `make entrega`. La [guía de defensa](documentos/DEFENSA.md) indica qué debe poder explicar cada integrante.

Los datos y el código también están disponibles en `codigo.zip`. Para consultar solo las entradas originales se incluye `datos.zip`.

## Reproducir desde un clon limpio

Entorno verificado sobre Linux con Python 3.14, GNU Make, TeX Live con pdfLaTeX y soporte de español, `latexmk` y `pdfinfo` de Poppler. En Debian o Ubuntu, los paquetes de compilación son `make latexmk texlive-latex-extra texlive-lang-spanish poppler-utils`. Usa un intérprete Python 3.14 disponible en tu sistema.

Desde la raíz del proyecto:

```bash
python3.14 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
make entrega
make comprobar
```

`make entrega` verifica las seis entradas mediante SHA-256, ejecuta las simulaciones y diagnósticos, genera las figuras en PDF y PNG, construye tablas desde los resultados y compila el informe. La regla verifica el límite de ocho páginas. `make comprobar` ejecuta las pruebas y prepara las cachés que requieren; en un clon limpio puede tardar unos minutos. Los comandos usan `.venv/bin/python` por defecto. Para usar otro entorno, pasa `PYTHON=/ruta/al/python` a Make.

| Producto | Ubicación |
| --- | --- |
| PDF para entregar | `resultados/reporte/reporte.pdf` |
| Diapositivas | `resultados/reporte/presentacion.pdf` |
| Diez figuras PDF y PNG | `resultados/figuras/` |
| Tablas JSON, CSV y resúmenes Markdown | `resultados/tablas/` |
| Estados procesados regenerables | `datos/procesados/estados.npz` |
| Hipótesis y decisiones | [Hoja de hipótesis](documentos/HIPOTESIS.md) |
| Registro de IA | [Bitácora 4D](documentos/BITACORA_4D.md) |
| Exposición del equipo | [Guía de defensa](documentos/DEFENSA.md) |

El reporte y las diapositivas identifican al equipo y el correo de contacto. La bitácora reserva espacios para las interacciones previas, la rotación de roles y la revisión manual que confirme el equipo.

## Animación interactiva

```bash
.venv/bin/python -m pip install -r requirements-vispy.txt
make animacion
./abrir_animaciones.sh
```

La primera preparación recalcula las trayectorias ajustadas al lanzamiento. Para una prueba automatizada de interfaz sin pantalla en Linux con Mesa/EGL usa `./abrir_animaciones.sh --comprobar --sin-pantalla`. La reproducción de figuras y del informe usa un backend de gráficos sin interfaz. [Controles de la animación](documentos/ANIMACIONES.md).

## Estructura versionada

```text
.
├── Makefile
├── README.md
├── main.py
├── animaciones.py
├── abrir_animaciones.sh
├── requirements.txt
├── requirements-vispy.txt
├── datos/
│   ├── SHA256SUMS
│   └── horizons/       # cuatro series diarias y dos series horarias
├── documentos/        # reporte.tex, hipótesis, bitácora, defensa y notas
├── herramientas/      # tablas de lectura y materiales para LaTeX
├── src/               # modelos, validación, gráficos e interfaz
└── tests/             # pruebas de física y trayectorias
```

El [inventario de auditoría](documentos/AUDITORIA.md) indica los archivos locales excluidos. La [formulación física](documentos/FISICA.md) y la [validación](documentos/VALIDACION.md) explican las ecuaciones y los límites de los resultados. [Procedencia de las efemérides](documentos/FUENTES.md).
