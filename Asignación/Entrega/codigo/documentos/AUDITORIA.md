# Inventario y alcance de la entrega

| Ruta | Decisión | Razón |
| --- | --- | --- |
| `src/`, `main.py`, `animaciones.py`, `abrir_animaciones.sh` | Versionar | Lectura de datos, modelos, validación, figuras y simulación. |
| `datos/horizons/*.txt`, `datos/SHA256SUMS` | Versionar | Seis fuentes originales y verificación de integridad. La pareja horaria aporta resolución cerca de Júpiter y las cuatro series diarias cubren el viaje. |
| `herramientas/`, `tests/`, `Makefile`, `requirements*.txt` | Versionar | Regeneración de tablas y figuras, dependencias y pruebas. |
| `documentos/reporte.tex`, `documentos/presentacion.tex`, `HIPOTESIS.md`, `BITACORA_4D.md`, `DEFENSA.md` | Versionar | Fuentes y entregables del curso. |
| `documentos/FISICA.md`, `VALIDACION.md`, `FUENTES.md`, `ANIMACIONES.md` | Versionar | Derivaciones, verificación, trazabilidad y uso de la interfaz. |
| `ModeloReporteUIS.pdf` | Referencia local excluida | Plantilla del profesor. El informe se compila con su propia fuente. |
| `originales/`, `documentos/MANIFIESTO_ORIGINALES.json` | Archivo local excluido | Versiones iniciales usadas para trazabilidad histórica. Ningún comando de entrega las necesita. |
| `datos/procesados/` | Candidato a limpieza local | El NPZ se genera de nuevo. Los tres CSV de compatibilidad solo servían al prototipo archivado. |
| `resultados/` | Candidato a limpieza local | Tablas, figuras, PDF final y cachés se regeneran. El PDF final se adjunta a la entrega del curso. |
| `documentos/material_reporte.*`, `informe_voyager_saturno.*`, `distancias_saturno_comparacion.pdf`, `RESUMEN_RESULTADOS.md`, `SATURNO.md` | Candidatos a archivar | Borradores y salidas previas. La entrega actual usa `reporte.tex` y escribe productos en `resultados/`. |
| `.venv/`, `__pycache__/`, `.pytest_cache/`, `.ipynb_checkpoints/`, auxiliares LaTeX | Candidatos a limpieza local | Productos específicos de una máquina o compilación. |
| `.agents/`, `.aws/`, `.codex/`, `.git/` | Administración local | El flujo de entrega los deja fuera de la copia de fuentes. |

Los candidatos están excluidos de Git y se dejaron en el directorio de trabajo para evitar pérdida de información. La verificación de entrega parte de una copia con únicamente archivos versionables.
