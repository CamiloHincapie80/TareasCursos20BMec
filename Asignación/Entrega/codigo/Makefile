PYTHON ?= .venv/bin/python
export MPLCONFIGDIR := $(CURDIR)/resultados/.matplotlib
export LC_ALL := C.UTF-8
export SOURCE_DATE_EPOCH := 1794873600

.PHONY: entrega resultados reporte comprobar animacion

entrega: reporte
	latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=resultados/reporte documentos/presentacion.tex

resultados:
	sha256sum -c datos/SHA256SUMS
	"$(PYTHON)" main.py todo
	"$(PYTHON)" herramientas/resumen.py
	"$(PYTHON)" herramientas/reporte.py

reporte: resultados
	mkdir -p resultados/reporte
	latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=resultados/reporte documentos/reporte.tex
	pdfinfo resultados/reporte/reporte.pdf | awk '/^Pages:/ { seen=1; if ($$2 > 8) exit 1; ok=1 } END { if (!seen || !ok) exit 1 }'

comprobar:
	"$(PYTHON)" -m unittest discover -s tests -v

animacion:
	"$(PYTHON)" animaciones.py --preparar
