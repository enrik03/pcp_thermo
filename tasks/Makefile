# PCP Thermodynamics – Task Automation
# Usage:
#   make figures   – generate all figures in figures/
#   make test      – run unit tests with pytest
#   make clean     – remove generated figures

PYTHON   ?= python
PYTEST   ?= pytest
OUTDIR   ?= figures
DPI      ?= 150
FORMAT   ?= png,pdf

.PHONY: figures test clean help

help:
	@echo "Available targets:"
	@echo "  make figures   – generate PCP pressure figures"
	@echo "  make test      – run unit tests"
	@echo "  make clean     – remove generated figures"

figures:
	@echo "Generating figures → $(OUTDIR)/"
	MPLBACKEND=Agg $(PYTHON) ../pcp_thermo_figures_corrected.py \
		--outdir $(OUTDIR) \
		--dpi $(DPI) \
		--format $(FORMAT)

test:
	@echo "Running tests…"
	$(PYTEST) ../tests/ -v

clean:
	@echo "Removing $(OUTDIR)/"
	rm -rf $(OUTDIR)
