SHELL := /bin/bash

ROOT := $(abspath ..)
PYTHON ?= python3
DATA_ROOT ?= $(CURDIR)/data/nstx
CASE ?= 141711Y01
CASES ?= 141711Y01 141711Y02
CASE_DIR ?= $(DATA_ROOT)/$(CASE)
CDF ?= $(CASE_DIR)/$(CASE).CDF
TRDAT ?= $(CASE_DIR)/$(CASE)TR.DAT
PLOT ?= $(CURDIR)/q-$(CASE).png

F90 ?= gfortran
FFLAGS ?= -cpp -g -O0 -fPIE
CPPFLAGS += -I$(IMAS_PREFIX)/include/gfortran
LDLIBS += -L$(IMAS_PREFIX)/lib -limas-gfortran -limas

.DEFAULT_GOAL := help

.PHONY: help quick verify verify-all verify-dd inspect plot-ids fortran clean

help:
	@echo "Diagnostic and example commands"
	@echo
	@echo "  make quick       Run checks that do not scan the full CDF"
	@echo "  make verify      Run source-side checks for CASE"
	@echo "  make verify-all  Run source-side checks for both NSTX cases"
	@echo "  make verify-dd   Inspect DD 3.39.0 and DD 4.1.1 schemas"
	@echo "  make inspect     Write a q-profile plot to PLOT"
	@echo "  make plot-ids    Show the converted-IDS plotting example usage"
	@echo "  make fortran     Build build/verify using IMAS_PREFIX"

quick:
	$(PYTHON) -m py_compile $(wildcard *.py)
	$(PYTHON) verify_writer_api.py

verify: quick
	@test -f "$(CDF)" || { echo "CDF not found: $(CDF)" >&2; exit 1; }
	@test -f "$(TRDAT)" || { echo "TR.DAT not found: $(TRDAT)" >&2; exit 1; }
	$(PYTHON) verify_distribution_sources_namelist.py "$(TRDAT)"
	$(PYTHON) verify_cdf_extension.py "$(CDF)" "$(TRDAT)"
	$(PYTHON) verify_cocos17.py "$(CDF)" --dat "$(TRDAT)"
	$(PYTHON) verify_nbi_distribution_sources.py "$(CDF)"
	$(PYTHON) verify_sources_transport.py "$(CDF)"
	$(PYTHON) verify_next_ids_inventory.py "$(CDF)" "$(TRDAT)"
	$(PYTHON) verify_metadata_provenance_inputs.py "$(CDF)" "$(TRDAT)"

verify-all:
	@set -e; for case in $(CASES); do \
		echo "Verifying NSTX case $$case"; \
		$(MAKE) --no-print-directory verify CASE="$$case" PYTHON="$(PYTHON)"; \
	done

verify-dd:
	$(PYTHON) verify_dd_compat.py 3.39.0 4.1.1

inspect:
	$(PYTHON) inspect_cdf.py "$(CDF)" --output "$(PLOT)"

plot-ids:
	@echo "$(PYTHON) examples/plot_ids.py <runid> <database> [HDF5|MDSPLUS] [quantity]"

build:
	@mkdir -p $@

build/verify: verify.f90 | build
	$(F90) $(FFLAGS) $(CPPFLAGS) $< $(LDLIBS) -o $@

fortran: build/verify

clean:
	$(RM) -r __pycache__ build q-*.png
