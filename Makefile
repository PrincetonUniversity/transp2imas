SHELL := /bin/bash

PYTHON ?= python3
VENV ?= .venv
VENV_PYTHON := $(VENV)/bin/python
DATA_ROOT ?= $(CURDIR)/examples/data/nstx
CASE ?= 141711Y01
CASES ?= 141711Y01 141711Y02
CASE_DIR ?= $(DATA_ROOT)/$(CASE)
CDF ?= $(CASE_DIR)/$(CASE).CDF
TRDAT ?= $(CASE_DIR)/$(CASE)TR.DAT
BACKEND ?= HDF5
DD_VERSION ?= auto
OUTPUT_DIR ?= $(CURDIR)/build/imas
CONVERT_FLAGS ?= -new
PF_PASSIVE_GEOMETRY ?=
PF_ACTIVE_GEOMETRY ?=
PYTHON_FILES := $(shell find . -type f -name '*.py' -not -path './.venv/*' -not -path './build/*' | sort)

.DEFAULT_GOAL := help

.PHONY: help venv install convert check check-shebang check-literals check-exceptions test verify verify-all verify-examples verify-deep verify-waves verify-pf-passive verify-pf-active verify-isolver verify-core-profiles verify-equilibrium verify-dd audit inspect plot-ids fortran clean

help:
	@echo "TRANSP-to-IMAS developer commands"
	@echo
	@echo "  make venv         Create a local Python virtual environment"
	@echo "  make install      Install requirements into the local environment"
	@echo "  make convert      Convert CDF into OUTPUT_DIR (fresh HDF5 entry by default)"
	@echo "  make test         Verify and convert both NSTX regression cases"
	@echo "  make check        Enforce Python headers, compile-check Python, and run fast regressions"
	@echo "  make verify       Run source diagnostics for CASE"
	@echo "  make verify-all   Verify EVERY bundled CDF automatically (all physics families)"
	@echo "  make verify-examples  Alias for verify-all (10 current CDF examples)"
	@echo "  make verify-deep  Run additional legacy numeric diagnostics for every CDF"
	@echo "  make verify-waves Validate EC/IC waves mappings on ALL bundled CDF files"
	@echo "  make verify-core-profiles Check temperature peaking/time grids in ALL CDF files"
	@echo "  make verify-equilibrium   Check equilibrium TIME/TIME3, Fourier LCFS boundaries in ALL CDFs"
	@echo "  make verify-inputs Validate TR.DAT/model namelist provenance and XML"
	@echo "  make export-inputs Preview TR.DAT and model XML in build/input-xml"
	@echo "  make verify-pf-passive Test PS_n/geometry mapping without an IMAS runtime"
	@echo "  make verify-pf-active  Test KK/VV/VS coil/supply mappings without IMAS"
	@echo "  make verify-isolver    Check real ITER/NSTX ISolver PF-current mappings (local large fixtures)"
	@echo "  make verify-dd    Check supported IMAS Data Dictionary layouts"
	@echo "  make audit        Write CDF mapping-coverage reports for CASE"
	@echo "  make inspect      Plot the final six q profiles from CDF"
	@echo "  make plot-ids     Print usage for the converted-IDS plotting example"
	@echo "  make fortran      Build the optional Fortran Access Layer verifier"
	@echo "  make clean        Remove generated local build/example artifacts"
	@echo
	@echo "Overrides: CASE, CASES, DATA_ROOT, CDF, TRDAT, OUTPUT_DIR, BACKEND, DD_VERSION, PF_PASSIVE_GEOMETRY, PF_ACTIVE_GEOMETRY, PYTHON"

venv:
	$(PYTHON) -m venv "$(VENV)"
	$(VENV_PYTHON) -m pip install --upgrade pip setuptools wheel

install: venv
	$(VENV_PYTHON) -m pip install -r requirements.txt

convert:
	@test -f "$(CDF)" || { echo "CDF not found: $(CDF)" >&2; exit 1; }
	@mkdir -p "$(OUTPUT_DIR)"
	$(PYTHON) convert.py -i "$(CDF)" -b "$(BACKEND)" -u "$(OUTPUT_DIR)" \
		--dd-version "$(DD_VERSION)" $(if $(PF_PASSIVE_GEOMETRY),--pf-passive-geometry "$(PF_PASSIVE_GEOMETRY)") \
		$(if $(PF_ACTIVE_GEOMETRY),--pf-active-geometry "$(PF_ACTIVE_GEOMETRY)") $(CONVERT_FLAGS)

check-shebang:
	@bad=0; \
	for file in $(PYTHON_FILES); do \
		first=$$(head -n 1 "$$file"); \
		if [ "$$first" != "#!/usr/bin/python" ]; then \
			echo "Invalid Python header: $$file: $$first" >&2; \
			bad=1; \
		fi; \
	done; \
	if [ $$bad -ne 0 ]; then exit 1; fi
	@echo "Python header policy OK: $$(printf '%s\n' $(PYTHON_FILES) | wc -l) files"

check-literals:
	$(PYTHON) examples/verify_python_literals.py

check-exceptions:
	$(PYTHON) examples/verify_exception_policy.py

check: check-shebang check-literals check-exceptions
	$(PYTHON) -m compileall -q -x '(^|/)(\.venv|build)/' .
	$(MAKE) -C examples quick PYTHON="$(PYTHON)"

test: check verify-all
	@set -e; for case in $(CASES); do \
		echo "Converting NSTX case $$case"; \
		$(MAKE) --no-print-directory convert CASE="$$case"; \
	done

verify:
	$(MAKE) -C examples verify PYTHON="$(PYTHON)" CDF="$(CDF)" TRDAT="$(TRDAT)"

verify-all verify-examples:
	$(MAKE) -C examples verify-all PYTHON="$(PYTHON)"

verify-deep:
	$(MAKE) -C examples verify-deep PYTHON="$(PYTHON)"

verify-waves:
	$(MAKE) -C examples verify-waves PYTHON="$(PYTHON)"

verify-inputs:
	$(MAKE) -C examples verify-inputs PYTHON="$(PYTHON)"

export-inputs:
	$(PYTHON) examples/export_input_xml.py "$(TRDAT)" --runid "$(CASE)" --output-dir build/input-xml

verify-pf-passive:
	$(MAKE) -C examples verify-pf-passive PYTHON="$(PYTHON)"

verify-pf-active:
	$(MAKE) -C examples verify-pf-active PYTHON="$(PYTHON)"

verify-isolver:
	$(MAKE) -C examples verify-isolver PYTHON="$(PYTHON)"

verify-core-profiles:
	$(MAKE) -C examples verify-core-profiles PYTHON="$(PYTHON)"

verify-equilibrium:
	$(MAKE) -C examples verify-equilibrium PYTHON="$(PYTHON)"

verify-dd:
	$(MAKE) -C examples verify-dd PYTHON="$(PYTHON)"

audit:
	$(MAKE) -C examples audit PYTHON="$(PYTHON)" CDF="$(CDF)" CASE="$(CASE)"

inspect:
	$(MAKE) -C examples inspect PYTHON="$(PYTHON)" CDF="$(CDF)"

plot-ids:
	$(MAKE) -C examples plot-ids PYTHON="$(PYTHON)"

fortran:
	$(MAKE) -C examples fortran

clean:
	$(RM) -r __pycache__ build
	$(MAKE) -C examples clean
