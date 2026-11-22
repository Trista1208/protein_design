# Kemp eliminase pipeline - one command per milestone.
#
#   make m0        audit + controls        (CPU, ~2h)
#   make m1        RFD3 smoke + VRAM       (GPU, ~2h)
#   make m2        LigandMPNN + Boltz      (GPU, ~2h)
#   make m3        scoring + gates         (CPU, ~3h)
#   make m4        scale                   (GPU, overnight)  [gated on m0]
#   make m5        select final 5          (CPU, ~2h)
#   make all       m0 -> m5, gated
#   make verify    run every tool for real
#   make vram      measure peak VRAM per GPU tool
#   make test      unit tests
#   make clean     remove run outputs (keeps inputs + weights)
#
# Every target is idempotent and logs to logs/.

SHELL      := /bin/bash
ROOT       := $(shell pwd)
ENVS       := /home/mandalorian/anaconda3/envs
PY_MM      := $(ENVS)/mm/bin/python
PY_FOUNDRY := $(ENVS)/foundry/bin/python
LOGS       := $(ROOT)/logs
STAMP      := $(shell date +%Y%m%d_%H%M%S)

export KEMP_ROOT := $(ROOT)

.PHONY: help all m0 m1 m2 m3 m4 m5 verify vram test lint clean envs weights doctor

help:
	@sed -n '2,24p' $(MAKEFILE_LIST) | sed 's/^# \{0,1\}//'

# ---------------------------------------------------------------- milestones
m0:
	@mkdir -p $(LOGS)
	$(PY_MM) -m kemp m0 2>&1 | tee $(LOGS)/m0_$(STAMP).log

m1:
	@mkdir -p $(LOGS)
	$(PY_FOUNDRY) -m kemp m1 2>&1 | tee $(LOGS)/m1_$(STAMP).log

m2:
	@mkdir -p $(LOGS)
	$(PY_MM) -m kemp m2 2>&1 | tee $(LOGS)/m2_$(STAMP).log

m3:
	@mkdir -p $(LOGS)
	$(PY_MM) -m kemp m3 2>&1 | tee $(LOGS)/m3_$(STAMP).log

# m4 refuses unless m0's control pair ranked correctly. That refusal is in code.
m4:
	@mkdir -p $(LOGS)
	$(PY_FOUNDRY) -m kemp m4 2>&1 | tee $(LOGS)/m4_$(STAMP).log

m5:
	@mkdir -p $(LOGS)
	$(PY_MM) -m kemp m5 2>&1 | tee $(LOGS)/m5_$(STAMP).log

all:
	@mkdir -p $(LOGS)
	@set -o pipefail; $(PY_MM) -m kemp all 2>&1 | tee $(LOGS)/all_$(STAMP).log

# ---------------------------------------------------------------- checks
doctor:
	@echo "== GPU =="; nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader
	@echo "== envs =="; for e in mm foundry aizymes chai1 her2_boltz2 qm viz q6; do \
	   if [ -x $(ENVS)/$$e/bin/python ] || [ -d $(ENVS)/$$e ]; then \
	     printf "  %-14s OK\n" $$e; else printf "  %-14s MISSING\n" $$e; fi; done
	@echo "== opt =="; for d in rosetta orca openmpi-4.1.8 Q6; do \
	   if [ -d $$HOME/opt/$$d ]; then printf "  %-16s OK\n" $$d; \
	   else printf "  %-16s MISSING\n" $$d; fi; done
	@echo "== weights =="
	@test -f $$HOME/.foundry/checkpoints/rfd3_latest.ckpt && echo "  rfd3  OK" || echo "  rfd3  MISSING"
	@ls $$HOME/.foundry/checkpoints/rf3_*.ckpt >/dev/null 2>&1 && echo "  rf3   OK" || echo "  rf3   MISSING"
	@test -d $$HOME/.boltz && echo "  boltz OK" || echo "  boltz MISSING"
	@ls tools/LigandMPNN/model_params/*.pt >/dev/null 2>&1 && echo "  lmpnn OK" || echo "  lmpnn MISSING"

verify:
	@mkdir -p $(LOGS)
	bash scripts/verify_tools.sh 2>&1 | tee $(LOGS)/verify_$(STAMP).log

vram:
	@mkdir -p $(LOGS)
	bash scripts/measure_vram.sh 2>&1 | tee $(LOGS)/vram_$(STAMP).log

test:
	$(PY_MM) -m pytest tests/ -v

lint:
	$(PY_MM) -m compileall -q kemp && echo "kemp/ compiles clean"

# ---------------------------------------------------------------- housekeeping
envs:
	@mkdir -p envs
	@for e in mm foundry aizymes chai1 qm viz q6 her2_boltz2; do \
	   /home/mandalorian/anaconda3/bin/conda env export -n $$e > envs/$$e.yml 2>/dev/null \
	     && echo "exported envs/$$e.yml" || echo "skip $$e"; done
	@$(ENVS)/foundry/bin/pip freeze > envs/pip-foundry.txt 2>/dev/null || true
	@$(ENVS)/aizymes/bin/pip freeze > envs/pip-aizymes.txt 2>/dev/null || true
	@$(ENVS)/mm/bin/pip freeze      > envs/pip-mm.txt      2>/dev/null || true
	@$(ENVS)/chai1/bin/pip freeze   > envs/pip-chai1.txt   2>/dev/null || true

weights:
	$(ENVS)/foundry/bin/foundry install rfd3 rf3
	$(ENVS)/foundry/bin/foundry list-installed
	bash tools/LigandMPNN/get_model_params.sh tools/LigandMPNN/model_params

clean:
	rm -rf outputs/runs outputs/verify outputs/verify_gen* outputs/vram outputs/state.json
	@echo "removed run outputs; inputs, tools and weights kept"
