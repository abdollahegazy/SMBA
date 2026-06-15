# Build NAMD systems for the 3 docking sources, per protein.
#
# Output layout (one shared toppar per protein, used by all 3 sims):
#   data/predictions/<id>/simulations/
#       toppar/                  <- CHARMM36 base + lig.rtf / lig.prm (one set)
#       boltz/                   <- boltz-predicted complex
#       smba/                    <- best vina pose complex
#       smba_af3_pocket/         <- best vina pose complex (af3 pocket redock)
#
# Only proteins with a docking_af3_pocket dir are simulated (Ki/Kd rows whose
# affinity converts to kcal/mol). Each protein has one ligand.
#
# One CGenFF param set per ligand (ligand_topos_params/{lig}/l{NN}.rtf+.prm) serves all
# three sims -- generate_structures.py names every pose to that single rtf. A sim
# is skipped if its ligand isn't parameterized yet. Run generate_structures.py
# (builds the complexes incl. the boltz one) before make.

.PHONY: all list clean-done

DATA          := ../../data/predictions
TOPPAR_BASE   := /serviceberry/tank/abdolla/SMBA/docking_methods_comparisons/vina_boltz/parameter_testing/toppar_from_af3_vina
LIG_TOPO_BASE := ligand_topos_params
VMD           := vmd

# --- discovery -------------------------------------------------------------
# proteins with an af3-pocket redock (the ones we simulate)
AF3_IDS := $(notdir $(patsubst %/docking_af3_pocket,%,$(wildcard $(DATA)/*/docking_af3_pocket)))

# single ligand id for a protein (subdir under docking/docking)
lig = $(notdir $(patsubst %/,%,$(wildcard $(DATA)/$1/docking/docking/*/)))

# complex pdb paths (all produced by generate_structures.py)   $1=id  $2=ligand
boltz_pdb = $(DATA)/$1/boltz/$2/$2_best_complex.pdb
smba_pdb  = $(DATA)/$1/docking/docking/$2/$2_best_complex.pdb
af3_pdb   = $(DATA)/$1/docking_af3_pocket/docking_af3_pocket/$2/$2_best_complex.pdb

# --- per-(protein, sim) rule generation ------------------------------------
# $1=id  $2=sim  $3=complex_pdb  $4=ligand
define sim_rule
ALL_DONE += $(DATA)/$1/simulations/$2/.done
$(DATA)/$1/simulations/$2/.done: ID  := $1
$(DATA)/$1/simulations/$2/.done: SIM := $2
$(DATA)/$1/simulations/$2/.done: LIG := $4
$(DATA)/$1/simulations/$2/.done: $3
	$$(build_one)
endef

# recipe (relies on target-specific ID/SIM/LIG; $< = complex pdb)
define build_one
@OUT="$(DATA)/$(ID)/simulations/$(SIM)"; \
TOP="$(DATA)/$(ID)/simulations/toppar"; \
LRTF=$$(ls $(LIG_TOPO_BASE)/$(LIG)/*.rtf 2>/dev/null | head -1); \
LPRM=$$(ls $(LIG_TOPO_BASE)/$(LIG)/*.prm 2>/dev/null | head -1); \
if [ -z "$$LRTF" ] || [ -z "$$LPRM" ]; then \
	echo "[skip] $(ID)/$(SIM): no ligand_topos params for $(LIG)"; \
	exit 0; \
fi; \
mkdir -p "$$OUT" "$$TOP"; \
if [ ! -f "$$TOP/.base_copied" ]; then \
	cp -r "$(TOPPAR_BASE)/." "$$TOP/" && touch "$$TOP/.base_copied"; \
fi; \
cp "$$LRTF" "$$TOP/lig.rtf"; \
cp "$$LPRM" "$$TOP/lig.prm"; \
ln -sfn ../toppar "$$OUT/toppar"; \
sed 's|__LIGAND_PARAMS__|parameters      toppar/lig.prm|' src/system_eq.namd.template  > "$$OUT/system_eq.namd"; \
sed 's|__LIGAND_PARAMS__|parameters      toppar/lig.prm|' src/system_run.namd.template > "$$OUT/system_run.namd"; \
sed 's/defaultjob/eq_$(ID)_$(SIM)/'  src/eq.sh  > "$$OUT/eq.sh"; \
sed 's/defaultjob/run_$(ID)_$(SIM)/' src/run.sh > "$$OUT/run.sh"; \
$(VMD) -dispdev text -e build_system.tcl -args "$<" "$$OUT" "$(LIG)"; \
if [ -f "$$OUT/system.psf" ] && [ -f "$$OUT/system.pdb" ] && [ -f "$$OUT/system.xsc" ]; then \
	touch "$$OUT/.done"; \
	echo "[ok] $(ID)/$(SIM)  ligand=$(LIG)"; \
else \
	echo "[FAIL] $(ID)/$(SIM): build_system.tcl produced no system.{psf,pdb,xsc}"; \
	rm -f "$$OUT"/PROTEIN.pdb "$$OUT"/LIG.pdb "$$OUT"/COORDS.psf "$$OUT"/COORDS.pdb; \
	exit 1; \
fi
endef

$(foreach id,$(AF3_IDS),\
  $(eval LG := $(call lig,$(id)))\
  $(if $(LG),\
    $(eval $(call sim_rule,$(id),boltz,$(call boltz_pdb,$(id),$(LG)),$(LG)))\
    $(eval $(call sim_rule,$(id),smba,$(call smba_pdb,$(id),$(LG)),$(LG)))\
    $(eval $(call sim_rule,$(id),smba_af3_pocket,$(call af3_pdb,$(id),$(LG)),$(LG)))))

all: $(ALL_DONE)

# show what will build vs. skip (ligand parameterized?)
list:
	@printf "%-12s %-10s %-16s %s\n" PROTEIN LIGAND SIM PARAMS
	@$(foreach id,$(AF3_IDS),$(eval LG := $(call lig,$(id)))$(if $(LG),\
	  $(foreach S,boltz smba smba_af3_pocket,\
	    printf "%-12s %-10s %-16s %s\n" "$(id)" "$(LG)" "$(S)" \
	      "$(if $(wildcard $(LIG_TOPO_BASE)/$(LG)/*.rtf),ok,MISSING)"; )))
	@echo
	@echo "proteins (af3_pocket): $(words $(AF3_IDS))"

clean-done:
	rm -f $(ALL_DONE)
