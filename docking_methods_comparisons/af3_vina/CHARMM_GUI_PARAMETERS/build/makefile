# Target structure: ./complexes/species/proteinid/ligandid/MD_sim (or AF3_sim)
# Source complex.pdb from:
#   MD_sim:  /fir/tank/abdolla/docking/fir_dock_comp/MD_complexes_c3/species/proteinid/ligandid/complex.pdb
#   AF3_sim: /fir/tank/abdolla/docking/fir_dock_comp/AF_complexes/species/proteinid/ligandid/complex.pdb

.PHONY: test all md_sims af3_sims

LIGANDS := 30999 6085 6309
TOPPAR_BASE := ../toppar_from_af3_vina


test_6085: ../complexes/DouglasFir/c6f797/6085/MD_sim ../complexes/DouglasFir/c6f797/6085/AF3_sim
test_6309: ../complexes/Arabidopsis/a0a7g2f626/6309/MD_sim ../complexes/Arabidopsis/a0a7g2f626/6309/AF3_sim
test_30999: ../complexes/Arabidopsis/a0a178wen7/30999/MD_sim ../complexes/Arabidopsis/a0a178wen7/30999/AF3_sim



# Find existing ligand dirs and strip the ../../complexes prefix
EXISTING_LIGAND_DIRS := $(patsubst ../../complexes/%,%,$(shell find ../../complexes -mindepth 3 -maxdepth 3 -type d))
# Generate targets in ../complexes/ 
MD_TARGETS := $(addprefix ../complexes/,$(addsuffix /MD_sim,$(EXISTING_LIGAND_DIRS)))
AF3_TARGETS := $(addprefix ../complexes/,$(addsuffix /AF3_sim,$(EXISTING_LIGAND_DIRS)))

md_sims: $(MD_TARGETS)
af3_sims: $(AF3_TARGETS)
all: md_sims af3_sims

define get_ligand_prefix
$(shell echo $(1) | cut -c1-2)
endef


#    $(eval LIGAND_PREFIX := $(call get_ligand_prefix,$(LIGAND)))

# Build recipe for MD simulations
define build_md_sim
    $(eval PARTS := $(subst /, ,$*))
    $(eval SPECIES := $(word 1,$(PARTS)))
    $(eval PROTEIN := $(word 2,$(PARTS)))
    $(eval LIGAND := $(word 3,$(PARTS)))
    $(eval SIM_TYPE := $(word 4,$(PARTS)))
    $(eval SOURCE_PDB := /fir/tank/abdolla/docking/fir_dock_comp/MD_complexes_c3/$(SPECIES)/$(PROTEIN)/$(LIGAND)/complex.pdb)
    $(eval LIGAND_PARAM_DIR := $(if $(filter 30999,$(LIGAND)),../toppar/$(LIGAND)_MD/lig,../toppar/$(LIGAND)/lig))
    mkdir -p $@/toppar
    cp -r $(TOPPAR_BASE)/* $@/toppar/
    cp $(LIGAND_PARAM_DIR)/lig.prm $@/toppar/
    cp $(LIGAND_PARAM_DIR)/lig.rtf $@/toppar/
    cp $(SOURCE_PDB) $@/complex.pdb
    @sed 's/defaultjob/md_$(PROTEIN)_$(LIGAND)/' ../src/eq.sh > $@/eq.sh
    @sed 's/defaultjob/md_$(PROTEIN)_$(LIGAND)/' ../src/run.sh > $@/run.sh
    vmd -dispdev text -e build_system.tcl -args $(SOURCE_PDB) $@ $(LIGAND) MD MD_sim

    cp ../src/system_eq.namd.template $@/system_eq.namd
	cp ../src/system_run.namd.template $@/system_run.namd
    
    @sed 's/defaultjob/eq_$(PROTEIN)_$(LIGAND)_$(DOCKING_ORIGIN)_pdb/' ../src/eq.sh > $@/eq.sh
	@sed 's/defaultjob/run_$(PROTEIN)_$(LIGAND)_$(DOCKING_ORIGIN)_pdb/' ../src/run.sh > $@/run.sh
	
    @echo "Built MD simulation: $@ with params from $(LIGAND_PARAM_DIR)"
endef
# Build recipe for AF3 simulations
define build_af3_sim
    $(eval PARTS := $(subst /, ,$*))
    $(eval SPECIES := $(word 1,$(PARTS)))
    $(eval PROTEIN := $(word 2,$(PARTS)))
    $(eval LIGAND := $(word 3,$(PARTS)))
    $(eval SIM_TYPE := $(word 4,$(PARTS)))
    $(eval SOURCE_PDB := /fir/tank/abdolla/docking/fir_dock_comp/AF_complexes/$(SPECIES)/$(PROTEIN)/$(LIGAND)/complex.pdb)
    $(eval LIGAND_PARAM_DIR := $(if $(filter 30999,$(LIGAND)),../toppar/$(LIGAND)_AF3/lig,../toppar/$(LIGAND)/lig))
    mkdir -p $@/toppar
    cp -r $(TOPPAR_BASE)/* $@/toppar/
    cp $(LIGAND_PARAM_DIR)/lig.prm $@/toppar/
    cp $(LIGAND_PARAM_DIR)/lig.rtf $@/toppar/
    cp $(SOURCE_PDB) $@/complex.pdb
    @sed 's/defaultjob/af3_$(PROTEIN)_$(LIGAND)/' ../src/eq.sh > $@/eq.sh
    @sed 's/defaultjob/af3_$(PROTEIN)_$(LIGAND)/' ../src/run.sh > $@/run.sh
    vmd -dispdev text -e build_system.tcl -args $(SOURCE_PDB) $@ $(LIGAND) MD MD_sim


    cp ../src/system_eq.namd.template $@/system_eq.namd
	cp ../src/system_run.namd.template $@/system_run.namd
    
    @sed 's/defaultjob/eq_$(PROTEIN)_$(LIGAND)_$(DOCKING_ORIGIN)_pdb/' ../src/eq.sh > $@/eq.sh
	@sed 's/defaultjob/run_$(PROTEIN)_$(LIGAND)_$(DOCKING_ORIGIN)_pdb/' ../src/run.sh > $@/run.sh
	
    @echo "Built AF3 simulation: $@ with params from $(LIGAND_PARAM_DIR)"
endef

# Rules for MD_sim targets
../complexes/%/MD_sim:
	$(build_md_sim)

# Rules for AF3_sim targets  
../complexes/%/AF3_sim:
	$(build_af3_sim)