# Get arguments from command line
set complex_pdb [lindex $argv 0]    ;# input: complex.pdb path
set output_dir [lindex $argv 1]     ;# output: build directory
set ligand_id [lindex $argv 2]      ;# e.g., "7930"
set sim_type [lindex $argv 4]     ;# "pdb_params" or "mol2_params"

puts "\n=== Building system (build_system.tcl)==="
puts "  Complex: $complex_pdb"
puts "  Output: $output_dir"
puts "  Ligand: $ligand_id"
puts "  Sim type: $sim_type"

proc solvate_and_ionize {psf_file pdb_file save_dir} {
    puts "\n\tRunning solvate_and_ionize"

    solvate $psf_file $pdb_file -o "${save_dir}/SOLVATED" -minmax {{-60 -60 -60} {60 60 60}}
    puts "\n\tSystem solvated."

    autoionize -psf "${save_dir}/SOLVATED.psf" -pdb "${save_dir}/SOLVATED.pdb" -sc 0.15 -o "${save_dir}/system"
    puts "\n\tSystem ionized."
}


package require psfgen
package require solvate
package require autoionize
package require pbctools


psfcontext reset
resetpsf


topology ${output_dir}/toppar/top_all36_prot.rtf
topology ${output_dir}/toppar/top_all36_cgenff.rtf

pdbalias residue HIS HSD
pdbalias atom ILE CD1 CD

# Set up temp files
set temp_PROTEIN "${output_dir}/PROTEIN.pdb"
set temp_LIG "${output_dir}/LIG.pdb"
set coordspsf "${output_dir}/COORDS.psf"
set coordspdb "${output_dir}/COORDS.pdb"

# Load complex
mol new $complex_pdb type pdb

# Process protein
set selProtein [atomselect top "protein"]
$selProtein set segname PROTEIN
$selProtein writepdb $temp_PROTEIN

segment PROTEIN { pdb $temp_PROTEIN }
coordpdb $temp_PROTEIN PROTEIN
regenerate angles dihedrals
guesscoord

# Process ligand
set selLIG [atomselect top "not (protein or water or ions)"]
if {[$selLIG num] <= 0} {
    puts "ERROR: No ligand found"
    exit 1
}

topology ${output_dir}/toppar/lig.rtf

# remake to capital bcos im used to that in PDB file
set shortid "L[string range $ligand_id 0 1]"

$selLIG set segname LIG
$selLIG set resname $shortid
$selLIG writepdb $temp_LIG

segment LIG {
    pdb $temp_LIG
    # i am telling psf gen not to generate dihedrals or angles. ligand 6085 got non-existtent dihedreal bc of this so thats th eonly one skipped
    auto none
}
coordpdb $temp_LIG LIG
regenerate angles dihedrals
guesscoord

# Write initial psf/pdb
writepsf $coordspsf
writepdb $coordspdb

# Solvate and ionize
mol new $coordspsf type psf
mol new $coordspdb type pdb
solvate_and_ionize $coordspsf $coordspdb $output_dir

# Write XSC file
mol new "${output_dir}/system.psf" type psf
mol new "${output_dir}/system.pdb" type pdb
pbc writexst "${output_dir}/system.xsc"

# Cleanup temp files
file delete -force $temp_PROTEIN $temp_LIG $coordspsf $coordspdb
file delete -force "${output_dir}/SOLVATED.pdb" "${output_dir}/SOLVATED.psf"
file delete -force "${output_dir}/SOLVATED.log"
mol delete all

puts "(COMPLETION) Finished building system"
exit


