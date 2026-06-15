set complex_pdb [lindex $argv 0]    ;# input protein-ligand complex pdb path
set output_dir [lindex $argv 1]     ;# output dir: data/predictions/id/simulations/{boltz,smba,smba_af3_pocket}/  (toppar shared one level up)
set ligand_id [lindex $argv 2]      ;# e.g., "10194105"

# Ligand atoms are already named to match toppar/lig.rtf (generate_structures.py),
# so no per-source renaming is needed. The resname comes from the rtf's RESI.
proc rtf_resi {rtf} {
    set f [open $rtf r]
    while {[gets $f line] >= 0} {
        if {[regexp {^RESI\s+(\S+)} $line -> name]} { close $f; return $name }
    }
    close $f
    error "no RESI line in $rtf"
}

puts "\n=== Building system (build_system.tcl)==="
puts "  Complex: $complex_pdb"
puts "  Output: $output_dir"
puts "  Ligand: $ligand_id"


# Helper function for atom renaming
proc fix_boltz_atom_names {selection} {

    # Ensure selection is valid
    if {[catch {$selection num}]} {
        error "Invalid atom selection passed to fix_boltz_atom_names"
    }

    # Get indices and elements
    set atom_indices [$selection get index]
    set atom_elements [$selection get element]

    # Initialize counters
    set counters [dict create O 1 N 1 C 1 H 1]

    # Loop over each atom
    for {set i 0} {$i < [llength $atom_indices]} {incr i} {
        set idx [lindex $atom_indices $i]
        set elem [lindex $atom_elements $i]

        if {[dict exists $counters $elem]} {
            set count [dict get $counters $elem]
            set newname "${elem}$count"
            dict incr counters $elem
        }

        set a [atomselect top "index $idx"]
        $a set name $newname
        $a delete
    }
}

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

set temp_PROTEIN "${output_dir}/PROTEIN.pdb"
set temp_LIG "${output_dir}/LIG.pdb"
set coordspsf "${output_dir}/COORDS.psf"
set coordspdb "${output_dir}/COORDS.pdb"

# Load complex
mol new $complex_pdb type pdb

# Process protein (chain P in all complex sources). chop off any UNK residues
set selProtein [atomselect top "chain P and not resname UNK"]
$selProtein set segname PROTEIN
$selProtein writepdb $temp_PROTEIN

segment PROTEIN { pdb $temp_PROTEIN }
coordpdb $temp_PROTEIN PROTEIN
regenerate angles dihedrals
guesscoord

# Process ligand (chain L in all complex sources), heavy atoms only
set selLIG [atomselect top "chain L and not element H"]
if {[$selLIG num] <= 0} {
    puts "ERROR: No ligand found"
    exit 1
}

topology ${output_dir}/toppar/lig.rtf

# resname must match the rtf's RESI exactly (the complex ligand is already named)
set shortid [rtf_resi ${output_dir}/toppar/lig.rtf]

$selLIG set segname LIG
$selLIG set resname $shortid
$selLIG writepdb $temp_LIG
segment LIG {
    pdb $temp_LIG
    # i am telling psf gen not to generate dihedrals or angles. ligand 6085 got some random ahh non-existtent dihedreal bc of this so thats th eonly one skipped
    auto none
}
coordpdb $temp_LIG LIG
if {$ligand_id ne "6085"} {
    regenerate dihedrals
}
regenerate angles
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
file delete -force -- [file normalize $temp_PROTEIN] [file normalize $temp_LIG] [file normalize $coordspsf] [file normalize $coordspdb]
file delete -force -- [file normalize "${output_dir}/SOLVATED.pdb"] [file normalize "${output_dir}/SOLVATED.psf"]
file delete -force -- [file normalize "${output_dir}/SOLVATED.log"]
mol delete all

puts "Finished building system"
exit


