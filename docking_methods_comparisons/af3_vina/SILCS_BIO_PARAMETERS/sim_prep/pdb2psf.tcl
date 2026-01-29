
# Load required packages
package require psfgen
package require solvate
package require autoionize
package require pbctools

proc solvate_and_ionize {psf_file pdb_file output_prefix} {
    # Solvate the system

    puts "running solvate_and_ionize"
    solvate $psf_file $pdb_file -o "${output_prefix}_solvated" -minmax {{-60 -60 -60} {60 60 60}}
    puts "System solvated: ${output_prefix}_solvated"

    # Ionize the system
    autoionize -psf "${output_prefix}_solvated.psf" -pdb "${output_prefix}_solvated.pdb" -sc 0.15 -o "${output_prefix}_system"
    puts "System ionized: ${output_prefix}_system"
}

# Load topology files - this may need to change depending on the system, but I think it should be fine
topology ../toppar/top_all36_prot.rtf
topology ../toppar/top_all36_cgenff.rtf
topology ../toppar/toppar_water_ions.str

# Alias for residues not recognized by CHARMM
#This is unlikely to be correct. You need to 
pdbalias residue HIS HSD ;# Alias HIS to HSD (neutral form with delta protonation)
pdbalias atom ILE CD1 CD ;# Example alias for atoms, adjust if required
# pdbalias residue UNL LIG


# loop through all complexes recursively
set complex_dirs [glob -nocomplain -type d ../complexes/*/*/*]

foreach dir $complex_dirs {
    puts "\nProcessing directory: $dir"

    set cif_files [glob -nocomplain ${dir}/MD_complex.pdb ${dir}/AF3_complex.pdb]

    foreach cif_file $cif_files {
        puts "\n\tProcessing CIF: $cif_file"


        set base [file rootname [file tail $cif_file]]
        set output_prefix "${dir}/${base}_sim"

        set checker "${output_prefix}_system.psf"

        if {[file exists $checker]} {
            puts "✅ Skipping $cif_file — Already processed: ${checker}"
            mol delete all
            continue
        }
        # puts $output_prefixå

        mol new $cif_file type pdb
        set selA [atomselect top "protein"]
        $selA set segname PA
        
        $selA writepdb chainA.pdb

        segment PA {
            pdb chainA.pdb
        }


        coordpdb chainA.pdb PA
        guesscoord

        # Extract base name of CIF file (e.g., AF3_complex -> AF3_complex_sim)


        #water/ions probably not necessary
        set selB [atomselect top "not (protein or water or ions)"] 
        $selB set segname LIG

        set ligand_id [file tail $dir] 
        set lig_param_file "../toppar/${ligand_id}.str"
        
        if {![file exists $lig_param_file]} {
            puts "⚠️  Skipping $cif_file — Missing topology file: $lig_param_file"
            resetpsf
            mol delete all
            continue
        }

        set shortid [string range $ligand_id 0 3]


        topology $lig_param_file
        pdbalias residue UNL $shortid
        pdbalias residue LIG $shortid
        #some file had ISP as the ligand residue??  Human/a0a5c2g238/6309/mdcmplex. might be sus
        pdbalias residue ISP $shortid
        

        $selB writepdb chainB.pdb

        segment LIG {
            pdb chainB.pdb
        }


        coordpdb chainB.pdb LIG
        guesscoord
        
        # break;  
        writepsf "${output_prefix}.psf"
        writepdb "${output_prefix}.pdb"



        mol new "${output_prefix}.psf" type psf
        mol addfile "${output_prefix}.pdb" type pdb waitfor all

        set wat [atomselect top "water"]
        puts "PSF+PDB has [llength [$wat get index]] water atoms"


        solvate_and_ionize "${output_prefix}.psf" "${output_prefix}.pdb" "${output_prefix}"


        mol new "${output_prefix}_system.psf" type psf
        mol addfile "${output_prefix}_system.pdb" type pdb waitfor all

        pbc writexst "${output_prefix}_system.xsc"
        resetpsf
        mol delete all

    }

}
