package require topotools

set protein_file [lindex $argv 0]
set ligand_file [lindex $argv 1]
set out [lindex $argv 2]

mol new $protein_file
mol new $ligand_file type pdb

set prot [atomselect 0 "all"]
set lig [atomselect 1 "all"]
$lig set chain "L"


set combined [::TopoTools::selections2mol [list $prot $lig]]
animate write pdb $out $combined

exit