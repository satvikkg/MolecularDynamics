## Interactive script for gromacs system generation

## Requirements
# MAESTRO for protein and ligand preperation
# GROMACS for simulation
# AMBER for generating ligand parametric files
# ACPYPE for ligand parameter generation

## Initial setup
# Prepare the complex in schrodinger and export both preotein and ligand as pdb files. Remove the CONNECT lines only in the ligand.pdb file

#conda activate md

## Preparing protein and ligand
#initial setup
echo Enter the ligand ID:
read ligname
echo "Enter the net charge of ligand: "
read charge
echo "Enter time to simulate in ns: "
read nstime
pstime=$(($nstime*1000))
intertime=$(($pstime*1000))
steps=$(($intertime/2))
echo $steps

# preparing protein
gmx pdb2gmx -f protein.pdb -ff oplsaa -water spc -o protein.gro -ignh
wait
echo PROTEIN CONVERSION COMPLETE SUCCESSFULLY
# -f is the input protein file in pdb format -ff is the forcefield -water is the water model -o is the output after conversion -ignh ignores the hydrogens atom coordinate files

## Making copies of topology file and gro file for complex generation
cp topol.top complex.top
cp protein.gro complex.gro
wait

## Preparing ligand parametric files
#echo "Enter the net charge of ligand: "
#read charge
grep  HETATM ligand.pdb > ligand-intermediate.pdb
rm ligand.pdb
mv ligand-intermediate.pdb ligand.pdb

antechamber -i ligand.pdb -fi pdb -o ligand.mol2 -fo mol2 -c bcc -nc $charge
wait

acpype -di ligand.mol2 -c bcc -n $charge
wait

cd ligand.acpype
cp ligand_GMX.gro ligand_GMX.top ligand_GMX.itp ../
cd ../
mv ligand_GMX.gro ligand.gro
mv ligand_GMX.itp ligand.itp
mv ligand_GMX.top ligand.top
echo LIGAND PARAMETERS GENERATED SUCCESSFULLY

## Generating parametric files for the complex
#Modifying the complex.top and complex.gro files
#complex.top file
cat complex.top | sed '/forcefield\.itp\"/a\
#include "ligand.itp"
' >| complex2.top

echo "Ligand   1" >> complex2.top
mv complex2.top complex.top
wait

#complex.gro file
natomsp=$(sed '2q;d' protein.gro)
natomsl=$(sed '2q;d' ligand.gro)
totalatoms=$(($natomsp + $natomsl))
echo "The total number of atoms from protein and ligands is: $totalatoms"
sed -i "2 s/.*/ $totalatoms/" complex.gro

# Modifying .mdp files
#echo Enter the ligand ID:
#read ligname
sed -i "s/ligname/$ligname/g" nvt.mdp
sed -i "s/ligname/$ligname/g" npt.mdp
sed -i "s/ligname/$ligname/g" md.mdp
sed -i "s/simulationtime/$steps/g" md.mdp

head -n -1 complex.gro > tmp-complex.gro
grep $ligname ligand.gro >> tmp-complex.gro
tail -n 1 complex.gro >> tmp-complex.gro
mv tmp-complex.gro complex.gro
rm tmp-complex.gro

echo SYSTEM BUILT SUCCESSFULLY

## Buliding the system
gmx editconf -f complex.gro -o newbox.gro -bt cubic -d 1
wait
gmx solvate -cp newbox.gro -cs spc216.gro -p complex.top -o solv.gro
wait
gmx grompp -f ions.mdp -c solv.gro -p complex.top -o ions.tpr
wait
echo SOL | gmx genion -s ions.tpr -o solv_ions.gro -p complex.top -pname NA -nname CL -neutral
wait

gmx grompp -f em.mdp -c solv_ions.gro -p complex.top -o em.tpr
wait
gmx mdrun -v -deffnm em -nb gpu -gpu_id 0

gmx make_ndx -f em.gro -o index.ndx <<EOF
1 | 13
q
EOF

gmx grompp -f nvt.mdp -c em.gro -r em.gro -p complex.top -n index.ndx -o nvt.tpr
wait
gmx mdrun -v -deffnm nvt -nb gpu -gpu_id 0
wait

gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p complex.top -n index.ndx -o npt.tpr
wait
gmx mdrun -v -deffnm npt -nb gpu -gpu_id 0
wait

gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p complex.top -n index.ndx -o md-run-1.tpr
wait
gmx mdrun -v -deffnm md-run-1 -nb gpu -gpu_id 0
