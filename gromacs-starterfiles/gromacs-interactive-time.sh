###############################################################################
#                   GROMACS SYSTEM GENERATION – INTERACTIVE SCRIPT            #
###############################################################################
#
# Author: Satvik Kotha
# Version: 2.0
# Citation: 
#
# Requirements:
#   • MAESTRO (Schrödinger) – for protein and ligand preparation
#   • GROMACS – for MD simulation (install GPU version; requires NVIDIA drivers
#     and CUDA toolkit)
#   • AMBERTOOLS + ACPYPE – for generating ligand/protein parameter files
#
# Notes:
#   • Install AMBERTOOLS and ACPYPE inside a dedicated conda environment
#   • Do NOT install GROMACS inside the conda environment
#
# Recommended Conda Environment Setup:
#   1. conda create -n md python=3.9
#   2. conda activate md
#   3. conda install -c conda-forge ambertools
#   4. pip install acpype
#
# Initial Preparation and Running:
#   • Prepare the protein–ligand complex in Maestro
#   • Export the protein and ligand separately as PDB files
#   • Remove all CONNECT lines from the ligand PDB file (protein PDB is unchanged)
#   • Note the charge of the ligand and the ligandname 
#     The ligandname is a 3-letter word in the ligand.pdb file eg LIG or UNK
#   • Place the ligand.pdb and protein.pdb in the same folder as the script 
#     and mdp files
#   • Run this script using sh gromacs-interactive-time.sh
#
###############################################################################

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
gmx pdb2gmx -f protein.pdb -ff amber99sb-ildn -water spc -o protein.gro -ignh
wait
echo PROTEIN CONVERSION COMPLETE SUCCESSFULLY

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

## Simulation
# Stage: EM
gmx grompp -f em.mdp -c solv_ions.gro -p complex.top -o em.tpr
wait
gmx mdrun -v -deffnm em -nb gpu -gpu_id 0
gmx make_ndx -f em.gro -o index.ndx <<EOF
1 | 13
q
EOF
# Stage: NVT
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p complex.top -n index.ndx -o nvt.tpr
wait
gmx mdrun -v -deffnm nvt -nb gpu -gpu_id 0
wait
# Stage: NPT
gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p complex.top -n index.ndx -o npt.tpr
wait
gmx mdrun -v -deffnm npt -nb gpu -gpu_id 0
wait
# Stage: Production MD
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p complex.top -n index.ndx -o md-run-1.tpr
wait
gmx mdrun -v -deffnm md-run-1 -nb gpu -gpu_id 0
