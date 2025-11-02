from pdbfixer import PDBFixer
from openmm.app import PDBFile
import os, time, subprocess, json
from cosolvkit.analysis import Report
from modSimulationResume import run_simulation
import shutil
import pymol

# Input parameters
pdb_name = "7ANQ_chainA-C-terminal-minimized.pdb"
md_duration = 300
base_path = os.getcwd()

# RESUME CONTROL - Set to True to skip completed simulations
RESUME_MODE = True
# Manually skip specific cosolvents if needed
SKIP_LIST = []

# Fix protein structure
fixed_protein = "fixedProtein.pdb"

# First fix to remove heterogens
fixer = PDBFixer(filename=pdb_name)
fixer.findMissingResidues()
fixer.findNonstandardResidues()
fixer.replaceNonstandardResidues()
fixer.removeHeterogens(True)
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(7.4)
PDBFile.writeFile(fixer.topology, fixer.positions, open("amber.pdb", 'w'))

# Second fix to prepare final structure
fixer = PDBFixer(filename=pdb_name)
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(7.4)
PDBFile.writeFile(fixer.topology, fixer.positions, open(fixed_protein, 'w'))

# Count residues function
def count_residues(pdb_file):
    residues = set()
    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith(('ATOM', 'HETATM')):
                residue_name = line[17:20].strip()
                residue_number = line[22:26].strip()
                residue_id = (residue_name, residue_number)
                residues.add(residue_id)
    return len(residues)

# Get residue count
num_residues = count_residues(fixed_protein)
print(f"Number of unique residues: {num_residues}")

# Function to check if simulation is complete
def is_simulation_complete(job_name):
    """Check if a simulation has been completed"""
    # Check multiple indicators of completion
    checks = [
        os.path.join(base_path, job_name, job_name, "analysis", "report.html"),
        os.path.join(base_path, job_name, job_name, "analysis", f"density_visualization.pse"),
        os.path.join(base_path, job_name, job_name, "clean.xtc"),
    ]
    return any(os.path.exists(check) for check in checks)

# Define all available cosolvents
all_cosolvents = [
    {
        "name": "6-terrazine",
        "concentration": 0.5,
        "smiles": "CNc1nncnn1",
        "resname": "SAT"
    },
    # {
    #     "name": "HingeBinder2",
    #     "concentration": 0.5,
    #     "smiles": "CNC(=O)c1ncc(nn1)NC",
    #     "resname": "HII"
    # },
    # {
    #     "name": "PyridineNCH3",
    #     "concentration": 0.5,
    #     "smiles": "CC(C)CN(C)c1ncccn1",
    #     "resname": "PYC"
    # },
    # Add other cosolvents here...
]

# Write forcefields.json
forcefields_data = {
    "OPENMM": [
        "amber14-all.xml", 
        "amber14/tip3pfb.xml"
    ],
    "AMBER": [
        "amber14-all.xml", 
        "amber14/tip3pfb.xml"
    ],
    "GROMACS": [
        "amber14-all.xml", 
        "amber14/tip3pfb.xml"
    ],
    "CHARMM": [
        "charmm36.xml",
        "charmm36/water.xml"
    ],
    "small_molecules": [
        "espaloma"
    ]
}

with open('forcefields.json', 'w') as json_file:
    json.dump(forcefields_data, json_file, indent=4)

# Function to run simulation for a single cosolvent
def run_cosolvent_simulation(cosolvent, repResname):
    cosolvent_name = cosolvent["name"]
    job_name = f"{cosolvent_name}-withRepulsion-{md_duration}ns"
    
    # Check if we should skip this simulation
    if RESUME_MODE and is_simulation_complete(job_name):
        print(f"\n{'='*80}\n✓ Skipping {cosolvent_name} - already completed\n{'='*80}")
        return "skipped"
    
    if cosolvent_name in SKIP_LIST:
        print(f"\n{'='*80}\n✗ Skipping {cosolvent_name} - in SKIP_LIST\n{'='*80}")
        return "skipped"
    
    print(f"\n\n{'='*80}\nRunning simulation for {cosolvent_name}\n{'='*80}")
    
    try:
        # Create a directory for this cosolvent
        if not os.path.exists(job_name):
            os.makedirs(job_name)
        
        # Change to cosolvent directory
        os.chdir(job_name)
        
        # Copy the fixed protein file
        shutil.copy(f"{base_path}/{fixed_protein}", fixed_protein)
        
        # Create cosolvents.json with just this cosolvent
        with open('cosolvents.json', 'w') as json_file:
            json.dump([cosolvent], json_file, indent=4)
        
        # Copy forcefields.json
        shutil.copy(f"{base_path}/forcefields.json", "forcefields.json")
        
        # Create config.json
        config_data = {
          "cosolvents": "cosolvents.json",
          "forcefields": "forcefields.json",
          "md_format": "openmm",
          "receptor": True,
          "protein_path": fixed_protein,
          "clean_protein": True,
          "keep_heterogens": False,
          "variants_d": {},
          "add_repulsive": False,
          "repulsive_residues": [repResname],
          "epsilon": 0.05,
          "sigma": 9.0,
          "solvent_smiles": "H2O",
          "solvent_copies": None,
          "membrane": False,
          "lipid_type": "POPC",
          "lipid_patch_path": None,
          "cosolvent_placement": 0,
          "waters_to_keep": [],
          "radius": None,
          "output": job_name,
          "run_cosolvent_system": True,
          "run_md": False
        }
        
        with open('config.json', 'w') as json_file:
            json.dump(config_data, json_file, indent=4)
        
        # Create cosolvent system (skip if already exists)
        if not os.path.exists(f"{job_name}/system.xml"):
            subprocess.run(["create_cosolvent_system", "-c", "config.json"])
        else:
            print("System files already exist, skipping system creation")
        
        # Run MD simulation
        print(f"Running MD simulation for {cosolvent_name}")
        start = time.time()
        run_simulation(
            simulation_format='OPENMM',
            results_path=job_name,
            topology=None,
            positions=None,
            pdb='system.pdb',
            system='system.xml',
            membrane_protein=False,
            traj_write_freq=25000,
            time_step=0.004,
            warming_steps=100000,
            simulation_steps=md_duration*250000.0,
            seed=42,
            resume=True  # This is default, so you don't even need to specify it
        )
        print(f"Simulation finished after {(time.time() - start)/60:.2f} min.")
        
        # Create cpptraj commands
        cpptraj_commands = f"""parm {os.getcwd()}/{job_name}/system.prmtop
            trajin {os.getcwd()}/{job_name}/trajectory.dcd
            center :1-{num_residues}@CA
            image
            reference {os.getcwd()}/{job_name}/system.pdb [myref]
            rms ref [myref] :1-{num_residues}@CA out protein.rmsd
            trajout {os.getcwd()}/{job_name}/clean.xtc
            run
            quit
            """
        
        # Create results directory
        results_dir = 'results'
        if not os.path.exists(results_dir):
            os.makedirs(results_dir)
        
        # Write cpptraj input
        with open(f'{results_dir}/centreProteinCpptraj.inp', 'w') as file:
            file.write(cpptraj_commands)
        
        # Run cpptraj
        os.chdir(results_dir)
        command = ['cpptraj', '-i', 'centreProteinCpptraj.inp']
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(result.stdout.decode())
        print(result.stderr.decode())
        os.chdir('..')
        
        # Analyze results
        out_path = f'{os.getcwd()}/{job_name}/analysis'
        if not os.path.exists(out_path):
            os.makedirs(out_path)
        
        log_file = f'{os.getcwd()}/{job_name}/statistics.csv'
        traj_file = f'{os.getcwd()}/{job_name}/clean.xtc'
        top_file = f'{os.getcwd()}/{job_name}/system.prmtop'
        cosolvents_file = f'{os.getcwd()}/cosolvents.json'
        
        # Generate analysis reports
        report = Report(log_file, traj_file, top_file, cosolvents_file)
        report.generate_report(out_path=out_path)
        report.generate_density_maps(out_path=out_path, analysis_selection_string="")
        
        # Extract resname
        with open(cosolvents_file, 'r') as file:
            cosolvents_data = json.load(file)
        
        # Generate density files for visualization
        resnames = [cs['resname'] for cs in cosolvents_data]
        denFiles_list = [f'{out_path}/map_density_{resname}.dx' for resname in resnames]
        
        # Clear PyMOL session before generating new reports
        pymol.cmd.reinitialize()
        
        report.generate_pymol_reports(
            topology=top_file,
            trajectory=traj_file,
            density_files=denFiles_list,
            selection_string="",
            out_path=out_path
        )
        
        # Return to base path
        os.chdir(base_path)
        print(f"\n✓ Completed simulation and analysis for {cosolvent_name}\n")
        return "completed"
        
    except Exception as e:
        print(f"\n✗ ERROR in {cosolvent_name}: {str(e)}\n")
        os.chdir(base_path)
        return "failed"

# Run simulations for each cosolvent with status tracking
results_summary = {}
for cosolvent in all_cosolvents:
    resname = cosolvent['resname']
    status = run_cosolvent_simulation(cosolvent, resname)
    results_summary[cosolvent['name']] = status

# Print summary
print("\n" + "="*80)
print("SUMMARY OF ALL SIMULATIONS")
print("="*80)
for name, status in results_summary.items():
    symbol = "✓" if status == "completed" else "○" if status == "skipped" else "✗"
    print(f"{symbol} {name}: {status}")
print("="*80)
