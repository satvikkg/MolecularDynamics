import os, shutil

simulation_duration = 1

if simulation_duration < 5:
   interval_duration = 4.8
else:
   interval_duration = simulation_duration

schrodinger_path = '/opt/schrodinger2023-1'

mae_folder_path = '/'
output_path = '/'

if os.path.exists(output_path) != True:
    print(f'Output Path Does Not Exist. Creating new path at {output_path}.')
    os.mkdir(output_path)

desmond_queu_batch = f'{output_path}/desmond_queu_batch.sh'

for mae_file in os.listdir(mae_folder_path):
    mae_name = os.path.splitext(mae_file)[0]
    desmond_folder = f'{output_path}/{mae_name}'
    desmond_setup = f'{desmond_folder}/desmond_setup_{mae_name}'
    desmond_md_job = f'{desmond_folder}/desmond_md_job_{mae_name}'

    if os.path.exists(desmond_folder) != True:
        os.mkdir(desmond_folder)
    else:
        print(f'{desmond_folder}: Exists')
    
    if os.path.exists(desmond_setup) != True:
        os.mkdir(desmond_setup)
    else:
        print(f'{desmond_setup}: Exists')

    if os.path.exists(desmond_md_job) != True:
        os.mkdir(desmond_md_job)
    else:
        print(f'{desmond_md_job}: Exists')

# Creating files to run desmond system builder:
    '''
    Files required to run desmond_system_builder:
    1) .mae file of interest.
    2) .msj file containing the options for system builder.
    3) .sh file to run the desmond_system_builder job.
    '''
    # Copying the .mae/maegz file:
    shutil.copyfile(f'{mae_folder_path}/{mae_file}', f'{desmond_setup}/{mae_file}')

    # Creating the .msj file:
    desmond_setup_msj = [
        'task {task =  "desmond:auto"}\n'
        '\n'
        'build_geometry {add_counterion = {ion = Na number = neutralize_system} box = {shape = orthorhombic size = [8.0 8.0 8.0 ] size_type = buffer} override_forcefield = S-OPLS rezero_system = false salt = {concentration = 0.15 negative_ion = Cl positive_ion = Na} solvent = TIP3P}\n'
        '\n'
        'assign_forcefield {forcefield = S-OPLS water = TIP3P}\n'
    ]

    desmond_setup_msj_file = f'{desmond_setup}/{mae_name}.msj'
    f = open(desmond_setup_msj_file, 'w')
    f.writelines(desmond_setup_msj)
    f.close()

# Creating files to run desmond molecular dynamics:
    '''
    1) .cms file from the system setup.
    2) .msj file containing the options for molecular dynamics.
    3) .cfg file containing the configuratuion options for molecular dynamics.
    4) .sh file to run the desmond_md_job.
    '''

    # Creating the .cfg file:
    desmond_md_job_cfg = [
    '''
        annealing = false
        backend = {
        }
        bigger_rclone = false
        box = ?
        checkpt = {
        first = 0.0
        interval = 240.06
        name = "$JOBNAME.cpt"
        write_last_step = true
        }
        cpu = 1
        cutoff_radius = 9.0
        dipole_moment = false
        ebias_force = false
        elapsed_time = 0.0
        energy_group = false
        eneseq = {
        first = 0.0
        interval = 1.2
        name = "$JOBNAME$[_replica$REPLICA$].ene"
        }
        ensemble = {
        barostat = {
            tau = 2.0
        }
        class = NPT
        method = MTK
        thermostat = {
            tau = 1.0
        }
        }
        gaussian_force = false
        glue = solute
        maeff_output = {
        center_atoms = solute
        first = 0.0
        interval = 120.0
        name = "$JOBNAME$[_replica$REPLICA$]-out.cms"
        periodicfix = true
        trjdir = "$JOBNAME$[_replica$REPLICA$]_trj"
        }
        meta = false
        meta_file = ?
        msd = false
        pressure = [1.01325 isotropic ]
        pressure_tensor = false
        randomize_velocity = {
        first = 0.0
        interval = inf
        seed = 2007
        temperature = "@*.temperature"
        }
        restrain = none
        restraints = {
        existing = ignore
        new = []
        }
        rnemd = false
        simbox = {
        first = 0.0
        interval = 1.2
        name = "$JOBNAME$[_replica$REPLICA$]_simbox.dat"
        }
        spatial_temperature = false
        surface_tension = 0.0
        taper = false
        temperature = [
        [300.0 0 ]
        ]
    '''
            f'time = {simulation_duration}'
    '''
        timestep = [0.002 0.002 0.006 ]
        trajectory = {
        center = []
        first = 0.0
        format = dtr
        frames_per_file = 250
    '''
            f'interval = {interval_duration}'
    '''
        name = "$JOBNAME$[_replica$REPLICA$]_trj"
        periodicfix = true
        write_last_vel = false
        write_velocity = false
        }
        wall_force = false
    '''
    ]

    desmond_md_job_cfg_file = f'{desmond_md_job}/desmond.cfg'
    f = open(desmond_md_job_cfg_file, 'w')
    f.writelines(desmond_md_job_cfg)
    f.close()

    # Creating the .msj file:
    desmond_md_job_msj = [
        '''
            task {
            task = "desmond:auto"
            set_family = {
                desmond = {
                    checkpt.write_last_step = no
                }
            }
            }

            simulate {
            title       = "Brownian Dynamics NVT, T = 10 K, small timesteps, and restraints on solute heavy atoms, 100ps"
            annealing   = off
            time        = 100
            timestep    = [0.001 0.001 0.003 ]
            temperature = 10.0
            ensemble = {
                class = "NVT"
                method = "Brownie"
                brownie = {
                    delta_max = 0.1
                }
            }
            restraints.new = [{
                name = posre_harm
                atoms = solute_heavy_atom
                force_constants = 50.0
            }]
            }

            simulate {
            title       = "NVT, T = 10 K, small timesteps, and restraints on solute heavy atoms, 12ps"
            annealing   = off
            time        = 12
            timestep    = [0.001 0.001 0.003]
            temperature = 10.0
            restraints.new = [{
                name = posre_harm
                atoms = solute_heavy_atom
                force_constants = 50.0
            }]
            ensemble = {
                class  = NVT
                method = Langevin
                thermostat.tau = 0.1
            }

            randomize_velocity.interval = 1.0
            eneseq.interval             = 0.3
            trajectory.center           = []
            }

            simulate {
            title       = "NPT, T = 10 K, and restraints on solute heavy atoms, 12ps"
            annealing   = off
            time        = 12
            temperature = 10.0
            restraints.existing = retain
            ensemble    = {
                class  = NPT
                method = Langevin
                thermostat.tau = 0.1
                barostat  .tau = 50.0
            }

            randomize_velocity.interval = 1.0
            eneseq.interval             = 0.3
            trajectory.center           = []
            }

            simulate {
            title       = "NPT and restraints on solute heavy atoms, 12ps"
            effect_if   = [["@*.*.annealing"] 'annealing = off temperature = "@*.*.temperature[0][0]"']
            time        = 12
            restraints.existing = retain
            ensemble    = {
                class  = NPT
                method = Langevin
                thermostat.tau = 0.1
                barostat  .tau = 50.0
            }

            randomize_velocity.interval = 1.0
            eneseq.interval             = 0.3
            trajectory.center           = []
            }

            simulate {
            title       = "NPT and no restraints, 24ps"
            effect_if   = [["@*.*.annealing"] 'annealing = off temperature = "@*.*.temperature[0][0]"']
            time        = 24
            ensemble    = {
                class  = NPT
                method = Langevin
                thermostat.tau = 0.1
                barostat  .tau = 2.0
            }

            eneseq.interval   = 0.3
            trajectory.center = solute
            }

            simulate {
            cfg_file = "desmond.cfg"
            jobname  = "$MAINJOBNAME"
            dir      = "."
            compress = ""
            }
        '''
    ]

    desmond_md_job_msj_file = f'{desmond_md_job}/desmond.msj'
    f = open(desmond_md_job_msj_file, 'w')
    f.writelines(desmond_md_job_msj)
    f.close()

# Creating the desmond_queu_batch.sh file:
    f = open(desmond_queu_batch, 'a')
    f.write(f'cd {desmond_setup}\n')
    f.write(f'echo "Preparing system for: {mae_name}"\n')
    f.write(f'"{schrodinger_path}/utilities/multisim" -JOBNAME desmond_setup_{mae_name} -m {desmond_setup_msj_file} {mae_file} -o desmond_setup_{mae_name}-out.cms -HOST localhost -TMPLAUNCHDIR -WAIT\n')
    f.write(f'echo "System setup complete for: {mae_name}"\n')
    f.write(f'cp {desmond_setup}/desmond_setup_{mae_name}-out.cms {desmond_md_job}/desmond_setup_{mae_name}-out.cms\n')
    f.write(f'cd {desmond_md_job}\n')
    f.write(f'echo "Starting MD simulation: {mae_name} for {simulation_duration}ns."\n')
    f.write(f'''"{schrodinger_path}/utilities/multisim" -JOBNAME desmond_md_job_{mae_name} -HOST localhost -maxjob 1 -cpu 1 -m desmond.msj -c desmond.cfg -description 'Molecular Dynamics' desmond_setup_{mae_name}-out.cms -mode umbrella -o desmond_md_job_{mae_name}-md-out.cms -lic DESMOND_GPGPU:16 -WAIT\n''')
    f.write(f'echo "MD simulation complete: {mae_name}"\n')
    f.write('\n')
    f.close()