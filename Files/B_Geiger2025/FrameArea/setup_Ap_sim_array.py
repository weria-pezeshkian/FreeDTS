import os
import random
import numpy as np
import shutil

# Define your arrays
# Frames/States to be simulated at their constant frame area A_p:
TSI = np.array([1,6,15,22,40,55,60,68,71,76,79,81,83,86,91,94,96,98,101,103,106,108,110,112,114,117,120,121,124,126,130,133,136,138,141,143,145,148,150]) 
#for extra data at transition points: fourth to twelvth TSI w. more replica:
#TSI = np.array([22,40,55,60,68,71,76,79,81]) 
gamma_TSI=np.array([0.3768])# For area constraint
kappa=np.array([5,10,20]) # bending rigidities
MCsteps=2000000 #simulation steps
nrep=10 #number of replica to run
current_directory = os.getcwd()

# Loop over the combinations of arrays
for i in range(0,len(TSI)):
    for j in range(0,len(kappa)):
    # Create directory:
        directory_name = "{}TSI_Kappa{}".format(TSI[i],kappa[j])
        os.makedirs(directory_name)

    # Create rep directories
        for rep in range(1, nrep+1):
            rep_directory = os.path.join(directory_name, "rep{}".format(rep))
            os.makedirs(rep_directory)

            # Generate a random seed number
            seed = random.randint(1, 100000)

            # Create and write to input.dts
            with open(os.path.join(rep_directory, "input.dts"), "w") as input_file:
                input_file.write("Integrator_Type = MC_Simulation\n\
Min_Max_Lenghts = 1 3\n\
MinfaceAngle = -0.5\n\
Temperature = 1 0\n\
Box_Centering_F = 0\n\
Set_Steps = 1 {}\n\
EnergyMethod = FreeDTS1.0_FF\n\
Kappa = {}  0 0\n\
Edge_Parameters = 5 0 0\n\
VertexArea = 0 0.7 0 0\n\
TimeSeriesData_Period = 100\n\
VertexPositionIntegrator = MetropolisAlgorithmOpenMP 1 1 0.05\n\
TotalAreaCoupling = HarmonicPotential 1000 {}\n\
AlexanderMove = MetropolisAlgorithmOpenMP 1\n\
VisualizationFormat = VTUFileFormat VTU_F 5000\n\
NonbinaryTrajectory = TSI TrajTSI 5000\n\
Restart_Period = 10000".format(MCsteps,kappa[j],gamma_TSI[0]))
                
                
            output_tsi_source ="output{}.tsi".format(TSI[i])
            output_tsi_dest = os.path.join(rep_directory, "output{}.tsi".format(TSI[i]))
            shutil.copy(output_tsi_source, output_tsi_dest)

            # Create and write to run.sh, run file for an HPC environment; modify paths accordingly!
            with open(os.path.join(rep_directory, "run.sh"), "w") as run_file:
                run_file.write("#!/bin/bash\n")
                run_file.write("/lustre/astro/beatrice/FreeDTSv2.0-master/DTS -in input.dts -top output{}.tsi -seed {}\n".format(TSI[i],seed))

            os.chmod(os.path.join(rep_directory, "run.sh"), 0o755)
            with open(os.path.join(current_directory, "directorySplit.txt"), "a") as dir_file:
                dir_file.write("{}\n".format(rep_directory))

print("Script execution completed.")
