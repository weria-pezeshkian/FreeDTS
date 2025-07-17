import os
import random
import numpy as np
import shutil

# Define your arrays
core=np.array([5.049,6.016, 6.999, 7.999,8.496,8.694, 8.999,9.164]) #core radius at that frame
TSI=np.array([5,112,240,390,472,506,560,590]) #names of frames to equilibrate
kappa=np.array([10.0]) #membrane bending rigidity
MCsteps=3000000 #Monte Carlo steps
nrep=1 #number of replica per simulation
current_directory = os.getcwd()

# Loop over the combinations of arrays
for i in range(0,len(kappa)):
    for j in range(0,len(TSI)):
        # Create directory
        directory_name = "Kappa{}TSI{}".format(kappa[i],TSI[j])
        os.makedirs(directory_name)

        # Create 3 rep directories
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
VolumeCoupling = SecondOrder 0  0.0000  0.15\n\
TotalAreaCoupling = HarmonicPotential 1000 0.37\n\
AlexanderMove = MetropolisAlgorithmOpenMP 1\n\
Boundary = EllipsoidalCore  {}  1  1  1\n\
VisualizationFormat = VTUFileFormat VTU_F 2000\n\
NonbinaryTrajectory = TSI TrajTSI 2000\n\
Restart_Period = 1000\n\
INCLUSION\n\
Define 2 Inclusions\n\
SRotation Type   K    KG  KP  KL  C0    C0P C0L  lambda   lkg   lkn    cn0\n\
0      Pro1      10   0   10   5   0.0    1   -1\n\
0      Pro2      20   0   0   0   -0.4   0   0\n\
Inclusion-Inclusion-Int\n\
1    1       1  2   0.0     0.0".format(MCsteps,kappa[i],core[j]))
                
                
                
            output_tsi_source ="dts{}.tsi".format(TSI[j])
            output_tsi_dest = os.path.join(rep_directory, "dts{}.tsi".format(TSI[j]))
            shutil.copy(output_tsi_source, output_tsi_dest)
                

            # Create and write to run.sh, a run file for an HPC environment; modify the paths accordingly!
            with open(os.path.join(rep_directory, "run.sh"), "w") as run_file:
                run_file.write("#!/bin/bash\n")
                run_file.write("/lustre/astro/beatrice/FreeDTSv2.0-master/DTS -in input.dts -top dts{}.tsi -seed {}\n".format(TSI[j],seed))

            os.chmod(os.path.join(rep_directory, "run.sh"), 0o755)
            with open(os.path.join(current_directory, "directory.txt"), "a") as dir_file:
                dir_file.write("{}\n".format(rep_directory))

print("Script execution completed.")
