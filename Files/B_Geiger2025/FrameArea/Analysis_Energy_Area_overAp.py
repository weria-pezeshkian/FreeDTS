#To analyse bending energy & surface area, save them as .npy array files, and give quick overview plots
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

# Frames/States that were simulated at their constant frame area A_p:
TSI = np.array([1,6,15,22,40,55,60,68,71,76,79,81,83,86,91,94,96,98,101,103,106,108,110,112,114,117,120,121,124,126,130,133,136,138,141,143,145,148,150])
kappa=np.array([5,10,20]) #bending rigidities
nrep=10 #number of replica ran
ALLenergies=[] #to be filled with TSI lists, each contains 3 kappa list, of which each contains ten tuples [E,errorE]
ALLareas=[] #to store surface area the same
for i in range(0,len(TSI)):
    TSIenergies,TSIareas=[],[]
    print("TSI{}".format(TSI[i]))
    for j in range(0,len(kappa)):
        Kappaenergies,Kappaareas=[],[]
        print("Kappa{}".format(kappa[j]))
        for k in range(1,11):
            file_path="{}TSI_Kappa{}/rep{}/dts-en.xvg".format(TSI[i],kappa[j],k)
            OneRepenergies,OneRepareas=[],[]
            with open(file_path, 'r') as file:
                data=file.readlines()
                for line in data[10000:]:
                    if line.startswith(' ##'):# Skip lines that start with ##
                        continue
                    columns = line.split()# Split the line into columns
                    OneRepenergies.append(float(columns[1]))# Append energy values to a list to then calc. the average & error
                    OneRepareas.append(float(columns[2]))
            OneRepenergies,Onerepareas=np.array(OneRepenergies),np.array(OneRepareas)
            n_blocks = 20 #for block averaging
            block_size=500 #len(OneRepenergies)/n_blocks
            blocksE,blocksA= [OneRepenergies[i*block_size:(i+1)*block_size] for i in range(n_blocks)],[OneRepareas[i*block_size:(i+1)*block_size] for i in range(n_blocks)]
            block_averagesE,block_averagesA = np.array([np.mean(block) for block in blocksE]),np.array([np.mean(block) for block in blocksA])
            overall_averageE,overall_averageA = np.mean(block_averagesE),np.mean(block_averagesA)
            block_varianceE,block_varianceA = np.var(block_averagesE, ddof=1),np.var(block_averagesA, ddof=1)  # Sample variance
            error_estimateE,error_estimateA = np.sqrt(block_varianceE/n_blocks),np.sqrt(block_varianceA/n_blocks)
            #E_rep= average of the block averages
            E_rep,A_rep=overall_averageE,overall_averageA
            #error=sqrt(variance(block average variances))
            E_rep_error,A_rep_error=error_estimateE,error_estimateA 
            Kappaenergies.append([E_rep,E_rep_error])
            Kappaareas.append([A_rep, A_rep_error])
    
        TSIenergies.append(Kappaenergies)
        TSIareas.append(Kappaareas)
    ALLenergies.append(TSIenergies)
    ALLareas.append(TSIareas)

#saving as numpy files to easily load again
ALLenergies,ALLareas = np.array(ALLenergies), np.array(ALLareas)
np.save('ALLenergiesGen1AconstALL-BlockAvg-size500.npy', ALLenergies)
np.save('ALLareasGen1AconstALL-BlockAvg-size500.npy',ALLareas)

# Load ALLenergies from the .npy file to confirm it's saved and loaded correctly
#ALLenergies_loaded = np.load('ALLenergies.npy')

#quick plotting for an overview 
#(E_B over box side length, not scaled by surface area A_0 yet, but should be done when making good plots!):
boxsizes=[33.95,34.11,34.42,34.64,34.83,34.92,35.01,35.12,35.23,35.34,35.45,35.56,35.69,35.79,35.98,36.09,36.24,36.37,36.51,
36.64,36.80,36.93,37.09,37.25,37.38,37.63,37.83,38.04,38.25,38.57,38.89,39.19,39.47,39.74,40.08,40.31,40.58,40.95,41.24]

#just E_B
fig, ax = plt.subplots()
colors=["r","g","b"] #for kappa=5,10,20
for j in range(3):
    for i in range(10):
        ax.errorbar(boxsizes, np.array(ALLenergies)[:, j, i, 0], yerr=np.array(ALLenergies)[:, j, i, 1], fmt="o",markersize=2,capsize=3,linewidth=1,color=colors[j], alpha=0.5)
ax.set_xlabel(r'$boxsize$ in $[l_{DTS}]$',fontsize=14)
ax.set_ylabel(r'$E_{B}$',fontsize=14)
ax.set_title('Bending Energy over Boxsize,const. area,\n red:$\kappa=5k_BT$, green:$\kappa=10k_BT$, blue:$\kappa=20k_BT$')
ax.grid(True)
ax.minorticks_on()
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.grid(which='minor', linestyle=':', linewidth='0.5', color='gray')
plt.savefig('BlockAveragingEB_over_boxGen1AconstALL.png', format='png', dpi=200)


#with stretching energy due to area constraint
gamma_TSI=0.3768
N_T=3368 #No. of triangles, can be found in tsi files too
A0=N_T*(np.sqrt(3)/4)*(1+2*gamma_TSI) #target surface area of the constraint
# stretching energy will be: E_A=(1000/N_T)*(A-A0)**2 with errorE_A=2*(1000/N_T)*(A-A0)*errorA
fig, ax = plt.subplots()
colors=["r","g","b"] #for kappa=5,10,20
for j in range(3):
    for i in range(10):
        ax.errorbar(boxsizes,np.array(ALLenergies)[:, j, i, 0]+(1000/N_T)*(np.array(ALLareas)[:, j, i, 0]-A0)**2, yerr=np.array(ALLenergies)[:, j, i, 1]+(2000/N_T)*(np.array(ALLareas)[:, j, i, 0]-A0)*np.array(ALLareas)[:, j, i, 1], fmt="o",markersize=2,capsize=3,linewidth=1,color=colors[j], alpha=0.5)

ax.set_xlabel(r'$boxsize$ in $[l_{DTS}]$',fontsize=14)
ax.set_ylabel(r'$E_{total}=E_B+E_A$',fontsize=14)
ax.set_title('Energy over Boxsize,const. area,\n red:$\kappa=5k_BT$, green:$\kappa=10k_BT$, blue:$\kappa=20k_BT$')
ax.grid(True)
ax.minorticks_on()
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.grid(which='minor', linestyle=':', linewidth='0.5', color='gray')
plt.savefig('BlockAveragingEB_stretch_over_boxGen1AconstALL.png', format='png', dpi=200)

#just stretching
fig, ax = plt.subplots()
for j in range(3):
    for i in range(10):
        ax.errorbar(boxsizes,(1000/N_T)*(np.array(ALLareas)[:, j, i, 0]-A0)**2, yerr=(2000/N_T)*(np.array(ALLareas)[:, j, i, 0]-A0)*np.array(ALLareas)[:, j, i, 1], fmt="o",markersize=2,capsize=3,linewidth=1,color=colors[j], alpha=0.5)

ax.set_xlabel(r'$boxsize$ in $[l_{DTS}]$',fontsize=14)
ax.set_ylabel(r'$E_A$',fontsize=14)
ax.set_title('Stretching Energy over Boxsize,const. area,\n red:$\kappa=5k_BT$, green:$\kappa=10k_BT$, blue:$\kappa=20k_BT$')
ax.grid(True)
ax.minorticks_on()
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.grid(which='minor', linestyle=':', linewidth='0.5', color='gray')
plt.savefig('BlockAveraging_Estretch_over_boxGen1AconstALL.png', format='png', dpi=200)