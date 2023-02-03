# FreeDTS

FreeDTS is software to perform computational research on biomembranes at messocpic length-scale. In FreeDTS, a membrane is represented by a dynamically triangulated surface equipped with vertex-based inclusions to integrate the effects of integral and peripheral membrane proteins. Several algorithms are included into the software that allow for simulation of framed membrane with constant tension, vesicles with various fixed volume or constant pressure difference, confined membranes into the fixed region of the space, constant fixed global curvature and application for external forces on regions of the membrane. In addition, the software allows one to turn off the shape evolution of the membrane and only explore inclusions organization. This allows to take realistic membrane shapes obtained from Cryo-ET and obtain heterogeneous organization of biomolecules which can be backmapped to finer simulations models. 

Weria Pezeshkian (weria.pezeshkian@gmail.com; weria.pezeshkian@nbi.ku.dk),
Niels Bohr International Academy, 
Niels Bohr Institute, 
University of Copenhagen, 
Blegdamsvej 17, 2100 Copenhagen, DENMARK


<img width="884" alt="Screenshot 2023-02-03 at 13 50 32" src="https://user-images.githubusercontent.com/47776510/216607774-2a3e9391-5a6d-43ce-998f-40960f38c273.png">


### Compiling the source code
To compile the OPENDTS source code just execute “compile.sh” script.

./compile.sh

This will create three binary files: DTS, CNV, and GEN. The DTS is for running simulations; GEN is for generating triangulated files (TS) and CNV is for converting different file formats.

## For More see User_Manual_Tutorials.pdf
