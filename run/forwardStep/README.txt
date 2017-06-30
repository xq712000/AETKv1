solver: allSpeedUnsteadyFoam_dualtime

Parallel to accelerate the computing process and the corresponding order is bellow:
decomposePar
mpirun -np 4 allSpeedUnsteadyFoam_dualtime -parallel | tee log
reconstructPar



