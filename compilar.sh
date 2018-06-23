
echo "Compilando..."
gfortran -o mein param.f90 utils.f90 linalg.f90 mesh.f90 iterativo.f90 rwg.f90 mom.f90 mlfma.f90 rcs.f90 main.f90 -fcheck=all
