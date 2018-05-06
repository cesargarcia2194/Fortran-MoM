#!/bin/bash

if [ -z "$1" ]
then
    echo "No especifico la prueba a compilar"

else
	echo $1".f90"

	echo "Compilando "$1

	gfortran -o $1 constants.f90 aus.f90 linalg.f90 mesh.f90 param.f90 utils.f90 iterativo.f90 modRWG.f90 mom.f90 mlfma.f90 $1.f90 -fcheck=all

	echo "Listo, ejecute "$1

fi

