# allSpeedUnsteadyFoam README


[TOC]


# SUMMARY
This is the code of this paper: [**Implementation of density-based solver for all speeds in the framework of OpenFOAM**](http://www.sciencedirect.com/science/article/pii/S0010465514002136)

The original version downloaded from [HERE](http://cpc.cs.qub.ac.uk/summaries/AETK_v1_0.html) can be debug or run based on OpenFOAM-2.0.1.

This README.md was modified from original README.txt to markdown format by Di Cheng.

My intention to to try to transplant it to OpenFOAM 4.1

---

# License

   allSpeedUnsteadyFoam is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or (at your
   option) any later version.

-------------------------------------------------------------------------------

# System requirements

   allSpeedUnsteadyFoam is developed and tested on Linux (ubuntu 10.04
   or later). Memory and disk space requirements depend on the mesh size of the
   specific problem. No additional requirements are recomended for computer that
   support a generic home linux.

-------------------------------------------------------------------------------

# OpenFOAM installation

- Skip this step if you have a functionally version of OpenFOAM-2.0.1 already installed on your pc.
  - Download the source code http://www.openfoam.org/archive/2.0.1/download/source.htm
  - Unpacking the Sources
  - The user should choose a directory location to unpack these files, which will become the installation directory of OpenFOAM. I established a folder named 'OpenFOAM' in the directionary '/home/username/' and put the source pack in this folder' /home/username/OpenFOAM'. 
  - After the installation directory is chosen (and, if necessary, created), simply copy the 2 source pack files into the directory and unpack using `tar xzf <filename>`, e.g. from the installation directory:

```shell
tar xzf OpenFOAM-2.0.1.gtgz
tar xzf ThirdParty-2.0.1.gtgz
# or right click on the icon and choose 'Extract Here'
# The files unpack to produce directories OpenFOAM-2.0.1 and ThirdParty-2.0.1.
```

## System Requirements: Ubuntu 10.04 (lucid)

- several packages need to be installed. This is a command line that will install most of the necessary (and some possibly unnecessary) packages for building OpenFOAM and ParaView:

```shell
sudo apt-get install binutils-dev flex bison git-core build-essential python-dev libreadline5-dev 
wget zlib1g-dev cmake libpng12-dev libxt-dev libxi-dev libxrender-dev libxrandr-dev libxcursor-dev \
libxinerama-dev libfreetype6-dev libfontconfig1-dev libglib2.0-dev freeglut3-dev libqt4-dev qt4-dev-tools
```

## User Configuration

In order to use the installed OpenFOAM package, complete the following 

Open the `.bashrc` file in the user's home directory in an editor, e.g. by typing in a terminal window (note the dot) 

```shell
gedit ~/.bashrc
```

At the bottom of that file, add the following line and save the file 

```shell
export FOAM_INST_DIR=/home/username/OpenFOAM
foamDotFile=$FOAM_INST_DIR/OpenFOAM-1.7.0/etc/bashrc
[ -f $foamDotFile ] && . $foamDotFile"
```
typing at the terminal
```shell
. $HOME/.bashrc
```

## Building the Sources

```shell
cd $WM_PROJECT_DIR
./Allwmake
```
Of course, Paraview should be Compiled according the corresponding instruction 
(http://www.openfoam.org/archive/2.0.1/download/source.html). 

---
# allSpeedUnsteadFoam and allSpeedUnsteadyFoam_dualtime installation

allSpeedUnsteadFoam is a collect of additional library for OpenFOAM.

Two solvers are included, **allSpeedUnsteadFoam** and **allSpeedUnsteadyFoam_dualtime**. The dual time step scheme is implemented in the second solver **allSpeedUnsteadyFoam_dualtime** which is based on the **allSpeedUnsteadFoam**.

## Extracting the source
Extract the package in temporary directory(i.e.  /home/username/test/) by typing in a terminal
```shell
tar -zxvf allSpeedUnsteadFoam.tar.gz
#or right click on the icon and choose 'Extract Here')
```
The tree-view of the files should look like this
```
	src
	|
	|
	|____allSpeedUnsteadyFoam
	|		|
	|       |____allSpeedUnsteadyFoam.C    <   main file   >
	|		|
	|       |____createFields.H 	<Create fields for the fluid domain   >
	|		|
	|       |____readTimeControls.H		<Read the control parameters   >
	|		|
	|       |____readMultiStage.H		<Read the control parameters of Runge-Kutta scheme   >
	|		|
	|       |____readFieldBounds.H		< Read field bounds   >
	|		|
	|       |____solveUnsteadyFluid.H		<Solve governing equations of the fluid   >
	|		|
	|       |____readFieldBounds.H		< Read field bounds   >
	|		|
	|       |____Make
	|				|
	|       		|____files     < Compiling Information about name of main file >
	|				|
	|       		|____options  <Compiling information about location of header files and basic library >
	|
	|____allSpeedUnsteadyFoam_dualtime
	|		|
	|       |____allSpeedUnsteadyFoam_dualtime.C    <   main file   >
	|		|
	|       |____createFields.H 	<Create fields for the fluid domain   >
	|		|
	|       |____solveUnsteadyFluid.H		<Solve governing equations of the fluid   >
	|		|
	|       |____createDualTimeSteppingFields.H		<   Create fields for dual time scheme    >
	|		|
	|       |____updateDualTimeSteppingFields.H		<   Update fields for dual time scheme    >
	|		|
	|       |____Make
	|				|
	|       		|____files     < Compiling Information about name of main file >
	|				|
	|       		|____options  <Compiling information about location of header files and basic library >
	|
	|____timeStepping
	|		|
	|       |____EulerLocalDdtScheme
	|		|		|
	|       |		|____EulerLocalDdtScheme.H  <Header file about definition of euler time step >
	|		|		|
	|       |		|____EulerLocalDdtScheme.C  < Concrete definition of euler time step >
	|		|		|
	|       |		|____EulerLocalDdtSchemes.C	 < Macro definition >
	|		|
	|       |____localTimeStep
	|		|		|
	|       |		|____localTimeStep.H	<Header file about definition of local time step(sub cycling) >
	|		|		|
	|       |		|____localTimeStep.C	< Concrete definition of local time step(sub cycling) >
	|		|
	|       |____temperatureDirectedInletOutletVelocity
	|		|		|
	|       |		|____temperatureDirectedInletOutletVelocityFvPatchVectorField.H	
	|       |		|			< Header file about definition of boundary temperatureDirectedInletOutletVelocity>
	|		|		|
	|       |		|____temperatureDirectedInletOutletVelocityFvPatchVectorField.C	
	|       |					< Concrete definition of boundary temperatureDirectedInletOutletVelocity>
	|		|
	| 		|____temperatureDirectedInletVelocity
	|		|		|
	|       |		|____temperatureDirectedInletVelocityFvPatchVectorField.H	
	|       |		|			< Header file about definition of boundary  temperatureDirectedInletVelocity>
	|		|		|
	|       |		|____temperatureDirectedInletVelocityFvPatchVectorField.C	
	|		|					< Concrete definition of boundary  temperatureDirectedInletVelocity>
	|		|
	| 		|____isentropicTotalTemperature		
	|		|		|
	|       |		|____isentropicTotalTemperatureFvPatchScalarField.H		
	|       |		|			< Header file about definition of boundary isentropicTotalTemperature>
	|		|		|
	|       |		|____isentropicTotalTemperatureFvPatchScalarField.C		
	|       |					< Concrete definition of boundary isentropicTotalTemperature>
	|		|
	| 		|____MRF
	|		|		|
	|       |		|____MRFZone.H  <Header file about definition of MRF (multi reference frame) zone >
	|		|		|
	|       |		|____MRFZone.C  <Concrete definition of definition of MRF (multi reference frame) zone >
	|		|		|
	|       |		|____MRFZoneTemplates.C <templates definition of MRFZone >
	|		|		|
	|       |		|____MRFZones.H < Header file of Container class for a set of MRF Zones>
	|		|		|
	|       |		|____MRFZones.C < Concrete definition of Container class for a set of MRF Zones>
	|
	|____godunovFlux
	|		|
	| 		|____godunovFlux.H <Header file about definition of face flux corresponding to Godunov scheme>
	|		|
	| 		|____godunovFlux.C<Concrete definition of face flux corresponding to Godunov scheme>
	|		|
	| 		|____AUSMplusPreFlux
	|		|
	| 		|____AUSMplusPreFlux.H <Header file about definition of preconditioned AUSM+ scheme>
	|		|
	| 		|____AUSMplusPreFlux.C<Concrete definition of preconditioned AUSM+ scheme>
	|		|
	| 		|____SlopeLimiter <folder about definition of limited function definition, including minmod, 
	|		|
	| 		|____vanleer, etc., and definition of multi-dimension limited function>
	|		|
	| 		|____Make
	|				|
	|       		|____files     <  Compiling Information about name of main file   >
	|				|
	|       		|____options  < Compiling information about location of header files and basic library  >
	|
	|____inv
	|		|
	| 		|____1RINV.H <Header file of definition of matrix inversion>
	|		|
	| 		|____1RINV.C <Concrete definition of matrix inversion>
	|		|
	| 		|___Make
	|				|
	|       		|____files     < Compiling Information about name of main file   >
	|				|
	|       		|____options  < Compiling information about location of header files and basic library  >
	run
	|
	|___cavity(test case - allSpeedUnsteadyFoam)
	|
	|___bump_Minf0675(test case - allSpeedUnsteadyFoam)	
	|
	|___forwardStep(test case -  allSpeedUnsteadyFoam_dualtime)
```



### Compiling the source

```shell
#update and compile `finiteVolume` library.
export temp=/home/username/mysolver
export OF=/home/username/OpenFOAM/OpenFOAM-2.0.1
cp $temp/timeStepping/MRF/* $OF/src/finiteVolume/cfdTools/general/MRF/
cd $OF/src/finiteVolume
wmake libso

cd $temp/allSpeedUnsteadyFoam/src/timeStepping
wmake libso

cd $temp/allSpeedUnsteadyFoam/src/godunov
wmake libso

cd $temp/allSpeedUnsteadyFoam/src/inv
wmake libso

cd $temp/allSpeedUnsteadyFoam/src/allSpeedUnsteadyFoam
wmake

cd $temp/allSpeedUnsteadyFoam/src/allSpeedUnsteadyFoam_dualtime
wmake
```
-------------------------------------------------------------------------------
#  case testing	

```shell
#Enter the directionary "temp/allSpeedUnsteadyFoam" in a terminal
cd $temp/allSpeedUnsteadFoam/run/cavity
#run the solver "allSpeedUnsteadyFoam" in a single processor, typing in the terminal
allSpeedUnsteadyFoam

#run the solver "allSpeedUnsteadyFoam" in two processor parallel, typing in the terminal

decomposePar
mpirun -np 2 allSpeedSteadyFoam - parallel
```
