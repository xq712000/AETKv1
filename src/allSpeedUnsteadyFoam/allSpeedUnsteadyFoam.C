/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    allSpeedUnsteadyFoam

Description
    Density-based compressible multi-stage  flow solver with preconditioner 
	on the basis densityBasedTurbo solver deceloped by Oliver Borm.
    The multi-stage coefficients and the number of multi-stages are runtime
    selectable. In order to evaluate the convective terms, a Riemann solver is
    utilized.

    
    TODO: Improve the operation of preconditioning matrix. 
		  Implicit method should be developed based on the multi-block lib tool.


    References:
    Jameson, Schmidt and Turkel, "Numerical Solution of the Euler Equations by

	J. Blazek, "Computational Fluid Dynamics: Principles and Applications",
	2nd ed, Elsevier, 2005.


	This is just a version that might work.  For estension or maintainability 
	of this program, the operation of preconditioning matrix should be improved 
	futher according to the high level code within OpenFOAM.

Author
    Chun Shen is copyright owner of the additional code, modified 
	on the basis of densityBasedTurbo solver. Chun Shen are are responsible 
	for these further modifications. 

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
//#include "MRFZone.H"
#include "localEulerDdtScheme.H"
#include "dynamicFvMesh.H"
#include "godunovFlux.H"

//#include "localTimeStep.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createDynamicFvMesh.H"
#   include "createFields.H"

	scalar residual=0;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
#       include "readTimeControls.H"

// this definition of compressibleCourantNo is not correct
// #       include "compressibleCourantNo.H"

// adjusting of physical deltaT is no longer needed with local-time stepping
// #       include "setDeltaT.H"

#       include "readMultiStage.H"
#       include "readFieldBounds.H"

/*
        if (adjustTimeStep)
        {
            localTimeStep.update(maxCo,adjustTimeStep);

            numberSubCycles = 1;
        }
*/

        runTime++;

        // rotate the mesh about the axis
        mesh.update();

#       include "solveUnsteadyFluid.H"
#		include "convergenceCheck.H"
        runTime.write();

        Info<< "\n    ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n" << endl;
    }

    Info<< "\n end \n";

    return(0);
}


// ************************************************************************* //
