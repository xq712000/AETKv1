 /*---------------------------------------------------------------------------*\
   =========                 |
   \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
    \\    /   O peration     |
     \\  /    A nd           | Copyright (C) 2004-2010 OpenCFD Ltd.
      \\/     M anipulation  |
 -------------------------------------------------------------------------------
 License
     This file is part of OpenFOAM.
 
     OpenFOAM is free software: you can redistribute it and/or modify it
     under the terms of the GNU General Public License as published by
     the Free Software Foundation, either version 3 of the License, or
     (at your option) any later version.
 
     OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
     ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
     FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
     for more details.
 
     You should have received a copy of the GNU General Public License
     along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
 
 Global
     readTimeControls
 
 Description
     Read the control parameters used by setDeltaT
 
 \*---------------------------------------------------------------------------*/
/* 
 const bool adjustTimeStep =
     runTime.controlDict().lookupOrDefault("adjustTimeStep", false);
 
 scalar maxCo =
     runTime.controlDict().lookupOrDefault<scalar>("maxCo", 1.0);
*/ 
// scalar maxDeltaT =
//     runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);
//add fixed local time step
/*
dimensionedScalar  deltaT_Local =
     runTime.controlDict().lookup("deltaT_Local");

scalar  deltaT_Local =
     runTime.controlDict().lookupOrDefault<scalar>("deltaT_Local",0.00000001);
*/
// fixed reference velocity 
/*
dimensionedScalar  velocity_Inlet =
     runTime.controlDict().lookup("velocity_Inlet");
*/
scalar  velocity_Inlet =
     runTime.controlDict().lookupOrDefault<scalar>("velocity_Inlet",1000000);
// scalar convergence_Res =
//     runTime.controlDict().lookupOrDefault<scalar>("convergence_Res", 1e-4);

 
 // ************************************************************************* //
