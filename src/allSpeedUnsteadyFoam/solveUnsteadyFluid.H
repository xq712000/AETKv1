{
// add


        // update the boundary absolute velocity
//        mrfZones.correctBoundaryVelocity(U);

        // set the pseudo face velocity
//        Godunov.dotX() = mrfZones.faceU();
/*
        // Compute cartesian velocity components from cylindrical ones
        coordsX = mesh.C().component(vector::X);
        coordsY = mesh.C().component(vector::Y);

        rhoU.replace(vector::X,(rhoCr*coordsX-rhoCu*coordsY)/radius);
        rhoU.replace(vector::Y,(rhoCr*coordsY+rhoCu*coordsX)/radius);
*/
/*
        rhoUOld.replace(vector::X,(rhoCrOld*coordsX-rhoCuOld*coordsY)/radius);
        rhoUOld.replace(vector::Y,(rhoCrOld*coordsY+rhoCuOld*coordsX)/radius);

        rhoUOldOld.replace(vector::X,(rhoCrOldOld*coordsX-rhoCuOldOld*coordsY)/radius);
*/
        // activate sub-cycle loop
        TimeState subCycleTimeState = runTime.subCycle(numberSubCycles);

        // get access to physical deltaT in ddt scheme!
        physDeltaT[1] = subCycleTimeState.deltaT().value();
        physDeltaT[2] = subCycleTimeState.deltaT0().value();

        // Begin sub-cycling - PseudoTime Integration
        // adjustTimeStep and numberSubCycles > 1 does not make any sense
        for (label j=0; j < numberSubCycles; j++)
        {
// add
			scalar res1;
			scalar res2;
//			scalar total_Res;
//			TOL=0;
			res1=0;
			res2=0;
//			total_Res=0;


	volScalarField Ux=U.component(0);
	volScalarField Uy=U.component(1);
	volScalarField Uz=U.component(2);	
	volScalarField U2=magSqr(U);
  	const volScalarField R_ = ( thermo.Cp() - thermo.Cv() );
	const volScalarField kappa(thermo.Cp()/thermo.Cv());



           // Update subCycle Time
            runTime++;

            // update local time step sizes just once each iteration for
            // all ddt schemes

//            if (!adjustTimeStep)
//            {
//                localTimeStep.update(maxCo,adjustTimeStep,deltaT_Local);
//            }



            // Plot Pseudo Time here, so that one can analyse the residuls with
            // pyFoam*
            Info<< "\n Time = " << runTime.value() << nl << endl;

            // Begin outer Loop for Runge - Kutta
            forAll (beta,i)
            {
                // Update primitive boundaries
                p.correctBoundaryConditions();
                U.correctBoundaryConditions();
                h.correctBoundaryConditions();

                // solve the approximate Riemann problem for this time step
                // reconstruct numerical fluxes at faces in a proper way
                Godunov.update(secondOrder);

                // get access to multi-stage coefficient in ddt scheme
                physDeltaT[0] = beta[i];

                // \f$ \mu \left( \nabla \vec{U} + \left( \nabla \vec{U}
                // \right)^T - \frac{2}{nD} \left( \nabla \bullet \vec{U}
                // \right) I \right) \f$
                // nD = 3 in three dimensional space
                // for two- and one-dimensional computations this
                // implementation is wrong
                // is equal to -turbulence->devRhoReff(), but we do not need to
                // evaluate gradU twice
                const volSymmTensorField tau
                (
                    "tau",
                    -turbulence->devRhoReff()
                    -((2.0/3.0)*I)*rho*turbulence->k()
                );

                volScalarField rhoFlux = -fvc::div(Godunov.rhoFlux());

                volVectorField rhoUFlux = -fvc::div(Godunov.rhoUFlux())
                    + fvc::div(tau);

                volScalarField rhoEFlux = -fvc::div(Godunov.rhoEFlux())
                    // Viscous heating with
                    + fvc::div( tau & U )
                    // Fourier law with static enthalpy
                    // with alphaEff() - the effective turbulent
                    // thermal diffusivity.
                    + fvc::laplacian(turbulence->alphaEff(), h)
                    // Molecular Diffusion and Turbulent Transport closure
                    // Wilcox (2006): Chapter 5.4.3
                    // should be better used DkEff(F1) instead of muEff(), but
                    // this function is not virtual, now it is assumed that
                    // \sigma_k = 5/3 is hard coded
                    + fvc::laplacian
                      (
                          (turbulence->mu()+0.6*turbulence->mut()),
                          turbulence->k()
                      )
                    ;

// add
		volScalarField rhoUFlux_x=rhoUFlux.component(0);
		volScalarField rhoUFlux_y=rhoUFlux.component(1);
		volScalarField rhoUFlux_z=rhoUFlux.component(2);


// Refer to section 2.2 Preconditioning Method in the papaer and Ref. Blazek P 234, eqn.9.67 and eqn.9.69.

// the preconditioning parameter sita
// Refer to Ref Blazek, Computational Fluid Dynamics: Principles and Applications, Elsevier, 2005. P334-P335
// the temporary parameter fai
// Refer to Ref Blazek, Computational Fluid Dynamics: Principles and Applications, Elsevier, 2005. P33, eqn.9.67
// the temporary parameter a2
// Refer to Ref Blazek, Computational Fluid Dynamics: Principles and Applications, Elsevier, 2005. P33, eqn.9.67
// the temporary parameter a3
// Refer to Ref Blazek, Computational Fluid Dynamics: Principles and Applications, Elsevier, 2005. P33, eqn.9.67

	dimensionedScalar velocity_Inlet_=dimensionedScalar("velocity_Inlet_", dimVelocity,velocity_Inlet);	
	Ur2=min
		(
			kappa*R_*T,
			max
				(
					U2, scalar(0.25)*sqr(velocity_Inlet_)
				)
		);

	tmp<volScalarField> sita =scalar(1.0)/Ur2 
							- scalar(1.0)/(kappa*R_*T) 
							+ scalar(1.0)/(R_*T);
	
	tmp<volScalarField> fai = scalar(1.0)/(kappa*R_*T*sita()-(kappa-1));	

	tmp<volScalarField> a2	= (kappa-1)*fai();
	
	tmp<volScalarField> a3 = (kappa-1)*fai()*T/rho;


// the elements of the inverse preconditioning matrix, refer to Ref. Blazek P334, eqn. 9.66
////////////////////////////////////////////////	
		scalarField taoinv11=a2().internalField()*(kappa.internalField()*R_.internalField()*T.internalField()/(kappa.internalField()-1)-(h.internalField()-0.5*U2.internalField()));

		scalarField taoinv12=-a2().internalField()*Ux.internalField();

		scalarField taoinv13=-a2().internalField()*Uy.internalField();

		scalarField taoinv14=-a2().internalField()*Uz.internalField();

		scalarField taoinv15=a2().internalField();


		scalarField taoinv21=-Ux.internalField()/rho.internalField();

		scalarField taoinv22=1.0/rho.internalField();

		scalarField taoinv23(0);

		scalarField taoinv24(0);

		scalarField taoinv25(0);


		scalarField taoinv31=-Uy.internalField()/rho.internalField();

		scalarField taoinv32(0);

		scalarField taoinv33=1.0/rho.internalField();

		scalarField taoinv34(0);

		scalarField taoinv35(0);


		scalarField taoinv41=-Uz.internalField()/rho.internalField();

		scalarField taoinv42(0);

		scalarField taoinv43(0);

		scalarField taoinv44=1.0/rho.internalField();

		scalarField taoinv45(0);


		scalarField taoinv51=a3().internalField()*(1-sita().internalField()*(h.internalField()-0.5*U2.internalField()));

		scalarField taoinv52=-sita().internalField()*a3().internalField()*Ux.internalField();

		scalarField taoinv53=-sita().internalField()*a3().internalField()*Uy.internalField();

		scalarField taoinv54=-sita().internalField()*a3().internalField()*Uz.internalField();

		scalarField taoinv55=sita().internalField()*a3().internalField();
////////////////////////////////////////////////	


// Refer to section 2.3.1., discrete form of source term of, eqn. (4). 
////////////////////////////////////////////////		
		rhoFlux2.internalField()=taoinv11*rhoFlux.internalField()+taoinv12*rhoUFlux_x.internalField()+taoinv13*rhoUFlux_y.internalField()+taoinv14*rhoUFlux_z.internalField()+taoinv15*rhoEFlux.internalField();

		forAll(rhoUFlux2.internalField(),cellI)
		{

			rhoUFlux2.internalField()[cellI].x()=taoinv21[cellI]*rhoFlux.internalField()[cellI]+taoinv22[cellI]*rhoUFlux.internalField()[cellI].x();
			rhoUFlux2.internalField()[cellI].y()=taoinv31[cellI]*rhoFlux.internalField()[cellI]+taoinv33[cellI]*rhoUFlux.internalField()[cellI].y();
			rhoUFlux2.internalField()[cellI].z()=taoinv41[cellI]*rhoFlux.internalField()[cellI]+taoinv44[cellI]*rhoUFlux.internalField()[cellI].z();
		}	

		rhoEFlux2.internalField()=taoinv51*rhoFlux.internalField()+taoinv52*rhoUFlux_x.internalField()+taoinv53*rhoUFlux_y.internalField()+taoinv54*rhoUFlux_z.internalField()+taoinv55*rhoEFlux.internalField();


		volScalarField Cp_rhoEFlux2=thermo.Cp()*rhoEFlux2;
////////////////////////////////////////////////	
                // Add source terms to the numerical fluxes due to the rotation
//                mrfZones.addCoriolis(rho,U,rhoUFlux);

                // Pseudo time integration
                solve
                (

                    fvm::ddt(p) == rhoFlux2

                );

                // time integration
                solve
                (

                    fvm::ddt(U) == rhoUFlux2
                );

                // time integration
                solve
                (

					fvm::ddt(h) == Cp_rhoEFlux2
                );

                // bound density
//                boundMinMax(rho,rhoMin,rhoMax);

                // bound rhoE
//                boundMinMax(rhoE,rhoEMin,rhoEMax);

                // Compute internal field of U
// delete
//                U.dimensionedInternalField() =
//                    rhoU.dimensionedInternalField()
//                   /rho.dimensionedInternalField();

                // Update boundary conditions of U
                U.correctBoundaryConditions();

                // Bound the velocity
                volScalarField magU = mag(U);

                if (max(magU) > UMax)
                {
                    Info<< "bounding " << U.name()
                        << " max: " << max(magU).value()
                        << endl;

                    volScalarField Ulimiter = pos(magU - UMax)*UMax
                        /(magU + smallU) + neg(magU - UMax);
                    Ulimiter.max(scalar(0));
                    Ulimiter.min(scalar(1));

                    U *= Ulimiter;
                    U.correctBoundaryConditions();
                }


                // store tmp fields in order to prevent computing twice
                const volScalarField kappa(thermo.Cp()/thermo.Cv());

//                 // bound static energy
// //                 volScalarField e = rhoE/rho - 0.5*magSqr(U);
//                 volScalarField e = rhoE/rho - 0.5*magSqr(U) - TKE;
// //                 e = rhoE/rho - 0.5*magSqr(U) - turbulence->k();
//                 e.correctBoundaryConditions();
//                 boundMinMax(e,eMin,eMax);

                // Update static enthalpy:
                // The turbulent kinetic energy k is part of the total energy E
                // Therefore it needs to be subtracted from E in order to get
                // the static enthalpy h
//                 h = kappa*e;

                // correct boundary conditions of static enthalpy
                h.correctBoundaryConditions();

                // bound enthalpy
//                boundMinMax(h,hMin,hMax);

                // compute complete field of p
//                 p = (kappa - 1.0)*rho*e;


                // correct boundary conditions of p
                p.correctBoundaryConditions();

                // bound pressure
//                boundMinMax(p,pMin,pMax);

                // correct thermo physical properties
                // therefore one needs correct p and e fields


		
                thermo.correct();

                // Update boundary field of rho
				rho=thermo.rho();
                rho.boundaryField() = thermo.rho()().boundaryField();

                // Update boundary field of rhoU
// add 2013.7.3
				rhoU = rho*U;

                rhoU.boundaryField() = rho.boundaryField()*U.boundaryField();

                // Update boundary field of rhoE

// add 2013.7.3
				rhoE = rho*(h + 0.5*magSqr(U) + turbulence->k()) - p;
                rhoE.boundaryField() =
                    rho.boundaryField()*
                    (
                        0.5*magSqr(U.boundaryField())
                      + turbulence->k()().boundaryField()
                    )
                  + p.boundaryField()/(kappa.boundaryField()-1.0);

                // needed for turbulence and CoEuler ddt scheme
                // and maybe elsewhere;
                // phi is already realtive to the mesh motion, if we are using
                // the ALE approach
                phi = Godunov.rhoFlux();
//                 mrfZones.relativeFlux(fvc::interpolate(rho),phi);

                // Convergence check
                // Averaged L2-Norm of fluxes
                scalar L2NormRho = max(Foam::sqrt(sum(sqr(rhoFlux.internalField()))
                    /mesh.nCells()),SMALL);
//                scalar LInfNormRho = max(mag(rhoFlux.internalField()));
//                 scalar L2NormRhoU = max(Foam::sqrt(sum(magSqr(rhoUFlux.internalField()))
//                     /mesh.nCells()),SMALL);
//                 scalar L2NormRhoE = max(Foam::sqrt(sum(sqr(rhoEFlux.internalField()))
//                     /mesh.nCells()),SMALL);

                // Averaged L2-Norm of fluxes
                Info<< "rho L2 Residual: "
//			<< residual << endl;
                    << Foam::log10(L2NormRho)  << endl
//                    << "rho LInf Residual: "
//                    << Foam::log10(LInfNormRho) << endl
//                     << "rhoU residual: "
//                     << Foam::log10(L2NormRhoU) << endl
//                     << "rhoE residual: "
//                     << Foam::log10(L2NormRhoE) << endl
                    << endl;

            // End outer Loop for Runge - Kutta
            }
            // Update relative velocity
//            mrfZones.relativeVelocity(U,URel);

            //Update turbulence after the multi-stage time integration
            turbulence->correct();



        }
        // End dual-time stepping
        runTime.endSubCycle();

        // Compute cylindrical ones from cartesian velocity components
//        rhoCr = (rhoU.component(vector::X)*coordsX+rhoU.component(vector::Y)*coordsY)/radius;
//        rhoCu = (rhoU.component(vector::Y)*coordsX-rhoU.component(vector::X)*coordsY)/radius;

//#       include "updateDualTimeSteppingFields.H"

}
