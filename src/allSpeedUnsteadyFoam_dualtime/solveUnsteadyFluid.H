{
// add 2013.7.4
	scalar flag;
	flag=0; 
//	scalar res;
	scalar flag1;
	flag1=0;

    const volScalarField& rhoOld =
        mesh.lookupObject<volScalarField>("rhoOld");
    const volScalarField& rhoOldOld =
        mesh.lookupObject<volScalarField>("rhoOldOld");

    const volVectorField& rhoUOld =
        mesh.lookupObject<volVectorField>("rhoUOld");
    const volVectorField& rhoUOldOld =
        mesh.lookupObject<volVectorField>("rhoUOldOld");

    const volScalarField& rhoEOld =
        mesh.lookupObject<volScalarField>("rhoEOld");
    const volScalarField& rhoEOldOld =
        mesh.lookupObject<volScalarField>("rhoEOldOld");



	volScalarField TMean=T;
	volScalarField TMean_old=TMean;
	volScalarField TMean_old2=TMean;
	
//	volScalarField& rhoOld = scalarFieldList[0];
//	volScalarField& rhoOldOld = scalarFieldList[0+scalarFieldNames.size()];

//     label rhoUOldFieldI = 0;
//     volVectorField& rhoUOld    = vectorFieldList[rhoUOldFieldI];
//     volVectorField& rhoUOldOld = vectorFieldList[rhoUOldFieldI+vectorFieldNames.size()];

//	volScalarField& rhoEOld = scalarFieldList[1];
//	volScalarField& rhoEOldOld = scalarFieldList[1+scalarFieldNames.size()];

        // update the boundary absolute velocity
//        mrfZones.correctBoundaryVelocity(U);

        // set the pseudo face velocity
//        Godunov.dotX() = mrfZones.faceU();
/*
        // Compute cartesian velocity components from cylindrical ones
        coordsX = mesh.C().component(vector::X);
        coordsY = mesh.C().component(vector::Y);
*/


        // activate sub-cycle loop
        TimeState subCycleTimeState = runTime.subCycle(numberSubCycles);

        // get access to physical deltaT in ddt scheme!
        physDeltaT[1] = subCycleTimeState.deltaT().value();
        physDeltaT[2] = subCycleTimeState.deltaT0().value();

        // Begin sub-cycling - PseudoTime Integration
        // adjustTimeStep and numberSubCycles > 1 does not make any sense
//        for (label j=0; j < numberSubCycles; j++)
        do 
        {
// // time average residual 

			scalar res1;
			scalar res2;
			scalar res3;


// add 2013.7.6
	flag1++;
	volScalarField T_old=T;
	TMean_old=TMean;
	if((runTime.timeIndex()%50)==1)
	{
		TMean_old2=TMean;
	}

	volScalarField Ux=U.component(0);
	volScalarField Uy=U.component(1);
	volScalarField Uz=U.component(2);	
	volScalarField U2=magSqr(U);
   	const volScalarField R_ = ( thermo.Cp() - thermo.Cv() );
	const volScalarField kappa(thermo.Cp()/thermo.Cv());




           // Update subCycle Time
            runTime++;
/*
            // update local time step sizes just once each iteration for
            // all ddt schemes
            if (!adjustTimeStep)
            {
//		localTimeStep.update(maxCo,adjustTimeStep);
                localTimeStep.update(maxCo,adjustTimeStep,deltaT_Local);
            }
*/

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

// the elements of the  preconditioning matrix, refer to Ref. Blazek P334, eqn. 9.66
////////////////////////////////////////////////	

		scalarField tao11=sita().internalField();

		scalarField tao12=tmp_zero.internalField()*rho.internalField()/p.internalField()/tmp_velocity.internalField();

		scalarField tao13=tmp_zero.internalField()*rho.internalField()/p.internalField()/tmp_velocity.internalField();

		scalarField tao14=tmp_zero.internalField()*rho.internalField()/p.internalField()/tmp_velocity.internalField();

		scalarField tao15=-rho.internalField()/T.internalField();


		scalarField tao21=sita().internalField()*Ux.internalField();

		scalarField tao22=rho.internalField();

		scalarField tao23=tmp_zero.internalField()*rho.internalField();

		scalarField tao24=tmp_zero.internalField()*rho.internalField();

		scalarField tao25=-rho.internalField()*Ux.internalField()/T.internalField();


		scalarField tao31=sita().internalField()*Uy.internalField();

		scalarField tao32=tmp_zero.internalField()*rho.internalField();

		scalarField tao33=rho.internalField();

		scalarField tao34=tmp_zero.internalField()*rho.internalField();


		scalarField tao35=-rho.internalField()*Uy.internalField()/T.internalField();


		scalarField tao41=sita().internalField()*Uz.internalField();

		scalarField tao42=tmp_zero.internalField()*rho.internalField();

		scalarField tao43=tmp_zero.internalField()*rho.internalField();

		scalarField tao44=rho.internalField();

		scalarField tao45=-rho.internalField()*Uz.internalField()/T.internalField();


		scalarField tao51=sita().internalField()*(h.internalField()+0.5*U2.internalField())-1;

		scalarField tao52=rho.internalField()*Ux.internalField();

		scalarField tao53=rho.internalField()*Uy.internalField();

		scalarField tao54=rho.internalField()*Uz.internalField();

		scalarField tao55=-rho.internalField()*U2.internalField()*0.5/T.internalField();

////////////////////////////////////////////////

// the elements of the  tranforrmation matrix, refer to Ref. Blazek P333, eqn. 9.63
////////////////////////////////////////////////

		scalarField P11=rho.internalField()/p.internalField();

		scalarField P12=tmp_zero.internalField()*rho.internalField()/p.internalField()/tmp_velocity.internalField();

		scalarField P13=tmp_zero.internalField()*rho.internalField()/p.internalField()/tmp_velocity.internalField();

		scalarField P14=tmp_zero.internalField()*rho.internalField()/p.internalField()/tmp_velocity.internalField();

		scalarField P15=-rho.internalField()/T.internalField();


		scalarField P21=rho.internalField()/p.internalField()*Ux.internalField();

		scalarField P22=rho.internalField();

		scalarField P23=tmp_zero.internalField()*rho.internalField();

		scalarField P24=tmp_zero.internalField()*rho.internalField();

		scalarField P25=-rho.internalField()*Ux.internalField()/T.internalField();


		scalarField P31=rho.internalField()/p.internalField()*Uy.internalField();

		scalarField P32=tmp_zero.internalField()*rho.internalField();

		scalarField P33=rho.internalField();

		scalarField P34=tmp_zero.internalField()*rho.internalField();

		scalarField P35=-rho.internalField()*Uy.internalField()/T.internalField();


		scalarField P41=rho.internalField()/p.internalField()*Uz.internalField();

		scalarField P42=tmp_zero.internalField()*rho.internalField();

		scalarField P43=tmp_zero.internalField()*rho.internalField();

		scalarField P44=rho.internalField();

		scalarField P45=-rho.internalField()*Uz.internalField()/T.internalField();


		scalarField P51=rhoE.internalField()/p.internalField();

		scalarField P52=rho.internalField()*Ux.internalField();

		scalarField P53=rho.internalField()*Uy.internalField();

		scalarField P54=rho.internalField()*Uz.internalField();

		scalarField P55=-rho.internalField()*U2.internalField()*0.5/T.internalField();
////////////////////////////////////////////////

// the elements of the matrix(tao+3*(pseudo time)/2(physical time)*P), refer to section 2.4.Dual Time eqn.20, in the paper
////////////////////////////////////////////////

		scalarField A11=tao11+3/2*runTime.deltaTValue()/physDeltaT[1]*P11;

		scalarField A12=tao12+3/2*runTime.deltaTValue()/physDeltaT[1]*P12;

		scalarField A13=tao13+3/2*runTime.deltaTValue()/physDeltaT[1]*P13;

		scalarField A14=tao14+3/2*runTime.deltaTValue()/physDeltaT[1]*P14;
	
		scalarField A15=tao15+3/2*runTime.deltaTValue()/physDeltaT[1]*P15;

		scalarField A21=tao21+3/2*runTime.deltaTValue()/physDeltaT[1]*P21;

		scalarField A22=tao22+3/2*runTime.deltaTValue()/physDeltaT[1]*P22;

		scalarField A23=tao23+3/2*runTime.deltaTValue()/physDeltaT[1]*P23;

		scalarField A24=tao24+3/2*runTime.deltaTValue()/physDeltaT[1]*P24;	
	
		scalarField A25=tao25+3/2*runTime.deltaTValue()/physDeltaT[1]*P25;	

		scalarField A31=tao31+3/2*runTime.deltaTValue()/physDeltaT[1]*P31;

		scalarField A32=tao32+3/2*runTime.deltaTValue()/physDeltaT[1]*P32;

		scalarField A33=tao33+3/2*runTime.deltaTValue()/physDeltaT[1]*P33;

		scalarField A34=tao34+3/2*runTime.deltaTValue()/physDeltaT[1]*P34;	
	
		scalarField A35=tao35+3/2*runTime.deltaTValue()/physDeltaT[1]*P35;	

		scalarField A41=tao41+3/2*runTime.deltaTValue()/physDeltaT[1]*P41;

		scalarField A42=tao42+3/2*runTime.deltaTValue()/physDeltaT[1]*P42;

		scalarField A43=tao43+3/2*runTime.deltaTValue()/physDeltaT[1]*P43;

		scalarField A44=tao44+3/2*runTime.deltaTValue()/physDeltaT[1]*P44;	
	
		scalarField A45=tao45+3/2*runTime.deltaTValue()/physDeltaT[1]*P45;	

		scalarField A51=tao51+3/2*runTime.deltaTValue()/physDeltaT[1]*P51;

		scalarField A52=tao52+3/2*runTime.deltaTValue()/physDeltaT[1]*P52;

		scalarField A53=tao53+3/2*runTime.deltaTValue()/physDeltaT[1]*P53;

		scalarField A54=tao54+3/2*runTime.deltaTValue()/physDeltaT[1]*P54;	
	
		scalarField A55=tao55+3/2*runTime.deltaTValue()/physDeltaT[1]*P55;
		
		scalarField Ainv11=A11;	

		scalarField Ainv12=A12;

		scalarField Ainv13=A13;	

		scalarField Ainv14=A14;

		scalarField Ainv15=A15;	

		scalarField Ainv21=A21;	

		scalarField Ainv22=A22;

		scalarField Ainv23=A23;	

		scalarField Ainv24=A24;

		scalarField Ainv25=A25;	

		scalarField Ainv31=A31;	

		scalarField Ainv32=A32;

		scalarField Ainv33=A33;	

		scalarField Ainv34=A34;

		scalarField Ainv35=A35;	

		scalarField Ainv41=A41;	

		scalarField Ainv42=A42;

		scalarField Ainv43=A43;	

		scalarField Ainv44=A44;

		scalarField Ainv45=A45;	

		scalarField Ainv51=A51;	

		scalarField Ainv52=A52;

		scalarField Ainv53=A53;	

		scalarField Ainv54=A54;

		scalarField Ainv55=A55;	
////////////////////////////////////////////////
		forAll(A55, cellI)
		{
			SquareMatrix<scalar> hmm(5);
			

			hmm[0][0] = A11[cellI];
			hmm[0][1] = A12[cellI];
			hmm[0][2] = A13[cellI];
			hmm[0][3] = A14[cellI];
			hmm[0][4] = A15[cellI];

			hmm[1][0] = A21[cellI];
			hmm[1][1] = A22[cellI];
			hmm[1][2] = A23[cellI];
			hmm[1][3] = A24[cellI];
			hmm[1][4] = A25[cellI];

			hmm[2][0] = A31[cellI];
			hmm[2][1] = A32[cellI];
			hmm[2][2] = A33[cellI];
			hmm[2][3] = A34[cellI];
			hmm[2][4] = A35[cellI];

			hmm[3][0] = A41[cellI];
			hmm[3][1] = A42[cellI];
			hmm[3][2] = A43[cellI];
			hmm[3][3] = A44[cellI];
			hmm[3][4] = A45[cellI];

			hmm[4][0] = A51[cellI];
			hmm[4][1] = A52[cellI];
			hmm[4][2] = A53[cellI];
			hmm[4][3] = A54[cellI];
			hmm[4][4] = A55[cellI];




			SquareMatrix<scalar> hmminv(5);

			hmminv[0][0] = Ainv11[cellI];
			hmminv[0][1] = Ainv12[cellI];
			hmminv[0][2] = Ainv13[cellI];
			hmminv[0][3] = Ainv14[cellI];
			hmminv[0][4] = Ainv15[cellI];

			hmminv[1][0] = Ainv21[cellI];
			hmminv[1][1] = Ainv22[cellI];
			hmminv[1][2] = Ainv23[cellI];
			hmminv[1][3] = Ainv24[cellI];
			hmminv[1][4] = Ainv25[cellI];

			hmminv[2][0] = Ainv31[cellI];
			hmminv[2][1] = Ainv32[cellI];
			hmminv[2][2] = Ainv33[cellI];
			hmminv[2][3] = Ainv34[cellI];
			hmminv[2][4] = Ainv35[cellI];

			hmminv[3][0] = Ainv41[cellI];
			hmminv[3][1] = Ainv42[cellI];
			hmminv[3][2] = Ainv43[cellI];
			hmminv[3][3] = Ainv44[cellI];
			hmminv[3][4] = Ainv45[cellI];

			hmminv[4][0] = Ainv51[cellI];
			hmminv[4][1] = Ainv52[cellI];
			hmminv[4][2] = Ainv53[cellI];
			hmminv[4][3] = Ainv54[cellI];
			hmminv[4][4] = Ainv55[cellI];

			rinv  c(5);
			c.input_matrix( hmm);         
			c.inv ();
			c.output_matrix (hmminv);

			Ainv11[cellI] = hmminv[0][0];
			Ainv12[cellI] = hmminv[0][1];
			Ainv13[cellI] = hmminv[0][2];
			Ainv14[cellI] = hmminv[0][3];
			Ainv15[cellI] = hmminv[0][4];

			Ainv21[cellI] = hmminv[1][0];
			Ainv22[cellI] = hmminv[1][1];
			Ainv23[cellI] = hmminv[1][2];
			Ainv24[cellI] = hmminv[1][3];
			Ainv25[cellI] = hmminv[1][4];

			Ainv31[cellI] = hmminv[2][0];
			Ainv32[cellI] = hmminv[2][1];
			Ainv33[cellI] = hmminv[2][2];
			Ainv34[cellI] = hmminv[2][3];
			Ainv35[cellI] = hmminv[2][4];

			Ainv41[cellI] = hmminv[3][0];
			Ainv42[cellI] = hmminv[3][1];
			Ainv43[cellI] = hmminv[3][2];
			Ainv44[cellI] = hmminv[3][3];
			Ainv45[cellI] = hmminv[3][4];

			Ainv51[cellI] = hmminv[4][0];
			Ainv52[cellI] = hmminv[4][1];
			Ainv53[cellI] = hmminv[4][2];
			Ainv54[cellI] = hmminv[4][3];
			Ainv55[cellI] = hmminv[4][4];



		}

		

// // Refer to section 2.3.4., discrete form of source term of eqn. (20).
////////////////////////////////////////////////
		rhoFlux2.internalField()=rhoFlux.internalField()-(3*rho.oldTime().internalField()-4*rhoOld.internalField()+rhoOldOld.internalField())/2/physDeltaT[1];

		rhoUFlux2.internalField()=rhoUFlux.internalField()-(3*rhoU.oldTime().internalField()-4*rhoUOld.internalField()+rhoUOldOld.internalField())/2/physDeltaT[1];
	
		rhoEFlux2.internalField()=rhoEFlux.internalField()-(3*rhoE.oldTime().internalField()-4*rhoEOld.internalField()+rhoEOldOld.internalField())/2/physDeltaT[1];	
////////////////////////////////////////////////

// add 2013.7.4
		volScalarField rhoUFlux2_x=rhoUFlux2.component(0);
		volScalarField rhoUFlux2_y=rhoUFlux2.component(1);
		volScalarField rhoUFlux2_z=rhoUFlux2.component(2);

// // Refer to section 2.3.4., final discrete form of source term of eqn. (20) of which is multiplied by the inverse matrix of (tao+3*(pseudo time)/2(physical time)*P).

////////////////////////////////////////////////
		rhoFlux3.internalField()=Ainv11*rhoFlux2.internalField()+Ainv12*rhoUFlux2_x.internalField()+Ainv13*rhoUFlux2_y.internalField()+Ainv14*rhoUFlux2_z.internalField()+Ainv15*rhoEFlux2.internalField();

		forAll(rhoUFlux3.internalField(),cellI)
		{

			rhoUFlux3.internalField()[cellI].x()=Ainv21[cellI]*rhoFlux2.internalField()[cellI]+Ainv22[cellI]*rhoUFlux2.internalField()[cellI].x()+Ainv23[cellI]*rhoUFlux2.internalField()[cellI].y()+Ainv24[cellI]*rhoUFlux2.internalField()[cellI].z()+Ainv25[cellI]*rhoEFlux2.internalField()[cellI];

			rhoUFlux3.internalField()[cellI].y()=Ainv31[cellI]*rhoFlux2.internalField()[cellI]+Ainv32[cellI]*rhoUFlux2.internalField()[cellI].x()+Ainv33[cellI]*rhoUFlux2.internalField()[cellI].y()+Ainv34[cellI]*rhoUFlux2.internalField()[cellI].z()+Ainv35[cellI]*rhoEFlux2.internalField()[cellI];

			rhoUFlux3.internalField()[cellI].z()=Ainv41[cellI]*rhoFlux2.internalField()[cellI]+Ainv42[cellI]*rhoUFlux2.internalField()[cellI].x()+Ainv43[cellI]*rhoUFlux2.internalField()[cellI].y()+Ainv44[cellI]*rhoUFlux2.internalField()[cellI].z()+Ainv45[cellI]*rhoEFlux2.internalField()[cellI];


		}	

		rhoEFlux3.internalField()=Ainv51*rhoFlux2.internalField()+Ainv52*rhoUFlux2_x.internalField()+Ainv53*rhoUFlux2_y.internalField()+Ainv54*rhoUFlux2_z.internalField()+Ainv55*rhoEFlux2.internalField();


		volScalarField Cp_rhoEFlux3=thermo.Cp()*rhoEFlux3;
////////////////////////////////////////////////
                // Add source terms to the numerical fluxes due to the rotation
//                mrfZones.addCoriolis(rho,U,rhoUFlux);

                // Pseudo time integration
                solve
                (
//                    fvm::ddt(rho) == rhoFlux
//                    fvm::ddt(p) == rhoFlux2
                    fvm::ddt(p) == rhoFlux3
                );

                // time integration
                solve
                (
//                    fvm::ddt(rhoU) == rhoUFlux
//                    fvm::ddt(U) == rhoUFlux2
                    fvm::ddt(U) == rhoUFlux3
                );

                // time integration
                solve
                (
//                    fvm::ddt(rhoE) == rhoEFlux
//					fvm::ddt(T) == rhoEFlux2
//					fvm::ddt(h) == Cp_rhoEFlux2
					fvm::ddt(h) == Cp_rhoEFlux3
                );

                // bound density
//                boundMinMax(rho,rhoMin,rhoMax);

                // bound rhoE
//                boundMinMax(rhoE,rhoEMin,rhoEMax);

                // Compute internal field of U


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
//                scalar L2NormRho = max(Foam::sqrt(sum(sqr(rhoFlux.internalField()))
                scalar L2NormRho = max(Foam::sqrt(gSum(sqr(rhoFlux.internalField()))
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


		Info << "res1=" << res1<< endl;
		Info << "res2=" << res2<< endl;
		Info << "res3=" << res3<< endl;
		Info << "flag= " << flag << endl;
		Info << "flag1= " << flag1 << endl;
            // End outer Loop for Runge - Kutta
            }
            // Update relative velocity
//            mrfZones.relativeVelocity(U,URel);

            //Update turbulence after the multi-stage time integration
            turbulence->correct();

//		if(Foam::log10(L2NormRho) <=0.0 || res_T <0.0001)

// time average residual 
		res1=gMax(mag(T.internalField()-T_old.internalField()));
		TMean=(TMean*(flag1-1)+T)/flag1;

		res2=gMax(mag(TMean.internalField()-TMean_old.internalField()));

		if((runTime.timeIndex()%50)==1)
		{
			res3=gMax(mag(TMean.internalField()-TMean_old2.internalField()));	
		}



// add 2013.7.11
//	runTime.controlDict().readOpt() = IOobject::MUST_READ_IF_MODIFIED;

		 scalar convergence_Res2;	
{
		    IOdictionary controlCoeffDict
		    (
		        IOobject
		        (
		            "controlDict",
		            runTime.system(),
		            mesh,
		            IOobject::MUST_READ,
		            IOobject::NO_WRITE
		        )
		    );

		convergence_Res2 =
			controlCoeffDict.lookupOrDefault<scalar>("convergence_Res", 1e-4);

		}


	Info << "convergence_Res2=" <<  convergence_Res2 << endl;

			if(res2 < convergence_Res2 )
			{
				flag=1;
				break;
			}
		reduce(flag,sumOp<scalar>());
		Info << "flag= " << flag << endl;

        }while(flag==0);

        runTime.endSubCycle();



#       include "updateDualTimeSteppingFields.H"

}
