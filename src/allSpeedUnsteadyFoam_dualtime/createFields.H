    Info<< "Reading thermophysical properties\n" << endl;

    autoPtr<psiThermo> pThermo
    (
        psiThermo::New(mesh)
    );
    psiThermo& thermo = pThermo();

    volScalarField& h = thermo.he();
//     volScalarField& e = thermo.e();

    volScalarField& p = thermo.p();
    const volScalarField& T = thermo.T();
//	volScalarField& T = thermo.T();
//     const volScalarField& psi = thermo.psi();

    Info<< "Allocating field rho\n" << endl;
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo.rho()
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
    mesh
    );
// The square of the reference velocity Ur
// Refer to Ref Blazek, Computational Fluid Dynamics: Principles and Applications, Elsevier, 2005. P335
    Info<< "Reading field U2\n" << endl;
    volScalarField Ur2
    (
        IOobject
        (
            "Ur2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
	dimensionedScalar("Ur2", dimVelocity*dimVelocity, 1.0)
    );
/*
// the preconditioning parameter sita
// Refer to Ref Blazek, Computational Fluid Dynamics: Principles and Applications, Elsevier, 2005. P334-P335
    Info<< "Reading field sita\n" << endl;
    volScalarField sita
    (
        IOobject
        (
            "sita",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
	dimensionedScalar("sita", pow(dimVelocity,-2), 1.0)

    );
// the temporary parameter fai
// Refer to Ref Blazek, Computational Fluid Dynamics: Principles and Applications, Elsevier, 2005. P33, eqn.9.67
    Info<< "Reading field fai\n" << endl;
    volScalarField fai
    (
        IOobject
        (
            "fai",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
	dimensionedScalar("fai", dimless, 1.0)

    );
// the temporary parameter a2
// Refer to Ref Blazek, Computational Fluid Dynamics: Principles and Applications, Elsevier, 2005. P33, eqn.9.67
    Info<< "Reading field a2\n" << endl;
    volScalarField a2
    (
        IOobject
        (
            "a2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
	dimensionedScalar("a2", dimless, 1.0)

    );
// the temporary parameter a3
// Refer to Ref Blazek, Computational Fluid Dynamics: Principles and Applications, Elsevier, 2005. P33, eqn.9.67
    Info<< "Reading field a3\n" << endl;
    volScalarField a3
    (
        IOobject
        (
            "a3",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
	dimensionedScalar("a3", dimTemperature/dimDensity, 1.0)

    );
*/
    Info<< "Reading field rhoUFlux2\n" << endl;
    volVectorField rhoUFlux2
    (
        IOobject
        (
            "rhoUFlux2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
	dimensionedVector("rhoUFlux2", dimensionSet(1, -2, -2, 0, 0, 0, 0)/dimDensity, vector(1.0, 1.0, 1.0))
    );

    Info<< "Reading field rhoFlux2\n" << endl;
    volScalarField rhoFlux2
    (
        IOobject
        (
            "rhoFlux2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
	dimensionedScalar("rhoFlux2", dimensionSet(1, -1, -3, 0, 0, 0, 0), 1.0)
    );

    Info<< "Reading field rhoEFlux2\n" << endl;
    volScalarField rhoEFlux2
    (
        IOobject
        (
            "rhoEFlux2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
	dimensionedScalar("rhoEFlux2", dimensionSet(0, 0, -1, 1, 0, 0, 0), 1.0)
    );

// add 2013.7.4
    Info<< "Reading field rhoUFlux3\n" << endl;
    volVectorField rhoUFlux3
    (
        IOobject
        (
            "rhoUFlux3",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
    mesh,
	dimensionedVector("rhoUFlux3", dimensionSet(1, -2, -2, 0, 0, 0, 0)/dimDensity, vector(1.0, 1.0, 1.0))
    );

    Info<< "Reading field rhoFlux3\n" << endl;
    volScalarField rhoFlux3
    (
        IOobject
        (
            "rhoFlux3",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
    mesh,
	dimensionedScalar("rhoFlux3", dimensionSet(1, -1, -3, 0, 0, 0, 0), 1.0)
    );

    Info<< "Reading field rhoEFlux3\n" << endl;
    volScalarField rhoEFlux3
    (
        IOobject
        (
            "rhoEFlux3",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
	dimensionedScalar("rhoEFlux3", dimensionSet(0, 0, -1, 1, 0, 0, 0), 1.0)
    );

//    Info<< "Reading field p2\n" << endl;
//    volScalarField p2
//    (
//        IOobject
//        (
//            "p2",
//            runTime.timeName(),
//            mesh,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE
//        ),
//    mesh,
//	dimensionedScalar("p2", dimPressure, 648.1)
//    );


#   include "compressibleCreatePhi.H"
/*
    Info<< "Creating MRF model\n" << endl;
    MRFZones mrfZones(mesh);
*/
/*
    volVectorField URel
    (
        IOobject
        (
            "URel",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U
    );
    URel.correctBoundaryConditions();
    mrfZones.relativeVelocity(U,URel);
*/
    Info<< "Creating turbulence model\n" << endl;
    // all turbulence models which use the skew of this velocity gradient will
    // fail for rotating meshs, like SpalartAllmaras or kOmegaSST_LowRe
    // Otherwise if using URel, then you will become problems at Rotor-Stator
    // Interfaces, as the grad(Urel) will be computed wrong there
    autoPtr<compressible::turbulenceModel> turbulence
    (
        compressible::turbulenceModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );

    Info<< "Create Riemann solver\n" << endl;
//   godunovFlux<hllcALEFlux> Godunov(p, U, T, thermo, turbulence());
//     godunovFlux<AUSMplusALEFlux> Godunov(p, U, T, thermo, turbulence());
//     godunovFlux<AUSMplusUpALEFlux> Godunov(p, U, T, thermo, turbulence());     
//     godunovFlux<roeALEFlux> Godunov(p, U, T, thermo, turbulence());
//     godunovFlux<rusanovALEFlux> Godunov(p, U, T, thermo, turbulence());
//     godunovFlux<AUSMALEFlux> Godunov(p, U, T, thermo, turbulence());

// add 2013.7.18
	godunovFlux<AUSMplusPreFlux> Godunov(p, U, T, thermo, turbulence());

    Info<< "Allocating field rhoU\n" << endl;
    volVectorField rhoU
    (
        IOobject
        (
            "rhoU",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho*U
    );

    Info<< "Allocating field rhoE\n" << endl;
    volScalarField rhoE
    (
        IOobject
        (
            "rhoE",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho*(h + 0.5*magSqr(U) + turbulence->k()) - p
//         rho*(e + 0.5*magSqr(U) + turbulence->k())
    );

/*
    volScalarField coordsX(mesh.C().component(vector::X));
    volScalarField coordsY(mesh.C().component(vector::Y));

    dimensionedScalar smallRadius("smallRadius", dimLength, SMALL);
    const volScalarField radius = max(sqrt(sqr(coordsX)+sqr(coordsY)),smallRadius);

    volScalarField rhoCr
    (
        IOobject
        (
            "rhoCr",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (rhoU.component(vector::X)*coordsX+rhoU.component(vector::Y)*coordsY)/radius
    );

    volScalarField rhoCu
    (
        IOobject
        (
            "rhoCu",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (rhoU.component(vector::Y)*coordsX-rhoU.component(vector::X)*coordsY)/radius
    );
*/
//    localTimeStep localTimeStep(URel, thermo, turbulence());

    Info<< "Allocating physDeltaT list for RK and Dual-Time Stepping\n" << endl;
    IOField<scalar> physDeltaT
    (
        IOobject
        (
            "physDeltaT",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        3
    );

    label numberSubCycles = 1;
    label numberRungeKuttaStages = 1;
    scalarList beta(numberRungeKuttaStages, 1.0);
    bool secondOrder = true;

#   include "readMultiStage.H"
#   include "createDualTimeSteppingFields.H"

 volScalarField TOL
 (
     IOobject
     (
         "TOL",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
     ),
     mesh,
     dimensionedScalar("TOL", dimless, 0)
 );
/*
 volScalarField TOL2
 (
     IOobject
     (
         "TOL2",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
     ),
     mesh,
     dimensionedScalar("TOL2", dimless, 0)
 );
*/
// add 2013.7.4
    volScalarField tmp_zero
    (
        IOobject
        (
            "tmp_zero",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
	dimensionedScalar("tmp_zero", dimless, 0.0)

    );

    volScalarField tmp_velocity
    (
        IOobject
        (
            "tmp_velocity",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh,
	dimensionedScalar("tmp_velocity", dimVelocity, 1.0)

    );

// add 2013.7.6
//	IOdictionary pLocs
//	(
//		IOobject
//		(
//			"probeLocations",
//			runTime.constant(),
//			mesh,
//			IOobject::MUST_READ,
//			IOobject::NO_WRITE
//		)
//	);

//	const pointField& probeLocations(pLocs.lookup("probeLocations"));
//	labelList probeCells(probeLocations.size(),-1);
