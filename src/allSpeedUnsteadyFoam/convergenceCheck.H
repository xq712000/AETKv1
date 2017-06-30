                volScalarField rhoFlux = -fvc::div(Godunov.rhoFlux());

//                volVectorField rhoUFlux = -fvc::div(Godunov.rhoUFlux())
//                    + fvc::div(tau);

                scalar L2NormRho = max(Foam::sqrt(sum(sqr(rhoFlux.internalField()))
                    /mesh.nCells()),SMALL);
//                scalar LInfNormRho = max(mag(rhoFlux.internalField()));


//add
		residual=Foam::log10(L2NormRho);
                // Averaged L2-Norm of fluxes
                Info<< "rho L2 Residual: "
                    << Foam::log10(L2NormRho)  << endl
		<< endl;
