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
Class
    Foam::fv::EulerLocalDdtScheme
Author
    Oliver Borm  All rights reserved.
\*---------------------------------------------------------------------------*/
//#pragma once
#include "EulerLocalDdtScheme.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
EulerLocalDdtScheme<Type>::fvcDdt
(
    const dimensioned<Type>& dt
)
{
    const objectRegistry& registry = this->mesh();

    // get access to the scalar beta[i]
    const scalarField& beta =
        registry.lookupObject<scalarField>(deltaTName_);

    volScalarField rDeltaT =
        1.0/(beta[0]*registry.lookupObject<volScalarField>(deltaTauName_));

    IOobject ddtIOobject
    (
        "ddt("+dt.name()+')',
        mesh().time().timeName(),
        mesh()
    );

    if (mesh().moving())
    {
        tmp<GeometricField<Type, fvPatchField, volMesh> > tdtdt
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                dimensioned<Type>
                (
                    "0",
                    dt.dimensions()/dimTime,
                    pTraits<Type>::zero
                )
            )
        );

        tdtdt.ref().ref() =
            rDeltaT.internalField()*dt.value()*(1.0 - mesh().V0()/mesh().V());

        return tdtdt;
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                dimensioned<Type>
                (
                    "0",
                    dt.dimensions()/dimTime,
                    pTraits<Type>::zero
                ),
                calculatedFvPatchField<Type>::typeName
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
EulerLocalDdtScheme<Type>::fvcDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const objectRegistry& registry = this->mesh();

    // get access to the scalar beta[i]
    const scalarField& beta =
        registry.lookupObject<scalarField>(deltaTName_);

    volScalarField rDeltaT =
        1.0/(beta[0]*registry.lookupObject<volScalarField>(deltaTauName_));

    IOobject ddtIOobject
    (
        "ddt("+vf.name()+')',
        mesh().time().timeName(),
        mesh()
    );

    if (mesh().moving())
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*vf.dimensions(),
                rDeltaT.internalField()*
                (
                    vf.internalField()
                  - vf.oldTime().internalField()*mesh().V0()/mesh().V()
                ),
                rDeltaT.boundaryField()*
                (
                    vf.boundaryField() - vf.oldTime().boundaryField()
                )
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*(vf - vf.oldTime())
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
EulerLocalDdtScheme<Type>::fvcDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const objectRegistry& registry = this->mesh();

    // get access to the scalar beta[i]
    const scalarField& beta =
        registry.lookupObject<scalarField>(deltaTName_);

    volScalarField rDeltaT =
        1.0/(beta[0]*registry.lookupObject<volScalarField>(deltaTauName_));

    IOobject ddtIOobject
    (
        "ddt("+rho.name()+','+vf.name()+')',
        mesh().time().timeName(),
        mesh()
    );

    if (mesh().moving())
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                rDeltaT.internalField()*rho.value()*
                (
                    vf.internalField()
                  - vf.oldTime().internalField()*mesh().V0()/mesh().V()
                ),
                rDeltaT.boundaryField()*rho.value()*
                (
                    vf.boundaryField() - vf.oldTime().boundaryField()
                )
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*rho*(vf - vf.oldTime())
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
EulerLocalDdtScheme<Type>::fvcDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const objectRegistry& registry = this->mesh();

    // get access to the scalar beta[i]
    const scalarField& beta =
        registry.lookupObject<scalarField>(deltaTName_);

    volScalarField rDeltaT =
        1.0/(beta[0]*registry.lookupObject<volScalarField>(deltaTauName_));

    IOobject ddtIOobject
    (
        "ddt("+rho.name()+','+vf.name()+')',
        mesh().time().timeName(),
        mesh()
    );

    if (mesh().moving())
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                rDeltaT.internalField()*
                (
                    rho.internalField()*vf.internalField()
                  - rho.oldTime().internalField()
                   *vf.oldTime().internalField()*mesh().V0()/mesh().V()
                ),
                rDeltaT.boundaryField()*
                (
                    rho.boundaryField()*vf.boundaryField()
                  - rho.oldTime().boundaryField()
                   *vf.oldTime().boundaryField()
                )
            )
        );
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh> >
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*(rho*vf - rho.oldTime()*vf.oldTime())
            )
        );
    }
}


template<class Type>
tmp<fvMatrix<Type> >
EulerLocalDdtScheme<Type>::fvmDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const objectRegistry& registry = this->mesh();

    // get access to the scalar beta[i]
    const scalarField& beta =
        registry.lookupObject<scalarField>(deltaTName_);

    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            vf.dimensions()*dimVol/dimTime
        )
    );

    fvMatrix<Type>& fvm = tfvm.ref();


//    scalarField rDeltaT =
//        1.0/(beta[0]*registry.lookupObject<volScalarField>(deltaTauName_).internalField());

    scalar rDeltaT =
        1.0/(beta[0]*mesh().time().deltaTValue());

// 	scalar rDeltaT = 1.0/mesh().time().deltaTValue();

    fvm.diag() = rDeltaT*mesh().V();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT*vf.oldTime().internalField()*mesh().V0();
    }
    else
    {
        fvm.source() = rDeltaT*vf.oldTime().internalField()*mesh().V();
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type> >
EulerLocalDdtScheme<Type>::fvmDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const objectRegistry& registry = this->mesh();

    // get access to the scalar beta[i]
    const scalarField& beta =
        registry.lookupObject<scalarField>(deltaTName_);

    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    scalarField rDeltaT =
        1.0/(beta[0]*registry.lookupObject<volScalarField>(deltaTauName_).internalField());

    fvm.diag() = rDeltaT*rho.value()*mesh().V();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT
            *rho.value()*vf.oldTime().internalField()*mesh().V0();
    }
    else
    {
        fvm.source() = rDeltaT
            *rho.value()*vf.oldTime().internalField()*mesh().V();
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type> >
EulerLocalDdtScheme<Type>::fvmDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const objectRegistry& registry = this->mesh();

    // get access to the scalar beta[i]
    const scalarField& beta =
        registry.lookupObject<scalarField>(deltaTName_);

    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    scalarField rDeltaT =
        1.0/(beta[0]*registry.lookupObject<volScalarField>(deltaTauName_).internalField());

    fvm.diag() = rDeltaT*rho.internalField()*mesh().V();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT
            *rho.oldTime().internalField()
            *vf.oldTime().internalField()*mesh().V0();
    }
    else
    {
        fvm.source() = rDeltaT
            *rho.oldTime().internalField()
            *vf.oldTime().internalField()*mesh().V();
    }

    return tfvm;
}


template<class Type>
tmp<typename EulerLocalDdtScheme<Type>::fluxFieldType>
EulerLocalDdtScheme<Type>::fvcDdtPhiCorr
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
)
{
    IOobject ddtIOobject
    (
        "ddtCorr("+ U.name() + ',' + phi.name() + ')',
        mesh().time().timeName(),
        mesh()
    );

    if (mesh().moving())
    {
        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                ddtIOobject,
                mesh(),
                dimensioned<typename flux<Type>::type>
                (
                    "0",
                    phi.dimensions()/dimTime,
                    pTraits<typename flux<Type>::type>::zero
                )
            )
        );
    }
    else
    {
        const objectRegistry& registry = this->mesh();

        // get access to the scalar beta[i]
        const scalarField& beta =
            registry.lookupObject<scalarField>(deltaTName_);

        volScalarField rDeltaT =
            1.0/(beta[0]*registry.lookupObject<volScalarField>(deltaTauName_));

        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                ddtIOobject,
                this->fvcDdtPhiCoeff(U.oldTime(), phi.oldTime())*
                (
                    fvc::interpolate(rDeltaT)*phi.oldTime()
                  - (fvc::interpolate(rDeltaT*U.oldTime()) & mesh().Sf())
                )
            )
        );
    }
}


template<class Type>
tmp<typename EulerLocalDdtScheme<Type>::fluxFieldType>
EulerLocalDdtScheme<Type>::fvcDdtPhiCorr
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
)
{
    IOobject ddtIOobject
    (
        "ddtCorr(" +rho.name() + ',' + U.name() + ',' + phi.name() + ')',
        mesh().time().timeName(),
        mesh()
    );

    if (mesh().moving())
    {
        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                ddtIOobject,
                mesh(),
                dimensioned<typename flux<Type>::type>
                (
                    "0",
                    phi.dimensions()/dimTime,
                    pTraits<typename flux<Type>::type>::zero
                )
            )
        );
    }
    else
    {
        const objectRegistry& registry = this->mesh();

        // get access to the scalar beta[i]
        const scalarField& beta =
            registry.lookupObject<scalarField>(deltaTName_);

        volScalarField rDeltaT =
            1.0/(beta[0]*registry.lookupObject<volScalarField>(deltaTauName_));

        if
        (
            U.dimensions() == dimVelocity
         && phi.dimensions() == dimVelocity*dimArea
        )
        {
            return tmp<fluxFieldType>
            (
                new fluxFieldType
                (
                    ddtIOobject,
                    this->fvcDdtPhiCoeff(U.oldTime(), phi.oldTime())
                   *(
                        fvc::interpolate(rDeltaT*rho.oldTime())*phi.oldTime()
                      - (fvc::interpolate(rDeltaT*rho.oldTime()*U.oldTime())
                      & mesh().Sf())
                    )
                )
            );
        }
        else if
        (
            U.dimensions() == dimVelocity
         && phi.dimensions() == dimDensity*dimVelocity*dimArea
        )
        {
            return tmp<fluxFieldType>
            (
                new fluxFieldType
                (
                    ddtIOobject,
                    this->fvcDdtPhiCoeff
                    (
                        U.oldTime(),
                        phi.oldTime()/fvc::interpolate(rho.oldTime())
                    )
                   *(
                        fvc::interpolate(rDeltaT*rho.oldTime())
                       *phi.oldTime()/fvc::interpolate(rho.oldTime())
                      - (
                            fvc::interpolate
                            (
                                rDeltaT*rho.oldTime()*U.oldTime()
                            ) & mesh().Sf()
                        )
                    )
                )
            );
        }
        else if
        (
            U.dimensions() == dimDensity*dimVelocity
         && phi.dimensions() == dimDensity*dimVelocity*dimArea
        )
        {
            return tmp<fluxFieldType>
            (
                new fluxFieldType
                (
                    ddtIOobject,
                    this->fvcDdtPhiCoeff(rho.oldTime(), U.oldTime(), phi.oldTime())
                   *(
                        fvc::interpolate(rDeltaT)*phi.oldTime()
                      - (
                            fvc::interpolate(rDeltaT*U.oldTime())&mesh().Sf()
                        )
                    )
                )
            );
        }
        else
        {
            FatalErrorIn
            (
                "EulerLocalDdtScheme<Type>::fvcDdtPhiCorr"
            )   << "dimensions of phi are not correct"
                << abort(FatalError);

            return fluxFieldType::null();
        }
    }
}


template<class Type>
tmp<surfaceScalarField> EulerLocalDdtScheme<Type>::meshPhi
(
    const GeometricField<Type, fvPatchField, volMesh>&
)
{
    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                "meshPhi",
                mesh().time().timeName(),
                mesh()
            ),
            mesh(),
            dimensionedScalar("0", dimVolume/dimTime, 0.0)
        )
    );
}



template<class Type>
tmp<typename EulerLocalDdtScheme<Type>::fluxFieldType>
EulerLocalDdtScheme<Type>::fvcDdtUfCorr
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
)
{
//      const surfaceScalarField rDeltaT(fvc::interpolate(localRDeltaT()));
        const objectRegistry& registry = this->mesh();
        // get access to the scalar beta[i]
        const scalarField& beta =
            registry.lookupObject<scalarField>(deltaTName_);
        volScalarField rDeltaT =
            1.0/(beta[0]*registry.lookupObject<volScalarField>(deltaTauName_));

    fluxFieldType phiUf0(mesh().Sf() & Uf.oldTime());
    fluxFieldType phiCorr
    (
        phiUf0 - fvc::dotInterpolate(mesh().Sf(), U.oldTime())
    );

    return tmp<fluxFieldType>
    (
        new fluxFieldType
        (
            IOobject
            (
                "ddtCorr(" + U.name() + ',' + Uf.name() + ')',
                mesh().time().timeName(),
                mesh()
            ),
            this->fvcDdtPhiCoeff(U.oldTime(), phiUf0, phiCorr)
           *rDeltaT*phiCorr 
        )
    );
}


template<class Type>
tmp<typename EulerLocalDdtScheme<Type>::fluxFieldType>
EulerLocalDdtScheme<Type>::fvcDdtUfCorr
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
)
{
    if
    (
        U.dimensions() == dimVelocity
     && Uf.dimensions() == dimDensity*dimVelocity
    )
    {
//      const surfaceScalarField rDeltaT(fvc::interpolate(localRDeltaT()));
        const objectRegistry& registry = this->mesh();
        // get access to the scalar beta[i]
        const scalarField& beta =
            registry.lookupObject<scalarField>(deltaTName_);
        volScalarField rDeltaT =
            1.0/(beta[0]*registry.lookupObject<volScalarField>(deltaTauName_));
			
        GeometricField<Type, fvPatchField, volMesh> rhoU0
        (
            rho.oldTime()*U.oldTime()
        );

        fluxFieldType phiUf0(mesh().Sf() & Uf.oldTime());
        fluxFieldType phiCorr(phiUf0 - fvc::dotInterpolate(mesh().Sf(), rhoU0));

        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                IOobject
                (
                    "ddtCorr("
                  + rho.name() + ',' + U.name() + ',' + Uf.name() + ')',
                    mesh().time().timeName(),
                    mesh()
                ),
                this->fvcDdtPhiCoeff(rhoU0, phiUf0, phiCorr)*rDeltaT*phiCorr
            )
        );
    }
    else if
    (
        U.dimensions() == dimDensity*dimVelocity
     && Uf.dimensions() == dimDensity*dimVelocity
    )
    {
        return fvcDdtUfCorr(U, Uf);
    }
    else
    {
        FatalErrorInFunction
            << "dimensions of Uf are not correct"
            << abort(FatalError);

        return fluxFieldType::null();
    }
}

template<class Type>
  tmp<GeometricField<Type, fvPatchField, volMesh>>
    EulerLocalDdtScheme<Type>::fvcDdt
      (
      const volScalarField& alpha,
      const volScalarField& rho,
      const GeometricField<Type, fvPatchField, volMesh>& vf
    )
    {
      NotImplemented;
      return tmp<GeometricField<Type, fvPatchField, volMesh>>();//return null tmp
    }

template<class Type>
  tmp<fvMatrix<Type>>
  EulerLocalDdtScheme<Type>::fvmDdt
    (
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
  )
  {
    NotImplemented;
    return tmp<fvMatrix<Type>>();//return null tmp
  }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
