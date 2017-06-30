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
    AUSMplusPreFlux

Author
    Oliver Borm  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "AUSMplusPreFlux.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //+


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::AUSMplusPreFlux::evaluateFlux
(
    scalar& rhoFlux,
    vector& rhoUFlux,
    scalar& rhoEFlux,
    const scalar pLeft,
    const scalar pRight,
    const vector ULeft,
    const vector URight,
// //
//     const scalar TLeft,
//     const scalar TRight,
// //
    const scalar rhoLeft,
    const scalar rhoRight,
// //
    const scalar kLeft,
    const scalar kRight,

	const scalar UrefLeft2,
	const scalar UrefRight2,
// //
//     const scalar RLeft,
//     const scalar RRight,
//     const scalar CvLeft,
//     const scalar CvRight,
// //
    const scalar kappaLeft,
    const scalar kappaRight,
// //
    const vector Sf,
    const scalar magSf,
    const vector dotX,
    const scalar Konstant
) const
{
        const scalar alpha = 3.0/16.0;
        const scalar beta = 1.0/8.0;

        // bounding variables
        const scalar rhoMin = SMALL;

        // normal vector
        const vector normalVector = Sf/magSf;

        // speed of sound, for left and right side, assuming perfect gas
        const scalar aLeft  = Foam::sqrt(max(SMALL,kappaLeft *pLeft /max(rhoLeft, rhoMin)));
        const scalar aRight = Foam::sqrt(max(SMALL,kappaRight*pRight/max(rhoRight,rhoMin)));

        // DensityTotalEnergy
        const scalar rhoHLeft =
            pLeft *kappaLeft /(kappaLeft -1.0)+rhoLeft *(0.5*magSqr(ULeft) +kLeft);
        const scalar rhoHRight =
            pRight*kappaRight/(kappaRight-1.0)+rhoRight*(0.5*magSqr(URight)+kRight);

        // Step 1: decode left and right:

        // Compute conservative variables assuming perfect gas law

        // DensityVelocity
        const vector rhoULeft  = rhoLeft *ULeft;
        const vector rhoURight = rhoRight*URight;

        // compute qLeft and qRight (q_{l,r} = U_{l,r} \bullet n)
        const scalar qLeft  = ((ULeft  - dotX) & normalVector);
        const scalar qRight = ((URight - dotX) & normalVector);

        const scalar aTilde = 0.5*(aLeft+aRight);

        const scalar MaRelLeft  = qLeft /aTilde;
        const scalar MaRelRight = qRight/aTilde;

        const scalar magMaRelLeft  = mag(MaRelLeft);
        const scalar magMaRelRight = mag(MaRelRight);

// add 2013.7.18
	const scalar MarefLeft2=UrefLeft2/sqr(aLeft);
	const scalar MarefRight2=UrefRight2/sqr(aRight);
//	const scalar Maref2Tilde=0.5*(MarefLeft2+MarefRight2); // Mref(1/2)2; Tilde ~ 1/2, interface
	const scalar Maref2Tilde=sqr(0.5*(sqrt(MarefLeft2)+sqrt(MarefRight2)));
	const scalar Ma2Tilde2=sqr(0.5*(MaRelLeft+MaRelRight));
	const scalar fTilde=sqrt(sqr(1-Maref2Tilde)* Ma2Tilde2+4*Maref2Tilde)/(1+Maref2Tilde);	

        const scalar aTilde_pie = fTilde*aTilde;
        const scalar MaRelLeft_pie  = qLeft /aTilde_pie;
        const scalar MaRelRight_pie = qRight/aTilde_pie;

        const scalar MaRelLeft_  = 0.5*((1+Maref2Tilde)*MaRelLeft_pie+(1-Maref2Tilde)*MaRelRight_pie);
        const scalar MaRelRight_ = 0.5*((1+Maref2Tilde)*MaRelRight_pie+(1-Maref2Tilde)*MaRelLeft_pie);

        const scalar magMaRelLeft_  = mag(MaRelLeft_);
        const scalar magMaRelRight_ = mag(MaRelRight_);
	
////////////////////////////////////////////

        const scalar Ma1PlusLeft   = 0.5*(MaRelLeft +magMaRelLeft );
        const scalar Ma1MinusRight = 0.5*(MaRelRight-magMaRelRight);

        const scalar Ma2PlusLeft   =  0.25*sqr(MaRelLeft +1.0);
        const scalar Ma2PlusRight  =  0.25*sqr(MaRelRight+1.0);
        const scalar Ma2MinusLeft  = -0.25*sqr(MaRelLeft -1.0);
        const scalar Ma2MinusRight = -0.25*sqr(MaRelRight-1.0);
// add 2013.11.29
//        const scalar Ma4BetaPlusLeft   = ((magMaRelLeft  >= 1.0) ? Ma1PlusLeft   : (Ma2PlusLeft  *(1.0-16.0*beta*Ma2MinusLeft)));
//        const scalar Ma4BetaMinusRight = ((magMaRelRight >= 1.0) ? Ma1MinusRight : (Ma2MinusRight*(1.0+16.0*beta*Ma2PlusRight)));

        const scalar P5alphaPlusLeft   = ((magMaRelLeft  >= 1.0) ?
            (Ma1PlusLeft/MaRelLeft)    : (Ma2PlusLeft  *(( 2.0-MaRelLeft) -16.0*alpha*MaRelLeft *Ma2MinusLeft )));
        const scalar P5alphaMinusRight = ((magMaRelRight >= 1.0) ?
            (Ma1MinusRight/MaRelRight) : (Ma2MinusRight*((-2.0-MaRelRight)+16.0*alpha*MaRelRight*Ma2PlusRight)));
// add 2013.11.29
//        const scalar MaRelTilde = Ma4BetaPlusLeft + Ma4BetaMinusRight;
        const scalar pTilde = pLeft*P5alphaPlusLeft + pRight*P5alphaMinusRight;

// add 2013.11.29
//        const scalar URelTilde = MaRelTilde*aTilde;
//        const scalar magURelTilde = mag(MaRelTilde)*aTilde;
        // There is a typo in Luo et. al, J. Comp. Physics 194 (2004), Chap 4.2 Eq. 4.8
        // refer to the origial Paper from Liou, J. Comp. Physics 129 (1996), Chap4, Eq. 42
//        rhoFlux  = (0.5*(URelTilde*(rhoLeft +rhoRight) -magURelTilde*(rhoRight -rhoLeft)))*magSf;
//        rhoUFlux = (0.5*(URelTilde*(rhoULeft+rhoURight)-magURelTilde*(rhoURight-rhoULeft))+pTilde*normalVector)*magSf;
//        rhoEFlux = (0.5*(URelTilde*(rhoHLeft+rhoHRight)-magURelTilde*(rhoHRight-rhoHLeft))+pTilde*(dotX & normalVector))*magSf;

// add 2013.7.18
//	const scalar rhoFlux_y  = (0.5*(URelTilde*(rhoLeft +rhoRight) -magURelTilde*(rhoRight -rhoLeft)))*magSf;

//	const scalar MaTilde_pos=0.5*(MaRelTilde+mag(MaRelTilde));
//	const scalar MaTilde_neg=0.5*(MaRelTilde-mag(MaRelTilde));
//	const scalar rhoFlux_AUSM = aTilde*(rhoLeft*MaTilde_pos+rhoRight*MaTilde_neg);  // flux at interface (19)

//

        const scalar Ma1PlusLeft_   = 0.5*(MaRelLeft_ +magMaRelLeft_ );
        const scalar Ma1MinusRight_ = 0.5*(MaRelRight_-magMaRelRight_);

        const scalar Ma2PlusLeft_   =  0.25*sqr(MaRelLeft_ +1.0);
        const scalar Ma2PlusRight_  =  0.25*sqr(MaRelRight_+1.0);
        const scalar Ma2MinusLeft_  = -0.25*sqr(MaRelLeft_ -1.0);
        const scalar Ma2MinusRight_ = -0.25*sqr(MaRelRight_-1.0);

        const scalar Ma4BetaPlusLeft_   = ((magMaRelLeft_  >= 1.0) ? Ma1PlusLeft_   : (Ma2PlusLeft_  *(1.0-16.0*beta*Ma2MinusLeft_)));
        const scalar Ma4BetaMinusRight_ = ((magMaRelRight_ >= 1.0) ? Ma1MinusRight_ : (Ma2MinusRight_*(1.0+16.0*beta*Ma2PlusRight_)));

        const scalar P5alphaPlusLeft_   = ((magMaRelLeft_  >= 1.0) ?
            (Ma1PlusLeft_/MaRelLeft_)    : (Ma2PlusLeft_  *(( 2.0-MaRelLeft_) -16.0*alpha*MaRelLeft_ *Ma2MinusLeft_ )));
        const scalar P5alphaMinusRight_ = ((magMaRelRight_ >= 1.0) ?
            (Ma1MinusRight_/MaRelRight_) : (Ma2MinusRight_*((-2.0-MaRelRight_)+16.0*alpha*MaRelRight_*Ma2PlusRight_)));

        const scalar MaRelTilde_ = Ma4BetaPlusLeft_ + Ma4BetaMinusRight_;
        const scalar pTilde_ = pLeft*P5alphaPlusLeft_ + pRight*P5alphaMinusRight_;

	const scalar MaTilde_pos_=0.5*(MaRelTilde_+mag(MaRelTilde_));
	const scalar MaTilde_neg_=0.5*(MaRelTilde_-mag(MaRelTilde_));

	const scalar rhoFlux_AUSM = aTilde_pie*(rhoLeft*MaTilde_pos_+rhoRight*MaTilde_neg_);  // flux at interface (19)

	const scalar rhoFlux_ = rhoFlux_AUSM + aTilde_pie*(1/Maref2Tilde-1)*(pLeft-pRight)/(pLeft/rhoLeft+pRight/rhoRight)*(Ma4BetaPlusLeft_-Ma1PlusLeft_-Ma4BetaMinusRight_+Ma1MinusRight_);		// (42)

	const scalar rhoFlux_pos = 0.5*(rhoFlux_+mag(rhoFlux_));
	const scalar rhoFlux_neg = 0.5*(rhoFlux_-mag(rhoFlux_));

	

	rhoFlux = rhoFlux_*magSf;
	rhoUFlux = magSf*(rhoFlux_pos*ULeft+rhoFlux_neg*URight+pTilde_*normalVector);
//	rhoEFlux = magSf*(rhoFlux_pos*rhoHLeft+rhoFlux_neg*rhoHRight+pTilde*(dotX & normalVector)) ;
	rhoEFlux = magSf*(rhoFlux_pos*rhoHLeft/rhoLeft+rhoFlux_neg*rhoHRight/rhoRight+pTilde*(dotX & normalVector)) ;


//        rhoFlux  = (0.5*(URelTilde*(rhoLeft +rhoRight) -magURelTilde*(rhoRight -rhoLeft)))*magSf;
//        rhoUFlux = (0.5*(URelTilde*(rhoULeft+rhoURight)-magURelTilde*(rhoURight-rhoULeft))+pTilde*normalVector)*magSf;
//        rhoEFlux = (0.5*(URelTilde*(rhoHLeft+rhoHRight)-magURelTilde*(rhoHRight-rhoHLeft))+pTilde*(dotX & normalVector))*magSf;


}

// ************************************************************************* //
