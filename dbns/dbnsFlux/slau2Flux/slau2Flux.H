/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

Class
    slau2Flux

Description
    Implementation of SLAU2 approximate Riemann solver scheme

SourceFiles
    buiFlux.H

\*---------------------------------------------------------------------------*/

#ifndef slau2Flux_H
#define slau2Flux_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------*\
                               Class slau2Flux Declaration
    \*---------------------------------------------------------------------------*/

    class slau2Flux
    {
    public:
        // Constructors

        //- Construct null
        inline slau2Flux()
        {
        }

        // Destructor - default

        // Member Functions

        inline void evaluateFlux(
            scalar &rhoFlux,
            vector &rhoUFlux,
            scalar &rhoEFlux,
            const scalar &pLeft,
            const scalar &pRight,
            const vector &ULeft,
            const vector &URight,
            const scalar &TLeft,
            const scalar &TRight,
            const scalar &RLeft,
            const scalar &RRight,
            const scalar &CvLeft,
            const scalar &CvRight,
            const vector &Sf,
            const scalar &magSf,
            const scalar &buiEps) const
        {
            // Step 1: decode rho left and right:
            scalar rhoLeft = pLeft / (RLeft * TLeft);
            scalar rhoRight = pRight / (RRight * TRight);
            const scalar rhoBar = 0.5 * (rhoLeft + rhoRight);

            // Decode left and right total energy:
            scalar eLeft = CvLeft * TLeft + 0.5 * magSqr(ULeft);
            scalar eRight = CvRight * TRight + 0.5 * magSqr(URight);

            // Adiabatic exponent is constant for ideal gas but if Cp=Cp(T)
            // it must be computed for each cell and evaluated at each face
            // through reconstruction
            const scalar kappaLeft = (CvLeft + RLeft) / CvLeft;
            const scalar kappaRight = (CvRight + RRight) / CvRight;

            // normal vector
            vector normalVector = Sf / magSf;

            // Compute left and right contravariant velocities:
            const scalar contrVLeft = (ULeft & normalVector);
            const scalar contrVRight = (URight & normalVector);
            const scalar contrVBar =
                (rhoLeft * mag(contrVLeft) + rhoRight * mag(contrVRight)) / (rhoLeft + rhoRight);

            // Compute left and right total enthalpies:
            const scalar hLeft = eLeft + pLeft / rhoLeft;
            const scalar hRight = eRight + pRight / rhoRight;

            // Speed of sound
            const scalar cLeft =
                sqrt(max((kappaLeft - 1) * (hLeft - 0.5 * magSqr(ULeft)), SMALL));
            const scalar cRight =
                sqrt(max((kappaRight - 1) * (hRight - 0.5 * magSqr(URight)), SMALL));
            const scalar cBar = 0.5 * (cLeft + cRight);

            const scalar MaLeft = contrVLeft / cBar;
            const scalar MaRight = contrVRight / cBar;
            const scalar g = -max(min(MaLeft, 0.0), -1.0) * min(max(MaRight, 0.0), 1.0);

            const scalar MaArc =
                min(1.0, 1.0 / cBar * sqrt(max(0.5 * (magSqr(ULeft) + magSqr(URight)), SMALL)));
            const scalar chi = sqr(1.0 - MaArc);

            const scalar mDot =
                0.5 * (1.0 - g) *
                    (rhoLeft * contrVLeft + rhoRight * contrVRight - contrVBar * (rhoRight - rhoLeft)) -
                0.5 * chi / cBar * (pRight - pLeft);

            scalar betaLeft = 0.0;
            scalar betaRight = 0.0;

            if (mag(MaLeft) < 1)
            {
                betaLeft = 0.25 * (2.0 - MaLeft) * sqr(MaLeft + 1.0);
            }
            else
            {
                betaLeft = 0.5 * (1.0 + sign(MaLeft));
            }

            if (mag(MaRight) < 1)
            {
                betaRight = 0.25 * (2.0 + MaRight) * sqr(MaRight - 1.0);
            }
            else
            {
                betaRight = 0.5 * (1.0 + sign(-MaRight));
            }

            const scalar pTilde = 0.5 * (pLeft + pRight) + 0.5 * (betaLeft - betaRight) * (pLeft - pRight)
                                  //+ (1.0 - chi)*(betaLeft + betaRight -1.0)*0.5*(pLeft + pRight);
                                  + sqrt(max(0.5 * (magSqr(ULeft) + magSqr(URight)), SMALL)) * (betaLeft + betaRight - 1.0) * rhoBar * cBar;

            if (mDot > 0)
            {
                rhoFlux = mDot * magSf;
                rhoUFlux = (mDot * ULeft + normalVector * pTilde) * magSf;
                rhoEFlux = (mDot * hLeft) * magSf;
            }
            else
            {
                rhoFlux = mDot * magSf;
                rhoUFlux = (mDot * URight + normalVector * pTilde) * magSf;
                rhoEFlux = (mDot * hRight) * magSf;
            }
        }
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
