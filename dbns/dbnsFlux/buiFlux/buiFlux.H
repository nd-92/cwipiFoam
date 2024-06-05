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

#ifndef buiFlux_H
#define buiFlux_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------*\
                               Class slau2Flux Declaration
    \*---------------------------------------------------------------------------*/

    class buiFlux
    {
    public:
        // Constructors

        //- Construct null
        inline buiFlux()
        {
        }

        // Destructor - default

        // Member Functions

        //- Update fluxes based on current left and right states
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
            const scalar rhoLeft = pLeft / (RLeft * TLeft);
            const scalar rhoRight = pRight / (RRight * TRight);

            // Adiabatic exponent is constant for ideal gas but if Cp=Cp(T)
            // it must be computed for each cell and evaluated at each face
            // through reconstruction
            const scalar kappaLeft = (CvLeft + RLeft) / CvLeft;
            const scalar kappaRight = (CvRight + RRight) / CvRight;

            // normal vector
            const vector normalVector = Sf / magSf;

            // Compute left and right total enthalpies:
            const scalar hLeft = (CvLeft * TLeft + 0.5 * magSqr(ULeft)) + pLeft / rhoLeft;
            const scalar hRight = (CvRight * TRight + 0.5 * magSqr(URight)) + pRight / rhoRight;

            // Step 2: compute Roe averged quantities for face:
            const scalar rhoTilde = sqrt(max(rhoLeft * rhoRight, SMALL));

            // Some temporary variables:
            const scalar rhoLeftSqrt = sqrt(max(rhoLeft, SMALL));
            const scalar rhoRightSqrt = sqrt(max(rhoRight, SMALL));

            const scalar wLeft = rhoLeftSqrt / (rhoLeftSqrt + rhoRightSqrt);
            const scalar wRight = 1 - wLeft;

            const vector UTilde = ULeft * wLeft + URight * wRight;
            const scalar hTilde = hLeft * wLeft + hRight * wRight;
            const scalar qTildeSquare = magSqr(UTilde);
            const scalar kappaTilde = kappaLeft * wLeft + kappaRight * wRight;

            // Speed of sound
            const scalar cTilde = sqrt(max((kappaTilde - 1) * (hTilde - 0.5 * qTildeSquare), SMALL));

            // Roe averaged contravariant velocity
            const scalar contrVTilde = (UTilde & normalVector);

            // Step 3: compute primitive differences:
            const scalar deltaP = pRight - pLeft;
            // const scalar deltaRho = rhoRight - rhoLeft;
            const vector deltaU = URight - ULeft;
            const scalar deltaContrV = (deltaU & normalVector);

            // Step 4: compute wave strengths:

            // Roe and Pike - formulation
            const scalar r1 = (deltaP - rhoTilde * cTilde * deltaContrV) / (2.0 * sqr(cTilde));
            const scalar r2 = (rhoRight - rhoLeft) - deltaP / sqr(cTilde);
            const scalar r3 = (deltaP + rhoTilde * cTilde * deltaContrV) / (2.0 * sqr(cTilde));

            scalar lambda1 = mag(contrVTilde - cTilde);
            scalar lambda2 = mag(contrVTilde);
            scalar lambda3 = mag(contrVTilde + cTilde);

            const scalar UL = ULeft & normalVector;
            const scalar UR = URight & normalVector;
            const scalar cLeft = sqrt(max((kappaLeft - 1) * (hLeft - 0.5 * magSqr(ULeft)), SMALL));
            const scalar cRight = sqrt(max((kappaRight - 1) * (hRight - 0.5 * magSqr(URight)), SMALL));

            // First eigenvalue: U - c
            scalar eps = 2 * max(0, (UR - cRight) - (UL - cLeft));
            if (lambda1 < eps)
            {
                lambda1 = (sqr(lambda1) + sqr(eps)) / (2.0 * eps);
            }

            // Second eigenvalue: U
            eps = 2 * max(0, UR - UL);
            if (lambda2 < eps)
            {
                lambda2 = (sqr(lambda2) + sqr(eps)) / (2.0 * eps);
            }

            // Third eigenvalue: U + c
            eps = 2 * max(0, (UR + cRight) - (UL + cLeft));
            if (lambda3 < eps)
            {
                lambda3 = (sqr(lambda3) + sqr(eps)) / (2.0 * eps);
            }

            const scalar fluxLeft11 = rhoLeft * (ULeft & normalVector);
            const scalar fluxRight11 = rhoRight * (URight & normalVector);

            rhoFlux = (0.5 * ((fluxLeft11) + (fluxRight11) - (buiEps * ((lambda1 * r1) + (lambda2 * r2) + (lambda3 * r3))))) * magSf;
            rhoUFlux = (0.5 * ((ULeft * fluxLeft11 + normalVector * pLeft) + (URight * fluxRight11 + normalVector * pRight) - (buiEps * ((lambda1 * r1 * (UTilde - cTilde * normalVector)) + (lambda2 * (r2 * UTilde + rhoTilde * (deltaU - deltaContrV * normalVector))) + (lambda3 * r3 * (UTilde + cTilde * normalVector)))))) * magSf;
            rhoEFlux = (0.5 * ((hLeft * fluxLeft11) + (hRight * fluxRight11) - (buiEps * ((lambda1 * r1 * (hTilde - cTilde * contrVTilde)) + (lambda2 * (r2 * (0.5 * qTildeSquare) + rhoTilde * ((UTilde & deltaU) - contrVTilde * deltaContrV))) + (lambda3 * r3 * (hTilde + cTilde * contrVTilde)))))) * magSf;
        }
        // {
        //     // Step 1: decode rho left and right:
        //     const scalar rhoLeft = pLeft / (RLeft * TLeft);
        //     const scalar rhoRight = pRight / (RRight * TRight);

        //     // Decode left and right total energy:
        //     const scalar eLeft = CvLeft * TLeft + 0.5 * magSqr(ULeft);
        //     const scalar eRight = CvRight * TRight + 0.5 * magSqr(URight);

        //     // Adiabatic exponent is constant for ideal gas but if Cp=Cp(T)
        //     // it must be computed for each cell and evaluated at each face
        //     // through reconstruction
        //     const scalar kappaLeft = (CvLeft + RLeft) / CvLeft;
        //     const scalar kappaRight = (CvRight + RRight) / CvRight;

        //     // normal vector
        //     const vector normalVector = Sf / magSf;

        //     // Compute left and right contravariant velocities:
        //     const scalar contrVLeft = (ULeft & normalVector);
        //     const scalar contrVRight = (URight & normalVector);

        //     // Compute left and right total enthalpies:
        //     const scalar hLeft = eLeft + pLeft / rhoLeft;
        //     const scalar hRight = eRight + pRight / rhoRight;

        //     // Step 2: compute Roe averged quantities for face:
        //     const scalar rhoTilde = sqrt(max(rhoLeft * rhoRight, SMALL));

        //     // Some temporary variables:
        //     const scalar rhoLeftSqrt = sqrt(max(rhoLeft, SMALL));
        //     const scalar rhoRightSqrt = sqrt(max(rhoRight, SMALL));

        //     const scalar wLeft = rhoLeftSqrt / (rhoLeftSqrt + rhoRightSqrt);
        //     const scalar wRight = 1 - wLeft;

        //     const vector UTilde = ULeft * wLeft + URight * wRight;
        //     const scalar hTilde = hLeft * wLeft + hRight * wRight;
        //     const scalar qTildeSquare = magSqr(UTilde);
        //     const scalar kappaTilde = kappaLeft * wLeft + kappaRight * wRight;

        //     // Speed of sound
        //     const scalar cTilde = sqrt(max((kappaTilde - 1) * (hTilde - 0.5 * qTildeSquare), SMALL));

        //     // Roe averaged contravariant velocity
        //     const scalar contrVTilde = (UTilde & normalVector);

        //     // Step 3: compute primitive differences:
        //     const scalar deltaP = pRight - pLeft;
        //     const scalar deltaRho = rhoRight - rhoLeft;
        //     const vector deltaU = URight - ULeft;
        //     const scalar deltaContrV = (deltaU & normalVector);

        //     // Step 4: compute wave strengths:

        //     // Roe and Pike - formulation
        //     const scalar r1 = (deltaP - rhoTilde * cTilde * deltaContrV) / (2.0 * sqr(cTilde));
        //     const scalar r2 = deltaRho - deltaP / sqr(cTilde);
        //     const scalar r3 = (deltaP + rhoTilde * cTilde * deltaContrV) / (2.0 * sqr(cTilde));

        //     // Step 5: compute l vectors

        //     // rho row:
        //     const scalar l1rho = 1;
        //     const scalar l2rho = 1;
        //     const scalar l3rho = 0;
        //     const scalar l4rho = 1;

        //     // first U column
        //     const vector l1U = UTilde - cTilde * normalVector;

        //     // second U column
        //     const vector l2U = UTilde;

        //     // third U column
        //     const vector l3U = deltaU - deltaContrV * normalVector;

        //     // fourth U column
        //     const vector l4U = UTilde + cTilde * normalVector;

        //     // E row
        //     const scalar l1e = hTilde - cTilde * contrVTilde;
        //     const scalar l2e = 0.5 * qTildeSquare;
        //     const scalar l3e = (UTilde & deltaU) - contrVTilde * deltaContrV;
        //     const scalar l4e = hTilde + cTilde * contrVTilde;

        //     // Step 6: compute eigenvalues

        //     // derived from algebra by hand, only for Euler equation usefull
        //     scalar lambda1 = mag(contrVTilde - cTilde);
        //     scalar lambda2 = mag(contrVTilde);
        //     scalar lambda3 = mag(contrVTilde + cTilde);

        //     // Step 7: check for Harten entropy correction

        //     //     const scalar eps = 0.1*cTilde; //adjustable parameter

        //     //     if (lambda1 < eps || lambda2 < eps || lambda3 < eps)
        //     //     {
        //     //         lambda1 = (sqr(lambda1) + sqr(eps))/(2.0*eps);
        //     //         lambda2 = (sqr(lambda2) + sqr(eps))/(2.0*eps);
        //     //         lambda3 = (sqr(lambda3) + sqr(eps))/(2.0*eps);
        //     //     }

        //     // Step 7a: Alternative entropy correction: Felipe Portela, 9/Oct/2013

        //     const scalar UL = ULeft & normalVector;
        //     const scalar UR = URight & normalVector;
        //     const scalar cLeft = sqrt(max((kappaLeft - 1) * (hLeft - 0.5 * magSqr(ULeft)), SMALL));

        //     const scalar cRight = sqrt(max((kappaRight - 1) * (hRight - 0.5 * magSqr(URight)), SMALL));

        //     // First eigenvalue: U - c
        //     scalar eps = 2 * max(0, (UR - cRight) - (UL - cLeft));
        //     if (lambda1 < eps)
        //     {
        //         lambda1 = (sqr(lambda1) + sqr(eps)) / (2.0 * eps);
        //     }

        //     // Second eigenvalue: U
        //     eps = 2 * max(0, UR - UL);
        //     if (lambda2 < eps)
        //     {
        //         lambda2 = (sqr(lambda2) + sqr(eps)) / (2.0 * eps);
        //     }

        //     // Third eigenvalue: U + c
        //     eps = 2 * max(0, (UR + cRight) - (UL + cLeft));
        //     if (lambda3 < eps)
        //     {
        //         lambda3 = (sqr(lambda3) + sqr(eps)) / (2.0 * eps);
        //     }

        //     // Step 8: Compute flux differences

        //     // Components of deltaF1
        //     const scalar diffF11 = lambda1 * r1 * l1rho;
        //     const vector diffF124 = lambda1 * r1 * l1U;
        //     const scalar diffF15 = lambda1 * r1 * l1e;

        //     // Components of deltaF2
        //     const scalar diffF21 = lambda2 * (r2 * l2rho + rhoTilde * l3rho);
        //     const vector diffF224 = lambda2 * (r2 * l2U + rhoTilde * l3U);
        //     const scalar diffF25 = lambda2 * (r2 * l2e + rhoTilde * l3e);

        //     // Components of deltaF3
        //     const scalar diffF31 = lambda3 * r3 * l4rho;
        //     const vector diffF324 = lambda3 * r3 * l4U;
        //     const scalar diffF35 = lambda3 * r3 * l4e;

        //     // Step 9: compute left and right fluxes

        //     // Left flux 5-vector
        //     const scalar fluxLeft11 = rhoLeft * contrVLeft;
        //     const vector fluxLeft124 = ULeft * fluxLeft11 + normalVector * pLeft;
        //     const scalar fluxLeft15 = hLeft * fluxLeft11;

        //     // Right flux 5-vector
        //     const scalar fluxRight11 = rhoRight * contrVRight;
        //     const vector fluxRight124 = URight * fluxRight11 + normalVector * pRight;
        //     const scalar fluxRight15 = hRight * fluxRight11;

        //     // Step 10: compute face flux 5-vector
        //     rhoFlux = (0.5 * (fluxLeft11 + fluxRight11 - (buiEps * (diffF11 + diffF21 + diffF31)))) * magSf;
        //     rhoUFlux = (0.5 * (fluxLeft124 + fluxRight124 - (buiEps * (diffF124 + diffF224 + diffF324)))) * magSf;
        //     rhoEFlux = (0.5 * (fluxLeft15 + fluxRight15 - (buiEps * (diffF15 + diffF25 + diffF35)))) * magSf;
        // }
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
