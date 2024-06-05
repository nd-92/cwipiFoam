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
    acousticCourantNo

Description
    Implementation of the acoustic Courant number

SourceFiles
    acousticCourantNo.H

\*---------------------------------------------------------------------------*/

#ifndef acousticCourantNo_H
#define acousticCourantNo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------*\
                            Class acousticCourantNo Declaration
    \*---------------------------------------------------------------------------*/

    template <class Thermo>
    class acousticCourantNo
    {
    public:
        // Constructor
        inline acousticCourantNo(
            const volVectorField &U,
            const fvMesh &mesh,
            const Thermo &thermo,
            Time &runTime)
            : U_(U),
              mesh_(mesh),
              thermo_(thermo),
              runTime_(runTime),
              amaxSf_(0.5 * (mag((fvc::interpolate(U) & mesh.Sf())) + (fvc::interpolate(sqrt(thermo.gamma() * (thermo.Cp() - thermo.Cv()) * thermo.T())) * mesh.magSf()))),
              sumAmaxSf_(fvc::surfaceSum(amaxSf_)().primitiveField()),
              CoNum_(0),
              meanCoNum_(0)
        {
            // Only need to evaluate this once at construction time
            // Avoid evaluation every time update is called
            if (mesh_.nInternalFaces() == false)
            {
                throw std::invalid_argument("No internal mesh faces.");
            }
        };

        // Destructor
        inline ~acousticCourantNo(){};

        // Calculate and print acoustic Courant number
        inline void calculate()
        {
            amaxSf_ = 0.5 * (mag((fvc::interpolate(U_) & mesh_.Sf())) + (fvc::interpolate(sqrt(thermo_.gamma() * (thermo_.Cp() - thermo_.Cv()) * thermo_.T())) * mesh_.magSf()));
            sumAmaxSf_ = fvc::surfaceSum(amaxSf_)().primitiveField();
            CoNum_ = 0.5 * gMax(sumAmaxSf_ / mesh_.V().field()) * runTime_.deltaTValue();
            meanCoNum_ = 0.5 * (gSum(sumAmaxSf_) / gSum(mesh_.V().field())) * runTime_.deltaTValue();

            // Output
            Info << "Mean and max acoustic Courant Numbers = " << meanCoNum_ << " " << CoNum_ << endl;
        };

        // Update the time step
        inline void setDeltaT(
            const scalar maxDeltaT,
            const scalar maxCo,
            const bool adjustTimeStep)
        {
            if (adjustTimeStep)
            {
                const scalar maxDeltaTFact = maxCo / (CoNum_ + SMALL);
                const scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1 * maxDeltaTFact), 1.2);

                runTime_.setDeltaT(
                    min(
                        deltaTFact * runTime_.deltaTValue(),
                        maxDeltaT));

                Info << "deltaT = " << runTime_.deltaTValue() << endl;
            }
        }

    private:
        // Const access to velocity field
        const volVectorField &U_;

        // Const access to mesh
        const fvMesh &mesh_;

        // Const access to thermo
        const Thermo &thermo_;

        // Non-const access to runTime
        Time &runTime_;

        surfaceScalarField amaxSf_;
        scalarField sumAmaxSf_;

        // Max Courant number
        scalar CoNum_;

        // Mean Courant number
        scalar meanCoNum_;
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
