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

\*---------------------------------------------------------------------------*/

#ifndef __CWIPI_FIELDS_H
#define __CWIPI_FIELDS_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    class cwipiFields
    {
    public:
        // Default constructor
        inline cwipiFields(
            const fvMesh &mesh,
            const Time &runTime,
            const volVectorField &U,
            const psiThermo &thermo)
            : U_(U),
              thermo_(thermo),
              L_(                                                                                                                                                                                 //
                  IOobject(                                                                                                                                                                       //
                      "L",                                                                                                                                                                        //
                      runTime.timeName(),                                                                                                                                                         //
                      mesh,                                                                                                                                                                       //
                      IOobject::NO_READ,                                                                                                                                                          //
                      IOobject::NO_WRITE),                                                                                                                                                        //
                  fvc::curl(U_) ^ U_),                                                                                                                                                            //
              s_(                                                                                                                                                                                 //
                  IOobject(                                                                                                                                                                       //
                      "s",                                                                                                                                                                        //
                      runTime.timeName(),                                                                                                                                                         //
                      mesh,                                                                                                                                                                       //
                      IOobject::NO_READ,                                                                                                                                                          //
                      IOobject::NO_WRITE),                                                                                                                                                        //
                  (thermo_.Cp() * log(thermo_.T() / Foam::constant::standard::Tstd)) - ((Foam::constant::physicoChemical::R / thermo_.W()) * log(thermo_.p() / Foam::constant::standard::Pstd))), //
              c_(                                                                                                                                                                                 //
                  IOobject(                                                                                                                                                                       //
                      "c",                                                                                                                                                                        //
                      runTime.timeName(),                                                                                                                                                         //
                      mesh,                                                                                                                                                                       //
                      IOobject::NO_READ,                                                                                                                                                          //
                      IOobject::NO_WRITE),                                                                                                                                                        //
                  Foam::sqrt(thermo_.gamma() * thermo_.p() / thermo_.rho()))
        {
            Info << "CWIPI source fields allocated." << endl;
        };

        // Default destructor
        inline ~cwipiFields(){};

        // Update Lamb vector, entropy and speed of sound
        inline void update()
        {
            // Update Lamb vector
            L_ = fvc::curl(U_) ^ U_;

            // Update entropy
            s_ = (thermo_.Cp() *                                       //
                  log(thermo_.T() / Foam::constant::standard::Tstd)) - //
                 ((Foam::constant::physicoChemical::R / thermo_.W()) * //
                  log(thermo_.p() / Foam::constant::standard::Pstd));  //

            // Update speed of sound
            c_ = Foam::sqrt(thermo_.gamma() * thermo_.p() / thermo_.rho());
            // c_ = sqrt(thermo_.Cp() / thermo_.Cv() / thermo_.psi());
        };

        // Access to Lamb vector by reference
        inline volVectorField const &L() const
        {
            return L_;
        };

        // Access to entropy by reference
        inline volScalarField const &s() const
        {
            return s_;
        };

        // Access to speed of sound by reference
        inline volScalarField const &c() const
        {
            return c_;
        };

    private:
        // Const reference to velocity
        const volVectorField &U_;

        // Const reference to thermo
        const psiThermo &thermo_;

        // Lamb vector
        volVectorField L_;

        // Entropy
        volScalarField s_;

        // Speed of sound
        volScalarField c_;
    };

}

#endif

// ************************************************************************* //