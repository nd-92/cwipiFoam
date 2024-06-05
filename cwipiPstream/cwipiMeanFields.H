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

#ifndef __CWIPI_MEANFIELDS_H
#define __CWIPI_MEANFIELDS_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    class cwipiMeanFields
    {
    public:
        // Default constructor
        inline cwipiMeanFields(
            const fvMesh &mesh,
            const Foam::Time &runTime)
            : UMean_(IOobject("UMean", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh),     // Time-averaged velocity
              rhoMean_(IOobject("rhoMean", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh), // Time-averaged density
              LMean_(IOobject("LMean", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh),     // Time-averaged Lamb vector
              sMean_(IOobject("sMean", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh),     // Time-averaged entropy
              cMean_(IOobject("cMean", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh),     // Time-averaged speed of sound
              TMean_(IOobject("TMean", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh){};

        // Default destructor
        inline ~cwipiMeanFields(){};

        // Access to speed of sound by reference
        inline volVectorField const &UMean() const
        {
            return UMean_;
        };

        // Access to speed of sound by reference
        inline volScalarField const &rhoMean() const
        {
            return rhoMean_;
        };

        // Access to speed of sound by reference
        inline volVectorField const &LMean() const
        {
            return LMean_;
        };

        // Access to speed of sound by reference
        inline volScalarField const &sMean() const
        {
            return sMean_;
        };

        // Access to speed of sound by reference
        inline volScalarField const &cMean() const
        {
            return cMean_;
        };

        // Access to speed of sound by reference
        inline volScalarField const &TMean() const
        {
            return TMean_;
        };

    private:
        // Time-averaged velocity
        const volVectorField UMean_;

        // Time-averaged density
        const volScalarField rhoMean_;

        // Time-averaged Lamb vector
        const volVectorField LMean_;

        // Time-averaged entropy
        const volScalarField sMean_;

        // Time-averaged speed of sound
        const volScalarField cMean_;

        // Time-averaged temperature
        const volScalarField TMean_;
    };

}

#endif

// ************************************************************************* //