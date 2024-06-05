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
    acousticSponge

Description
    Implementation of an absorption-based acoustic sponge

SourceFiles
    acousticSponge.H

\*---------------------------------------------------------------------------*/

#ifndef acousticSponge_H
#define acousticSponge_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    /*---------------------------------------------------------------------------*\
                               Class acousticSponge Declaration
    \*---------------------------------------------------------------------------*/

    // class acousticSponge
    // {
    // public:
    //     // Constructor
    //     inline acousticSponge(
    //         const fvMesh &mesh)
    //     {
    //         IOdictionary spongeProperties(
    //             IOobject(
    //                 "spongeProperties",
    //                 mesh.time().constant(),
    //                 mesh,
    //                 IOobject::MUST_READ,
    //                 IOobject::NO_WRITE));

    //         const scalar myScalar = readScalar(spongeProperties.lookup("myScalar"));
    //         Info << myScalar << endl;
    //     };

    //     // Get reference values for damping
    //     // const dimensionedScalar T_inf(
    //     //     "T_inf",
    //     //     thermo.T().dimensions(),
    //     //     readScalar(mesh.solutionDict().subDict("freestreamProperties").lookup("T")));
    //     // const dimensionedScalar p_inf(
    //     //     "p_inf",
    //     //     thermo.p().dimensions(),
    //     //     readScalar(mesh.solutionDict().subDict("freestreamProperties").lookup("p")));
    //     // const dimensionedScalar rho_inf(
    //     //     "rho_inf",
    //     //     rho.dimensions(),
    //     //     readScalar(mesh.solutionDict().subDict("freestreamProperties").lookup("rho")));
    //     // const dimensionedVector U_inf(
    //     //     "U_inf",
    //     //     U.dimensions(),
    //     //     vector(mesh.solutionDict().subDict("freestreamProperties").lookup("U")));
    //     // const dimensionedScalar e_inf(
    //     //     "h_inf",
    //     //     average(thermo.he(p_inf + (0 * thermo.p()), T_inf + (0 * thermo.T()))));

    //     // Destructor
    //     inline ~acousticSponge(){};

    // private:
    //     // volScalarField &f_;

    //     // const dimensionedScalar fInf_;

    //     // const volScalarField &acousticDamping_;
    // };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //