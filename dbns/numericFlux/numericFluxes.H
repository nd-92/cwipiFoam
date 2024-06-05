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

#include "makeBasicNumericFlux.H"
#include "buiFlux.H"
#include "slauFlux.H"
#include "slau2Flux.H"
#include "firstOrderLimiter.H"
#include "BarthJespersenLimiter.H"
#include "VenkatakrishnanLimiter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /* * * * * * * * * * * * * * * Private Static Data * * * * * * * * * * * * * */

#define makeBasicNumericFluxForAllLimiters(Flux)        \
    makeBasicNumericFlux(Flux, firstOrderLimiter);      \
    makeBasicNumericFlux(Flux, BarthJespersenLimiter);  \
    makeBasicNumericFlux(Flux, VenkatakrishnanLimiter); \
    makeBasicNumericFlux(Flux, WangLimiter);            \
    makeBasicNumericFlux(Flux, MichalakGoochLimiter);

    makeBasicNumericFluxForAllLimiters(buiFlux);
    makeBasicNumericFluxForAllLimiters(slauFlux);
    makeBasicNumericFluxForAllLimiters(slau2Flux);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //