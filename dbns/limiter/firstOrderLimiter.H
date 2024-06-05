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
    firstOrderLimiter

Description
    First order limiter: all second order terms are removed

Author
    Hrvoje Jasak

\*---------------------------------------------------------------------------*/

#ifndef firstOrderLimiter_H
#define firstOrderLimiter_H

#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------*\
                          Class firstOrderLimiter Declaration
    \*---------------------------------------------------------------------------*/

    class firstOrderLimiter
    {
    public:
        // Constructor

        //- Construct null
        inline firstOrderLimiter()
        {
        }

        // Destructor - default

        // Member functions

        //- Set scalar limiter value
        inline void limiter(
            scalar &lim,
            const scalar &cellVolume,
            const scalar &deltaOneMax,
            const scalar &deltaOneMin,
            const scalar &deltaTwo)
        {
            lim = 0;
        }

        //- Set Type limiter
        template <class Type>
        inline void limiter(
            Type &lim,
            const scalar &cellVolume,
            const Type &deltaOneMax,
            const Type &deltaOneMin,
            const Type &extrapolate)
        {
            lim = pTraits<Type>::zero;
        }
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
