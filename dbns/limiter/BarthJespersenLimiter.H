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
    BarthJespersenLimiter

Description
    Barth-Jespersen limiter

Author
    Aleksandar Jemcov
    Rewrite by Hrvoje Jasak

\*---------------------------------------------------------------------------*/

#ifndef BarthJespersenLimiter_H
#define BarthJespersenLimiter_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------*\
                        Class BarthJespersenLimiter Declaration
    \*---------------------------------------------------------------------------*/

    class BarthJespersenLimiter
    {
    public:
        // Constructors

        //- Construct null
        inline BarthJespersenLimiter()
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
            const scalar &extrapolate)
        {
            if (mag(extrapolate) < SMALL)
            {
                return;
            }
            else if (extrapolate - deltaOneMax > SMALL)
            {
                lim = min(lim, deltaOneMax / extrapolate);
            }
            else if (extrapolate - deltaOneMin < -SMALL)
            {
                lim = min(lim, deltaOneMin / extrapolate);
            }
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
            for (direction cmpt = 0; cmpt < Type::nComponents; cmpt++)
            {
                limiter(
                    lim.component(cmpt),
                    cellVolume,
                    deltaOneMax.component(cmpt),
                    deltaOneMin.component(cmpt),
                    extrapolate.component(cmpt));
            }
        }
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
