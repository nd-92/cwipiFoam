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
    rungeKutta

Description
    Implementation of constexpr Runge-Kutta time integration coefficients

SourceFiles
    rungeKutta.H

\*---------------------------------------------------------------------------*/

#ifndef rungeKutta_H
#define rungeKutta_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------*\
                               Class rungeKutta Declaration
    \*---------------------------------------------------------------------------*/

// Macro for Runge-Kutta time integration
#define rkLoop(i) for (uint8_t i = 0; i < 4; i++)

    class rungeKutta
    {
    public:
        // Constructor
        constexpr rungeKutta(){};

        // Destructor
        // constexpr ~rungeKutta(){};

        // Return constexpr Runge-Kutta coefficient
        constexpr const std::array<scalar, 4> &beta() const
        {
            return beta_;
        }

    private:
        // Runge-Kutta coefficients
        static constexpr const std::array<scalar, 4> beta_ = {0.1100, 0.2766, 0.500, 1.000};
        // static constexpr const scalar beta_[4] = {0.1100, 0.2766, 0.500, 1.000};
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
