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

Application
    rhoGodunovTurbFoam

Description
    Density-based compressible transient solver using the flux difference
    splitting scheme of Bui. Primarily designed for LES.

\*---------------------------------------------------------------------------*/

#include "cwipiPstream.H"
#include "fvCFD.H"
#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
#include "numericFlux.H"
#include "pointMesh.H"
#include "pointFields.H"
#include "volPointInterpolation.H"
#include "acousticCourantNo.H"
#include "shockSensor.H"
#include "numericFluxes.H"
#include "cwipiFields.H"

int main(int argc, char *argv[])
{

    argList::addNote(
        "Density-based compressible transient solver using the flux difference"
        "splitting scheme of Bui. Primarily designed for LES.");

#include "postProcess.H"
#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"
#include "createFields.H"
#include "createTimeControls.H"

    // Validate turbulence
    turbulence->validate();

    // Runge-Kutta coefficient
    constexpr const std::array<scalar, 4> beta = {0.1100, 0.2766, 0.5, 1};

    // Create CWIPI coupling
    cwipiPstream coupling(
        runTime,
        mesh,
        thermo,
        U);

    Info << "Starting time loop" << endl;

    while (runTime.run())
    {
        // Send sources at correct time step
        coupling.send();

        // Execute main solver loop
#include "cwipiRhoGodunovTurbFoam.H"

        // Do I/O
        runTime.write();

        // Update time step of coupling
        coupling.updateTime();

        // Print execution time
        runTime.printExecutionTime(Info);
    }

    Info << "End" << endl;

    return 0;
}

// ************************************************************************* //
