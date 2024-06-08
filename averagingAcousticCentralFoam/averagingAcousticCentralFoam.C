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
    cwipiAcousticCentralFoam

Group
    grpCompressibleSolvers

Description
    Density-based compressible flow solver based on
    central-upwind schemes of Kurganov and Tadmor with
    support for mesh-motion and topology changes.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
#include "fixedRhoFvPatchScalarField.H"
#include "directionInterpolate.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "cwipiFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote(
        "Density-based compressible flow solver based on"
        " central-upwind schemes of Kurganov and Tadmor with"
        " support for mesh-motion and topology changes.");

#define NO_CONTROL
#include "postProcess.H"
#include "addCheckCaseOptions.H"
#include "setRootCaseLists.H"
#include "createTime.H"
#include "createMesh.H"
#include "createFields.H"
#include "createFieldRefs.H"
#include "createTimeControls.H"

    // Create CWIPI fields
    cwipiFields couplingFields(
        mesh,
        runTime,
        U,
        thermo);

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "readFluxScheme.H"

    const dimensionedScalar v_zero(dimVolume / dimTime, Zero);

    // Courant numbers used to adjust the time-step
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;

    Info << endl;
    Info << "Starting time loop" << endl;
    Info << endl;

    while (runTime.run())
    {
        // Update CWIPI fields
        couplingFields.update();

        // Execute main body of solver
#include "averagingAcousticCentralFoam.H"

        // Do I/O
        runTime.write();

        // Print execution time
        runTime.printExecutionTime(Info);
    }

    Info << "End" << endl;

    return 0;
}

// ************************************************************************* //
