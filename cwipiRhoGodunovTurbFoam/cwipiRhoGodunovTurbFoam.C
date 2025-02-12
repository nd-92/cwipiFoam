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

    // Do the coupling if defined
    // Must not be averaging simultaneously
    if (runTime.controlDict().lookupOrDefault("cwipiSwitch", false) && !(runTime.controlDict().lookupOrDefault("cwipiAveraging", false)))
    {
        // Create CWIPI coupling
        cwipiPstream coupling(runTime, mesh, thermo, U);

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
    }

    // Do the averaging if defined
    // Must not be coupling simultaneously
    if (runTime.controlDict().lookupOrDefault("cwipiAveraging", false) && !(runTime.controlDict().lookupOrDefault("cwipiSwitch", false)))
    {
        // Create CWIPI fields
        cwipiFields couplingFields(mesh, runTime, U, thermo);

        Info << "Starting time loop" << endl;

        while (runTime.run())
        {
            // Update CWIPI fields
            couplingFields.update();

            // Execute main solver loop
#include "averagingRhoGodunovTurbFoam.H"

            // Write runtime output
            runTime.write();
            runTime.printExecutionTime(Info);
        }
    }

    // Otherwise neither are defined, so run as normal and compute the source terms
    if (!(runTime.controlDict().lookupOrDefault("cwipiSwitch", false)) && !(runTime.controlDict().lookupOrDefault("cwipiAveraging", false)))
    {
        cwipiFields couplingFields(mesh, runTime, U, thermo);
        const cwipiMeanFields baseFlow(mesh, runTime);

        volScalarField DsDt(
            IOobject(
                "entropyMonopole",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE),
            baseFlow.cSqMean() * (baseFlow.rhoMean() / thermo.Cp()) * (fvc::ddt(couplingFields.s()) + (baseFlow.UMean() & fvc::grad((couplingFields.s() - baseFlow.sMean())))));

        volVectorField TGrads(
            IOobject(
                "entropyGradient",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE),
            ((thermo.T() - baseFlow.TMean()) * fvc::grad(baseFlow.sMean())) - ((couplingFields.s() - baseFlow.sMean()) * fvc::grad(baseFlow.TMean())));

        volVectorField LPrime(
            IOobject(
                "LPrime",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE),
            couplingFields.L() - baseFlow.LMean());

        const volScalarField pMean(
            IOobject(
                "pMean",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE),
            mesh);

        Info << "Starting time loop" << endl;

        while (runTime.run())
        {
            // Execute main solver loop
#include "rhoGodunovTurbFoam.H"

            // Update CWIPI fields
            couplingFields.update();

            DsDt = baseFlow.cSqMean() * (baseFlow.rhoMean() / thermo.Cp()) * (fvc::ddt(couplingFields.s()) + (baseFlow.UMean() & fvc::grad((couplingFields.s() - baseFlow.sMean()))));
            TGrads = ((thermo.T() - baseFlow.TMean()) * fvc::grad(baseFlow.sMean())) - ((couplingFields.s() - baseFlow.sMean()) * fvc::grad(baseFlow.TMean()));
            LPrime = couplingFields.L() - baseFlow.LMean();

            // Write runtime output
            runTime.write();
            runTime.printExecutionTime(Info);
        }
    }

    Info << "End" << endl;

    return 0;
}

// ************************************************************************* //
