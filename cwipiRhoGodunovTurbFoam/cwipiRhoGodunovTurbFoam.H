#include "readTimeControls.H"
Co.calculate();
Co.setDeltaT(maxDeltaT, maxCo, adjustTimeStep);

runTime++;
Info << "\nTime = " << runTime.value() << "\n";

for (size_t i = 0; i < 4; i++)
{
    //  Solve the approximate Riemann problem for this time step
    dbnsFlux.computeFlux();

    // Time integration
    solve                                  //
        (                                  //
            1.0 / beta[i] * fvm::ddt(rho)  //
            + fvc::div(dbnsFlux.rhoFlux()) //
        );                                 //

    solve                                        //
        (                                        //
            1.0 / beta[i] * fvm::ddt(rhoU)       //
            + fvc::div(dbnsFlux.rhoUFlux())      //
            + fvc::div(turbulence->devRhoReff()) //
        );                                       //

    solve                                                        //
        (                                                        //
            1.0 / beta[i] * fvm::ddt(rhoE)                       //
            + fvc::div(dbnsFlux.rhoEFlux())                      //
            + fvc::div(turbulence->devRhoReff() & rhoU / rho)    //
            - fvc::laplacian(turbulence->alphaEff(), rhoE / rho) //
        );                                                       //

    // Compute U
    U = rhoU / rho;
    U.correctBoundaryConditions();
    rhoU.boundaryFieldRef() == rho.boundaryField() * U.boundaryField();

    // Calculate internal energy
    e = ((rhoE / rho) - (0.5 * magSqr(U)));
    e.correctBoundaryConditions();

    // Correct thermodynamics
    thermo.correct();
    rhoE.boundaryFieldRef() == rho.boundaryField() * (e.boundaryField() + 0.5 * magSqr(U.boundaryField()));

    // Compute p from rho
    thermo.p().ref() = rho() / thermo.psi();
    thermo.p().correctBoundaryConditions();
    // rho.boundaryFieldRef() == thermo.psi().boundaryField() * thermo.p().boundaryField();

    // Apply damping
    thermo.T() = (T_inf * acousticDamping) + (thermo.T() * noDamping);
    thermo.T().correctBoundaryConditions();
    thermo.p() = (p_inf * acousticDamping) + (thermo.p() * noDamping);
    thermo.p().correctBoundaryConditions();
    rho = (rho_inf * acousticDamping) + (rho * noDamping);
    rho.correctBoundaryConditions();
    U = (U_inf * acousticDamping) + (U * noDamping);
    U.correctBoundaryConditions();
}

// Update turbulence fields
turbulence->correct();

// Update the shock sensor
shockDetector.update();
