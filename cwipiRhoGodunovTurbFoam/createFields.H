// Thermo fields
Info << "Reading thermophysical properties" << endl;
autoPtr<psiThermo> pThermo(psiThermo::New(mesh));
psiThermo &thermo = pThermo();
volScalarField &e = thermo.he();

// Momentum fields
Info << "Reading fields rho, U, rho * U, rho * E" << endl;
volVectorField U(
    IOobject(
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE),
    mesh);
volScalarField rho(
    IOobject(
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    thermo.rho());
volVectorField rhoU(
    IOobject(
        "rhoU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    thermo.rho() * U);
volScalarField rhoE(
    IOobject(
        "rhoE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    thermo.rho() * (e + 0.5 * magSqr(U)));

// Density gradient-based shock sensor
shockSensor<psiThermo> shockDetector(U, mesh, thermo, runTime);

// Numeric flux
autoPtr<basicNumericFlux> dbnsFluxPtr = basicNumericFlux::New(
    thermo.p(),
    U,
    thermo.T(),
    shockDetector.upwindingFactor(),
    thermo);
basicNumericFlux &dbnsFlux = dbnsFluxPtr();
const surfaceScalarField &phi = dbnsFlux.rhoFlux();

// Turbulence fields
Info << "Creating turbulence model" << endl;
autoPtr<compressible::turbulenceModel> turbulence(
    compressible::turbulenceModel::New(
        rho,
        U,
        phi,
        thermo));

// Acoustic damping
const volScalarField acousticDamping(
    IOobject(
        "acousticDamping",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE),
    mesh);
const volScalarField noDamping = 1 - acousticDamping;

// LES delta
// volScalarField LESdelta(
//     IOobject(
//         "LESdelta",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE),
//     mesh,
//     dimensionedScalar("LESdelta", dimLength, 0));

// if (mesh.objectRegistry::foundObject<volScalarField>("delta") == true)
// {
//     const volScalarField &delta = mesh.objectRegistry::lookupObject<volScalarField>("delta");

//     LESdelta.internalFieldRef() = delta.internalField();
//     LESdelta.boundaryFieldRef() = delta.boundaryField();
// }

// Acoustic Courant number CFL control
acousticCourantNo<psiThermo> Co(U, mesh, thermo, runTime);

// Get reference values for damping
const dimensionedScalar T_inf(
    "T_inf",
    thermo.T().dimensions(),
    readScalar(mesh.solutionDict().subDict("freestreamProperties").lookup("T")));
const dimensionedScalar p_inf(
    "p_inf",
    thermo.p().dimensions(),
    readScalar(mesh.solutionDict().subDict("freestreamProperties").lookup("p")));
const dimensionedScalar rho_inf(
    "rho_inf",
    rho.dimensions(),
    readScalar(mesh.solutionDict().subDict("freestreamProperties").lookup("rho")));
const dimensionedVector U_inf(
    "U_inf",
    U.dimensions(),
    vector(mesh.solutionDict().subDict("freestreamProperties").lookup("U")));
