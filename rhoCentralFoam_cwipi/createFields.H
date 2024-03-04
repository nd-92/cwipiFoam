#include "createRDeltaT.H"

Info << "Reading thermophysical properties\n"
     << endl;

autoPtr<psiThermo> pThermo(
    psiThermo::New(mesh));
psiThermo &thermo = pThermo();

volScalarField &e = thermo.he();

Info << "Reading field U\n"
     << endl;
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

const surfaceScalarField pos(
    IOobject(
        "pos",
        runTime.timeName(),
        mesh),
    mesh,
    dimensionedScalar("pos", dimless, 1.0));

const surfaceScalarField neg(
    IOobject(
        "neg",
        runTime.timeName(),
        mesh),
    mesh,
    dimensionedScalar("neg", dimless, -1.0));

surfaceScalarField phi("phi", fvc::flux(rhoU));

Info << "Creating turbulence model\n"
     << endl;
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
