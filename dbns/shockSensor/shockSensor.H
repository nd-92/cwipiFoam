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
    shockSensor

Description
    Implementation of a density gradient based shock sensor

SourceFiles
    shockSensor.H

\*---------------------------------------------------------------------------*/

#ifndef shockSensor_H
#define shockSensor_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    /*---------------------------------------------------------------------------*\
                               Class shockSensor Declaration
    \*---------------------------------------------------------------------------*/

    template <class Thermo>
    class shockSensor
    {
    public:
        // Constructor
        inline shockSensor(
            const volVectorField &U,
            const fvMesh &mesh,
            const Thermo &thermo,
            const Time &runTime)
            : U_(U),
              mesh_(mesh),
              thermo_(thermo),
              runTime_(runTime),
              rhoGradRef_(dimensionedScalar(
                  "refVal",
                  dimMass / (dimLength * dimLength * dimLength * dimLength),
                  readScalar(mesh.solutionDict().subDict("freestreamProperties").lookup("rhoGrad")))),
              shockSensor_(
                  IOobject(
                      "buiEps",
                      runTime.timeName(),
                      mesh,
                      IOobject::NO_READ,
                      IOobject::NO_WRITE),
                  mag(fvc::grad(thermo_.rho())) / rhoGradRef_),
              USmall_(dimensionedScalar("USmall", dimLength / dimTime, SMALL)),
              cellScales_(cellScalesI(mesh))
        {
            shockSensor_ = erf(erf(mag(fvc::grad(thermo_.rho()) & (U_ / (mag(U_) + USmall_))) / rhoGradRef_) + cellScales_);
            shockSensor_ = shockSensor_ - min(shockSensor_);
            shockSensor_ = shockSensor_ / max(shockSensor_);
        };

        // Destructor
        inline ~shockSensor(){};

        // Update the shock sensor
        inline void update()
        {
            shockSensor_ = erf(erf(mag(fvc::grad(thermo_.rho()) & (U_ / (mag(U_) + USmall_))) / rhoGradRef_) + cellScales_);
            shockSensor_ = shockSensor_ - min(shockSensor_);
            shockSensor_ = shockSensor_ / max(shockSensor_);
        }

        // Public const access to the shock sensor
        [[nodiscard]] inline const volScalarField &upwindingFactor() const
        {
            return shockSensor_;
        }

    private:
        // Const access to velocity field
        const volVectorField &U_;

        // Const access to mesh
        const fvMesh &mesh_;

        // Const access to thermo
        const Thermo &thermo_;

        // Const access to runTime
        const Time &runTime_;

        // Density gradient clip value
        const dimensionedScalar rhoGradRef_;

        // Shock sensor
        volScalarField shockSensor_;

        // Small amount of velocity to help with direction vector
        const dimensionedScalar USmall_;

        // Cell scales
        const volScalarField cellScales_;

        // Initialise the cell dimensions as a volScalarField
        inline const volScalarField cellScalesI(const fvMesh &mesh) const
        {
            volScalarField cellScales_(
                IOobject(
                    "cellScales",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE),
                mesh,
                dimensionedScalar("one", dimless, 1.0));
            // cellScales_.primitiveFieldRef() = cbrt(this->mesh().V());
            const cellList &cells = mesh.cells();
            const vectorField &cellC = mesh.cellCentres();
            const vectorField &faceC = mesh.faceCentres();
            const vectorField faceN(mesh.faceAreas() / mag(mesh.faceAreas()));
            scalarField hmax(cells.size());
            forAll(cells, celli)
            {
                scalar deltaMaxTmp = 0.0;
                const labelList &cFaces = cells[celli];
                const point &cc = cellC[celli];

                forAll(cFaces, cFacei)
                {
                    label facei = cFaces[cFacei];
                    const point &fc = faceC[facei];
                    const vector &n = faceN[facei];

                    scalar tmp = mag(n & (fc - cc));
                    if (tmp > deltaMaxTmp)
                    {
                        deltaMaxTmp = tmp;
                    }
                }

                hmax[celli] = deltaMaxTmp;
            }

            cellScales_.primitiveFieldRef() = hmax;

            const dimensionedScalar cellScaleRef_ = dimensionedScalar(
                "cellScaleRef",
                dimless,
                readScalar(mesh.solutionDict().subDict("freestreamProperties").lookup("cellDim")));

            cellScales_ = erf(cellScales_ / cellScaleRef_);

            return cellScales_;
        }
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //