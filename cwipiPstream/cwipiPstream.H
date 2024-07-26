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

\*---------------------------------------------------------------------------*/

#ifndef __CWIPI_PSTREAM_H
#define __CWIPI_PSTREAM_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fvCFD.H"
#include "volPointInterpolation.H"
#include "volFields.H"
#include "psiThermo.H"
#include "cwipiSwitch.H"
#include "cwipiFields.H"
#include "cwipiMeanFields.H"
#include "fixedValueFvPatchFields.H"

namespace Foam
{
    // Define application names
    constexpr const char *LOCAL_APPLICATION_NAME = "cwipiFoam";
    constexpr const char *DISTANT_APPLICATION_NAME = "FOAM_APE";

    class cwipiPstream
    {
    public:
        // Default constructor
        cwipiPstream(
            const Time &runTime,
            const fvMesh &mesh,
            const psiThermo &thermo,
            const volVectorField &U);

        // Default destructor
        ~cwipiPstream();

        // Send fields
        void send();

        // Advance time step
        void updateTime();

        // Public access to speed of sound
        inline volScalarField const &c() const
        {
            return fields_.c();
        }

        // Evaluate send condition
        [[nodiscard]] inline constexpr const Foam::Switch sendNow() const
        {
            return (cwipiTimeStep == cwipiStep);
        };

    private:
        // Return the names of the source fields
        [[nodiscard]] inline const char *cwipiSourceFieldNames(const size_t dim)
        {
            switch (dim)
            {
            case 2:
                Info << "Coupling enabled with 2 physical dimensions" << endl;
                return "F_0_u,F_0_v,F_0_p";
                break;
            case 3:
                Info << "Coupling enabled with 3 physical dimensions" << endl;
                return "F_0_u,F_0_v,F_0_w,F_0_p";
                break;
            default:
                throw std::invalid_argument("Variable cwipiDim should be 2 or 3.");
                break;
            }
        };

        [[nodiscard]] inline const volScalarField initialiseSmoothingCoefficient(
            const fvMesh &__mesh,
            const Time &__runTime) const
        {
            // Initialise the field to zero
            volScalarField DT(
                IOobject(
                    "DT",
                    __runTime.timeName(),
                    __mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE),
                __mesh,
                dimensionedScalar(
                    "zero",
                    dimensionSet(0, 2, -1, 0, 0),
                    0));

            // Scale it via volume ^ (2/3) * dt
            forAll(__mesh.cells(), i)
            {
                DT[i] = pow(__mesh.V()[i], 2 / 3) * __runTime.deltaTValue();
            }

            // Multiply by a dimensionless scale factor
            return DT * dimensionedScalar("zero", dimless, __runTime.controlDict().lookupOrDefault("smoothingScale", 1));
        };

        // Smoothen source fields
        inline void smoothenSources()
        {
            // Smoothen the continuity sources
            pressureSmoother_ == F_p_;
            pressureSmoother_.correctBoundaryConditions();
            pressureSmoother_.oldTime() = F_p_;
            pressureSmoother_.oldTime().correctBoundaryConditions();
            solve(fvm::ddt(pressureSmoother_) - fvm::laplacian(DT_, pressureSmoother_));
            F_p_ = pressureSmoother_;
            F_p_.correctBoundaryConditions();

            // Smoothen the momentum sources
            velocitySmoother_ == F_u_;
            velocitySmoother_.correctBoundaryConditions();
            velocitySmoother_.oldTime() = F_u_;
            velocitySmoother_.oldTime().correctBoundaryConditions();
            solve(fvm::ddt(velocitySmoother_) - fvm::laplacian(DT_, velocitySmoother_));
            F_u_ = velocitySmoother_;
            F_u_.correctBoundaryConditions();
        }

        // Update source fields to latest time step
        inline void updateSources()
        {
            // Compute continuity equation source
            // F_p = (rhoMean / c_p) * (DMean/Dt(sPrime))
            // ddt(s - sMean) = ddt(s) since ddt(sMean) = 0
            F_p_ = entropyDerivativeSwitch * (sqr(baseFlow_.cMean()) * (baseFlow_.rhoMean() / thermo_.Cp()) * (fvc::ddt(fields_.s()) + (baseFlow_.UMean() & fvc::grad((fields_.s() - baseFlow_.sMean()))))) * sourceDamping_;

            // Compute momentum equation source
            // F_u = (TPrime * grad(sMean)) - (sPrime * grad(TMean)) - LPrime;
            F_u_ = ((entropyGradientSwitch * (((thermo_.T() - baseFlow_.TMean()) * fvc::grad(baseFlow_.sMean())) - ((fields_.s() - baseFlow_.sMean()) * fvc::grad(baseFlow_.TMean())))) - (lambVectorSwitch * (fields_.L() - baseFlow_.LMean()))) * sourceDamping_;

            if (smoothenSourcesSwitch)
            {
                smoothenSources();
                F_0_p_ = pInterp_.interpolate(pressureSmoother_);
                F_0_u_ = pInterp_.interpolate(velocitySmoother_);
            }
            else
            {
                // Pointwise interpolation of sources
                F_0_p_ = pInterp_.interpolate(F_p_);
                F_0_u_ = pInterp_.interpolate(F_u_);
            }

            // Assign fieldsToSend_ for either 2D or 3D case
            if (isThreeDimensional_)
            {
                forAll(mesh_.points(), i)
                {
                    fieldsToSend_[(4 * i) + 0] = F_0_u_[i].x();
                    fieldsToSend_[(4 * i) + 1] = F_0_u_[i].y();
                    fieldsToSend_[(4 * i) + 2] = F_0_u_[i].z();
                    fieldsToSend_[(4 * i) + 3] = F_0_p_[i];
                }
            }
            else
            {
                forAll(mesh_.points(), i)
                {
                    fieldsToSend_[(3 * i) + 0] = F_0_u_[i].x();
                    fieldsToSend_[(3 * i) + 1] = F_0_u_[i].y();
                    fieldsToSend_[(3 * i) + 2] = F_0_p_[i];
                }
            }
        };

        // Run time
        const Time &runTime_;

        // Mesh
        const fvMesh &mesh_;

        // Thermo
        const psiThermo &thermo_;

        // Pointwise interpolation
        const volPointInterpolation pInterp_;

        // Send tag and status
        int sendTag;
        int status;

        // Spatial dimension
        const size_t dim_;
        const bool isThreeDimensional_;

        // Parameters for coupling
        const scalar lambVectorSwitch;
        const scalar entropyGradientSwitch;
        const scalar entropyDerivativeSwitch;
        const label cwipiStep;
        label cwipiTimeStep;

        // Mean flow fields
        const cwipiMeanFields baseFlow_;

        // Cwipi fields object
        cwipiFields fields_;

        // Acoustic sources and damping
        const volScalarField sourceDamping_;
        volScalarField F_p_;
        volVectorField F_u_;

        // Point fields
        pointScalarField F_0_p_;
        pointVectorField F_0_u_;

        // List of argument names
        const char *sourceFieldNames_;

        // Mesh info for cwipi
        std::vector<scalar> pointCoords_;
        std::vector<int> connecIdx_;
        std::vector<int> connec_;

        std::vector<scalar> pointCoordsInitialise(const fvMesh &mesh) const
        {
            std::vector<scalar> pointCoords = std::vector<scalar>(3 * mesh.nPoints(), 0);
            // Create mesh connectivity list
            forAll(mesh.points(), i)
            {
                pointCoords[3 * i + 0] = mesh.points()[i].x();
                pointCoords[3 * i + 1] = mesh.points()[i].y();
                pointCoords[3 * i + 2] = mesh.points()[i].z();
            }

            return pointCoords;
        }

        std::vector<int> connecIdxInitialise(const fvMesh &mesh) const
        {
            std::vector<int> connecIdx = std::vector<int>(mesh.nCells() + 1, 0);

            connecIdx[0] = 0;
            forAll(mesh.cells(), i)
            {
                connecIdx[i + 1] = connecIdx[i] + 8;
            }

            return connecIdx;
        };

        std::vector<int> connecInitialise(const fvMesh &mesh) const
        {
            std::vector<int> connec = std::vector<int>(mesh.nCells() * 8, 0);

            forAll(mesh.cells(), i)
            {
                forAll(mesh.cellShapes()[i], j)
                {
                    connec[8 * i + j] = mesh.cellShapes()[i][j] + 1;
                }
            }

            return connec;
        }

        // Cwipi source array
        std::vector<scalar> fieldsToSend_;

        std::vector<scalar> fieldsToSendInitialise(
            const fvMesh &mesh,
            const pointScalarField &F_p,
            const pointVectorField &F_U,
            const size_t dim) const
        {
            std::vector<scalar> fieldsToSend = std::vector<scalar>((dim + 1) * mesh.nPoints(), 0);

            if (dim == 3)
            {
                forAll(mesh.points(), i)
                {
                    fieldsToSend[(4 * i) + 0] = F_U[i].x();
                    fieldsToSend[(4 * i) + 1] = F_U[i].y();
                    fieldsToSend[(4 * i) + 2] = F_U[i].z();
                    fieldsToSend[(4 * i) + 3] = F_p[i];
                }
            }
            else
            {
                forAll(mesh.points(), i)
                {
                    fieldsToSend[(3 * i) + 0] = F_U[i].x();
                    fieldsToSend[(3 * i) + 1] = F_U[i].y();
                    fieldsToSend[(3 * i) + 2] = F_p[i];
                }
            }

            return fieldsToSend;
        }

        // Laplacian smoothing coefficient
        const volScalarField DT_;
        volScalarField pressureSmoother_;
        volVectorField velocitySmoother_;
        const bool smoothenSourcesSwitch;
    };

}

#endif

// ************************************************************************* //
