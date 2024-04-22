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

#include "fvCFD.H"
#include "volPointInterpolation.H"
#include "volFields.H"
#include "psiThermo.H"
#include "cwipiSwitch.H"
#include "cwipiFields.H"
#include "cwipiMeanFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    class cwipiPstream
    {
    public:
        // Default constructor
        cwipiPstream(
            const Foam::Time &runTime,
            const fvMesh &mesh,
            const psiThermo &thermo,
            const cwipiFields &sourceFields);

        // Default destructor
        ~cwipiPstream();

        // Send fields
        void send();

        // Advance time step
        void updateTime();

        // Evaluate send condition
        inline constexpr const Foam::Switch sendNow() const
        {
            return (cwipiTimeStep == cwipiStep);
        };

    private:
        // Return the names of the source fields
        inline const char *cwipiSourceFieldNames(const uint8_t dim)
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
        // Smoothen source fields
        // inline void smoothenSources()
        // {
        //     // Create scalar smooth field from virgin scalar smooth field template
        //     volScalarField &pressureSmoother = pressureSmoother_;
        //     // pressureSmoother.dimensions().reset(F_p_.dimensions());
        //     pressureSmoother == F_p_;
        //     pressureSmoother.correctBoundaryConditions();
        //     // pressureSmoother.oldTime().dimensions().reset(F_p_.dimensions());
        //     pressureSmoother.oldTime() = F_p_;
        //     pressureSmoother.oldTime().correctBoundaryConditions();

        //     // Create vector smooth field from virgin vector smooth field template
        //     volVectorField &velocitySmoother = velocitySmoother_;
        //     // velocitySmoother.dimensions().reset(F_u_.dimensions());
        //     velocitySmoother == F_u_;
        //     velocitySmoother.correctBoundaryConditions();
        //     // velocitySmoother.oldTime().dimensions().reset(F_u_.dimensions());
        //     velocitySmoother.oldTime() = F_u_;
        //     velocitySmoother.oldTime().correctBoundaryConditions();

        //     // Laplacian smoothing
        //     solve(fvm::ddt(pressureSmoother) - fvm::laplacian(DT_, pressureSmoother));
        //     solve(fvm::ddt(velocitySmoother) - fvm::laplacian(DT_, velocitySmoother));

        //     // Assign smoothened values
        //     F_p_ = pressureSmoother;
        //     F_p_.correctBoundaryConditions();

        //     F_u_ = velocitySmoother;
        //     F_u_.correctBoundaryConditions();
        // }

        // Update source fields to latest time step
        inline void updateSources()
        {
            // Compute continuity equation source
            // F_p = (rhoMean / c_p) * (DMean/Dt(sPrime))
            // ddt(s - sMean) = ddt(s) since ddt(sMean) = 0
            F_p_ = entropyDerivativeSwitch * ((baseFlow_.rhoMean() / thermo_.Cp()) * (fvc::ddt(sourceFields_.s()) + (baseFlow_.UMean() & fvc::grad((sourceFields_.s() - baseFlow_.sMean()))))) * sourceDamping_;

            // Compute momentum equation source
            // F_u = (TPrime * grad(sMean)) - (sPrime * grad(TMean)) - LPrime;
            F_u_ = ((entropyGradientSwitch * ((thermo_.T() - baseFlow_.TMean()) * fvc::grad(baseFlow_.sMean())) - ((sourceFields_.s() - baseFlow_.sMean()) * fvc::grad(baseFlow_.TMean()))) - (lambVectorSwitch * (sourceFields_.L() - baseFlow_.LMean()))) * sourceDamping_;

            // smoothenSources();

            // Pointwise interpolation of sources
            F_0_p_ = pInterp_.interpolate(F_p_);
            F_0_u_ = pInterp_.interpolate(F_u_);

            // Assign fieldsToSend for either 2D or 3D case
            if (isThreeDimensional_)
            {
                forAll(mesh_.points(), i)
                {
                    fieldsToSend[(4 * i) + 0] = F_0_u_[i].x();
                    fieldsToSend[(4 * i) + 1] = F_0_u_[i].y();
                    fieldsToSend[(4 * i) + 2] = F_0_u_[i].z();
                    fieldsToSend[(4 * i) + 3] = F_0_p_[i];
                }
            }
            else
            {
                forAll(mesh_.points(), i)
                {
                    fieldsToSend[(3 * i) + 0] = F_0_u_[i].x();
                    fieldsToSend[(3 * i) + 1] = F_0_u_[i].y();
                    fieldsToSend[(3 * i) + 2] = F_0_p_[i];
                }
            }
        };

// Pointwise interpolation
        const volPointInterpolation pInterp_;

        // Send tag and status
        int sendTag;
        int status;

        // Spatial dimension
        const uint8_t dim_;
        const bool isThreeDimensional_;

        // Parameters for coupling
        const scalar lambVectorSwitch;
        const scalar entropyGradientSwitch;
        const scalar entropyDerivativeSwitch;
        const label cwipiStep;
        label cwipiTimeStep;

        // Mean flow fields
        const cwipiMeanFields baseFlow_;

        // Acoustic sources and damping
        const volScalarField sourceDamping_;
        volScalarField F_p_;
        volVectorField F_u_;

        // Run time
        const Time &runTime_;

        // Mesh
        const fvMesh &mesh_;

        // Thermo
        const psiThermo &thermo_;

        // Cwipi fields object
        const cwipiFields &sourceFields_;

        // Point fields
        pointScalarField F_0_p_;
        pointVectorField F_0_u_;

        // Laplacian smoothing length
        // const dimensionedScalar smoothingLength_;
        // const dimensionedScalar DT_;

        // Laplacian smoothing fields
        // volScalarField pressureSmoother_;
        // volVectorField velocitySmoother_;

        // List of argument names
        const char *sourceFieldNames_;

        // Mesh info for cwipi
        std::vector<scalar> pointCoords = std::vector<scalar>();
        std::vector<int> connecIdx = std::vector<int>();
        std::vector<int> connec = std::vector<int>();

        // Cwipi source array
        std::vector<scalar> fieldsToSend = std::vector<scalar>();
    };

}

// ************************************************************************* //
