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

#include "argList.H"
#include "fvMesh.H"
#include "fvCFD.H"
#include "volPointInterpolation.H"
#include "cwipiPstream.H"
#include <cwipi.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Constructor
    cwipiPstream::cwipiPstream(
        const Foam::Time &runTime,
        const fvMesh &mesh,
        const psiThermo &thermo,
        const volVectorField &U)
        : runTime_(runTime),                                                                                                                                                                                                                                                                                                                       // Run time
          mesh_(mesh),                                                                                                                                                                                                                                                                                                                             // Mesh
          thermo_(thermo),                                                                                                                                                                                                                                                                                                                         // Thermo model
          pInterp_(mesh),                                                                                                                                                                                                                                                                                                                          // Pointwise interpolation
          sendTag(0),                                                                                                                                                                                                                                                                                                                              // Set send tag to 0
          status(0),                                                                                                                                                                                                                                                                                                                               // Set status to 0
          dim_(static_cast<uint8_t>(readInt(runTime.controlDict().lookup("cwipiDim")))),                                                                                                                                                                                                                                                           // Get dimension
          isThreeDimensional_(static_cast<bool>(readInt(runTime.controlDict().lookup("cwipiDim")) - 2)),                                                                                                                                                                                                                                           // Get switch for 3d
          lambVectorSwitch(static_cast<Foam::scalar>(readBool(runTime.controlDict().lookup("cwipiLambVector")))),                                                                                                                                                                                                                                  // Cast lamb vector coefficient to scalar
          entropyGradientSwitch(static_cast<Foam::scalar>(readBool(runTime.controlDict().lookup("cwipiEntropy")))),                                                                                                                                                                                                                                // Cast entropy gradient coefficient to scalar
          entropyDerivativeSwitch(static_cast<Foam::scalar>(readBool(runTime.controlDict().lookup("cwipiDsDt")))),                                                                                                                                                                                                                                 // Cast entropy material derivative coefficient to scalar
          cwipiStep(readInt(runTime.controlDict().lookup("cwipiStep"))),                                                                                                                                                                                                                                                                           // Get time step
          cwipiTimeStep(readInt(runTime.controlDict().lookup("cwipiStep"))),                                                                                                                                                                                                                                                                       // Assign time step
          baseFlow_(cwipiMeanFields(mesh, runTime)),                                                                                                                                                                                                                                                                                               // Mean flow fields
          fields_(cwipiFields(mesh, runTime, U, thermo)),                                                                                                                                                                                                                                                                                          // Instantaneous flow fields
          sourceDamping_(IOobject("sourceDamping", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE), mesh, dimensionedScalar("one", dimless, 1)),                                                                                                                                                                        // Source damping coefficient
          F_p_(IOobject("F_p", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE), entropyDerivativeSwitch * ((baseFlow_.rhoMean() / thermo_.Cp()) * (fvc::ddt(fields_.s()) + (baseFlow_.UMean() & fvc::grad((fields_.s() - baseFlow_.sMean()))))) * sourceDamping_),                                                              // Continuity equation sources
          F_u_(IOobject("F_u", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE), ((entropyGradientSwitch * ((thermo_.T() - baseFlow_.TMean()) * fvc::grad(baseFlow_.sMean())) - ((fields_.s() - baseFlow_.sMean()) * fvc::grad(baseFlow_.TMean()))) - (lambVectorSwitch * (fields_.L() - baseFlow_.LMean()))) * sourceDamping_), // Momentum equation sources
          F_0_p_(pInterp_.interpolate(F_p_)),                                                                                                                                                                                                                                                                                                      // Pointwise interpolation of continuity equation sources
          F_0_u_(pInterp_.interpolate(F_u_)),                                                                                                                                                                                                                                                                                                      // Pointwise interpolation of momentum equation sources
          sourceFieldNames_(cwipiSourceFieldNames(static_cast<uint8_t>(readInt(runTime.controlDict().lookup("cwipiDim"))))),
          fieldsToSend_(std::vector<scalar>((dim_ + 1) * mesh_.nPoints(), 0))
    {
        // Resize mesh vectors to fit
        pointCoords.resize(3 * mesh_.nPoints());
        connecIdx.resize(mesh_.nCells() + 1);
        connec.resize(mesh_.nCells() * 8);

        // Create mesh connectivity list
        forAll(mesh_.points(), i)
        {
            pointCoords[3 * i + 0] = mesh_.points()[i].x();
            pointCoords[3 * i + 1] = mesh_.points()[i].y();
            pointCoords[3 * i + 2] = mesh_.points()[i].z();
        }
        connecIdx[0] = 0;
        forAll(mesh_.cells(), i)
        {
            connecIdx[i + 1] = connecIdx[i] + 8;
            forAll(mesh_.cellShapes()[i], j)
            {
                connec[8 * i + j] = mesh_.cellShapes()[i][j] + 1;
            }
        }

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

        // Add local control parameters
        cwipi_add_local_int_control_parameter(
            "nSendVars",
            dim_ + 1);
        cwipi_add_local_string_control_parameter(
            "sendFieldNames",
            sourceFieldNames_);
        cwipi_add_local_int_control_parameter(
            "receiveTag",
            sendTag);

        // Create coupling
        cwipi_create_coupling(
            "cwipiFoam",
            CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING,
            "FOAM_APE",
            3,
            1,
            CWIPI_STATIC_MESH,
            CWIPI_SOLVER_CELL_VERTEX,
            0,
            "EnSight Gold",
            "text");

        // Sync and dump
        cwipi_synchronize_control_parameter(
            "FOAM_APE");
        cwipi_dump_application_properties();

        // Get send tag
        sendTag = cwipi_get_distant_int_control_parameter(
            "FOAM_APE",
            "receiveTag");

        // Define mesh
        cwipi_define_mesh(
            "cwipiFoam",
            mesh_.nPoints(),
            mesh_.nCells(),
            pointCoords.data(),
            connecIdx.data(),
            connec.data());

        // Locate interpolation
        Info << "About to call cwipi_locate." << endl;
        Info << "This is a common point of failure due to a high oversampling rate defined in Nektar++." << endl;
        Info << "If the application crashes at this point, check the Oversample property of your coupling entry." << endl;
        cwipi_locate(
            "cwipiFoam");
        Info << "Interpolation located successfully." << endl;
    };

    // Destructor
    cwipiPstream::~cwipiPstream(){};

    // Send fields
    void cwipiPstream::send()
    {
        fields_.update();

        // Update source fields
        updateSources();

        if (sendNow())
        {
            // Send data
            cwipi_issend(
                "cwipiFoam",
                "ex1",
                sendTag,
                dim_ + 1,
                1,
                0,
                sourceFieldNames_,
                fieldsToSend_.data(),
                &status);

            // Handle exchange status
            switch (status)
            {
            case CWIPI_EXCHANGE_OK:
                // Info << "Exchange Ok" << endl;
                break;
            case CWIPI_EXCHANGE_BAD_RECEIVING:
                throw std::runtime_error("Error: Bad receive status.");
                // Info << "Bad receiving" << endl;
                break;
            default:
                throw std::runtime_error("Error: Undefined receive status.");
                // Info << "Error: bad exchange status" << endl;
                break;
            }
        }
    };

    // Advance time step
    void cwipiPstream::updateTime()
    {
        // Compare time step to step
        if (sendNow())
        {
            // Wait
            cwipi_wait_issend(
                "cwipiFoam",
                status);
        }

        // Use modulo operator and increment by 1, faster
        cwipiTimeStep = (cwipiTimeStep % cwipiStep) + 1;
    };

}

// ************************************************************************* //
