/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

#include "UOPstream.H"
#include "PstreamGlobals.H"
#include "profilingPstream.H"

#include <mpi.h>

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::UOPstream::bufferIPCsend()
{
    return UOPstream::write(
        commsType(),
        toProcNo_,
        sendBuf_.cdata(),
        sendBuf_.size(),
        tag_,
        comm_);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::UOPstream::write(
    const UPstream::commsTypes commsType,
    const int toProcNo,
    const char *buf,
    const std::streamsize bufSize,
    const int tag,
    const label communicator)
{
    if (debug)
    {
        Pout << "UOPstream::write : starting write to:" << toProcNo
             << " tag:" << tag
             << " comm:" << communicator << " size:" << label(bufSize)
             << " commType:" << UPstream::commsTypeNames[commsType]
             << Foam::endl;
    }
    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout << "UOPstream::write : starting write to:" << toProcNo
             << " tag:" << tag
             << " comm:" << communicator << " size:" << label(bufSize)
             << " commType:" << UPstream::commsTypeNames[commsType]
             << " warnComm:" << UPstream::warnComm
             << Foam::endl;
        error::printStack(Pout);
    }

    PstreamGlobals::checkCommunicator(communicator, toProcNo);

    bool failed = true;

    profilingPstream::beginTiming();

    if (commsType == commsTypes::blocking)
    {
        failed = MPI_Bsend(
            const_cast<char *>(buf),
            bufSize,
            MPI_BYTE,
            toProcNo,
            tag,
            PstreamGlobals::MPICommunicators_[communicator]);

        // Assume these are from scatters ...
        profilingPstream::addScatterTime();

        if (debug)
        {
            Pout << "UOPstream::write : finished write to:" << toProcNo
                 << " tag:" << tag << " size:" << label(bufSize)
                 << " commsType:" << UPstream::commsTypeNames[commsType]
                 << Foam::endl;
        }
    }
    else if (commsType == commsTypes::scheduled)
    {
        failed = MPI_Send(
            const_cast<char *>(buf),
            bufSize,
            MPI_BYTE,
            toProcNo,
            tag,
            PstreamGlobals::MPICommunicators_[communicator]);

        // Assume these are from scatters ...
        profilingPstream::addScatterTime();

        if (debug)
        {
            Pout << "UOPstream::write : finished write to:" << toProcNo
                 << " tag:" << tag << " size:" << label(bufSize)
                 << " commsType:" << UPstream::commsTypeNames[commsType]
                 << Foam::endl;
        }
    }
    else if (commsType == commsTypes::nonBlocking)
    {
        MPI_Request request;

        failed = MPI_Isend(
            const_cast<char *>(buf),
            bufSize,
            MPI_BYTE,
            toProcNo,
            tag,
            PstreamGlobals::MPICommunicators_[communicator],
            &request);

        profilingPstream::addWaitTime();

        if (debug)
        {
            Pout << "UOPstream::write : started write to:" << toProcNo
                 << " tag:" << tag << " size:" << label(bufSize)
                 << " commType:" << UPstream::commsTypeNames[commsType]
                 << " request:" << PstreamGlobals::outstandingRequests_.size()
                 << Foam::endl;
        }

        PstreamGlobals::outstandingRequests_.push_back(request);
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported communications type " << int(commsType)
            << Foam::abort(FatalError);
    }

    return !failed;
}

// ************************************************************************* //
