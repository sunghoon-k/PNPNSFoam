/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    pUCoupledFoam

Description
    Steady-state solver for incompressible, turbulent flow, with implicit
    coupling between pressure and velocity achieved by fvBlockMatrix.
    Turbulence is in this version solved using the existing turbulence
    structure.

Authors
    Klas Jareteg, Chalmers University of Technology,
    Vuko Vukcevic, FMENA Zagreb.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvBlockMatrix.H"
#include "singlePhaseTransportModel.H"
//#include "RASModel.H"
//#include "simpleControl.H"
#include "ElectrochemicalSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

//    simpleControl simple(mesh);
    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "createControls.H"
    #include "initConvergenceCheck.H"

    scalar phiInstant_ = 0;
    double timeEnd = runTime.endTime().value();
    double timeLap = 1.0;

    bool solveTransient(readBool(runTime.controlDict().lookup("PNPNStransient")));
    bool bounded = true;
    Info<< "\nStarting time loop\n" << endl;

    while(ECsystem.phiRun())
    {
        phiInstant_ = ECsystem.phiInstant();
        Info << "\n-------------- | psiE = " << phiInstant_<< " |--------------"<<endl;
        forAll(mesh.boundary(),patchI)
        {
            if(mesh.boundary()[patchI].name() == phiDyMBound)
            {
                scalarField& psiEb = psiE.boundaryField()[patchI];
                forAll(psiEb, faceI)
                {
                    psiEb[faceI] = phiInstant_; // psiE0: thermal voltage
                }
            }
        }
        
        runTime.setEndTime(timeEnd*timeLap);

        Info<< "\nStarting time loop\n" << endl;
        while (runTime.loop())
        {
            #include "readBlockSolverControls.H"
            #include "readFieldBounds.H"
            // #include "CourantNo.H" // 얘는 delta time 수정용
            
            maxResidual = 10;

            Info<< "Time = " << runTime.timeName() << nl << endl;

            scalar PNPIter = 0;

            Info<< "#############################################" << nl 
            << "Get initial value from PNPEqn" << nl << endl;            
            #include "PNPEqn.H"

            if(solveTransient)
            {
                while(PNPIter++ < nPNPIter and maxResidual > 1e-6) // for(int PNPIter = 0; PNPIter < nOuterIter; OuterIter++)
                {
                    Info <<"         Transient solver/PNP-NS Iteration      # " << PNPIter+1 <<endl;
                    #include "PNPNSEqn.H"
                }
            }
            else // solve steady-state
            {
                while(PNPIter++ < nPNPIter and maxResidual > 1e-6) // for(int PNPIter = 0; PNPIter < nOuterIter; OuterIter++)
                {
                    Info <<"         Steady solver/PNP-NS Iteration      # " << PNPIter+1 <<endl;
                    bounded = false;
                    #include "PNPNSEqn.H"
                    netCharge = (zPlus*cPlus + zMinus*cMinus);
                    bodyForce = - netCharge * fvc::grad(psiE);
                }
            }
            
            Info<< "maxResidual = " << maxResidual << nl
                << endl;
            
            // turbulence->correct();
            runTime.write();

            Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;

        } // Time loop closed

        timeLap+=1.0;

        ECsystem.changePhi();

    } // phiRun() loop closed


    Info<< "End\n" << endl;

    return 0;
}
