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
    #include "createControls.H"
    #include "initConvergenceCheck.H"

    scalar phiInstant_ = 0;
    double timeEnd = runTime.endTime().value();
    double timeLap = 1.0;

    bool solveTransient(readBool(runTime.controlDict().lookup("PNPNStransient")));
    bool reachedResidual = false;
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
                    psiEb[faceI] = phiInstant_; //  * psiE0.value(); // psiE0: thermal voltage
                }
            }
        }
        
        runTime.setEndTime(timeEnd*timeLap);

        Info<< "\nStarting time loop\n" << endl;
        while (runTime.loop())
        {
            #include "readBlockSolverControls.H"
            #include "readFieldBounds.H"

            maxResidual = 10;

            reachedResidual = false;
            Info<< "Time = " << runTime.timeName() << nl << endl;

            // scalar PNPIter = 0;
            // while(PNPIter++ < nPNPIter and maxResidual > 1e-6) // 
            for(int PNPIter = 0; PNPIter < nPNPIter; PNPIter++)
            {
                bool bounded = false;

                fvBlockMatrix<vector3> PNPEqn(PNP);

                #include "psiEEqn.H"
                #include "cEqn.H"
                #include "couplingTerms.H"

                volScalarField F1 = fvc::surfaceIntegrate(J1)
                maxResidual = cmptMax(PNPEqn.solve().initialResidual()); // residual = b - A*x_n
                
                // Retrieve solution
                PNPEqn.retrieveSolution(0, psiE.internalField());
                PNPEqn.retrieveSolution(1, cPlus.internalField());
                PNPEqn.retrieveSolution(2, cMinus.internalField());

                psiE.correctBoundaryConditions();
                cPlus.correctBoundaryConditions();
                cMinus.correctBoundaryConditions();

                #include "boundC.H"
            }
            
            Info<< "maxResidual = " << maxResidual << nl
                << endl;

/*
                // Segregated solver에서는 이 기준을 사용할 필요가 없다. 
                // 왜냐? 이미 fvSolution 의 기준으로 이미 tolerance를 통과했기 때문
                if (maxResidual < 1e-6 )
                {
                    reachedResidual = true;
                    Info<< "\nConvergence reached\n" << endl;
                    runTime.writeAndEnd();

                    break;
                }

*/
            
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
