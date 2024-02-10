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
#include "simpleControl.H"
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

    #include "NormalizedParameters.H"

    bool solveTransient(readBool(runTime.controlDict().lookup("PNPNStransient")));
    scalar phiInstant_ = 0;
    double timeEnd = runTime.endTime().value();
    double timeLap = 1.0;

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
                    psiEb[faceI] = phiInstant_;
                }
            }
        }
    
        runTime.setEndTime(timeEnd*timeLap);

        Info<< "\nStarting time loop\n" << endl;
        while (runTime.loop())
        {
            #include "readBlockSolverControls.H"
            #include "readFieldBounds.H"
            #include "CourantNo.H"

            Info<< "Time = " << runTime.timeName() << nl << endl;

            for (int i=0; i<nInIter; i++)
            {
                p.storePrevIter();

                // Initialize the Up block system (matrix, source and reference to Up)
                fvBlockMatrix<vector7> PNPNSEqn(PNPNS);

                // Assemble and insert momentum equation
                #include "UEqn.H"

                // Assemble and insert pressure equation
                #include "pEqn.H"

                // Assemble and insert pressure equation
                #include "psiEEqn.H"

                // Assemble and insert pressure equation
                #include "cEqn.H"

                // Assemble and insert coupling terms
                #include "couplingTerms.H"


                // Solve the block matrix
        //        maxResidual = cmptMax(PNPNSEqn.solve().initialResidual());
                maxResidual = cmptMax(PNPNSEqn.solve().finalResidual()); // residual = b - A*x_n

                // Retrieve solution
                PNPNSEqn.retrieveSolution(0, U.internalField());
                PNPNSEqn.retrieveSolution(3, p.internalField());
                PNPNSEqn.retrieveSolution(4, psiE.internalField());
                PNPNSEqn.retrieveSolution(5, cPlus.internalField());
                PNPNSEqn.retrieveSolution(6, cMinus.internalField());

                U.correctBoundaryConditions();
                p.correctBoundaryConditions();
                psiE.correctBoundaryConditions();
                cPlus.correctBoundaryConditions();
                cMinus.correctBoundaryConditions();

                phi = (fvc::interpolate(U) & mesh.Sf()) + pEqn.flux() + presSource;

                #include "continuityErrs.H"

                #include "boundPU.H"

                p.relax();            
            } // Inner loop closed

            // turbulence->correct();
            runTime.write();

            // 반복문을 돌다가 수렴하지 않으면 break
            #include "convergenceCheck.H"

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
