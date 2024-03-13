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
    scalar sqrDebL_nonD = 1/sqr(DebL_nonD.value());          

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
                    psiEb[faceI] = phiInstant_; // psiE0: thermal voltage
                }
            }
        }
        
        runTime.setEndTime(timeEnd*timeLap);

        scalar dcPlus(10);
        scalar PNPIter = 0;
        scalar dcPlusMin(readScalar(runTime.controlDict().lookup("C1C2convergence")));
        word dVstring(" >>> dV = ");
    
        Info<< "\nStarting time loop\n" << endl;
        while (runTime.loop())
        {            
            Info<< "Time = " << runTime.timeName() << nl << endl;
            #include "readFieldBounds.H"
            #include "readBlockSolverControls.H"
            maxResidual = 10;
            dcPlus=dcPlusMin+10;

            while(PNPIter < nPNPIter and dcPlus > dcPlusMin)
            {
                Info <<"         SteadyState solver/PNP-NS Iteration      # " << PNPIter <<endl;
                volScalarField cPlusold = cPlus;
                #include "PNPNSEqn.H"
                dcPlus = mag(gMax((cPlus.internalField() - cPlusold.internalField())));
                Info <<endl<<nl<<dVstring<<phiInstant_<<": C convergence: "<<dcPlus<<endl<<endl;
                PNPIter++;
            }
            
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
