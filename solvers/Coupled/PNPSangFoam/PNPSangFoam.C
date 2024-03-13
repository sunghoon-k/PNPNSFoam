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
                    psiEb[faceI] = phiInstant_ * psiE0.value(); // psiE0: thermal voltage
                }
            }
        }
        
        runTime.setEndTime(timeEnd*timeLap);

        Info<< "\nStarting time loop\n" << endl;
        while (runTime.loop())
        {            
            #include "readFieldBounds.H"
            #include "readBlockSolverControls.H"
            maxResidual = 10;

            reachedResidual = false;
            Info<< "Time = " << runTime.timeName() << nl << endl;

            /********************************************************************************************************/
            volScalarField cPlusT(cPlus);        
            volScalarField cMinusT(cMinus);        
            volScalarField psiET(psiE);        
            volScalarField cV ( IOobject ( "cV", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE ), mesh, dimensionedScalar("zero", dimensionSet(0, 3, 0, 0, 0, 0, 0), 0.0) );   

            scalarField& cV_ = cV.internalField();        

            forAll(cV_, k) {cV_[k] = mesh.V()[k];}        
            /********************************************************************************************************/

            for(int PNPIter = 0; PNPIter < nPNPIter; PNPIter++)
            {
                if(solveTransient)
                {
                    Info <<"         Transient solver/PNP-NS Iteration      # " << PNPIter+1 <<endl;
                }
                else
                {
                    Info <<"         Steady solver/PNP-NS Iteration      # " << PNPIter+1 <<endl;
                }

                //Info <<"before Sang" <<endl;
                #include "Sang.H"
                //Info <<"after Sang" <<endl;

                bool bounded = false;

                fvBlockMatrix<vector3> PNPEqn(PNP);

                #include "psiEEqn.H"
                #include "cEqn.H"
                //#include "couplingTerms.H"

                // Assign F to blockM's right-hand side            
                forAll(PNPEqn.source(),cellI)            
                {                
                    PNPEqn.source()[cellI](0) = -FcPlus[cellI];                
                    PNPEqn.source()[cellI](1) = -FcMinus[cellI];                
                    PNPEqn.source()[cellI](2) = -FpsiE[cellI];            
                }            

                maxResidual = cmptMax(PNPEqn.solve().initialResidual()); // residual = b - A*x_n
                
                // Retrieve solution
                PNPEqn.retrieveSolution(0, cPlus.internalField());
                PNPEqn.retrieveSolution(1, cMinus.internalField());
                PNPEqn.retrieveSolution(2, psiE.internalField());

                scalar lengthcPlus(0);
                scalar lengthcMinus(0);
                scalar lengthpsiE(0);

                // cPlus residual-length

                forAll(cPlus.internalField(),i)
                {
                    lengthcPlus += sqr(cPlus.internalField()[i]);
                }

                // cMinus residual-length

                forAll(cMinus.internalField(),i)
                {
                    lengthcMinus += sqr(cMinus.internalField()[i]);
                }

                // Phi residual-length

                forAll(psiE.internalField(),i)
                {
                    lengthpsiE += sqr(psiE.internalField()[i]);
                }

                Info <<" Iter # "<< PNPIter  << " - Residual: " <<"[" << "cPlus:" << Foam::sqrt(lengthcPlus)
                << "  " << "cMinus:" << Foam::sqrt(lengthcMinus) << "  " << "Phi:" << Foam::sqrt(lengthpsiE) << "]"  << endl;

                // Correct cPlus,cMinus,Phi

                cPlusT  += cPlus;
                cMinusT  += cMinus;
                psiET += psiE;
                // Assign cPlusT,cMinusT,PhiT to cPlus,cMinus,Phi
                cPlus  = cPlusT;
                cMinus  = cMinusT;
                psiE = psiET;

                psiE.correctBoundaryConditions();
                cPlus.correctBoundaryConditions();
                cMinus.correctBoundaryConditions();


                // Correct Boundary value according to J-Flux
                forAll(cPlus.boundaryField(),patchI)
                {
                    if(cPlus.boundaryField()[patchI].type() == "fixedIonicFlux")
                    {
                        scalarField& Cb1 = cPlus.boundaryField()[patchI];
                        const tmp<scalarField>&  Ci1 = cPlus.boundaryField()[patchI].patchInternalField();
                        scalarField& Phib1 = psiE.boundaryField()[patchI];
                        const tmp<scalarField>& Phii1 = psiE.boundaryField()[patchI].patchInternalField();
                        Cb1 = Ci1/(1.0 + zPlus.value()*(Phib1 - Phii1));
                    }
                }


                forAll(cMinus.boundaryField(),patchI)
                {
                    if(cMinus.boundaryField()[patchI].type() == "fixedIonicFlux")
                    {
                        scalarField& Cb2 = cMinus.boundaryField()[patchI];
                        const tmp<scalarField>& Ci2 = cMinus.boundaryField()[patchI].patchInternalField();
                        scalarField& Phib2 = psiE.boundaryField()[patchI];
                        const tmp<scalarField>&  Phii2 = psiE.boundaryField()[patchI].patchInternalField();
                        Cb2 = Ci2/(1.0 + zMinus.value()*(Phib2 - Phii2));
                    }
                }


                scalar residual;
                residual = Foam::sqrt(lengthcPlus + lengthcMinus + lengthpsiE);
                if (residual < convergenceCriterion)
                {
                    break;
                }
                #include "boundC.H"
                //#include "convergenceCheck.H"
                netCharge = e*(zPlus*cPlus + zMinus*cMinus);
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
