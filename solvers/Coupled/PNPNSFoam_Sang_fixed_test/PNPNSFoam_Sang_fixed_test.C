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

    bool solveTransient(readBool(runTime.controlDict().lookup("PNPNStransient")));
    bool solveNS(readBool(runTime.controlDict().lookup("solveNS")));
    bool debugMode(readBool(runTime.controlDict().lookup("debugMode")));

//    simpleControl simple(mesh);
    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "createControls.H"
//    #include "initConvergenceCheck.H"

    Info<< "\nCalculating scalar transport  \n" << endl;
    double timeEnd = runTime.endTime().value();
    double timeLap = 1.0;
    scalar phiInstant_ = ECsystem.phiInstant();

    Info << "phi start: " << ECsystem.phiStart() << endl;

    while(ECsystem.phiRun())
    {
        phiInstant_ = ECsystem.phiInstant();
        Info << "\n-------------- | Applied Voltage = " << phiInstant_<< " |--------------"<<endl;
        forAll(mesh.boundary(),patchI)
        {
            if(mesh.boundary()[patchI].name() == phiDyMBound)
            {
                scalarField& Phib = psiE.boundaryField()[patchI];
                forAll(Phib, faceI)
                {
                    Phib[faceI] = phiInstant_;
                }
            }
        }
        //setTime(value() + deltaT_, timeIndex_ + 1);

        Info << "delta time: " << runTime.deltaT() << endl;
        Info << "End time: " << runTime.endTime() << endl;
        Info << "current time: " << runTime.value() << endl;
        runTime.setEndTime(phiInstant_);// (timeEnd*timeLap);
        //Info << "time Index: " << runTime.timeIndex() << endl;
        runTime.setTime(runTime.endTime() - runTime.deltaT(), runTime.timeIndex());        
        Info << "time Index: " << runTime.timeIndex() << endl;
        Info << "current time: " << runTime.value() << endl;
        scalarList Is(nPNPNSIter,0);
        scalar ipnp_ns = 0,dC(10), dU(10);
        scalar dCmin(readScalar(runTime.controlDict().lookup("C1C2convergence")));
        word dVstring(" >>> dV = ");
        char phiInstant_str[6];
        gcvt(phiInstant_,6,phiInstant_str);

        bool reachedResidual = false;
        Info<< "\nStarting time loop\n" << endl;

        while (runTime.loop())
        {
            #include "readBlockSolverControls.H"
            #include "readFieldBounds.H"
            ipnp_ns = 0;dC=dCmin+10;

            maxResidual_PNP = 10;

            reachedResidual = false;
            Info<< "Time = " << runTime.timeName() << nl << endl;


            Info <<dVstring<<phiInstant_<<":         SteadyState solver"<<endl;
            while(ipnp_ns < nPNPNSIter and dC>dCmin)
            {
                Info <<"         SteadyState solver/PNP-NS Iteration      # " << ipnp_ns <<endl;
                volScalarField cPlusold = cPlus;
                volVectorField Uold = U;
                volScalarField psiEold = psiE;
//                Info<<endl<<"*** PNPEqn ***"<<endl;
                if(runTime.controlDict().lookupOrDefault<Switch>("solveNS",true))
                {
                    Info<<endl<<"*** PNPNS ***"<<endl;
                    //#include "PNPNSEqn.H"
                    //#include "PNPEqn.H"
                    //#include "NSNewtonEqn.H"
                } else 
                {
                    Info<<endl<<"****** only PNP ******"<<endl;
                    //#include "PNPEqn.H"
                    /*
                    */
                    // test laplacian(psiE) Newton raphson
                    volScalarField F = fvc::laplacian(psiE); // fvc::surfaceSum(fvc::snGrad(psiE));// fvc::laplacian(psiE);
                    volScalarField psiET = psiE;
                    for(int NewtonIter = 0; NewtonIter < nNewtonIterPNP; NewtonIter++)
                    {
                        fvBlockMatrix<vector3> PNPEqn(PNP);
                        fvBlockMatrix<vector2> psiEFieldEqn(psiEField);

                        volScalarField psiEF = psiE;
                        forAll(psiEF.boundaryField(), patchI)
                        {
                            scalarField& psiEFb = psiEF.boundaryField()[patchI];
                            //Info << psiEFb << endl;
                            psiEFb = 0;
                            //Info << psiEFb << endl;
                            //Info<< psiEF.boundaryField()[patchI] << endl;
                        }

                        fvScalarMatrix psiEEqn
                        (
                            fvm::laplacian(psiE) // == - fvc::laplacian(psiEF)
                        );
                        Info<< psiEEqn.source() << endl;
                        //psiEEqn += fvc::laplacian(psiEF);
                        //Info<< psiEEqn.source() << endl;

                        psiEFieldEqn.insertEquation(0, psiEEqn);
                        fvScalarMatrix cPlusEqn
                        (
                            fvm::laplacian(cPlus)//  == - fvc::laplacian(psiEF)
                        );
                        psiEFieldEqn.insertEquation(1, cPlusEqn);
                        //Info << fvc::laplacian
                        
                        Field<tensor2> d = psiEFieldEqn.diag().asSquare();
                        forAll(d, cellI)
                        {
                            Info<< "diag of psiE: "<< d[cellI](0, 0) << endl;
                        }
                        Info << "source of psiE"<<psiEFieldEqn.source() << endl;

                        //psiEEqn += fvc::laplacian(psiEF);
                        //psiEEqn.solve();
                        //psiEFieldEqn.insertEquation(1, psiEEqn);

                        volScalarField F1 = fvc::surfaceIntegrate(fvc::snGrad(psiE)); // fvc::laplacian(psiE);
                        // Assign F to blockM's right-hand side
                        forAll(psiEFieldEqn.source(),cellI)
                        {
                            psiEFieldEqn.source()[cellI](0) = -F1[cellI]*mesh.V()[cellI];
                        }

                        maxResidual_PNP = cmptMax(psiEFieldEqn.solve().initialResidual()); // residual = b - A*x_n
                        // Retrieve solution
                        psiEFieldEqn.retrieveSolution(0, psiE.internalField());
                        psiEFieldEqn.retrieveSolution(1, cPlus.internalField());

                        forAll(psiE.internalField(), cellI)
                        {
                            Info << "dpsiE's cell" << cellI << ": " << psiE.internalField()[cellI] << endl;
                        }
                        scalar lengthpsiE(0);

                        // residual-length
                        forAll(psiE.internalField(),i){ lengthpsiE     += sqr(psiE.internalField()[i]); }
                        psiET += psiE;
                        psiE = psiET;
                        psiE.correctBoundaryConditions();
                        forAll(psiE.internalField(), cellI)
                        {
                            Info << "psiE's cell" << cellI << ": " << psiE.internalField()[cellI] << endl;
                        }

                        Info <<" Iter # "<< NewtonIter  << " - Residual: " <<"[" << "psiE:" << Foam::sqrt(lengthpsiE) << "]"  << endl;

                        scalar residual = Foam::sqrt(lengthpsiE);
                        if (residual < convergenceCriterionPNP){ break; }


/*
                        volScalarField psiET
                        (
                            IOobject
                            (
                                "psiET",
                                psiE.instance(),
                                mesh,
                                IOobject::NO_READ,
                                IOobject::NO_WRITE
                            ),
                            mesh,
                            dimensioned<scalar>("0", psiE.dimensions(), pTraits<scalar>::zero),
                            zeroGradientFvPatchScalarField::typeName//zeroGradientFvPatchField<Type>::typeName
                        );
                        volScalarField dpsiE = psiE;
                        psiET.internalField() = psiE.internalField();
                        //dpsiE.internalField() = 0; // psiE.internalField();

                        psiET.internalField() = 10;

                        forAll(psiE.internalField(), cellI)
                        {
                            Info << "before psiE's cell" << cellI << ": " << psiE.internalField()[cellI] << endl;
                        }

                        //(fvc::laplacian(psiET)).internalField() *= mesh.V();
                        //F.internalField() *= mesh.V();
                        fvScalarMatrix psiEEqn
                        (
                            fvm::laplacian(psiE)
                        );
                        Info << "psiEEqn's source: " << psiEEqn.source() << endl;

                        //psiEEqn.source() = -F;
                        psiEEqn.solve();
                        //dpsiE += psiE;
                        //psiE = dpsiE; 
                        F = fvc::laplacian(psiE); // fvc::surfaceSum(fvc::snGrad(psiE));// fvc::laplacian(psiE);

                        forAll(F.internalField(), cellI)
                        {
                            Info << "F's cell" << cellI << ": " << F.internalField()[cellI] << endl;
                        }
                    
                        forAll(psiE.internalField(), cellI)
                        {
                            Info << "after psiE's cell" << cellI << ": " << psiE.internalField()[cellI] << endl;
                        }
                        Info << "psiEEqn's diag: " << psiEEqn.diag() << endl;
                        //Info << fvc::snGrad(psiE) << endl;
                        
                        psiE.correctBoundaryConditions();
                        scalar lengthpsiE(0);

                        // residual-length
                        forAll(dpsiE.internalField(),i){ lengthpsiE     += sqr(dpsiE.internalField()[i]); }

                        Info <<" Iter # "<< NewtonIter  << " - Residual: " <<"[" << "psiE:" << Foam::sqrt(lengthpsiE) << "]"  << endl;

                        scalar residual;
                        residual = Foam::sqrt(lengthpsiE);
                        if (residual < convergenceCriterionPNP){ break; }

*/
                    }
                }
                
                dC = mag(gMax((psiE.internalField() - psiEold.internalField())));// mag(gMax((cPlus.internalField() - cPlusold.internalField())));
                Info <<endl<<nl<<dVstring<<phiInstant_<<": C convergence: "<<dC<<endl<<endl;
                
                dU = mag(gMax((U.internalField() - Uold.internalField())));
                Info <<endl<<nl<<dVstring<<phiInstant_<<": U convergence: "<<dU<<endl<<endl;
                ipnp_ns++;
            }

            runTime.write();

            Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;

            if(max(maxResidual_PNP, maxResidual_NS) < convergenceCriterionTotal and !solveTransient)
            {
                Info<< "reached convergence criterion for steady solver: " << convergenceCriterionTotal << endl;
                runTime.writeAndEnd();
                Info<< "latestTime = " << runTime.timeName() << endl;
            }
        } // Time loop closed
        timeLap+=1.0;
        ECsystem.changePhi();
    }
    return 0;
}
