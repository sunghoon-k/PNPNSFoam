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
    bool pseudoTransient(readBool(runTime.controlDict().lookup("pseudoTransient")));

//    simpleControl simple(mesh);
    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "createControls.H"
//    #include "initConvergenceCheck.H"

    Info<< "\nCalculating scalar transport  \n" << endl;
    //double timeEnd = runTime.endTime().value();
    double timeLap = 1.0;
    scalar phiInstant_ = ECsystem.phiInstant();

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

        if(ECsystem.phiStart() != ECsystem.phiEnd())
        {
            runTime.setEndTime(phiInstant_);// (timeEnd*timeLap);
            runTime.setTime(runTime.endTime() - runTime.deltaT(), runTime.timeIndex());        
        }

        scalar dC(10), dU(0);
        scalar dCmin(readScalar(runTime.controlDict().lookup("C1C2convergence")));
        word dVstring(" >>> dV = ");
        //char phiInstant_str[6];
        //gcvt(phiInstant_,6,phiInstant_str);

        bool reachedResidual = false;
        Info<< "\nStarting time loop\n" << endl;

        while (runTime.loop())
        {
            #include "readBlockSolverControls.H"
            #include "readFieldBounds.H"
            //ipnp_ns = 0;
            dC=dCmin+10;

            maxResidual_PNP = 10;

            reachedResidual = false;
            Info<< "Time = " << runTime.timeName() << nl << endl;

            if(solveTransient)
            {
                Info <<dVstring<<phiInstant_<<":         TransientState solver"<<endl;
            }
            else if(pseudoTransient)
            {
                Info <<dVstring<<phiInstant_<<":         PseudoTransientState solver"<<endl;
            }   
            else
            {
                Info <<dVstring<<phiInstant_<<":         SteadyState solver"<<endl;
            }                 

            //while(ipnp_ns < nPNPNSIter and max(dC, dU)>dCmin)
            for(int PNPNSIter=0; PNPNSIter < nPNPNSIter; PNPNSIter++) 
            {
                if(solveTransient)
                {
                    Info <<"         Transient solver/PNP-NS Iteration      # " << PNPNSIter+1 <<endl;
                }
                else if(pseudoTransient)
                {
                    Info <<"         PseudoTransient solver/PNP-NS Iteration      # " << PNPNSIter+1 <<endl;
                }   
                else
                {
                    Info <<"         SteadyState solver/PNP-NS Iteration      # " << PNPNSIter+1 <<endl;
                }                 
                volScalarField cPlusold = cPlus;
                volScalarField cMinusold = cMinus;
                volVectorField Uold = U;

                volScalarField ddtcPlus = fvc::ddt(cPlus);
                volScalarField ddtcMinus = fvc::ddt(cMinus);
                volVectorField ddtU = fvc::ddt(U);

                if(solveNS)// (runTime.controlDict().lookupOrDefault<Switch>("solveNS",true))
                {
                    Info<<endl<<"*** PNPNS ***"<<endl;
                    #include "PNPNSEqn.H"
                } else 
                {
                    Info<<endl<<"****** only PNP ******"<<endl;
                    #include "PNPEqn.H"
                }
                
                dC = mag(gMax((cPlus.internalField() - cPlusold.internalField())));///psiEold.internalField()
                Info <<endl<<nl<<dVstring<<phiInstant_<<": C convergence: "<<dC<<endl<<endl;
                
                dU = mag(gMax((U.internalField() - Uold.internalField()))); // /Uold.internalField()
                Info <<endl<<nl<<dVstring<<phiInstant_<<": U convergence: "<<dU<<endl<<endl;
                //ipnp_ns++;
                if(max(dC, dU)<dCmin)
                {
                    reachedResidual = true;
                    break;
                }
            }

            Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;

            runTime.write();

            if(reachedResidual)
            {
                Info << "Reached!!!!!!!!!!!!!" << endl;
                break;
            }
        } // Time loop closed
        timeLap+=1.0;
        ECsystem.changePhi();
    }
    return 0;
}
