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

//    simpleControl simple(mesh);
    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "createControls.H"
//    #include "initConvergenceCheck.H"

    Info<< "\nCalculating scalar transport  \n" << endl;
    double timeEnd = runTime.endTime().value();
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

        runTime.setEndTime(timeEnd*timeLap);
        scalarList Is(nPNPNSIter,0);
        scalar ipnp_ns = 0,dC(10);
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
                Info<<endl<<"*** PNPEqn ***"<<endl;
                #include "PNPEqn.H"
                if(runTime.controlDict().lookupOrDefault<Switch>("solveNS",true))
                {
                Info<<endl<<"*** NS ***"<<endl;
                #include "NSEqn.H"
                } else Info<<endl<<"****** not solve NS ******"<<endl;
                dC = mag(gMax((cPlus.internalField() - cPlusold.internalField())));
                Info <<endl<<nl<<dVstring<<phiInstant_<<": C convergence: "<<dC<<endl<<endl;
                ipnp_ns++;

            }
            
            //Info<< "maxResidual = " << maxResidual << nl
            //    << endl;

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
