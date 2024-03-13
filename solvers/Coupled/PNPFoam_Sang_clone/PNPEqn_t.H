        volScalarField C1T(C1);
        volScalarField C2T(C2);
        volScalarField PhiT(Phi);

        volScalarField cV
        (
            IOobject
            (
                "cV",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimensionSet(0, -1, 0, 0, 0, 0, 0), 0.0)
        );
        scalarField& cV_ = cV.internalField();
        forAll(cV_, k) {cV_[k] = mesh.V()[k];}

        // Solve PNP Equations
        scalar nNewtonIteration(readScalar(runTime.controlDict().lookup("nNewtonIteration")));
        scalar convergenceCriterion(readScalar(mesh.solutionDict().subDict("solvers").subDict("blockPNP").lookup("convergenceCriterion")));
        dimensionedScalar dimL("dimL",dimensionSet(0,1,0,0,0,0,0), 1);
        dimensionedScalar dimU("dimPhi",dimensionSet(0,1,-1,0,0,0,0), 1);
        dimensionedScalar dimt("dimt",dimensionSet(0,-1,1,0,0,0,0), 1);
        int i; i=0;
        while((i++)<nNewtonIteration)
        {
            surfaceScalarField gradC1f = fvc::snGrad(C1);
            surfaceScalarField gradC2f = fvc::snGrad(C2);
            surfaceScalarField C1f = fvc::interpolate(C1);
            surfaceScalarField C2f = fvc::interpolate(C2);
            surfaceScalarField gradPhif = fvc::snGrad(Phi);
            surfaceScalarField magSf = mag(mesh.Sf());

            // Calculate J-Flux
            surfaceScalarField J1 = -DimlD1*(gradC1f + z1*C1f*gradPhif)*magSf + Pe*C1f*phi/(dimL*dimU);
            surfaceScalarField J2 = -DimlD2*(gradC2f + z2*C2f*gradPhif)*magSf + Pe*C2f*phi/(dimL*dimU);

            // Assign J-Flux to specific value
            forAll(C1.boundaryField(),patchI)
            {
                if(C1.boundaryField()[patchI].type() == "fixedIonicFlux")
                {
                    scalarField& j1f = J1.boundaryField()[patchI];
                    forAll(j1f, i)
                    {
                        j1f[i] = 0;
                    }
                }
            }

            forAll(C2.boundaryField(),patchI)
            {
                if(C2.boundaryField()[patchI].type() == "fixedIonicFlux")
                {
                    scalarField& j2f = J2.boundaryField()[patchI];
                    forAll(j2f, i)
                    {
                        j2f[i] = 0;
                    }
                }
            }

            // Calculate Right-Hand Side F
            scalar D = 1.0/Epsilon2;
            scalar iEpsilon = 1.0/Epsilon;
            volScalarField zC1Vol = z1*C1*cV*D/2.0;
            volScalarField zC2Vol = z2*C2*cV*D/2;
            volScalarField iEpsdC1dtVol = iEpsilon*fvc::ddt(C1)*cV;
            volScalarField iEpsdC2dtVol = iEpsilon*fvc::ddt(C2)*cV;

            volScalarField F1 = iEpsdC1dtVol*dimt*dimL + fvc::surfaceIntegrate(J1)*cV*sqr(dimL);
            volScalarField F2 = iEpsdC2dtVol*dimt*dimL + fvc::surfaceIntegrate(J2)*cV*sqr(dimL);
            volScalarField F3 = fvc::surfaceIntegrate(gradPhif*magSf)*cV*sqr(dimL) + zC1Vol + zC2Vol;
            /// Floating charged region:
            if(transportProperties.isDict("floatingChargeVolumes"))
            {
              dictionary fcVol=transportProperties.subDict("floatingChargeVolumes");
              wordList volumeList = fcVol.lookup("volumes");
              forAll(volumeList, ivol)
              {
                label volID = mesh.cellZones().findZoneID(volumeList[ivol]);
                if(volID>-1)
                {
                  const cellZone& ListCellsOfVol = mesh.cellZones()[volID];
                  scalar N = readScalar(fcVol.subDict(volumeList[ivol]).lookup("N"));
                  forAll(ListCellsOfVol, icells)
                  {
                    F3[ListCellsOfVol[icells]] += N*cV[ListCellsOfVol[icells]]*D/2;
                  }
                }
              }
            }
            ///
            // Jacobian-Matrix Assemble
            surfaceScalarField D1z1gradPhiSUf("D1z1gradPhiSUf", -DimlD1*z1*gradPhif*magSf + Pe*phi/(dimL*dimU));
            surfaceScalarField D2z2gradPhiSUf("D2z2gradPhiSUf", -DimlD2*z2*gradPhif*magSf + Pe*phi/(dimL*dimU));
            surfaceScalarField D1z1C1f("D1z1C1f", -DimlD1*z1*C1f);
            surfaceScalarField D2z2C2f("D2z2C2f", -DimlD2*z2*C2f);

            forAll(C1.boundaryField(),patchI)
            {
                if(C1.boundaryField()[patchI].type() == "fixedIonicFlux")
                {
                    scalarField& D1z1gradPhiSUf_ = D1z1gradPhiSUf.boundaryField()[patchI];
                    forAll(D1z1gradPhiSUf_, i) { D1z1gradPhiSUf_[i] = 0; }
                    scalarField& D1z1C1f_ = D1z1C1f.boundaryField()[patchI];
                    forAll(D1z1C1f_, i) { D1z1C1f_[i] = 0; }
                }
            }

            forAll(C2.boundaryField(),patchI)
            {
                if(C2.boundaryField()[patchI].type() == "fixedIonicFlux")
                {
                    scalarField& D2z2gradPhiSUf_ = D2z2gradPhiSUf.boundaryField()[patchI];
                    forAll(D2z2gradPhiSUf_, i) { D2z2gradPhiSUf_[i] = 0; }
                    scalarField& D2z2C2f_ = D2z2C2f.boundaryField()[patchI];
                    forAll(D2z2C2f_, i) { D2z2C2f_[i] = 0; }
                }
            }
            // f1 equation respect to C1
            dimensionedScalar two("two",dimensionSet(0,-2,1,0,0,0,0), 1);
            fvScalarMatrix f1C1Eqn
            (
                  two*iEpsilon*fvm::ddt(C1)
                + fvm::laplacian(-DimlD1,C1)
                + fvm::div(D1z1gradPhiSUf,C1)
            );

            f1C1Eqn.relax();
            // f1 equation respect to Phi
            fvScalarMatrix f1PhiEqn
            (
                fvm::laplacian(D1z1C1f,Phi)
            );

            f1PhiEqn.relax();
            // f2 equation respect to C2
            fvScalarMatrix f2C2Eqn
            (
                  two*iEpsilon*fvm::ddt(C2)
                + fvm::laplacian(-DimlD2,C2)
                + fvm::div(D2z2gradPhiSUf,C2)
            );

            f2C2Eqn.relax();
            // f2 equation respect to Phi
            fvScalarMatrix f2PhiEqn
            (
                fvm::laplacian(D2z2C2f,Phi)
            );

            f2PhiEqn.relax();
            // fPhi Equation
            fvScalarMatrix fPhiEqn
            (
                fvm::laplacian(Phi)
            );

            fPhiEqn.relax();
            // Prepare block system
            fvBlockMatrix<vector3> blockM(blockPNP);
            // Insert equations into block Matrix
            blockM.insertEquation(0, f1C1Eqn);
            blockM.insertEquation(1, f2C2Eqn);
            blockM.insertEquation(2, fPhiEqn);
            // Add off-diagonal coupling terms
            scalarField Phicoupling1(mesh.nCells(),z1.value()*D/2);
            scalarField Phicoupling2(mesh.nCells(),z2.value()*D/2);

            blockM.insertEquationCoupling(0, 2, f1PhiEqn);
            blockM.insertEquationCoupling(1, 2, f2PhiEqn);
            blockM.insertEquationCoupling(2, 0, Phicoupling1);
            blockM.insertEquationCoupling(2, 1, Phicoupling2);

            // Assign F to blockM's right-hand side

            forAll(blockM.source(),cellI)
            {
                blockM.source()[cellI](0) = -F1[cellI];
                blockM.source()[cellI](1) = -F2[cellI];
                blockM.source()[cellI](2) = -F3[cellI];
            }
            //- Block coupled solver call
            blockM.solve();
            // Retrieve solution
            blockM.retrieveSolution(0, C1.internalField());
            blockM.retrieveSolution(1, C2.internalField());
            blockM.retrieveSolution(2, Phi.internalField());
            //-Calculate and screen residual-length
            // Initialize residual-length

            scalar lengthC1(0);
            scalar lengthC2(0);
            scalar lengthPhi(0);

            // C1 residual-length

            forAll(C1.internalField(),i)
            {
                lengthC1 += sqr(C1.internalField()[i]);
            }

            // C2 residual-length

            forAll(C2.internalField(),i)
            {
                lengthC2 += sqr(C2.internalField()[i]);
            }

            // Phi residual-length

            forAll(Phi.internalField(),i)
            {
                lengthPhi += sqr(Phi.internalField()[i]);
            }

            Info <<" Iter # "<< i << " - Residual: " <<"[" << "C1:" << Foam::sqrt(lengthC1)
            << "  " << "C2:" << Foam::sqrt(lengthC2) << "  " << "Phi:" << Foam::sqrt(lengthPhi)
            << "]"  << endl;

            // Correct C1,C2,Phi

            C1T  += C1;
            C2T  += C2;
            PhiT += Phi;

            // Assign C1T,C2T,PhiT to C1,C2,Phi
            C1  = C1T;
            C2  = C2T;
            Phi = PhiT;

            // CorrectBoundaryCondition

            C1.correctBoundaryConditions();
            C2.correctBoundaryConditions();
            Phi.correctBoundaryConditions();

            // Correct Boundary value according to J-Flux
            forAll(C1.boundaryField(),patchI)
            {
                if(C1.boundaryField()[patchI].type() == "fixedIonicFlux")
                {
                    scalarField& Cb1 = C1.boundaryField()[patchI];
                    scalarField  Ci1 = C1.boundaryField()[patchI].patchInternalField();
                    scalarField& Phib1 = Phi.boundaryField()[patchI];
                    scalarField  Phii1 = Phi.boundaryField()[patchI].patchInternalField();
                    Cb1 = Ci1/(1.0 + z1.value()*(Phib1 - Phii1));
                }
            }

            forAll(C2.boundaryField(),patchI)
            {
                if(C2.boundaryField()[patchI].type() == "fixedIonicFlux")
                {
                    scalarField& Cb2 = C2.boundaryField()[patchI];
                    scalarField  Ci2 = C2.boundaryField()[patchI].patchInternalField();
                    scalarField& Phib2 = Phi.boundaryField()[patchI];
                    scalarField  Phii2 = Phi.boundaryField()[patchI].patchInternalField();
                    Cb2 = Ci2/(1.0 + z2.value()*(Phib2 - Phii2));
                }
            }
            scalar residual;
            residual = Foam::sqrt(lengthC1 + lengthC2 + lengthPhi);
            if (residual < convergenceCriterion)
            {
               break;
            }
        }

