
namespace Foam
{
	defineTypeNameAndDebug(petscSolver, 0);
}

Foam::petscSolver::petscSolver(int argc, char *argv[])
{
	PetscInitialize(&argc, &argv, NULL, NULL);
}

Foam::petscSolver::~petscSolver()
{
	PetscFinalize();
}

template<class Type>
void Foam::petscSolver::solve
(
    fvMatrix<Type>& ofMat,
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    GeometricField<Type, fvPatchField, volMesh>& T =
    const_cast< GeometricField<Type, fvPatchField, volMesh>& >
    (
        ofMat.psi()
    );

    typename Type::labelType validComponents
    (
        pow
        (
            T.mesh().solutionD(),
            pTraits<typename powProduct<Vector<label>, Type::rank>::type>::zero
        )
    );

    label n = ofMat.psi().mesh().C().size();

    Info<<"============= Solve fvMatrix using PETSC Solver ============="<<nl<<endl;
    Info<<"N cells = "<<n<<nl<<endl;

    Mat A;
    Vec rhs, sol;
    KSP ksp;
    PC  pc;


    MatCreate(PETSC_COMM_WORLD, &A);
    VecCreate(PETSC_COMM_WORLD, &rhs);
    VecCreate(PETSC_COMM_WORLD, &sol);
    KSPCreate(PETSC_COMM_WORLD, &ksp);

    for (direction cmpt=0; cmpt<Type::nComponents; cmpt++)
    {
        if (validComponents[cmpt] == -1) continue;
        else
        {
            Info<<" Solving for direction "<<cmpt<<nl<<endl;
            if (cmpt==0)
            {
                MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n,n);
                MatSetFromOptions(A);
                MatSetUp(A);
                VecSetSizes(rhs, PETSC_DECIDE, n);
                VecSetFromOptions(rhs);
                VecDuplicate(rhs, &sol);
            }
            else
            {
                MatZeroEntries(A);
                VecSet(rhs,0.);
            }

        //- Matrix coefficient
            //- Off diag elements
                const lduAddressing& offDiag = ofMat.lduAddr();
                label nFaces = offDiag.lowerAddr().size();

                label col, row;

                for (label facei=0; facei<nFaces; facei++)
                {
                    row = offDiag.upperAddr()[facei];
                    col = offDiag.lowerAddr()[facei];
                    MatSetValues(A,1,&row,1,&col,&ofMat.lower()[facei],INSERT_VALUES);
                    MatSetValues(A,1,&col,1,&row,&ofMat.upper()[facei],INSERT_VALUES);
                }

            //- Diag elements
                for (label celli=0; celli<n; celli++)
                {
                    row = celli;
                    MatSetValues(A,1,&row,1,&row,&ofMat.diag()[celli],INSERT_VALUES);
                }

        //- RHS vector
            scalarField source = ofMat.source().component(cmpt);
            source += vf.mesh().V()*vf.field().component(cmpt);

            for (label celli=0; celli<n; celli++)
            {
                row = celli;
                VecSetValues(rhs,1,&row,&source[celli],INSERT_VALUES);
            }

            MatAssemblyBegin(A, MAT_FLUSH_ASSEMBLY);
            MatAssemblyEnd(A, MAT_FLUSH_ASSEMBLY);
            VecAssemblyBegin(rhs);
            VecAssemblyEnd(rhs);

        // - Contribution from BCs
            forAll(T.boundaryField(), patchi)
            {
                scalarField bC = ofMat.boundaryCoeffs()[patchi].component(cmpt);
                scalarField iC = ofMat.internalCoeffs()[patchi].component(cmpt);
                const unallocLabelList& addr = ofMat.lduMatrix::lduAddr().patchAddr(patchi);

                if (!T.boundaryField()[patchi].coupled())
                {
                    forAll(addr, facei)
                    {
                        // Matrix of coefs - Diagonal
                        row = addr[facei];
                        MatSetValues(A,1,&row,1,&row,&iC[facei],ADD_VALUES);

                        // Source vector
                        VecSetValues(rhs,1,&row,&bC[facei],ADD_VALUES);
                    }
                }
            }

            MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
            VecAssemblyBegin(rhs);
            VecAssemblyEnd(rhs);
            Info<<"Converting fvMatrix to Mat : Ok-DONE"<<nl<<endl;
        //- KSP solver
            label maxKspIter=500;
            PetscScalar maxKspTol = 1e-06;
            PetscScalar maxKspRTol= 1e-10;
            KSPSetOperators(ksp, A, A, DIFFERENT_NONZERO_PATTERN);
            KSPSetFromOptions(ksp);
            KSPSetTolerances(ksp,maxKspRTol,maxKspTol,PETSC_DEFAULT,maxKspIter);
            KSPGetPC(ksp,&pc);
            PCSetType(pc,PCJACOBI);
            KSPSolve(ksp,rhs,sol);
            Info<<"Solving the system with maxKspTol = "<<maxKspTol<<", maxKspRTol = "<<maxKspRTol<<" : Ok"<<nl<<endl;

        //- Transfer solution
            scalarField transSol(n, 0.);
            const PetscScalar *xx;
            VecGetArrayRead(sol, &xx);
            for (int i=0; i<n; i++)
            {
                transSol[i] = xx[i];
            }
            VecRestoreArrayRead(sol,&xx);
            T.internalField().replace(cmpt, transSol);
            Info<<"Transfering solution : Ok"<<nl<<endl;
        }
    }

    //- Destroy
    MatDestroy(&A);
    KSPDestroy(&ksp);
    VecDestroy(&rhs);  VecDestroy(&sol);
    Info<<"Destroyed petsc data."<<nl<<endl;
}

template<class Type>
void Foam::petscSolver::solve
(
    fvMatrix<Type>& ofMat,
    tmp<GeometricField<Type, fvPatchField, volMesh> > tvf
)
{
    solve(ofMat, tvf());
}

template<class Type>
void Foam::petscSolver::solve
(
    tmp<fvMatrix<Type> > ofMat,
    tmp<GeometricField<Type, fvPatchField, volMesh> > tvf
)
{
    solve(ofMat(), tvf());
}

template<class Type>
void Foam::petscSolver::solve
(
    fvBlockMatrix<Type>& ofBlockMat
)
{
    GeometricField<Type, fvPatchField, volMesh>& T =
    const_cast< GeometricField<Type, fvPatchField, volMesh>& >
    (
        ofBlockMat.psi()
    );

    const direction nCmpts = pTraits<Type>::nComponents;
    label nCells = T.mesh().nCells();
    label matSizeGlobal = nCells * nCmpts;

    Info<<"============ Solve fvBlockMatrix using PETSC Solver ============"<<endl;
    Info<<"Block size = "<<nCmpts<<endl;
    Info<<"N cells = "<<nCells<<nl<<endl;

    Mat A;
    Vec rhs, sol, Ax;
    KSP ksp;
    PC  pc;

    MatCreate(PETSC_COMM_WORLD, &A);
    VecCreate(PETSC_COMM_WORLD, &rhs);
    VecCreate(PETSC_COMM_WORLD, &sol);

    MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, matSizeGlobal, matSizeGlobal);
    MatSetBlockSize(A, nCmpts);
    MatSetFromOptions(A);

    const unallocLabelList& owner = T.mesh().owner();
    const unallocLabelList& neibr = T.mesh().neighbour();

    labelList nnzblock(nCells, 1);
    forAll(owner, facei)
    {
        nnzblock[owner[facei]] += 1;
        nnzblock[neibr[facei]] += 1;
    }

    PetscInt nz = 0, nnz[matSizeGlobal];
    for (int r=0; r < matSizeGlobal; r++)
    {
        scalar k = r/nCmpts;
        label blockRowId = static_cast<int>(k);
        nnz[r] = nnzblock[blockRowId] * nCmpts;
    }

    MatSeqAIJSetPreallocation(A, nz, nnz);

    VecSetSizes(rhs, PETSC_DECIDE, matSizeGlobal);
    VecSetBlockSize(rhs, nCmpts);
    VecSetFromOptions(rhs);
    VecDuplicate(rhs, &sol);
    VecDuplicate(rhs, &Ax);

    scalar tTemp, dT;

    //- Start assembling matrix

    tTemp = T.mesh().time().elapsedCpuTime();

    //- 1. Off diag elements
        const lduAddressing& offDiag = ofBlockMat.lduAddr();
        label nFaces = offDiag.lowerAddr().size();
        label col, row;

        for (label facei=0; facei<nFaces; facei++)
        {
            row = offDiag.upperAddr()[facei];
            col = offDiag.lowerAddr()[facei];
            if (ofBlockMat.lower().activeType() == 3)
            {
                PetscScalar upperValue[ofBlockMat.upper().asSquare()[facei].size()];
                PetscScalar lowerValue[ofBlockMat.lower().asSquare()[facei].size()];
                forAll(ofBlockMat.lower().asSquare()[facei], valI)
                {
                    upperValue[valI] = ofBlockMat.upper().asSquare()[facei][valI];
                    lowerValue[valI] = ofBlockMat.lower().asSquare()[facei][valI];
                }
                MatSetValuesBlocked(A,1,&row,1,&col,lowerValue,INSERT_VALUES);
                MatSetValuesBlocked(A,1,&col,1,&row,upperValue,INSERT_VALUES);
            }
            else
            {
                Info<<"OFF-DIAG BLOCKS ARE NOT SQUARE !!!"<<nl<<endl;
            }
        }

    //- 2. Diag elements
        for (label celli=0; celli<nCells; celli++)
        {
            if (ofBlockMat.diag().activeType() == 3)
            {
                label blockSize = ofBlockMat.diag().asSquare()[celli].size();
                PetscScalar diagValue[blockSize];
                for (label valI=0; valI<blockSize; valI++)
                {
                    diagValue[valI] = ofBlockMat.diag().asSquare()[celli][valI];
                }

                MatSetValuesBlocked(A,1,&celli,1,&celli,diagValue,INSERT_VALUES);
            }
            else
            {
                Info<<"OFF-DIAG BLOCKS ARE NOT SQUARE !!!"<<nl<<endl;
            }
        }

    //- 3. RHS vector
        Field<Type> source = ofBlockMat.source();
        for (label celli=0; celli<nCells; celli++)
        {
            PetscScalar sourceValue[source.size()];
            forAll(source[valI], valI)
            {
                sourceValue[valI] = source[celli][valI];
            }
            VecSetValuesBlocked(rhs,1,&celli,sourceValue,INSERT_VALUES);
        }

    //- 4. Finalize matrix assembly
        MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

        VecAssemblyBegin(rhs);
        VecAssemblyEnd(rhs);

        dT = T.mesh().time().elapsedCpuTime() - tTemp;
        Info<<"Converting fvBlockMatrix to Mat : Ok (ExecutionTime = "<<dT<<" s)"<<nl<<endl;
        tTemp = T.mesh().time().elapsedCpuTime();

    //- KSP solver
        dictionary PetscSolver(T.mesh().solutionDict().subDict("PetscSolver"));
        scalar maxKspIter(readScalar(PetscSolver.lookup("maxKspIter")));
        scalar maxKspTol(readScalar(PetscSolver.lookup("maxKspTol")));
        scalar maxKspRTol(readScalar(PetscSolver.lookup("maxKspRTol")));
        PetscInt iter;
        PetscReal PetscSolverError;

        KSPCreate(PETSC_COMM_WORLD, &ksp);
        KSPSetOperators(ksp, A, A, DIFFERENT_NONZERO_PATTERN);
        KSPGetPC(ksp,&pc);
        PCSetType(pc,PCLU);
        KSPSetType(ksp,KSPBCGS);
        KSPSetTolerances(ksp,maxKspRTol,maxKspTol,PETSC_DEFAULT,maxKspIter);
        KSPSetFromOptions(ksp);
        KSPSolve(ksp,rhs,sol);
        KSPGetIterationNumber(ksp,&iter);
        KSPGetResidualNorm(ksp, &PetscSolverError);

        dT = T.mesh().time().elapsedCpuTime() - tTemp;
        Info<<"Solving the system"<<endl
            <<"   * MaxKspTol = "<<maxKspTol<<endl
            <<"   * MaxKspRTol = "<<maxKspRTol<<endl
            <<"   * Norm = "<<PetscSolverError<<". Iteration = "<<iter<<" (ExecutionTime = "<<dT<<" s)"<<nl<<endl;
        tTemp = T.mesh().time().elapsedCpuTime();

    //- Transfer solution

        //- Conver solution vector to scalarField
        scalarField transSol(matSizeGlobal, 0.);
        const PetscScalar *xx;
        VecGetArrayRead(sol, &xx);
        forAll(transSol, i)
        {
            transSol[i] = xx[i];
        }
        VecRestoreArrayRead(sol,&xx);

        //- Split converted solution
        List<scalarField> splitSolution(nCmpts);
        forAll(splitSolution, i)
        {
            splitSolution[i].setSize(nCells);
        }

        label idx=0;
        forAll(transSol, i)
        {
            label k = i % nCmpts;
            if(i != 0 and k == 0) idx++;

            splitSolution[k][idx] = transSol[i];
        }

        //- Replace components of T by splitted solution
        for (direction cmpt=0; cmpt<nCmpts; cmpt++)
        {
            T.internalField().replace(cmpt, splitSolution[cmpt]);
        }

        dT = T.mesh().time().elapsedCpuTime() - tTemp;
        Info<<"Transfering solution : Ok (ExecutionTime = "<<dT<<" s)"<<nl<<endl;

    //- Destroy
    MatDestroy(&A);
    KSPDestroy(&ksp);
    VecDestroy(&rhs);  VecDestroy(&sol); VecDestroy(&Ax);
    Info<<"Destroyed petsc data."<<endl;
    Info<<"================================================================"<<nl<<endl;

}
