{
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    /* Poisson Equation */

    // Add off-diagonal coupling terms
//    fvScalarMatrix Sp_fvScalarMatrix_cPlus(fvm::Sp((e*zPlus)/(eps0*epsr), cPlus));
//    fvScalarMatrix Sp_fvScalarMatrix_cMinus(fvm::Sp((e*zMinus)/(eps0*epsr), cMinus));
//    PNPEqn.insertEquationCoupling(0, 1, Sp_fvScalarMatrix_cPlus); // 4: psiE <- 5: cPlus
//    PNPEqn.insertEquationCoupling(0, 2, Sp_fvScalarMatrix_cMinus); // 4: psiE <- 6: cMinus


    scalarField PhicouplingPlus(mesh.nCells(),((e*zPlus)/(eps0*epsr)).value());
    scalarField PhicouplingMinus(mesh.nCells(),((e*zMinus)/(eps0*epsr)).value());
    //PNPEqn.insertEquationCoupling(0, 1, PhicouplingPlus); // 0: psiE <- 1: cPlus
    //PNPEqn.insertEquationCoupling(0, 2, PhicouplingMinus); // 0: psiE <- 2: cMinus

//    fvScalarMatrix interpolatecPlusDPluszPlus(fvm::laplacian(-(DPlus*zPlus)*fvc::interpolate(cPlus)/psiE0, psiE));
//    fvScalarMatrix interpolatecMinusDMinuszMinus(fvm::laplacian(-(DMinus*zMinus)*fvc::interpolate(cMinus)/psiE0, psiE));
    fvScalarMatrix interpolatecPlusDPluszPlus(fvm::laplacian(-(DPlus*zPlus)*(cPlus)/psiE0, psiE));
    fvScalarMatrix interpolatecMinusDMinuszMinus(fvm::laplacian(-(DMinus*zMinus)*(cMinus)/psiE0, psiE));
    //PNPEqn.insertEquationCoupling(1, 0, interpolatecPlusDPluszPlus); // 0: psiE -> 1: cPlus
    //PNPEqn.insertEquationCoupling(2, 0, interpolatecMinusDMinuszMinus); // 0: psiE -> 2: cMinus

}
