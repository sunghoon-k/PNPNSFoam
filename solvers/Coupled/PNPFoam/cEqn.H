{
  phiE = fvc::snGrad(psiE) * mesh.magSf();

  surfaceScalarField DPluszPlusphiE("DPluszPlusphiE", DPlus*zPlus*phiE/psiE0);
  surfaceScalarField DPluszPluscPlusf("DPluszPluscPlusf", DPlus*zPlus*cPlusf/psiE0);

  forAll(cPlus.boundaryField(),patchI)
  {
      if(cPlus.boundaryField()[patchI].type() == "fixedIonicFlux")
      {
          scalarField& DPluszPlusphiE_ = DPluszPlusphiE.boundaryField()[patchI];
          forAll(DPluszPlusphiE_, i) { DPluszPlusphiE_[i] = 0; }
          scalarField& DPluszPluscPlusf_ = DPluszPluscPlusf.boundaryField()[patchI];
          forAll(DPluszPluscPlusf_, i) { DPluszPluscPlusf_[i] = 0; }
      }
  }

  surfaceScalarField DMinuszMinusphiE("DMinuszMinusphiE", DMinus*zMinus*phiE/psiE0);
  surfaceScalarField DMinuszMinuscMinusf("DMinuszMinuscMinusf", DMinus*zMinus*cMinusf/psiE0);

  forAll(cMinus.boundaryField(),patchI)
  {
      if(cMinus.boundaryField()[patchI].type() == "fixedIonicFlux")
      {
          scalarField& DMinuszMinusphiE_ = DMinuszMinusphiE.boundaryField()[patchI];
          forAll(DMinuszMinusphiE_, i) { DMinuszMinusphiE_[i] = 0; }
          scalarField& DMinuszMinuscMinusf_ = DMinuszMinuscMinusf.boundaryField()[patchI];
          forAll(DMinuszMinuscMinusf_, i) { DMinuszMinuscMinusf_[i] = 0; }
      }
  }

  fvScalarMatrix cPlusEqn
  (
    - fvm::laplacian(DPlus, cPlus)
    - fvm::div(DPluszPlusphiE, cPlus)
    //- fvm::div(mobilityPlus * phiE, cPlus)
    //- fvc::laplacian(mobilityPlus * cPlus, psiE)
    //+ fvc::div(mobilityPlus * phiE, cPlus)
  );

  fvScalarMatrix cMinusEqn
  (
    - fvm::laplacian(DMinus, cMinus)
    - fvm::div(DMinuszMinusphiE, cMinus)
    //- fvm::div(mobilityMinus * phiE, cMinus)
    //- fvc::laplacian(mobilityMinus * cMinus, psiE)
    //+ fvc::div(mobilityPlus * phiE, cPlus)

  );

  if(solveTransient)
  {
    cPlusEqn += fvm::ddt(cPlus);
    cMinusEqn += fvm::ddt(cMinus);
  }

  cPlusEqn.relax();
  cMinusEqn.relax();

  PNPEqn.insertEquation(0, cPlusEqn);
  PNPEqn.insertEquation(1, cMinusEqn);

  fvScalarMatrix electroMigrationLaplacianPlus
  (
    fvm::laplacian(-DPluszPluscPlusf, psiE)
  );
  electroMigrationLaplacianPlus.relax();

  fvScalarMatrix electroMigrationLaplacianMinus
  (
    fvm::laplacian(-DMinuszMinuscMinusf, psiE)
  );
  electroMigrationLaplacianMinus.relax();

  PNPEqn.insertEquationCoupling(0, 2, electroMigrationLaplacianPlus); // 0: psiE -> 1: cPlus
  PNPEqn.insertEquationCoupling(1, 2, electroMigrationLaplacianMinus); // 0: psiE -> 2: cMinus

}