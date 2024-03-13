surfaceScalarField DPluszPlusphiE("DPluszPlusphiE", DPlus_nonD*zPlus*gradpsiEf*magSf);
surfaceScalarField DPluszPluscPlusf("DPluszPluscPlusf", DPlus_nonD*zPlus*cPlusf);

forAll(cPlus.boundaryField(),patchI)
{
    if(cPlus.boundaryField()[patchI].type() == "zeroIonicFlux")
    {
        scalarField& DPluszPlusphiE_ = DPluszPlusphiE.boundaryField()[patchI];
        forAll(DPluszPlusphiE_, i) { DPluszPlusphiE_[i] = 0; }
        scalarField& DPluszPluscPlusf_ = DPluszPluscPlusf.boundaryField()[patchI];
        forAll(DPluszPluscPlusf_, i) { DPluszPluscPlusf_[i] = 0; }
    }
}

surfaceScalarField DMinuszMinusphiE("DMinuszMinusphiE", DMinus_nonD*zMinus*gradpsiEf*magSf);
surfaceScalarField DMinuszMinuscMinusf("DMinuszMinuscMinusf", DMinus_nonD*zMinus*cMinusf);

forAll(cMinus.boundaryField(),patchI)
{
    if(cMinus.boundaryField()[patchI].type() == "zeroIonicFlux")
    {
        scalarField& DMinuszMinusphiE_ = DMinuszMinusphiE.boundaryField()[patchI];
        forAll(DMinuszMinusphiE_, i) { DMinuszMinusphiE_[i] = 0; }
        scalarField& DMinuszMinuscMinusf_ = DMinuszMinuscMinusf.boundaryField()[patchI];
        forAll(DMinuszMinuscMinusf_, i) { DMinuszMinuscMinusf_[i] = 0; }
    }
}

fvScalarMatrix cPlusEqn
(
  - fvm::laplacian(DPlus_nonD, cPlus)
  - fvm::div(DPluszPlusphiE, cPlus)
  //- fvm::div(mobilityPlus * phiE, cPlus)
  //- fvc::laplacian(DPluszPluscPlusf, psiE)
  //+ fvc::div(mobilityPlus * phiE, cPlus)
  //==
  //- FcPlus
);
cPlusEqn.relax();

fvScalarMatrix cMinusEqn
(
  - fvm::laplacian(DMinus_nonD, cMinus)
  - fvm::div(DMinuszMinusphiE, cMinus)
  //- fvm::div(mobilityMinus * phiE, cMinus)
  //- fvc::laplacian(DMinuszMinuscMinusf, psiE)
  //+ fvc::div(mobilityPlus * phiE, cPlus)
  //==
  //- FcMinus
);
cMinusEqn.relax();

if(solveTransient)
{
  cPlusEqn += fvm::ddt(cPlus);
  cMinusEqn += fvm::ddt(cMinus);
}

/*
  if (currentIter == nOuterIter -1)
  {
    cPlusEqn.relax(1);
  }
  else
  {
    cPlusEqn.relax();
  }
*/

