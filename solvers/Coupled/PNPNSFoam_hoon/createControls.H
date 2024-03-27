scalar maxResidual = 0;
scalar maxResidual_PNP = 0;

scalar nPNPNSIter(readScalar(runTime.controlDict().lookup("nPNPNSIter")));
scalar nNewtonIter(readScalar(runTime.controlDict().lookup("nNewtonIter")));
scalar nNewtonIterPNP(readScalar(runTime.controlDict().lookup("nNewtonIterPNP")));

scalar convergenceCriterion(readScalar(mesh.solutionDict().subDict("solvers").subDict("PNPNS").lookup("convergenceCriterion")));
scalar convergenceCriterionPNP(readScalar(mesh.solutionDict().subDict("solvers").subDict("PNP").lookup("convergenceCriterionPNP")));
