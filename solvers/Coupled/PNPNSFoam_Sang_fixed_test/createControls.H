scalar maxResidual_PNP = 0;
scalar maxResidual_NS = 0;

//int    nPNPNSIter = mesh.solutionDict().subDict("SIMPLE").lookupOrDefault<int>("nPNPIter", 3);
scalar nPNPNSIter(readScalar(runTime.controlDict().lookup("nPNPNSIter")));
scalar nNewtonIterNS(readScalar(runTime.controlDict().lookup("nNewtonIterNS")));
scalar nNewtonIterPNP(readScalar(runTime.controlDict().lookup("nNewtonIterPNP")));
scalar convergenceCriterionNS(readScalar(mesh.solutionDict().subDict("solvers").subDict("Up").lookup("convergenceCriterionNS")));
scalar convergenceCriterionPNP(readScalar(mesh.solutionDict().subDict("solvers").subDict("PNP").lookup("convergenceCriterionPNP")));
scalar convergenceCriterionTotal(readScalar(mesh.solutionDict().lookup("convergenceCriterionTotal")));

//scalar nNSPNPIter(readScalar(runTime.controlDict().lookup("nNSPNPIter")));
