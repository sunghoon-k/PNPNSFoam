// Check convergence
if (maxResidual > convergenceCriterion)
{
    Info<< "not reached convergence criterion: " << convergenceCriterion << endl;
    runTime.writeAndEnd();
    Info<< "latestTime = " << runTime.timeName() << endl;
}

