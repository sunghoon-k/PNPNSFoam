// Check convergence
if ((maxResidual < convergenceCriterion) && (sumLocalContErr < convergenceCriterion) && !bounded)
{
    Info<< "reached convergence criterion: " << convergenceCriterion << endl;
//    runTime.writeAndEnd();
    // Info<< "latestTime = " << runTime.timeName() << endl;
    //convergenceInnerLoop = true;
    break;
}

