// Check convergence
if ((maxResidual < convergenceNS))//  && !bounded)
{
    Info<< "reached convergence criterion: " << convergenceNS << endl;
//    runTime.writeAndEnd();
    // Info<< "latestTime = " << runTime.timeName() << endl;
    break;
}

