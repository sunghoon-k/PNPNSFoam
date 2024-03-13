{
/**************************************************************************/

    // Bound the concentration
    dimensionedScalar dimensioned_C1 = min(C1);    
    dimensionedScalar dimensioned_C2 = min(C2);

    if (dimensioned_C1.value() < 0) // when the concentration is negative
    {
        Info<< "C1: " << dimensioned_C1.value() << " "
            << ".  Bounding." << endl;

        //C1.max(C1Min); // max(C1, C1Min)
        //C1.correctBoundaryConditions();
        //bounded = true;
    }

    if (dimensioned_C2.value() < 0) // when the concentration is negative
    {
        Info<< "C2: " << dimensioned_C2.value() << " "
            << ".  Bounding." << endl;

        //C2.max(C2Min); // max(C1, C1Min)
        //C2.correctBoundaryConditions();
        //bounded = true;
    }

}
