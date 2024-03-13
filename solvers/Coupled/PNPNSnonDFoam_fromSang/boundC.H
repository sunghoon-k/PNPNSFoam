{
/**************************************************************************/

    // Bound the concentration
    dimensionedScalar dimensioned_cPlus = min(cPlus);    
    dimensionedScalar dimensioned_cMinus = min(cMinus);

    if (dimensioned_cPlus < cPlusMin) // when the concentration is negative
    {
        Info<< "cPlus: " << dimensioned_cPlus.value() << " "
            << ".  Bounding." << endl;

        //cPlus.max(cPlusMin); // max(cPlus, cPlusMin)
        //cPlus.correctBoundaryConditions();
        bounded = true;
    }

    if (dimensioned_cMinus < cMinusMin) // when the concentration is negative
    {
        Info<< "cMinus: " << dimensioned_cMinus.value() << " "
            << ".  Bounding." << endl;

        //cMinus.max(cMinusMin); // max(cPlus, cPlusMin)
        //cMinus.correctBoundaryConditions();
        bounded = true;
    }

}
