/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "slipEOFFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvCFD.H"
#include "ElectrochemicalSystem.H"
// #include "EDFEquation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::slipEOFFvPatchVectorField::
slipEOFFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    zetaPotential_(),
    epsr_(),
    mu_(),
    rho_()
//    elecM_()
{}

Foam::slipEOFFvPatchVectorField::
slipEOFFvPatchVectorField
(
    const slipEOFFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
    
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    zetaPotential_(ptf.zetaPotential_),
    epsr_(ptf.epsr_),
    mu_(ptf.mu_),
    rho_(ptf.rho_)
//    elecM_(ptf.elecM_)   
{}

Foam::slipEOFFvPatchVectorField::
slipEOFFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    zetaPotential_(readScalar(dict.lookup("zetaPotential"))),
    epsr_(readScalar(dict.lookup("epsr"))),
    mu_(readScalar(dict.lookup("mu"))),
    rho_(readScalar(dict.lookup("rho")))
//    elecM_(readScalar(dict.lookup("elecMobility")))
{
    fvPatchField<vector>::operator=
    (
        vectorField("value", dict, p.size())
    );
    
}
    
Foam::slipEOFFvPatchVectorField::
slipEOFFvPatchVectorField
(
    const slipEOFFvPatchVectorField& tppsf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(tppsf, iF),
    zetaPotential_(tppsf.zetaPotential_),
    epsr_(tppsf.epsr_),
    mu_(tppsf.mu_),
    rho_(tppsf.rho_)
//    elecM_(tppsf.elecM_) 
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::slipEOFFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    const volScalarField& psiE_ =
        db().lookupObject<volScalarField>("psiE");
        
    volVectorField Ef(-fvc::grad(psiE_));
    
    // independent constant
    dimensionedScalar eps0 = ElectrochemicalSystem::eps0;
    scalar zetaPotential = zetaPotential_;
    scalar epsr = epsr_;
    scalar mu = mu_;
    scalar rho = rho_;
    scalar nu = mu/rho;
            
    scalar elecM = (eps0.value()*epsr)*zetaPotential/nu;

    vectorField::operator=( elecM * Ef.boundaryField()[patch().index()] );
         
    fixedValueFvPatchVectorField::updateCoeffs();
}

 
void Foam::slipEOFFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("zetaPotential")
        << zetaPotential_ << token::END_STATEMENT << nl;
    os.writeKeyword("epsr")
        << epsr_ << token::END_STATEMENT << nl;
    os.writeKeyword("mu")
        << mu_ << token::END_STATEMENT << nl;
    os.writeKeyword("rho")
        << rho_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        slipEOFFvPatchVectorField
    );
}

// ************************************************************************* //
