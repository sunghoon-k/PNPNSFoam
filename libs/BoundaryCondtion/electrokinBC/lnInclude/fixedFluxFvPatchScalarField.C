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

#include "fixedFluxFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"
#include "ElectrochemicalSystem.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::fixedFluxFvPatchScalarField::fixedFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    cName_("c"),
    z_(1.0),
    T_(300)
{}

Foam::fixedFluxFvPatchScalarField::fixedFluxFvPatchScalarField
(
    const fixedFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    cName_(ptf.cName_),
    z_(ptf.z_),
    T_(ptf.T_)
{}

Foam::fixedFluxFvPatchScalarField::fixedFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),    
    cName_(dict.lookup("c")),
    z_(readScalar(dict.lookup("z"))),
    T_(readScalar(dict.lookup("T")))

{
    if (dict.found("gradient")) 
    {
        gradient() = scalarField("gradient", dict, p.size());
        fixedGradientFvPatchScalarField::updateCoeffs();
        fixedGradientFvPatchScalarField::evaluate();
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
        gradient() = 0.0;
    }
}

Foam::fixedFluxFvPatchScalarField::fixedFluxFvPatchScalarField
(
    const fixedFluxFvPatchScalarField& wbppsf
)
:
    fixedGradientFvPatchScalarField(wbppsf),
    cName_(wbppsf.cName_),
    z_(wbppsf.z_),
    T_(wbppsf.T_)

{}

Foam::fixedFluxFvPatchScalarField::fixedFluxFvPatchScalarField
(
    const fixedFluxFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(wbppsf, iF),
    cName_(wbppsf.cName_),
    z_(wbppsf.z_),
    T_(wbppsf.T_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    scalar kB = (ElectrochemicalSystem::kB).value();
    scalar e = (ElectrochemicalSystem::e).value();
    scalar psiE0 = kB*T_/e;

    const fvPatchField<scalar>& cp =
        patch().lookupPatchField<volScalarField, vector>(cName_);

    const volScalarField& psiE_ =
        db().lookupObject<volScalarField>("psiE");
        
    scalarField Epatch( psiE_.boundaryField()[patch().index()].snGrad() );

    gradient() =  - z_ * cp * Epatch/psiE0;

    fixedGradientFvPatchScalarField::updateCoeffs();

}


void Foam::fixedFluxFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntry("value",os);
    os.writeKeyword("c") << cName_ << token::END_STATEMENT << nl;
    os.writeKeyword("z") << z_ << token::END_STATEMENT << nl;
    os.writeKeyword("T") << T_ << token::END_STATEMENT << nl;

    gradient().writeEntry("gradient", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixedFluxFvPatchScalarField
    );
}

// ************************************************************************* //
