/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "zeroIonicFluxFvPatchScalarField.H"
// #include "EDFEquation.H"
#include "ElectrochemicalSystem.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zeroIonicFluxFvPatchScalarField::zeroIonicFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    zib_(),
    T_()
{}


Foam::zeroIonicFluxFvPatchScalarField::zeroIonicFluxFvPatchScalarField
(
    const zeroIonicFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    zib_(ptf.zib_),
    T_(ptf.T_)
{}


Foam::zeroIonicFluxFvPatchScalarField::zeroIonicFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    zib_(0),
    T_(0)
{

    if (dict.found("value") && dict.found("gradient"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
        gradient() = scalarField("gradient", dict, p.size());
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
        gradient() = 0.0;
    }

/*
    const dictionary& elecDict = db().lookupObject<IOdictionary>("electricProperties");
        
    PtrList<entry> specEntries_(elecDict.subDict("parameters").lookup("species"));
    
    forAll (specEntries_, specI)
    {    
       if ( specEntries_[specI].keyword() == this->dimensionedInternalField().name() )
        { 
          dimensionedScalar zid_(specEntries_[specI].dict().lookup("z"));
          zib_ = zid_.value();
          break;
        }
    } 
*/    

}


Foam::zeroIonicFluxFvPatchScalarField::zeroIonicFluxFvPatchScalarField
(
    const zeroIonicFluxFvPatchScalarField& wbppsf
)
:
    fixedGradientFvPatchScalarField(wbppsf),
    zib_(wbppsf.zib_),
    T_(wbppsf.T_)
{}


Foam::zeroIonicFluxFvPatchScalarField::zeroIonicFluxFvPatchScalarField
(
    const zeroIonicFluxFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(wbppsf, iF),
    zib_(wbppsf.zib_),
    T_(wbppsf.T_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::zeroIonicFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;       
    }

/*
    const dictionary& elecDict = db().lookupObject<IOdictionary>("electricProperties");
    dimensionedScalar T_(elecDict.subDict("parameters").lookup("T"));

    const volScalarField& phiE_ =
        db().lookupObject<volScalarField>("phiE");
*/    
    scalar e = (ElectrochemicalSystem::e).value();
    scalar kB = (ElectrochemicalSystem::kB).value();
    scalar zib = zib_;
    scalar T = T_;
    
    const volScalarField& psiE_ =
        db().lookupObject<volScalarField>("psiE");
        
    scalarField Epatch( (psiE_).boundaryField()[patch().index()].snGrad() );
           
    gradient() = -(*this) * e * zib * Epatch / (kB * T);      
     
    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::zeroIonicFluxFvPatchScalarField::evaluate(const Pstream::commsTypes)
{

/*
    const dictionary& elecDict = db().lookupObject<IOdictionary>("electricProperties");
        
    dimensionedScalar T_(elecDict.subDict("parameters").lookup("T"));
*/        
    
    const fvPatchField<scalar>& psib_ = patch().lookupPatchField<volScalarField, scalar>("psiE");
    
    scalarField deltaPsi = psib_ - psib_.patchInternalField();
     
    scalar e = (ElectrochemicalSystem::e).value();
    scalar kB = (ElectrochemicalSystem::kB).value();
    scalar zib = zib_;
    scalar T = T_;
    
    Field<Foam::scalar>::operator=
    (
        this->patchInternalField() * Foam::exp( -e * zib * deltaPsi/ (kB * T) ) 
    );
    
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    fvPatchField<Foam::scalar>::evaluate();
}

// Render the BC fully-explicit
 
Foam::tmp<Field<Foam::scalar> > Foam::zeroIonicFluxFvPatchScalarField::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    // return tmp<Field<Foam::scalar> >(new Field<Foam::scalar>(this->size(), pTraits<Foam::scalar>::zero));
    return (*this);
}

// Render the BC fully-explicit
 
Foam::tmp<Field<Foam::scalar> > Foam::zeroIonicFluxFvPatchScalarField::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    // return (*this)/this->patchInternalField();
    return tmp<Field<Foam::scalar> >(new Field<Foam::scalar>(this->size(), pTraits<Foam::scalar>::zero));
}

 
void Foam::zeroIonicFluxFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        zeroIonicFluxFvPatchScalarField
    );
}

// ************************************************************************* //
