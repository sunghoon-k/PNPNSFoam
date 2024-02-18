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

Class
    Foam::fixedFluxFvPatchScalarField

Description
    The BC are intended to work with electrokinFOAM as no flux boundary conditions
    for the ions subject to the Nernst-Poisson-Planck eqauations dynamics
    They are specified as

    {
        type            fixedFlux;     // BC type
        n               nMinus;        // name of the ion field that they will refer to
        charge          -1;            // singed charge of ion
        D               1.91e-09;      // diffusion coefficient of the ion
        mu              7.4379062e-08; // mobility coefficient od the ion
        value           calculated;
    }
    Written by Pawel Jan Zuk

SourceFiles
    fixedFluxFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedFluxFvPatchScalarFields_H
#define fixedFluxFvPatchScalarFields_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fixedFluxPressureFvPatch Declaration
\*---------------------------------------------------------------------------*/

class fixedFluxFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    // Private data

        //- Name of the - Electrolyte field
        word cName_;

        //- charge of ions
        scalar z_;
        scalar T_;


public:

    //- Runtime type information
    TypeName("fixedFlux");


    // Constructors

        //- Construct from patch and internal field
        fixedFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedFluxFvPatchScalarField onto
        //  a new patch
        fixedFluxFvPatchScalarField
        (
            const fixedFluxFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedFluxFvPatchScalarField
        (
            const fixedFluxFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new fixedFluxFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedFluxFvPatchScalarField
        (
            const fixedFluxFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new fixedFluxFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
