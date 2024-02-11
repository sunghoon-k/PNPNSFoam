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
    slipEOFFvPatchVectorField

Description
    Slip velocity following the Helmoltz-EOF theory.
    It considers: E=-grad(phiE). 
    
    Example:
    
      {
        type            slipEOF;
	elecMobility    4e-8;              // electroosmotic mobility (U/E)
	value           uniform (0 0 0);   // place-holder    
      }
      
    This BC is part of rheoTool.

SourceFiles
    slipEOFFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef slipEOFFvPatchVectorField_H
#define slipEOFFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class slipEOFFvPatch Declaration
\*---------------------------------------------------------------------------*/

class slipEOFFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{

    // Private data
      
       //- Electric mobility (=U/E)
       
        scalar zetaPotential_;
        scalar epsr_;
        scalar nu_;
        scalar rho_;

public:

    //- Runtime type information
    TypeName("slipEOF");


    // Constructors

        //- Construct from patch and internal field
        slipEOFFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        slipEOFFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given slipEOFFvPatchVectorField
        //  onto a new patch
        slipEOFFvPatchVectorField
        (
            const slipEOFFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
            
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new slipEOFFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        slipEOFFvPatchVectorField
        (
            const slipEOFFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new slipEOFFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

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
