/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Description

\*---------------------------------------------------------------------------*/

#include "writeFuns.H"
#include "interpolatePointToCell.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Store List in dest
template<class Type>
void Foam::writeFuns::insert
(
    const List<Type>& source,
    DynamicList<floatScalar>& dest
)
{
    forAll(source, i)
    {
        insert(source[i], dest);
    }
}


template<class Type>
void Foam::writeFuns::write
(
    std::ostream& os,
    const bool binary,
    const GeometricField<Type, fvPatchField, volMesh>& vvf,
    const vtkMesh& vMesh, bool vtu
)
{
    const fvMesh& mesh = vMesh.mesh();

    const labelList& superCells = vMesh.topo().superCells();

    label nValues = mesh.nCells() + superCells.size();

    if(vtu)
    {
      os<<"   <DataArray Name=\""<<vvf.name()<<"\" type=\"Float32\" format=\"ascii\" NumberOfComponents=\""<< pTraits<Type>::nComponents <<"\">"<<endl;
    }
    else
    os  << vvf.name() << ' ' << pTraits<Type>::nComponents << ' '
        << nValues << " float" << std::endl;

    DynamicList<floatScalar> fField(pTraits<Type>::nComponents*nValues);

    insert(vvf.internalField(), fField);

    forAll(superCells, superCellI)
    {
        label origCellI = superCells[superCellI];

        insert(vvf[origCellI], fField);
    }
    write(os, binary, fField);

    if(vtu)
    {
      os<<"   </DataArray>"<<endl;
    }
}


template<class Type>
void Foam::writeFuns::write
(
    std::ostream& os,
    const bool binary,
    const GeometricField<Type, pointPatchField, pointMesh>& pvf,
    const vtkMesh& vMesh, bool vtu
)
{
    const fvMesh& mesh = vMesh.mesh();
    const vtkTopo& topo = vMesh.topo();

    const labelList& addPointCellLabels = topo.addPointCellLabels();
    const label nTotPoints = mesh.nPoints() + addPointCellLabels.size();
    if(vtu)
    {
      os<<"   <DataArray Name=\""<<pvf.name()<<"\" type=\"Float32\" format=\"ascii\" NumberOfComponents=\""<< pTraits<Type>::nComponents <<"\">"<<endl;
    }
    else
    os  << pvf.name() << ' ' << pTraits<Type>::nComponents << ' '
        << nTotPoints << " float" << std::endl;

    DynamicList<floatScalar> fField(pTraits<Type>::nComponents*nTotPoints);

    insert(pvf, fField);

    forAll(addPointCellLabels, api)
    {
        label origCellI = addPointCellLabels[api];

        insert(interpolatePointToCell(pvf, origCellI), fField);
    }
    write(os, binary, fField);
    if(vtu)
    {
      os<<"   </DataArray>"<<endl;
    }
}


template<class Type>
void Foam::writeFuns::write
(
    std::ostream& os,
    const bool binary,
    const GeometricField<Type, fvPatchField, volMesh>& vvf,
    const GeometricField<Type, pointPatchField, pointMesh>& pvf,
    const vtkMesh& vMesh, bool vtu
)
{
    const fvMesh& mesh = vMesh.mesh();
    const vtkTopo& topo = vMesh.topo();

    const labelList& addPointCellLabels = topo.addPointCellLabels();
    const label nTotPoints = mesh.nPoints() + addPointCellLabels.size();

    os  << vvf.name() << ' ' << pTraits<Type>::nComponents << ' '
        << nTotPoints << " float" << std::endl;

    DynamicList<floatScalar> fField(pTraits<Type>::nComponents*nTotPoints);

    insert(pvf, fField);

    forAll(addPointCellLabels, api)
    {
        label origCellI = addPointCellLabels[api];

        insert(vvf[origCellI], fField);
    }
    write(os, binary, fField);
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::writeFuns::write
(
    std::ostream& os,
    const bool binary,
    const PtrList<GeometricField<Type, PatchField, GeoMesh> >& flds,
    const vtkMesh& vMesh, bool vtu
)
{
    forAll(flds, i)
    {
        write(os, binary, flds[i], vMesh,vtu);
    }
}


template<class Type>
void Foam::writeFuns::write
(
    std::ostream& os,
    const bool binary,
    const volPointInterpolation& pInterp,
    const PtrList<GeometricField<Type, fvPatchField, volMesh> >& flds,
    const vtkMesh& vMesh, bool vtu
)
{
    forAll(flds, i)
    {
        write(os, binary, flds[i], pInterp.interpolate(flds[i])(), vMesh,vtu);
    }
}


// ************************************************************************* //
