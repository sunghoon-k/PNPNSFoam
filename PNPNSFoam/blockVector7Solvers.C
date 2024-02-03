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

\*---------------------------------------------------------------------------*/

#include "vector7Field.H"
#include "tensor7Field.H"
#include "ExpandTensorN.H"
#include "ExpandTensorNField.H"
#include "blockLduMatrices.H"
#include "addToRunTimeSelectionTable.H"

#include "blockLduPrecons.H"
#include "BlockNoPrecon.H"
#include "blockDiagonalPrecons.H"
#include "blockGaussSeidelPrecons.H"
#include "BlockCholeskyPrecon.H"

#include "blockLduSmoothers.H"
#include "blockGaussSeidelSmoothers.H"
#include "BlockILUSmoother.H"

#include "blockLduSolvers.H"
#include "BlockDiagonalSolver.H"
#include "BlockBiCGStabSolver.H"
#include "BlockCGSolver.H"
#include "BlockGaussSeidelSolver.H"
#include "BlockGMRESSolver.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Preconditioners
typedef BlockLduPrecon<vector7> blockVector7Precon;
defineNamedTemplateTypeNameAndDebug(blockVector7Precon, 0);
defineTemplateRunTimeSelectionTable(blockVector7Precon, dictionary);

typedef BlockNoPrecon<vector7> blockNoPreconVector7;
makeBlockPrecon(blockVector7Precon, blockNoPreconVector7);

typedef BlockDiagonalPrecon<vector7> blockDiagonalPreconVector7;
makeBlockPrecon(blockVector7Precon, blockDiagonalPreconVector7);

typedef BlockGaussSeidelPrecon<vector7> blockGaussSeidelPreconVector7;
makeBlockPrecon(blockVector7Precon, blockGaussSeidelPreconVector7);

typedef BlockCholeskyPrecon<vector7> blockCholeskyPreconVector7;
makeBlockPrecon(blockVector7Precon, blockCholeskyPreconVector7);


// Smoothers
typedef BlockLduSmoother<vector7> blockVector7Smoother;
defineNamedTemplateTypeNameAndDebug(blockVector7Smoother, 0);
defineTemplateRunTimeSelectionTable(blockVector7Smoother, dictionary);

typedef BlockGaussSeidelSmoother<vector7> blockGaussSeidelSmootherVector7;
makeBlockSmoother(blockVector7Smoother, blockGaussSeidelSmootherVector7);

typedef BlockILUSmoother<vector7> blockILUSmootherVector7;
makeBlockSmoother(blockVector7Smoother, blockILUSmootherVector7);


// Solvers
typedef BlockLduSolver<vector7> blockVector7Solver;
defineNamedTemplateTypeNameAndDebug(blockVector7Solver, 0);
defineTemplateRunTimeSelectionTable
(
    blockVector7Solver,
    symMatrix
);

defineTemplateRunTimeSelectionTable
(
    blockVector7Solver,
    asymMatrix
);

typedef BlockDiagonalSolver<vector7> blockDiagonalSolverVector7;
defineNamedTemplateTypeNameAndDebug(blockDiagonalSolverVector7, 0);

typedef BlockBiCGStabSolver<vector7> blockBiCGStabSolverVector7;
makeBlockSolverTypeName(blockBiCGStabSolverVector7);
addSolverToBlockMatrix(Vector7, blockBiCGStabSolverVector7, symMatrix);
addSolverToBlockMatrix(Vector7, blockBiCGStabSolverVector7, asymMatrix);

typedef BlockCGSolver<vector7> blockCGSolverVector7;
makeBlockSolverTypeName(blockCGSolverVector7);
addSolverToBlockMatrix(Vector7, blockCGSolverVector7, symMatrix);

typedef BlockGaussSeidelSolver<vector7> blockGaussSeidelSolverVector7;
makeBlockSolverTypeName(blockGaussSeidelSolverVector7);
addSolverToBlockMatrix(Vector7, blockGaussSeidelSolverVector7, symMatrix);
addSolverToBlockMatrix(Vector7, blockGaussSeidelSolverVector7, asymMatrix);

typedef BlockGMRESSolver<vector7> blockGMRESSolverVector7;
makeBlockSolverTypeName(blockGMRESSolverVector7);
addSolverToBlockMatrix(Vector7, blockGMRESSolverVector7, symMatrix);
addSolverToBlockMatrix(Vector7, blockGMRESSolverVector7, asymMatrix);



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
