#ifndef petscSolver_H
#define petscSolver_H

#include <petsc.h>
#include "volMesh.H"
#include "fvScalarMatrix.H"
namespace Foam
{
	template<class Type, template<class> class PatchField, class GeoMesh>
	class GeometricField;

	template<class Type>
	class fvPatchField;

	template<class Type>
	class fvMatrix;

	template<class Type>
	class fvBlockMatrix;

	class petscSolver
	{

		public:

			ClassName("petscSolver");

			//- Constructor
			petscSolver(int argc, char *argv[]);
			
			//- Destructor
			~petscSolver();

			template<class Type>
			void solve
			(
				fvMatrix<Type>&,
				GeometricField<Type, fvPatchField, volMesh>&
			);

			template<class Type>
			void solve
			(
				fvMatrix<Type>&,
				tmp<GeometricField<Type, fvPatchField, volMesh> > 
			);

			template<class Type>
			void solve
			(
				tmp<fvMatrix<Type> >,
				tmp<GeometricField<Type, fvPatchField, volMesh> > 
			);
			
			template<class Type>
			void solve
			(
				fvBlockMatrix<Type>&
			);

	}; // End class


} // End namespace Foam

#include "petscSolverImp.H"

#endif
