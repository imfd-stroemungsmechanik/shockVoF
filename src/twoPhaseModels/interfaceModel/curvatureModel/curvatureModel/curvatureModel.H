/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 Martin Heinrich
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of shockVoF, based on OpenFOAM.

    shockVoF is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    shockVoF is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with shockVoF.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::curvatureModel

Description
    Abstract base class for runtime-selectable interface-curvature models
    of the continuum-surface-force (CSF) surface-tension treatment
    (Brackbill et al. 1992).

    The model owns the curvature field \c kappa_ [1/m], which a concrete
    model fills in correct().  A model that needs the interface normal
    takes the shared interfaceNormal.  interfaceModel owns the model and
    forms the face capillary coefficient sigma*kappa_f.

    The model is selected by \c curvatureModel in the \c surfaceTension
    sub-dictionary of \c phaseProperties, which also holds the entries of
    the selected model.

    References:
    \verbatim
        Brackbill, J. U., Kothe, D. B. & Zemach, C. (1992).
        A continuum method for modeling surface tension.
        Journal of Computational Physics, 100(2), 335-354.
    \endverbatim

SourceFiles
    curvatureModel.C
    curvatureModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef curvatureModel_H
#define curvatureModel_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "interfaceNormal.H"
#include "surfaceFields.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class curvatureModel Declaration
\*---------------------------------------------------------------------------*/

class curvatureModel
{
protected:

    // Protected Data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Reference to the phase-1 volume fraction
        const volScalarField& alpha1_;

        //- Interface normal shared with rho-THINC
        const interfaceNormal& normal_;

        //- Interface curvature [1/m]
        volScalarField kappa_;


public:

    //- Runtime type information
    TypeName("curvatureModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            curvatureModel,
            dictionary,
            (
                const fvMesh& mesh,
                const volScalarField& alpha1,
                const interfaceNormal& normal,
                const dictionary& dict
            ),
            (mesh, alpha1, normal, dict)
        );


    // Constructors

        //- Construct from mesh, phase-1 fraction, interface normal and the
        //  surfaceTension sub-dictionary
        curvatureModel
        (
            const fvMesh& mesh,
            const volScalarField& alpha1,
            const interfaceNormal& normal,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        curvatureModel(const curvatureModel&) = delete;


    // Selector

        //- Select from the surfaceTension sub-dictionary
        static autoPtr<curvatureModel> New
        (
            const fvMesh& mesh,
            const volScalarField& alpha1,
            const interfaceNormal& normal,
            const dictionary& dict
        );


    //- Destructor
    virtual ~curvatureModel();


    // Member Functions

        //- Recompute the curvature field from the current volume fraction
        virtual void correct() = 0;

        //- Interface curvature [1/m]
        const volScalarField& kappa() const
        {
            return kappa_;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const curvatureModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
