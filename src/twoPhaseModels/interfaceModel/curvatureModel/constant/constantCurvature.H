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
    Foam::curvatureModels::constant

Description
    Constant (analytically prescribed) interface curvature.

    Dictionary entry (in the \c surfaceTension sub-dictionary):
    \verbatim
    surfaceTension
    {
        sigma           0.072;

        curvatureModel  constant;
        kappa           2000;
    }
    \endverbatim

SourceFiles
    constantCurvature.C

\*---------------------------------------------------------------------------*/

#ifndef constantCurvature_H
#define constantCurvature_H

#include "curvatureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace curvatureModels
{

/*---------------------------------------------------------------------------*\
                           Class constant Declaration
\*---------------------------------------------------------------------------*/

class constant
:
    public curvatureModel
{
    // Private Data

        //- Prescribed constant curvature [1/m]
        dimensionedScalar kappaConst_;


public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from mesh, phase-1 fraction, interface normal and
        //  the surfaceTension sub-dictionary
        constant
        (
            const fvMesh& mesh,
            const volScalarField& alpha1,
            const interfaceNormal& normal,
            const dictionary& dict
        );


    //- Destructor
    virtual ~constant();


    // Member Functions

        //- Set the curvature field to the prescribed constant
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace curvatureModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
