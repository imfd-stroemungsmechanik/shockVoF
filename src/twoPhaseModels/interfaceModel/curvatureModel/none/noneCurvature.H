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
    Foam::curvatureModels::none

Description
    No surface tension: the curvature is left at zero, so the CSF force and the
    capillary flux terms all vanish.

    \verbatim
    surfaceTension
    {
        sigma           0;
        curvatureModel  none;
    }
    \endverbatim

    interfaceModel reads \c sigma for every model, so the
    entry is required here too.

SourceFiles
    noneCurvature.C

\*---------------------------------------------------------------------------*/

#ifndef noneCurvature_H
#define noneCurvature_H

#include "curvatureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace curvatureModels
{

/*---------------------------------------------------------------------------*\
                            Class none Declaration
\*---------------------------------------------------------------------------*/

class none
:
    public curvatureModel
{
public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from mesh, phase-1 fraction, interface normal and
        //  the surfaceTension sub-dictionary
        none
        (
            const fvMesh& mesh,
            const volScalarField& alpha1,
            const interfaceNormal& normal,
            const dictionary& dict
        );


    //- Destructor
    virtual ~none();


    // Member Functions

        //- No-op: curvature stays zero
        virtual void correct()
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace curvatureModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
