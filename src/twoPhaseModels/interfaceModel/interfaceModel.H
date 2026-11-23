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
    Foam::interfaceModel

Description
    Interface properties of the two-phase flow: the interface normal, the
    curvature, the face capillary coefficient sigma*kappa_f and the
    near-interface indicator.

    The CSF source and the Riemann fluxes both read the face coefficient
    sigma*kappa_f, so the capillary force stays balanced whichever curvature
    model is selected; the pressure relaxation uses the cell value
    sigma*kappa.  kappa_f is interpolated linearly from the cell curvature,
    independent of fvSchemes.

    Read from phaseProperties:
    \verbatim
    surfaceTension
    {
        sigma           0.072;
        curvatureModel  smoothed;
    }

    interfaceNormal
    {
        psiExp          0.1;
        nSmooth         0;
    }
    \endverbatim

SourceFiles
    interfaceModel.C

See also
    Foam::interfaceNormal
    Foam::curvatureModel

\*---------------------------------------------------------------------------*/

#ifndef interfaceModel_H
#define interfaceModel_H

#include "interfaceNormal.H"
#include "curvatureModel.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class interfaceModel Declaration
\*---------------------------------------------------------------------------*/

class interfaceModel
{
    // Private Data

        //- Reference to the phase-1 volume fraction
        const volScalarField& alpha1_;

        //- Surface-tension coefficient [N/m]
        dimensionedScalar sigma_;

        //- Interface normal
        interfaceNormal normal_;

        //- Interface-curvature model
        autoPtr<curvatureModel> curvature_;

        //- Face capillary coefficient sigma*kappa_f [Pa]
        surfaceScalarField sigmaKappaf_;


public:

    // Constructors

        //- Construct from the phase-1 volume fraction and the
        //  phaseProperties dictionary.  The curvature and sigmaKappaf are
        //  zero until the first correct()
        interfaceModel
        (
            const volScalarField& alpha1,
            const dictionary& phaseProperties
        );

        //- Disallow default bitwise copy construction
        interfaceModel
        (
            const interfaceModel&
        ) = delete;


    //- Destructor
    ~interfaceModel();


    // Member Functions

        //- Surface-tension coefficient [N/m]
        const dimensionedScalar& sigma() const
        {
            return sigma_;
        }

        //- Interface normal
        const interfaceNormal& normal() const
        {
            return normal_;
        }

        //- Interface curvature [1/m]
        const volScalarField& kappa() const
        {
            return curvature_->kappa();
        }

        //- Cell capillary pressure sigma*kappa [Pa]
        tmp<volScalarField> sigmaKappa() const
        {
            return sigma_*curvature_->kappa();
        }

        //- Face capillary coefficient sigma*kappa_f [Pa]
        const surfaceScalarField& sigmaKappaf() const
        {
            return sigmaKappaf_;
        }

        //- Indicator of the cells with 0.01 <= alpha1 <= 0.99 (1 there,
        //  0 elsewhere), as in OpenFOAM's interfaceProperties
        tmp<volScalarField> nearInterface() const;

        //- Update the normal, the curvature and sigmaKappaf after alpha1
        //  changed
        void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const interfaceModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
