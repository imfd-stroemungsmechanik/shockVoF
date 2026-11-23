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
    Foam::fluxScheme

Description
    Abstract base class for the inviscid face fluxes of the 6-equation
    two-phase model with the additional mixture total-energy equation.

    The phasic internal-energy fluxes are purely advective,

        F_k = alpha_k*rho_k*e_k*u,

    the pressure work alpha_k*p_k*div(u) of the phase pressure p_k being a
    cell source.  The mixture total-energy flux includes the pressure work,

        F_E = (rho*E + p)*u,   rho*E = sum_k alpha_k*rho_k*e_k + 0.5*rho*|u|^2.

SourceFiles
    fluxScheme.C
    fluxSchemeNew.C

\*---------------------------------------------------------------------------*/

#ifndef fluxScheme_H
#define fluxScheme_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class fluxScheme Declaration
\*---------------------------------------------------------------------------*/

class fluxScheme
{

public:

    //- Runtime type information
    TypeName("fluxScheme");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            fluxScheme,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& dict
            ),
            (mesh, dict)
        );


    // Constructors

        //- Construct null
        fluxScheme();


    //- Destructor
    virtual ~fluxScheme();


    // Selectors

        //- Select the type named in the fluxScheme dictionary
        static autoPtr<fluxScheme> New
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Select the given type, e.g., the low-order HLL scheme
        static autoPtr<fluxScheme> New
        (
            const word& type,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Member Functions

        //- Compute all face fluxes from left (pos/owner) and right
        //  (neg/neighbor) face states
        virtual void evaluateFlux
        (
            // Output fluxes
            scalar& phi,
            vector& phiUp,
            scalar& phiAlpha1Rho1,
            scalar& phiAlpha2Rho2,
            scalar& phiAlpha1,
            scalar& phiAlpha1RhoE1,
            scalar& phiAlpha2RhoE2,
            scalar& phiRhoE,
            scalar& phiv,
            scalar& amaxSf,

            // Output: contact-upwinded face volume fraction for the
            // surface-tension source
            scalar& alphaf,

            // Left (pos / owner) state
            const scalar& p_pos,
            const vector& U_pos,
            const scalar& rho_pos,
            const scalar& c_pos,
            const scalar& alpha1_pos,
            const scalar& alpha1Rho1_pos,
            const scalar& alpha2Rho2_pos,
            const scalar& alphaRhoE1_pos,
            const scalar& alphaRhoE2_pos,

            // Right (neg / neighbor) state
            const scalar& p_neg,
            const vector& U_neg,
            const scalar& rho_neg,
            const scalar& c_neg,
            const scalar& alpha1_neg,
            const scalar& alpha1Rho1_neg,
            const scalar& alpha2Rho2_neg,
            const scalar& alphaRhoE1_neg,
            const scalar& alphaRhoE2_neg,

            // Geometry
            const vector& Sf,
            const scalar& magSf,

            // Face capillary coefficient sigma*kappa_f [Pa] for the
            // well-balanced surface-tension treatment (Garrick et al.
            // 2017a); 0 disables it
            const scalar& sigmaKappaf
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
