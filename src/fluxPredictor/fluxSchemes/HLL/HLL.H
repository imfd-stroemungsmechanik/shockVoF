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
    Foam::HLL

Description
    HLL approximate Riemann solver for the 6-equation two-phase model.

    Evaluated with first-order states, it serves as the low-order flux
    for the flux-corrected-transport limiter. It includes a
    well-balanced treatment of surface tension (Garrick et al. 2017a)
    using a capillary-reduced pressure to lower spurious currents at
    static interfaces.

    References:
    \verbatim
        Harten, A., Lax, P. D. & van Leer, B. (1983).
        On upstream differencing and Godunov-type schemes for
        hyperbolic conservation laws.
        SIAM Review, 25(1), 35-61.

        Davis, S. F. (1988).
        Simplified second-order Godunov-type methods.
        SIAM Journal on Scientific and Statistical Computing, 9(3), 445-473.

        Garrick, D. P., Owkes, M. & Regele, J. D. (2017a).
        A finite-volume HLLC-based scheme for compressible interfacial
        flows with surface tension.
        Journal of Computational Physics, 339, 46-67.
    \endverbatim

SourceFiles
    HLL.C

\*---------------------------------------------------------------------------*/

#ifndef HLL_H
#define HLL_H

#include "fluxScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class HLL Declaration
\*---------------------------------------------------------------------------*/

class HLL
:
    public fluxScheme
{

public:

    //- Runtime type information
    TypeName("HLL");


    // Constructors

        //- Construct from components
        HLL
        (
            const fvMesh&,
            const dictionary&
        );


    //- Destructor
    virtual ~HLL() = default;


    // Member Functions

        //- Evaluate HLL fluxes for the 6-equation two-phase model
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

            // Output: HLL contact-upwinded face volume fraction for the
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

            // Face capillary coefficient sigma*kappa_f [Pa]
            const scalar& sigmaKappaf
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
