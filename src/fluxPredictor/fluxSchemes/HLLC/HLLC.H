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
    Foam::HLLC

Description
    HLLC approximate Riemann solver for the 6-equation two-phase model.

    The partial densities and phasic internal energies are advected with
    the contact wave; the mixture total energy uses the standard HLLC star
    state (Bharate & Mandal 2025).  The volume fraction is advected
    non-conservatively with the HLLC-consistent face velocity (Johnsen &
    Colonius 2006).  Surface tension enters as a capillary pressure jump in
    the contact-wave speed S* (Garrick et al. 2017a).

    A low-Mach correction (Chen et al. 2020) scales the acoustic part of
    the star-region pressures,

        p_d = phiL*phiR*(qR - qL)/(phiR - phiL),   phiK = rhoK*(SK - qK),

    with the face-normal velocities qK, by the Mach-number function

        f = M*sqrt(4 + (1 - M^2)^2)/(1 + M^2),
        M = min(1, max(max(|U_pos|, Ucut)/c_pos, max(|U_neg|, Ucut)/c_neg)).

    The pressure fix pf = (f - 1)*p_d enters the star-region momentum and
    total-energy fluxes as pf*(0, n, S*); the face states, the wave speeds
    and S* are those of plain HLLC, which is recovered for M -> 1.  The
    shock sensor and the shear viscosity of Chen et al. (2020) are not
    included.

    The face velocities in M are floored by the capillary velocity

        Ucut = sqrt(|sigma*kappa_f|/(rho1 + rho2))

    with the volume-fraction-weighted phase densities rho1 and rho2 of the
    two face states.

    \verbatim
    fluxScheme
    {
        type                HLLC;

        lowMachCorrection   true;
    }
    \endverbatim

    References:
    \verbatim
        Toro, E. F., Spruce, M. & Speares, W. (1994).
        Restoration of the contact surface in the HLL-Riemann solver.
        Shock Waves, 4(1), 25-34.

        Davis, S. F. (1988).
        Simplified second-order Godunov-type methods.
        SIAM Journal on Scientific and Statistical Computing, 9(3), 445-473.

        Johnsen, E. & Colonius, T. (2006).
        Implementation of WENO schemes in compressible multicomponent
        flow problems.
        Journal of Computational Physics, 219(2), 715-732.

        Chen, S.-S., Lin, B.-X., Li, Y.-S. & Yan, C. (2020).
        HLLC+: Low-Mach shock-stable HLLC-type Riemann solver for
        all-speed flows.
        SIAM Journal on Scientific Computing, 42(4), B921-B950.

        Saurel, R., Petitpas, F. & Berry, R. A. (2009).
        Simple and efficient relaxation methods for interfaces
        separating compressible fluids, cavitating flows and shocks
        in multiphase mixtures.
        Journal of Computational Physics, 228(5), 1678-1712.

        Garrick, D. P., Owkes, M. & Regele, J. D. (2017a).
        A finite-volume HLLC-based scheme for compressible interfacial
        flows with surface tension.
        Journal of Computational Physics, 339, 46-67.

        Bharate, G. & Mandal, J. C. (2025).
        Enhanced diffuse interface method for multiphase flow
        simulations across all Mach numbers.
        Journal of Computational Physics, 543, 114397.
    \endverbatim

SourceFiles
    HLLC.C

\*---------------------------------------------------------------------------*/

#ifndef HLLC_H
#define HLLC_H

#include "fluxScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class HLLC Declaration
\*---------------------------------------------------------------------------*/

class HLLC
:
    public fluxScheme
{
    // Private Data

        //- Apply the low-Mach correction (Chen et al. 2020)
        const bool lowMachCorrection_;


public:

    //- Runtime type information
    TypeName("HLLC");


    // Constructors

        //- Construct from components
        HLLC
        (
            const fvMesh&,
            const dictionary&
        );


    //- Destructor
    virtual ~HLLC() = default;


    // Member Functions

        //- Evaluate HLLC fluxes for the 6-equation two-phase model
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
            // surface-tension source (Garrick et al. 2017a)
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
