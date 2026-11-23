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
    Foam::phaseRelaxation

Description
    Infinite (instantaneous) pressure relaxation of two stiffened-gas phases
    in the phasic-internal-energy formulation (see fluxScheme.H).  The two
    phases are driven to the mechanical equilibrium p1 - p2 = pc, with the
    capillary pressure pc = sigma*kappa limited so that both decoded phase
    pressures stay above phaseFloors::p, in three steps:

        Step A  volume-fraction relaxation from the per-phase internal
                energies with the pre-relaxation interface pressure,
                solved in closed form via the quadratic of Bharate &
                Mandal (2026), the explicit interface-pressure variant
                of the stiffened-gas solution of Pelanti & Shyue (2014)
        Step B  Kapila mixture-EoS decode (Kapila et al. 2001) of the phase
                pressures p2 and p1 = p2 + pc from the conserved mixture
                energy with the relaxed volume fractions
        Step C  per-phase internal energy reinitialization from p1, p2,
                plus the mixture T, the heat-capacity-weighted mean of
                the phase temperatures

    update() fixes (alpha, alphaE, p, T) with the mixture pressure
    p = alpha1*p1 + alpha2*p2; the per-phase densities and sound speeds are
    not updated here - the mixture refreshes them from the phase pressures
    afterward.

    References:
    \verbatim
        Kapila, A. K., Menikoff, R., Bdzil, J. B., Son, S. F. &
        Stewart, D. S. (2001).
        Two-phase modeling of deflagration-to-detonation transition in
        granular materials: Reduced equations.
        Physics of Fluids, 13(10), 3002-3024.

        Saurel, R., Petitpas, F. & Berry, R. A. (2009).
        Simple and efficient relaxation methods for interfaces
        separating compressible fluids, cavitating flows and shocks
        in multiphase mixtures.
        Journal of Computational Physics, 228(5), 1678-1712.

        Pelanti, M. & Shyue, K.-M. (2014).
        A mixture-energy-consistent six-equation two-phase numerical
        model for fluids with interfaces, cavitation and evaporation
        waves.
        Journal of Computational Physics, 259, 331-357.

        Bharate, G. & Mandal, J. C. (2026).
        A robust all-Mach six-equation diffuse-interface method for
        multiphase flows with surface tension.
        arXiv:2609.10727.
    \endverbatim

SourceFiles
    phaseRelaxation.C
    phaseRelaxationI.H

\*---------------------------------------------------------------------------*/

#ifndef phaseRelaxation_H
#define phaseRelaxation_H

#include "volFields.H"
#include "phaseFloors.H"
#include "stiffenedGasEoS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration
class phaseModel;

/*---------------------------------------------------------------------------*\
                       Class phaseRelaxation Declaration
\*---------------------------------------------------------------------------*/

class phaseRelaxation
{
    // Private Data

        //- Reference to phase 1's stiffened-gas EoS
        const stiffenedGasEoS& sgEoS1_;

        //- Reference to phase 2's stiffened-gas EoS
        const stiffenedGasEoS& sgEoS2_;


    // Private Member Functions

        //- Pointwise relaxation kernel: Steps A, B, C and the mixture
        //  temperature for a single cell
        void relaxCell
        (
            // EoS parameters (consts)
            const scalar g1,
            const scalar g2,
            const scalar pInf1,
            const scalar pInf2,
            const scalar Cv1,
            const scalar Cv2,
            // Inputs (cell-local)
            const scalar rho,
            const scalar eMixture,
            const scalar alphaRho1,
            const scalar alphaRho2,
            // In/out
            scalar& pc,
            scalar& alpha1,
            scalar& alpha2,
            scalar& alphaRhoE1,
            scalar& alphaRhoE2,
            scalar& p,
            scalar& T
        ) const;


public:

    // Constructors

        //- Construct from the stiffened-gas EoS of the two phases
        phaseRelaxation
        (
            const stiffenedGasEoS& sgEoS1,
            const stiffenedGasEoS& sgEoS2
        );


    //- Destructor
    ~phaseRelaxation() = default;


    // Member Functions

        //- Execute the relaxation (Steps A-C) and update T; see the class
        //  description. eMixture is the conservative mixture specific
        //  internal energy rhoE/rho - 1/2|U|^2, pc the capillary pressure
        //  sigma*kappa, returned limited, and p the mixture pressure
        //  alpha1*p1 + alpha2*p2
        void update
        (
            phaseModel& phase1,
            phaseModel& phase2,
            const volScalarField& rho,
            const volScalarField& eMixture,
            volScalarField& pc,
            volScalarField& p,
            volScalarField& T
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const phaseRelaxation&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "phaseRelaxationI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
