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
    Foam::fluxPredictor

Description
    Inviscid flux predictor of the density-based two-phase VoF solver.

    update() evaluates the maximum signal speed and two flux sets, both
    with evaluateFluxes():
      - high-order: face states from the selected reconstructionScheme
        through the selected flux scheme (typically HLLC),
      - low-order: cell-center states, from a fixed oneDim upwind
        reconstruction without sharpening or positivity limiting,
        through HLL.
    blend() combines them per face with the flux-corrected-transport
    factor theta, flux = fluxLO + theta*(flux - fluxLO).

    On non-coupled patches both sets use the boundary state on both sides.

SourceFiles
    fluxPredictor.C

See also
    Foam::fluxScheme
    Foam::rhoThinc
    Foam::compressibleTwoPhaseMixture
    Foam::interfaceModel

\*---------------------------------------------------------------------------*/

#ifndef fluxPredictor_H
#define fluxPredictor_H

#include "volFields.H"
#include "surfaceFields.H"
#include "compressibleTwoPhaseMixture.H"
#include "interfaceModel.H"

#include "reconstructionScheme.H"
#include "fluxScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class fluxPredictor Declaration
\*---------------------------------------------------------------------------*/

class fluxPredictor
{
    // Private Data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Reference to the mixture
        compressibleTwoPhaseMixture& mixture_;

        //- Reference to the interface properties
        const interfaceModel& interface_;


        // Reference to primitive fields

            //- Phase fraction
            const volScalarField& alpha1_;

            //- Velocity
            const volVectorField& U_;

            //- Phase 1 density
            const volScalarField& rho1_;

            //- Phase 2 density
            const volScalarField& rho2_;


        //- Runtime-selected primitive face reconstruction
        autoPtr<reconstructionScheme> reconstruction_;

        //- First-order (upwind) face reconstruction of the low-order fluxes
        autoPtr<reconstructionScheme> reconstructionLowOrder_;


        //- High-order flux scheme
        autoPtr<fluxScheme> fluxScheme_;

        //- Low-order HLL flux scheme
        autoPtr<fluxScheme> fluxSchemeLowOrder_;

        //- Maximum wave speed * face area, of the high- and the low-order
        //  Riemann problems
        surfaceScalarField amaxSf_;


        // Face fluxes

            //- Mass flux
            surfaceScalarField phi_;

            //- Momentum flux
            surfaceVectorField phiUp_;

            //- Partial density flux for phase 1
            surfaceScalarField phiAlpha1Rho1_;

            //- Partial density flux for phase 2
            surfaceScalarField phiAlpha2Rho2_;

            //- Volume fraction advective flux
            surfaceScalarField phiAlpha1_;

            //- Per-phase internal energy flux for phase 1
            surfaceScalarField phiAlpha1RhoE1_;

            //- Per-phase internal energy flux for phase 2
            surfaceScalarField phiAlpha2RhoE2_;

            //- Mixture total energy flux (conservative)
            surfaceScalarField phiRhoE_;

            //- Volumetric flux
            surfaceScalarField phiv_;

            //- Face volume fraction for the surface-tension source
            surfaceScalarField alphaf_;


        // Low-order face fluxes (HLL + first-order reconstruction)

            surfaceScalarField phiLO_;
            surfaceVectorField phiUpLO_;
            surfaceScalarField phiAlpha1Rho1LO_;
            surfaceScalarField phiAlpha2Rho2LO_;
            surfaceScalarField phiAlpha1LO_;
            surfaceScalarField phiAlpha1RhoE1LO_;
            surfaceScalarField phiAlpha2RhoE2LO_;
            surfaceScalarField phiRhoELO_;
            surfaceScalarField phivLO_;
            surfaceScalarField alphafLO_;


        // Flux-corrected-transport diagnostics

            //- Cell-centered flux limiter for diagnostics
            autoPtr<volScalarField> fluxLimiterPtr_;

            //- Current time index
            label curTimeIndex_;


    // Private Member Functions

        //- Disallow default bitwise copy construction
        fluxPredictor(const fluxPredictor&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const fluxPredictor&) = delete;

        //- Evaluate the fluxes of the owner/neighbor face states with the
        //  given flux scheme, and raise amaxSf_ to its signal speed
        void evaluateFluxes
        (
            const fluxScheme& scheme,
            const surfaceScalarField& alpha1_pos,
            const surfaceScalarField& alpha1_neg,
            const surfaceScalarField& p1_pos,
            const surfaceScalarField& p1_neg,
            const surfaceScalarField& p2_pos,
            const surfaceScalarField& p2_neg,
            const surfaceVectorField& U_pos,
            const surfaceVectorField& U_neg,
            const surfaceScalarField& rho1_pos,
            const surfaceScalarField& rho1_neg,
            const surfaceScalarField& rho2_pos,
            const surfaceScalarField& rho2_neg,
            surfaceScalarField& phi,
            surfaceVectorField& phiUp,
            surfaceScalarField& phiAlpha1Rho1,
            surfaceScalarField& phiAlpha2Rho2,
            surfaceScalarField& phiAlpha1,
            surfaceScalarField& phiAlpha1RhoE1,
            surfaceScalarField& phiAlpha2RhoE2,
            surfaceScalarField& phiRhoE,
            surfaceScalarField& phiv,
            surfaceScalarField& alphaf
        );


public:

    //- Runtime type information (enables the "fluxPredictor" debug switch)
    ClassName("fluxPredictor");


    // Constructors

        //- Construct from components
        fluxPredictor
        (
            const volScalarField& p,
            const volVectorField& U,
            compressibleTwoPhaseMixture& mixture,
            const interfaceModel& interface
        );


    //- Destructor
    ~fluxPredictor() = default;


    // Member Functions

        // Return fluxes

            //- Return mass flux
            const surfaceScalarField& phi() const
            {
                return phi_;
            }

            //- Return momentum flux
            const surfaceVectorField& phiUp() const
            {
                return phiUp_;
            }

            //- Return density flux for phase 1
            const surfaceScalarField& phiAlpha1Rho1() const
            {
                return phiAlpha1Rho1_;
            }

            //- Return density flux for phase 2
            const surfaceScalarField& phiAlpha2Rho2() const
            {
                return phiAlpha2Rho2_;
            }

            //- Return volume fraction advective flux
            const surfaceScalarField& phiAlpha1() const
            {
                return phiAlpha1_;
            }

            //- Return per-phase internal energy flux for phase 1
            const surfaceScalarField& phiAlpha1RhoE1() const
            {
                return phiAlpha1RhoE1_;
            }

            //- Return per-phase internal energy flux for phase 2
            const surfaceScalarField& phiAlpha2RhoE2() const
            {
                return phiAlpha2RhoE2_;
            }

            //- Return mixture total energy flux
            const surfaceScalarField& phiRhoE() const
            {
                return phiRhoE_;
            }

            //- Return volumetric flux
            const surfaceScalarField& phiv() const
            {
                return phiv_;
            }

            //- Return the face volume fraction for the surface-tension source
            const surfaceScalarField& alphaf() const
            {
                return alphaf_;
            }


        // Return low-order fluxes

            //- Return low-order density flux for phase 1
            const surfaceScalarField& phiAlpha1Rho1LO() const
            {
                return phiAlpha1Rho1LO_;
            }

            //- Return low-order density flux for phase 2
            const surfaceScalarField& phiAlpha2Rho2LO() const
            {
                return phiAlpha2Rho2LO_;
            }

            //- Return low-order per-phase internal energy flux for phase 1
            const surfaceScalarField& phiAlpha1RhoE1LO() const
            {
                return phiAlpha1RhoE1LO_;
            }

            //- Return low-order per-phase internal energy flux for phase 2
            const surfaceScalarField& phiAlpha2RhoE2LO() const
            {
                return phiAlpha2RhoE2LO_;
            }

            //- Return low-order mixture total energy flux
            const surfaceScalarField& phiRhoELO() const
            {
                return phiRhoELO_;
            }

            //- Return low-order volumetric flux
            const surfaceScalarField& phivLO() const
            {
                return phivLO_;
            }


        //- Return the maximum signal speed * face area of both flux sets
        const surfaceScalarField& amaxSf() const
        {
            return amaxSf_;
        }

        //- Reconstruct the primitives and evaluate both the high-order
        //  (typically HLLC) and low-order (HLL, first-order) fluxes
        void update();

        //- Blend the high-order fluxes toward the low-order ones
        void blend(const surfaceScalarField& theta);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
