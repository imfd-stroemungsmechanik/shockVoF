/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2026 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2026 Martin Heinrich
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
    Foam::solvers::shockVoF

Description
    Solver module for the density-based, compressible two-phase VoF method
    using a 6-equation model with a redundant 7th equation for mixture energy.

    Uses explicit SSP Runge-Kutta time integration with a single-pass
    flux-corrected-transport bound-preserving limiter: the s-stage, second-
    order SSP(s,2) scheme, s = PIMPLE/nStages (default 2).  maxCo and
    maxAlphaCo bound the Courant numbers of the full time step.  Each of the
    s forward-Euler stages has the length deltaT/(s - 1), so both limits can
    be raised in proportion to s - 1 at the same per-stage Courant number.
    The module is intended to run with a single outer corrector
    (PIMPLE/nOuterCorrectors 1).

    Setting the PIMPLE \c flow switch to \c off freezes the flow state
    and advects only the interface. maxCo then bounds the convective
    Courant number of the volumetric flux instead of the signal-speed one.

    References:
    \verbatim
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

        Ketcheson, D. I. (2008).
        Highly efficient strong stability-preserving Runge-Kutta methods
        with low-storage implementations.
        SIAM Journal on Scientific Computing, 30(4), 2113-2136.

        Brackbill, J. U., Kothe, D. B. & Zemach, C. (1992).
        A continuum method for modeling surface tension.
        Journal of Computational Physics, 100(2), 335-354.

        Garrick, D. P., Owkes, M. & Regele, J. D. (2017a).
        A finite-volume HLLC-based scheme for compressible interfacial
        flows with surface tension.
        Journal of Computational Physics, 339, 46-67.
    \endverbatim

SourceFiles
    shockVoF.C
    fluxLimiter.C
    momentumPredictor.C
    solveEquations.C
    updatePrimitives.C

See also
    Foam::solvers::basicFluidSolver
    Foam::solvers::shockFluid

\*---------------------------------------------------------------------------*/

#ifndef shockVoF_H
#define shockVoF_H

#include "basicFluidSolver.H"
#include "compressibleTwoPhaseMixture.H"
#include "interfaceModel.H"
#include "fluxPredictor.H"
#include "twoPhaseMomentumTransportModel.H"
#include "twoPhaseThermophysicalTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class shockVoF Declaration
\*---------------------------------------------------------------------------*/

class shockVoF
:
    public basicFluidSolver
{

protected:

    // Kinematic properties

        //- Velocity field
        volVectorField U_;


    // Thermophysical properties

        //- The compressible two-phase mixture (owns p, T, phases,
        //  relaxation)
        compressibleTwoPhaseMixture mixture_;

        //- Reference to the mixture pressure field
        volScalarField& p_;

        //- Reference to the mixture density field
        volScalarField& rho_;


    // Conserved / primitive references into the mixture phases

        //- Phase-1 volume fraction
        volScalarField& alpha1_;

        //- Phase-2 volume fraction
        volScalarField& alpha2_;

        //- Phase-1 partial density (alpha1*rho1)
        volScalarField& alpha1Rho1_;

        //- Phase-2 partial density (alpha2*rho2)
        volScalarField& alpha2Rho2_;

        //- Phase-1 internal energy density (alpha1*rho1*e1)
        volScalarField& alpha1RhoE1_;

        //- Phase-2 internal energy density (alpha2*rho2*e2)
        volScalarField& alpha2RhoE2_;


    // Mixture conserved fields

        //- Momentum density
        volVectorField rhoU_;

        //- Mixture total energy density (7th equation, conservative)
        volScalarField rhoE_;


    // Interface properties

        //- Interface normal, curvature and capillary coefficient
        interfaceModel interface_;


    // Inviscid fluxes and transport models

        //- Inviscid flux predictor (reconstruction + Riemann solver)
        fluxPredictor flux_;

        //- Momentum transport model
        twoPhaseMomentumTransportModel momentumTransport_;

        //- Thermophysical transport model
        twoPhaseThermophysicalTransportModel thermophysicalTransport_;


    // Time-step control

        //- Interface (phase-fraction) Courant number
        scalar alphaCoNum_;


private:

    // Private Member Functions

        //- Number of stages s of the SSP(s,2) Runge-Kutta scheme
        label nRungeKuttaStages() const;

        //- Correct the cached signal-speed (convective with the flow
        //  frozen) and interface Courant numbers, both per time step
        void correctCoNum();

        //- Solve the 6-equation model and its redundant mixture-total-
        //  energy (7th) equation for one explicit (forward-Euler)
        //  Runge-Kutta sub-stage of length deltaT
        void solveEquations(const dimensionedScalar& deltaT);

        //- Flux-corrected-transport limiter for a forward-Euler sub-stage
        //  of length deltaT
        tmp<surfaceScalarField> fluxLimiter(const dimensionedScalar& deltaT);

        //- Conservative-to-primitive decode (incl. pressure relaxation)
        void updatePrimitives();


public:

    //- Runtime type information
    TypeName("shockVoF");


    // Constructors

        //- Construct from region mesh
        shockVoF(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        shockVoF(const shockVoF&) = delete;


    //- Destructor
    virtual ~shockVoF();


    // Member Functions

        //- Return the current maximum time-step for stable solution
        virtual scalar maxDeltaT() const;

        //- Called at the start of the time-step, before the PIMPLE loop
        virtual void preSolve();

        //- Called at the start of the PIMPLE loop to move the mesh
        virtual void moveMesh();

        //- Corrections that follow mesh motion
        virtual void motionCorrector();

        //- Called at the start of the PIMPLE loop.  Drives the coupled
        //  advance (momentumPredictor()) when the flow is frozen (flow off),
        //  so that the phase fraction is always advected
        virtual void prePredictor();

        //- Predict the momentum transport
        virtual void momentumTransportPredictor();

        //- Predict thermophysical transport
        virtual void thermophysicalTransportPredictor();

        //- Advance the conserved fields through the full SSP Runge-Kutta
        //  cycle with single-pass flux-corrected-transport limiting
        virtual void momentumPredictor();

        //- Not used: the energy equations are advanced inside
        //  momentumPredictor()
        virtual void thermophysicalPredictor();

        //- Not used: the pressure is recovered by relaxation inside
        //  momentumPredictor()
        virtual void pressureCorrector();

        //- Correct the momentum transport
        virtual void momentumTransportCorrector();

        //- Correct the thermophysical transport
        virtual void thermophysicalTransportCorrector();

        //- Called after the PIMPLE loop at the end of the time-step
        virtual void postSolve();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const shockVoF&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
