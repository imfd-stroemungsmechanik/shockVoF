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
    Foam::compressibleTwoPhaseMixture

Description
    Two-phase compressible mixture for the density-based solver.

    Phasic-internal-energy formulation with redundant mixture total
    energy: correct() takes the conservative mixture specific internal
    energy eMixture (from the 7th equation: rhoE/rho - 1/2|U|^2) and
    drives the pressure relaxation.

    Both phases use a stiffened-gas equation of state: stiffenedGas, or
    idealGas, its special case with pInf = 0.  primitivesToState() applies
    the same closure to convert primitives into the conserved states of the
    flux schemes.

SourceFiles
    compressibleTwoPhaseMixture.C
    compressibleTwoPhaseMixtureI.H

\*---------------------------------------------------------------------------*/

#ifndef compressibleTwoPhaseMixture_H
#define compressibleTwoPhaseMixture_H

#include "twoPhaseMixture.H"
#include "viscosity.H"
#include "phaseModel.H"
#include "stiffenedGasEoS.H"
#include "phaseRelaxation.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class compressibleTwoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class compressibleTwoPhaseMixture
:
    public twoPhaseMixture,
    public viscosity
{
    // Private Data

        // Shared primitives owned by the mixture

            //- Mixture pressure [Pa]
            volScalarField p_;

            //- Mixture temperature [K]
            volScalarField T_;

            //- Capillary pressure p1 - p2 [Pa]: sigma*kappa, limited so
            //  that both phase pressures stay above phaseFloors::p
            volScalarField pc_;

            //- Phase-1 pressure p1 = p + alpha2*pc [Pa]
            volScalarField p1_;

            //- Phase-2 pressure p2 = p - alpha1*pc [Pa]
            volScalarField p2_;


        // Phase models

            //- Phase 1 model
            phaseModel phase1_;

            //- Phase 2 model
            phaseModel phase2_;


        // Phase equations of state

            //- Stiffened-gas EoS of phase 1
            const stiffenedGasEoS& sgEoS1_;

            //- Stiffened-gas EoS of phase 2
            const stiffenedGasEoS& sgEoS2_;


        // Derived mixture fields

            //- Mixture density [kg/m^3]
            volScalarField rho_;

            //- Mixture dynamic viscosity [Pa s]
            volScalarField mu_;

            //- Mixture kinematic viscosity [m^2/s]
            volScalarField nu_;

            //- Frozen mixture sound speed [m/s]
            volScalarField c_;

            //- Mixture compressibility [s^2/m^2]
            volScalarField psi_;


        // Pressure relaxation

            //- Pressure relaxation operator
            phaseRelaxation relaxation_;


    // Private Member Functions

        //- Update mixing-rule transport and frozen mixture acoustics
        void updateMixtureTransportAndAcoustics();


public:

    // Constructors

        //- Construct from mesh
        explicit compressibleTwoPhaseMixture(const fvMesh& mesh);


    //- Destructor
    virtual ~compressibleTwoPhaseMixture();


    // Member Functions

        // Phase access

            //- Phase 1 model
            phaseModel& phase1()
            {
                return phase1_;
            }
            const phaseModel& phase1() const
            {
                return phase1_;
            }

            //- Phase 2 model
            phaseModel& phase2()
            {
                return phase2_;
            }
            const phaseModel& phase2() const
            {
                return phase2_;
            }


        // Access to shared primitive state

            //- Mixture pressure [Pa]
            volScalarField& p()
            {
                return p_;
            }
            const volScalarField& p() const
            {
                return p_;
            }

            //- Mixture temperature [K]
            volScalarField& T()
            {
                return T_;
            }
            const volScalarField& T() const
            {
                return T_;
            }

            //- Capillary pressure pc = p1 - p2 [Pa]
            const volScalarField& pc() const
            {
                return pc_;
            }

            //- Phase pressures p1 = p + alpha2*pc and p2 = p - alpha1*pc
            const volScalarField& p1() const
            {
                return p1_;
            }
            const volScalarField& p2() const
            {
                return p2_;
            }


        // Access to derived mixture fields

            //- Mixture density [kg/m^3]
            volScalarField& rho()
            {
                return rho_;
            }
            const volScalarField& rho() const
            {
                return rho_;
            }

            //- Mixture dynamic viscosity [Pa s]
            const volScalarField& mu() const
            {
                return mu_;
            }

            //- Frozen mixture sound speed [m/s]
            const volScalarField& c() const
            {
                return c_;
            }

            //- Mixture compressibility [s^2/m^2]
            const volScalarField& psi() const
            {
                return psi_;
            }


        // State conversion

            //- Pointwise conversion from the primitives (p1, p2, alpha1,
            //  rho1, rho2) with the phase pressures p1, p2 to the mixture
            //  density and frozen sound speed, the partial densities and
            //  the per-phase internal energy densities.  Evaluated for
            //  both face states of every face, for the high- and the
            //  low-order fluxes
            inline void primitivesToState
            (
                const scalar p1,
                const scalar p2,
                const scalar alpha1,
                const scalar rho1,
                const scalar rho2,
                scalar& rho,
                scalar& c,
                scalar& alphaRho1,
                scalar& alphaRho2,
                scalar& alphaRhoE1,
                scalar& alphaRhoE2
            ) const;


        // Viscosity interface

            //- Mixture kinematic viscosity [m^2/s]
            virtual tmp<volScalarField> nu() const;

            //- Mixture kinematic viscosity for a patch [m^2/s]
            virtual tmp<scalarField> nu(const label patchi) const;


        // Update

            //- Required by the twoPhaseMixture interface
            virtual void correct()
            {}

            //- Re-read the phaseProperties dictionary
            virtual bool read()
            {
                return twoPhaseMixture::read();
            }

            //- Impose alpha2 = 1 - alpha1 on the non-coupled patches and
            //  update the phase pressures from p, pc and the volume
            //  fractions.  Called by initialize() and correct(), and by the
            //  solver after it advances alpha1 with the flow frozen
            void correctPhasePressures();

            //- Mixture density from the transported partial densities,
            //  floored as the density with which the solver forms the
            //  velocity and the specific energy.  Called by the solver
            //  before correct()
            void correctRho();

            //- Seed the phase states from the phase pressures at the
            //  mixture pressure p and temperature T, with the capillary
            //  pressure sigmaKappa.  Phase densities read on a restart are
            //  kept
            void initialize(const volScalarField& sigmaKappa);

            //- Drive the full conservative-to-primitive decode, including
            //  the boundary state (correctBoundariesFromPrimitives()), with
            //  the capillary pressure sigmaKappa and the mixture density of
            //  correctRho()
            void correct
            (
                const volScalarField& eMixture,
                const volScalarField& sigmaKappa
            );

            //- Impose the non-coupled boundary values of the phase
            //  densities and conserved phase fields and the mixture density
            //  from the primitive variables.  Called by correct() after
            //  correctPhasePressures() and the cell decode
            void correctBoundariesFromPrimitives();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "compressibleTwoPhaseMixtureI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
