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
    Foam::phaseModel

Description
    Per-phase container (phasic-internal-energy formulation): owns the
    volume fraction, the conserved fields alphaRho and alphaE
    (= alpha*rho*e, no kinetic energy), the derived rho and c, and the
    phase EoS and transport models.

    correct(p, T) performs the conservative decode at every RK stage:
    rho_k = alphaRho_k/alpha_k where the phase is present, the EoS density
    rho_k(p, T) where it has vanished, then c_k from the EoS at
    (p, rho_k). Temperature-dependent transport properties are
    evaluated at the shared mixture temperature (single-temperature
    closure), so the phase temperatures themselves are not stored.

    The phase density rho.<phase> is written, and on a restart it is read
    back and initialize() re-seeds alphaRho and alphaE from it: the phases
    are not in thermal equilibrium, so in a mixed cell the phase density
    cannot be recovered from p and the mixture temperature.  On a cold start
    (no rho.<phase> file) initialize() sets the phase density from the EoS
    at p and T, in thermal equilibrium.

SourceFiles
    phaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef phaseModel_H
#define phaseModel_H

#include "phaseEquationOfState.H"
#include "phaseTransport.H"
#include "IOdictionary.H"
#include "viscosity.H"
#include "volFields.H"
#include "autoPtr.H"
#include "phaseFloors.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class phaseModel Declaration
\*---------------------------------------------------------------------------*/

class phaseModel
:
    public IOdictionary,
    public viscosity
{
    // Private Data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Phase name
        const word name_;

        //- Equation of state for this phase
        autoPtr<phaseEquationOfState> eos_;

        //- Transport model for this phase
        autoPtr<phaseTransport> transport_;


        // Primitive fields

            //- Volume fraction
            volScalarField& alpha_;


        // Conserved fields

            //- alpha*rho [kg/m^3]
            volScalarField alphaRho_;

            //- alpha*rho*e [J/m^3]
            volScalarField alphaE_;


        // Derived fields

            //- Phase density [kg/m^3]
            volScalarField rho_;

            //- Speed of sound in this phase [m/s]
            volScalarField c_;


        //- Was the phase density read (restart)?  initialize() then keeps it
        const bool rhoRead_;


    // Private Member Functions

        //- Disallow default bitwise copy construction
        phaseModel(const phaseModel&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const phaseModel&) = delete;


public:

    // Constructors

        //- Construct from components
        phaseModel
        (
            const fvMesh& mesh,
            const word& phaseName,
            volScalarField& alpha,
            const volScalarField& p,
            const volScalarField& T
        );


    //- Destructor
    ~phaseModel();


    // Member Functions

        // Access

            //- Phase name
            const word& name() const
            {
                return name_;
            }

            //- Const access to the equation of state
            const phaseEquationOfState& eos() const
            {
                return eos_();
            }

            //- Const access to the transport model
            const phaseTransport& transport() const
            {
                return transport_();
            }


        // Transport property access

            //- Dynamic viscosity [Pa s]
            const volScalarField& mu() const
            {
                return transport_().mu();
            }

            //- Thermal conductivity [W/(m K)]
            const volScalarField& kappa() const
            {
                return transport_().kappa();
            }


        // Primitive field access

            //- Volume fraction
            volScalarField& alpha()
            {
                return alpha_;
            }
            const volScalarField& alpha() const
            {
                return alpha_;
            }


        // Conserved field access

            //- Partial density alpha*rho [kg/m^3]
            volScalarField& alphaRho()
            {
                return alphaRho_;
            }
            const volScalarField& alphaRho() const
            {
                return alphaRho_;
            }

            //- Phasic internal energy density [J/m^3]
            //  alpha_k*rho_k*e_k (no kinetic energy)
            volScalarField& alphaE()
            {
                return alphaE_;
            }
            const volScalarField& alphaE() const
            {
                return alphaE_;
            }


        // Derived field access

            //- Speed of sound [m/s]
            const volScalarField& c() const
            {
                return c_;
            }

            //- Phase density [kg/m^3]
            volScalarField& rho()
            {
                return rho_;
            }
            const volScalarField& rho() const
            {
                return rho_;
            }


        // Viscosity interface (inherited from Foam::viscosity)

            //- Kinematic viscosity of the phase [m^2/s]
            virtual tmp<volScalarField> nu() const;

            //- Kinematic viscosity of the phase for a patch [m^2/s]
            virtual tmp<scalarField> nu(const label patchi) const;


        // Update

            //- Seed the phase density, sound speed and conserved fields at
            //  the phase pressure p, the phase density from the EoS at p and
            //  T unless it was read
            void initialize(const volScalarField& p, const volScalarField& T);

            //- Conservative decode of derived per-phase fields from the
            //  post-relaxation phase pressure and the conservatively
            //  transported partial density alphaRho.  T supplies the EoS
            //  density of a phase that has vanished, where alphaRho/alpha
            //  is 0/0
            void correct(const volScalarField& p, const volScalarField& T);

            //- Refresh the temperature-dependent transport properties
            //  (mu, kappa) at the supplied mixture temperature
            void correctTransport(const volScalarField& T);

            //- Correct the boundary conditions of the conservative
            //  fields based on the primitive variables
            void correctBoundariesFromPrimitives
            (
                const volScalarField& p,
                const volScalarField& T
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
