/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2024 OpenFOAM Foundation
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
    Foam::twoPhaseThermophysicalTransportModel

Description
    Thermophysical transport model for the density-based shockVoF
    solver.

    Provides the effective thermal conductivity of the mixture as a
    volume-fraction-weighted mixing rule with a turbulent contribution:

        kappaEff = sum_k  alpha_k * (kappa_k + rho_k * Cp_k * nut_k)

    i.e., with a turbulent Prandtl number of unity for each phase.  The
    turbulent kinematic viscosity nut_k is that of the phase-k model in the
    per-phase (twoPhaseTransport) mode and the mixture nut in the mixture
    mode.

SourceFiles
    twoPhaseThermophysicalTransportModel.C

See also
    Foam::twoPhaseMomentumTransportModel
    Foam::compressibleTwoPhaseMixture

\*---------------------------------------------------------------------------*/

#ifndef twoPhaseThermophysicalTransportModel_H
#define twoPhaseThermophysicalTransportModel_H

#include "compressibleTwoPhaseMixture.H"
#include "twoPhaseMomentumTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class twoPhaseThermophysicalTransportModel Declaration
\*---------------------------------------------------------------------------*/

class twoPhaseThermophysicalTransportModel
{
    // Private Data

        //- Reference to the two-phase mixture
        const compressibleTwoPhaseMixture& mixture_;

        //- Reference to the momentum transport model
        const twoPhaseMomentumTransportModel& momentumTransport_;


public:

    // Constructors

        //- Construct from components
        twoPhaseThermophysicalTransportModel
        (
            const compressibleTwoPhaseMixture& mixture,
            const twoPhaseMomentumTransportModel& momentumTransport
        );

        //- Disallow default bitwise copy construction
        twoPhaseThermophysicalTransportModel
        (
            const twoPhaseThermophysicalTransportModel&
        ) = delete;


    //- Destructor
    ~twoPhaseThermophysicalTransportModel()
    {}


    // Member Functions

        // Transport properties

            //- Effective thermal conductivity of the mixture [W/(m K)]
            tmp<volScalarField> kappaEff() const;


        // Predict / correct

            //- kappaEff is algebraic and requires no prediction step
            void predict();

            //- kappaEff is algebraic and requires no correction step
            void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=
        (
            const twoPhaseThermophysicalTransportModel&
        ) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
