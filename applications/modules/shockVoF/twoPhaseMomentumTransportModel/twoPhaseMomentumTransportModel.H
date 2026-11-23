/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2024 OpenFOAM Foundation
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
    Foam::twoPhaseMomentumTransportModel

Description
    Momentum transport model for the density-based shockVoF
    solver.

    Supports two modes selected via \c simulationType in
    constant/momentumTransport: per-phase Euler-Euler models
    (\c twoPhaseTransport) or a single mixture model (any other type, e.g.,
    laminar, RAS or LES).  Provides muEff, nut and the explicit viscous
    stress tensors tau and tauMC.

SourceFiles
    twoPhaseMomentumTransportModel.C

See also
    Foam::compressibleTwoPhaseMixture
    Foam::phaseModel

\*---------------------------------------------------------------------------*/

#ifndef twoPhaseMomentumTransportModel_H
#define twoPhaseMomentumTransportModel_H

#include "compressibleTwoPhaseMixture.H"
#include "compressibleMomentumTransportModel.H"
#include "phaseCompressibleMomentumTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration
class twoPhaseThermophysicalTransportModel;

/*---------------------------------------------------------------------------*\
               Class twoPhaseMomentumTransportModel Declaration
\*---------------------------------------------------------------------------*/

class twoPhaseMomentumTransportModel
{
    // Private Data

        //- Switch to select two-phase or mixture transport modeling
        Switch twoPhaseTransport_;

        //- Reference to the two-phase mixture
        const compressibleTwoPhaseMixture& mixture_;

        //- Mixture volumetric flux (from flux predictor)
        const surfaceScalarField& phiv_;

        //- Phase-1 mass flux (from flux predictor)
        const surfaceScalarField& alphaRhoPhi1_;

        //- Phase-2 mass flux (from flux predictor)
        const surfaceScalarField& alphaRhoPhi2_;

        //- Mixture transport model (constructed for mixture transport)
        autoPtr<compressible::momentumTransportModel>
            mixtureMomentumTransport_;

        //- Phase-1 transport model (constructed for two-phase transport)
        autoPtr<phaseCompressible::momentumTransportModel>
            momentumTransport1_;

        //- Phase-2 transport model (constructed for two-phase transport)
        autoPtr<phaseCompressible::momentumTransportModel>
            momentumTransport2_;

        //- Allow thermophysical transport model to access internals
        friend class twoPhaseThermophysicalTransportModel;


public:

    // Constructors

        //- Construct from components
        twoPhaseMomentumTransportModel
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const surfaceScalarField& phiv,
            const surfaceScalarField& alphaRhoPhi1,
            const surfaceScalarField& alphaRhoPhi2,
            const compressibleTwoPhaseMixture& mixture
        );

        //- Disallow default bitwise copy construction
        twoPhaseMomentumTransportModel
        (
            const twoPhaseMomentumTransportModel&
        ) = delete;


    //- Destructor
    ~twoPhaseMomentumTransportModel()
    {}


    // Member Functions

        // Access

            //- Is two-phase transport modeling enabled?
            bool twoPhaseTransport() const
            {
                return twoPhaseTransport_;
            }


        // Turbulent viscosity

            //- Return the turbulent kinematic viscosity (mixture or
            //  volume-fraction-weighted)
            tmp<volScalarField> nut() const;


        // Explicit methods (density-based formulation)

            //- Effective dynamic viscosity of the mixture
            tmp<volScalarField> muEff() const;

            //- Explicit viscous stress tensor [Pa]
            tmp<volSymmTensorField> tau
            (
                const volScalarField& muEff,
                const volTensorField& gradU
            ) const;

            //- Explicit cross-stress tensor [Pa]
            tmp<volTensorField> tauMC
            (
                const volScalarField& muEff,
                const volTensorField& gradU
            ) const;


        // Predict / correct

            //- Predict the phase or mixture transport models
            void predict();

            //- Correct the phase or mixture transport models
            void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const twoPhaseMomentumTransportModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
