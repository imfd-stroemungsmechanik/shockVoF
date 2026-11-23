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

\*---------------------------------------------------------------------------*/

#include "twoPhaseMomentumTransportModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseMomentumTransportModel::twoPhaseMomentumTransportModel
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const surfaceScalarField& phiv,
    const surfaceScalarField& alphaRhoPhi1,
    const surfaceScalarField& alphaRhoPhi2,
    const compressibleTwoPhaseMixture& mixture
)
:
    twoPhaseTransport_(false),
    mixture_(mixture),
    phiv_(phiv),
    alphaRhoPhi1_(alphaRhoPhi1),
    alphaRhoPhi2_(alphaRhoPhi2)
{
    // Read the momentum transport dictionary to determine transport type
    {
        IOdictionary momentumTransport
        (
            IOobject
            (
                momentumTransportModel::typeName,
                U.time().constant(),
                U.db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        word simulationType
        (
            momentumTransport.lookup("simulationType")
        );

        if (simulationType == "twoPhaseTransport")
        {
            twoPhaseTransport_ = true;
        }
    }

    if (twoPhaseTransport_)
    {
        Info<< "Selecting two-phase transport model" << endl;

        momentumTransport1_ =
            phaseCompressible::momentumTransportModel::New
            (
                mixture_.phase1().alpha(),
                mixture_.phase1().rho(),
                U,
                alphaRhoPhi1_,
                phiv_,
                mixture_.phase1()
            );

        momentumTransport2_ =
            phaseCompressible::momentumTransportModel::New
            (
                mixture_.phase2().alpha(),
                mixture_.phase2().rho(),
                U,
                alphaRhoPhi2_,
                phiv_,
                mixture_.phase2()
            );

        // Form nut from the initial state, as for the mixture model
        momentumTransport1_->validate();
        momentumTransport2_->validate();
    }
    else
    {
        // Mixture transport: single model for the entire mixture.
        // phi here is the mixture mass flux; mixture_ is-a viscosity.
        mixtureMomentumTransport_ =
            compressible::momentumTransportModel::New
            (
                rho,
                U,
                phi,
                mixture_
            );

        mixtureMomentumTransport_->validate();
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::twoPhaseMomentumTransportModel::nut() const
{
    if (twoPhaseTransport_)
    {
        return
        (
            mixture_.phase1().alpha()*momentumTransport1_->nut()
          + mixture_.phase2().alpha()*momentumTransport2_->nut()
        );
    }
    else
    {
        return mixtureMomentumTransport_->nut();
    }
}


Foam::tmp<Foam::volScalarField>
Foam::twoPhaseMomentumTransportModel::muEff() const
{
    if (twoPhaseTransport_)
    {
        return volScalarField::New
        (
            "muEff",
            mixture_.phase1().alpha()
           *(
                mixture_.phase1().mu()
              + mixture_.phase1().rho()*momentumTransport1_->nut()
            )
          + mixture_.phase2().alpha()
           *(
                mixture_.phase2().mu()
              + mixture_.phase2().rho()*momentumTransport2_->nut()
            )
        );
    }
    else
    {
        return volScalarField::New
        (
            "muEff",
            mixture_.mu() + mixture_.rho()*mixtureMomentumTransport_->nut()
        );
    }
}


Foam::tmp<Foam::volSymmTensorField>
Foam::twoPhaseMomentumTransportModel::tau
(
    const volScalarField& muEff,
    const volTensorField& gradU
) const
{
    return volSymmTensorField::New
    (
        "tau",
        muEff*dev(twoSymm(gradU))
    );
}


Foam::tmp<Foam::volTensorField>
Foam::twoPhaseMomentumTransportModel::tauMC
(
    const volScalarField& muEff,
    const volTensorField& gradU
) const
{
    return volTensorField::New
    (
        "tauMC",
        muEff*dev2(T(gradU))
    );
}


void Foam::twoPhaseMomentumTransportModel::predict()
{
    if (twoPhaseTransport_)
    {
        momentumTransport1_->predict();
        momentumTransport2_->predict();
    }
    else
    {
        mixtureMomentumTransport_->predict();
    }
}


void Foam::twoPhaseMomentumTransportModel::correct()
{
    if (twoPhaseTransport_)
    {
        momentumTransport1_->correct();
        momentumTransport2_->correct();
    }
    else
    {
        mixtureMomentumTransport_->correct();
    }
}


// ************************************************************************* //
