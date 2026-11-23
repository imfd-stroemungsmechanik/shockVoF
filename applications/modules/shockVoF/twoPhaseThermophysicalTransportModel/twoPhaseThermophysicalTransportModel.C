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

\*---------------------------------------------------------------------------*/

#include "twoPhaseThermophysicalTransportModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseThermophysicalTransportModel::
twoPhaseThermophysicalTransportModel
(
    const compressibleTwoPhaseMixture& mixture,
    const twoPhaseMomentumTransportModel& momentumTransport
)
:
    mixture_(mixture),
    momentumTransport_(momentumTransport)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::twoPhaseThermophysicalTransportModel::kappaEff() const
{
    if (momentumTransport_.twoPhaseTransport())
    {
        // Per-phase turbulent viscosities (friend access)
        const tmp<volScalarField> tnut1
        (
            momentumTransport_.momentumTransport1_->nut()
        );
        const tmp<volScalarField> tnut2
        (
            momentumTransport_.momentumTransport2_->nut()
        );

        return volScalarField::New
        (
            "kappaEff",
            mixture_.phase1().alpha()
           *(
                mixture_.phase1().kappa()
              + mixture_.phase1().rho()
               *mixture_.phase1().eos().Cp()
               *tnut1()
            )
          + mixture_.phase2().alpha()
           *(
                mixture_.phase2().kappa()
              + mixture_.phase2().rho()
               *mixture_.phase2().eos().Cp()
               *tnut2()
            )
        );
    }
    else
    {
        // Mixture mode: single turbulence model, single nut
        const tmp<volScalarField> tnut(momentumTransport_.nut());

        return volScalarField::New
        (
            "kappaEff",
            mixture_.phase1().alpha()
           *(
                mixture_.phase1().kappa()
              + mixture_.phase1().rho()
               *mixture_.phase1().eos().Cp()
               *tnut()
            )
          + mixture_.phase2().alpha()
           *(
                mixture_.phase2().kappa()
              + mixture_.phase2().rho()
               *mixture_.phase2().eos().Cp()
               *tnut()
            )
        );
    }
}


void Foam::twoPhaseThermophysicalTransportModel::predict()
{}


void Foam::twoPhaseThermophysicalTransportModel::correct()
{}


// ************************************************************************* //
