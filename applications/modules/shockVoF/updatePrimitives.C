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

Description
    Conservative-to-primitive decode after one Runge-Kutta stage.
    Phasic-internal-energy formulation with redundant mixture total energy.

\*---------------------------------------------------------------------------*/

#include "shockVoF.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::shockVoF::updatePrimitives()
{
    // Clamp volume fraction - always advected (see shockVoF.H)
    alpha1_.max(0);
    alpha1_.min(1);

    // Skip the primitive update if the flow is not solved.  The phase
    // pressures follow alpha1 through the capillary pressure
    if (!pimple.flow())
    {
        alpha1_.correctBoundaryConditions();
        alpha2_ = 1.0 - alpha1_;
        alpha2_.correctBoundaryConditions();
        mixture_.correctPhasePressures();

        return;
    }

    // With the flow solved, the pressure relaxation in mixture_.correct()
    // changes alpha1 and alpha2 again and evaluates their boundary
    // conditions after it
    alpha2_ = 1.0 - alpha1_;

    // Refresh mixture density and U
    mixture_.correctRho();

    U_ = rhoU_/rho_;
    U_.correctBoundaryConditions();

    // Conservative mixture specific internal energy from the 7th
    // equation.  This is the authoritative energy for pressure
    // decode at shocks.
    const volScalarField e
    (
        IOobject
        (
            "e",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        rhoE_/rho_ - 0.5*magSqr(U_)
    );

    // Full conservative-to-primitive decode
    mixture_.correct(e, interface_.sigmaKappa());

    // Reinitialize rhoE from the relaxed state
    rhoE_ == alpha1RhoE1_ + alpha2RhoE2_ + 0.5*rho_*magSqr(U_);
}


// ************************************************************************* //
