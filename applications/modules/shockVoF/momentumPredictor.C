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

\*---------------------------------------------------------------------------*/

#include "shockVoF.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::shockVoF::momentumPredictor()
{
    // Complete hyperbolic + viscous advance for the time step
    if (!pimple.firstIter())
    {
        return;
    }

    // s-stage, second-order SSP Runge-Kutta, SSP(s,2)
    const label nStages = nRungeKuttaStages();

    // Length of each forward-Euler sub-stage
    const dimensionedScalar deltaT(runTime.deltaT()/scalar(nStages - 1));

    for (label stage = 0; stage < nStages; ++stage)
    {
        Info<< "Runge-Kutta stage " << stage + 1 << "/" << nStages << endl;

        // Interface normal, curvature and capillary coefficient of this stage
        interface_.correct();

        // High-order and low-order inviscid fluxes
        flux_.update();

        // Face blending factor theta
        const surfaceScalarField theta(fluxLimiter(deltaT));

        // Apply the flux limiter to blend the high- and low-order fluxes
        flux_.blend(theta);

        // Explicit forward-Euler sub-stage update
        solveEquations(deltaT);

        // Final stage: convex blend with the time-level-n state
        if (stage == nStages - 1)
        {
            const scalar an = scalar(1)/nStages;
            const scalar as = scalar(nStages - 1)/nStages;

            alpha1_ = an*alpha1_.oldTime() + as*alpha1_;

            if (pimple.flow())
            {
                alpha1Rho1_  = an*alpha1Rho1_.oldTime()  + as*alpha1Rho1_;
                alpha2Rho2_  = an*alpha2Rho2_.oldTime()  + as*alpha2Rho2_;
                rhoU_        = an*rhoU_.oldTime()        + as*rhoU_;
                alpha1RhoE1_ = an*alpha1RhoE1_.oldTime() + as*alpha1RhoE1_;
                alpha2RhoE2_ = an*alpha2RhoE2_.oldTime() + as*alpha2RhoE2_;
                rhoE_        = an*rhoE_.oldTime()        + as*rhoE_;
            }
        }

        // Recompute primitive variables
        updatePrimitives();
    }
}


// ************************************************************************* //
