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
    Solve the 6-equation model and its redundant mixture-total-energy
    (7th) equation for one explicit Runge-Kutta sub-stage (phasic-
    internal-energy formulation).

\*---------------------------------------------------------------------------*/

#include "shockVoF.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcLaplacian.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::shockVoF::solveEquations
(
    const dimensionedScalar& deltaT
)
{
    // Divergences shared by the phase-fraction equation and the phasic
    // pressure-work sources
    const volScalarField divPhiv("divPhiv", fvc::div(flux_.phiv()));
    const volScalarField divPhiAlpha1
    (
        "divPhiAlpha1",
        fvc::div(flux_.phiAlpha1())
    );

    // The flow equations are advanced before the phase fraction, so that
    // every term of this forward-Euler stage is evaluated on the stage-start
    // state.
    if (pimple.flow())
    {
        // Phase density (mass continuity) equations
        {
            alpha1Rho1_ -=
                deltaT
              * (
                    fvc::div(flux_.phiAlpha1Rho1())
                );
            alpha2Rho2_ -=
                deltaT
              * (
                    fvc::div(flux_.phiAlpha2Rho2())
                );
        }


        // Continuum-surface-force (CSF) capillary force (Brackbill et al.
        // 1992) in the well-balanced form of Garrick et al. (2017a)
        const surfaceScalarField& sigmaKappaf = interface_.sigmaKappaf();

        const volVectorField fSigma
        (
            "fSigma",
            fvc::div(sigmaKappaf*mesh.Sf()*flux_.alphaf())
          - alpha1_*fvc::div(sigmaKappaf*mesh.Sf())
        );


        // Effective viscosity and velocity gradient, formed once for the
        // viscous terms of the momentum and the total-energy equation
        const volScalarField muEff(momentumTransport_.muEff());
        const volTensorField gradU(fvc::grad(U_));


        // Momentum equation
        {
            // Cross-stress tensor for momentum equation
            const volTensorField tauMC
            (
                "tauMC",
                momentumTransport_.tauMC(muEff, gradU)
            );

            rhoU_ -=
                deltaT
              * (
                    fvc::div(flux_.phiUp())
                  - fvc::laplacian(muEff, U_)
                  - fvc::div(tauMC)
                );

            rhoU_ += deltaT*fSigma;
        }


        // Per-phase internal energy equations with the pressure work of
        // the phase pressures p1 = p + alpha2*pc and p2 = p - alpha1*pc
        {
            alpha1RhoE1_ -=
                deltaT
              * (
                    fvc::div(flux_.phiAlpha1RhoE1())
                  + alpha1_*mixture_.p1()*divPhiv
                );

            alpha2RhoE2_ -=
                deltaT
              * (
                    fvc::div(flux_.phiAlpha2RhoE2())
                  + alpha2_*mixture_.p2()*divPhiv
                );
        }


        // Total energy equation
        {
            // Full stress tensor
            const volSymmTensorField tau
            (
                "tau",
                momentumTransport_.tau(muEff, gradU)
            );

            const surfaceScalarField tauDotU
            (
                "tauDotU",
                (fvc::interpolate(tau) & fvc::interpolate(U_)) & mesh.Sf()
            );

            rhoE_ -=
                deltaT
              * (
                    fvc::div(flux_.phiRhoE())
                  - fvc::div(tauDotU)
                  - fvc::laplacian
                    (
                        thermophysicalTransport_.kappaEff(),
                        mixture_.T()
                    )
                );

            // Surface-tension work source
            rhoE_ += deltaT*(U_ & fSigma);
        }
    }


    // Phase fraction equation
    {
        alpha1_ -=
            deltaT
          * (
                divPhiAlpha1
              - alpha1_*divPhiv
            );
    }
}


// ************************************************************************* //
