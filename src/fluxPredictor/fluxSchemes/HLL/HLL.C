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

#include "HLL.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(HLL, 0);
    addToRunTimeSelectionTable(fluxScheme, HLL, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::HLL::HLL(const fvMesh&, const dictionary&)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::HLL::evaluateFlux
(
    scalar& phi,
    vector& phiUp,
    scalar& phiAlpha1Rho1,
    scalar& phiAlpha2Rho2,
    scalar& phiAlpha1,
    scalar& phiAlpha1RhoE1,
    scalar& phiAlpha2RhoE2,
    scalar& phiRhoE,
    scalar& phiv,
    scalar& amaxSf,

    // Output: HLL contact-upwinded face volume fraction
    scalar& alphaf,

    // Left state
    const scalar& p_pos,
    const vector& U_pos,
    const scalar& rho_pos,
    const scalar& c_pos,
    const scalar& alpha1_pos,
    const scalar& alpha1Rho1_pos,
    const scalar& alpha2Rho2_pos,
    const scalar& alphaRhoE1_pos,
    const scalar& alphaRhoE2_pos,

    // Right state
    const scalar& p_neg,
    const vector& U_neg,
    const scalar& rho_neg,
    const scalar& c_neg,
    const scalar& alpha1_neg,
    const scalar& alpha1Rho1_neg,
    const scalar& alpha2Rho2_neg,
    const scalar& alphaRhoE1_neg,
    const scalar& alphaRhoE2_neg,

    // Geometry
    const vector& Sf,
    const scalar& magSf,

    // Face capillary coefficient sigma*kappa_f [Pa]
    const scalar& sigmaKappaf
) const
{
    const vector n = Sf/max(magSf, small);

    const scalar qL = U_pos & n;
    const scalar qR = U_neg & n;

    const vector rhoUL = rho_pos*U_pos;
    const vector rhoUR = rho_neg*U_neg;

    // Davis wave-speed estimates
    const scalar SL = min(qL - c_pos, qR - c_neg);
    const scalar SR = max(qL + c_pos, qR + c_neg);

    // Mixture total energy at the faces
    const scalar rhoEL =
        alphaRhoE1_pos + alphaRhoE2_pos + 0.5*rho_pos*magSqr(U_pos);
    const scalar rhoER =
        alphaRhoE1_neg + alphaRhoE2_neg + 0.5*rho_neg*magSqr(U_neg);

    // Capillary-reduced pressure
    const scalar PiL = p_pos - sigmaKappaf*alpha1_pos;
    const scalar PiR = p_neg - sigmaKappaf*alpha1_neg;

    if (SL >= 0)
    {
        // Supersonic from the left - left physical flux
        const scalar qA = qL*magSf;

        phi             = rho_pos*qA;
        phiUp           = rhoUL*qA + p_pos*Sf;
        phiAlpha1Rho1   = alpha1Rho1_pos*qA;
        phiAlpha2Rho2   = alpha2Rho2_pos*qA;
        phiAlpha1       = alpha1_pos*qA;
        phiAlpha1RhoE1  = alphaRhoE1_pos*qA;
        phiAlpha2RhoE2  = alphaRhoE2_pos*qA;
        phiRhoE         = (rhoEL + p_pos)*qA;
        phiv            = qL*magSf;
        alphaf          = alpha1_pos;
    }
    else if (SR <= 0)
    {
        // Supersonic from the right - right physical flux
        const scalar qA = qR*magSf;

        phi             = rho_neg*qA;
        phiUp           = rhoUR*qA + p_neg*Sf;
        phiAlpha1Rho1   = alpha1Rho1_neg*qA;
        phiAlpha2Rho2   = alpha2Rho2_neg*qA;
        phiAlpha1       = alpha1_neg*qA;
        phiAlpha1RhoE1  = alphaRhoE1_neg*qA;
        phiAlpha2RhoE2  = alphaRhoE2_neg*qA;
        phiRhoE         = (rhoER + p_neg)*qA;
        phiv            = qR*magSf;
        alphaf          = alpha1_neg;
    }
    else
    {
        // Subsonic - single HLL averaged region:
        //   F_HLL = (SR*F_L - SL*F_R + SL*SR*(U_R - U_L))/(SR - SL)
        const scalar invDS = 1.0/stabilise(SR - SL, small);

        // Conservative fluxes.  Every advected scalar X has physical flux
        // X*q, so F_HLL collapses to the compact form below.
        phi =
            (SR*rho_pos*qL - SL*rho_neg*qR + SL*SR*(rho_neg - rho_pos))
           *invDS*magSf;

        phiAlpha1Rho1 =
            (SR*alpha1Rho1_pos*qL - SL*alpha1Rho1_neg*qR
              + SL*SR*(alpha1Rho1_neg - alpha1Rho1_pos))*invDS*magSf;
        phiAlpha2Rho2 =
            (SR*alpha2Rho2_pos*qL - SL*alpha2Rho2_neg*qR
              + SL*SR*(alpha2Rho2_neg - alpha2Rho2_pos))*invDS*magSf;

        // Per-phase internal energy: purely advective (no pressure work)
        phiAlpha1RhoE1 =
            (SR*alphaRhoE1_pos*qL - SL*alphaRhoE1_neg*qR
              + SL*SR*(alphaRhoE1_neg - alphaRhoE1_pos))*invDS*magSf;
        phiAlpha2RhoE2 =
            (SR*alphaRhoE2_pos*qL - SL*alphaRhoE2_neg*qR
              + SL*SR*(alphaRhoE2_neg - alphaRhoE2_pos))*invDS*magSf;

        // Non-conservative volume-fraction advection
        const scalar rhoHLL =
            (SR*rho_neg - SL*rho_pos - (rho_neg*qR - rho_pos*qL))*invDS;
        const scalar rhoUnHLL =
            (SR*rho_neg*qR - SL*rho_pos*qL
              - ((rho_neg*qR*qR + PiR) - (rho_pos*qL*qL + PiL)))*invDS;
        const scalar uI = rhoUnHLL/stabilise(rhoHLL, small);

        phiv      = uI*magSf;
        alphaf    = (uI >= 0) ? alpha1_pos : alpha1_neg;
        phiAlpha1 = alphaf*uI*magSf;

        // Momentum
        const vector FmL = rhoUL*qL + PiL*n;
        const vector FmR = rhoUR*qR + PiR*n;
        phiUp =
            (
                (SR*FmL - SL*FmR + SL*SR*(rhoUR - rhoUL))*invDS
              + sigmaKappaf*alphaf*n
            )*magSf;

        // Mixture total energy; satisfies
        //
        //     phiRhoE(sigmaKappaf, p) = phiRhoE(0, Pi) + sigmaKappaf*phiAlpha1
        const scalar FEL = (rhoEL + PiL)*qL;
        const scalar FER = (rhoER + PiR)*qR;
        phiRhoE =
            (
                (SR*FEL - SL*FER + SL*SR*(rhoER - rhoEL))*invDS
              + sigmaKappaf*alphaf*uI
            )*magSf;
    }

    // Maximum signal speed for the CFL condition (same SL, SR as the flux)
    amaxSf = max(mag(SL), mag(SR))*magSf;
}


// ************************************************************************* //
