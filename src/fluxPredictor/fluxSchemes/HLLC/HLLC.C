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

#include "HLLC.H"
#include "addToRunTimeSelectionTable.H"
#include "Switch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(HLLC, 0);
    addToRunTimeSelectionTable(fluxScheme, HLLC, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::HLLC::HLLC(const fvMesh&, const dictionary& dict)
:
    lowMachCorrection_
    (
        dict.lookupOrDefault<Switch>("lowMachCorrection", true)
    )
{
    Info<< "    lowMachCorrection: " << Switch(lowMachCorrection_) << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::HLLC::evaluateFlux
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

    // Output: contact-upwinded face volume fraction for the
    // surface-tension source (Garrick et al. 2017a)
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

    // Davis wave-speed estimates from the sound speeds of the reconstructed
    // face states
    const scalar SL = min(qL - c_pos, qR - c_neg);
    const scalar SR = max(qL + c_pos, qR + c_neg);

    // Mixture total energy at the faces
    const scalar rhoEL =
        alphaRhoE1_pos + alphaRhoE2_pos + 0.5*rho_pos*magSqr(U_pos);
    const scalar rhoER =
        alphaRhoE1_neg + alphaRhoE2_neg + 0.5*rho_neg*magSqr(U_neg);

    // Contact-wave speed S*; denom = phiR - phiL with phiK = rhoK*(SK - qK)
    const scalar denom =
        stabilise(rho_neg*(SR - qR) - rho_pos*(SL - qL), small);

    const scalar SStar =
        (rho_neg*qR*(SR - qR) - rho_pos*qL*(SL - qL)
          + p_pos - p_neg
          + sigmaKappaf*(alpha1_neg - alpha1_pos))
      / denom;

    // Star-region pressures
    const scalar pStarL =
        rho_pos*(qL - SL)*(qL - SStar) + p_pos;
    const scalar pStarR =
        rho_neg*(qR - SR)*(qR - SStar) + p_neg;

    // Low-Mach pressure fix of the star-region momentum and total-energy
    // fluxes, pf*(0, n, S*) (Chen et al. 2020, Eq. 3.9)
    scalar pf = 0;

    if (lowMachCorrection_)
    {
        // Capillary cutoff velocity from the alpha-weighted face phase
        // densities
        scalar Ucut = 0;

        const scalar alpha1f = alpha1_pos + alpha1_neg;
        const scalar alpha2f = (1 - alpha1_pos) + (1 - alpha1_neg);

        if
        (
            mag(sigmaKappaf) > 0
         && alpha1f > small
         && alpha2f > small
        )
        {
            const scalar rho1f = (alpha1Rho1_pos + alpha1Rho1_neg)/alpha1f;
            const scalar rho2f = (alpha2Rho2_pos + alpha2Rho2_neg)/alpha2f;

            Ucut = sqrt(mag(sigmaKappaf)/(rho1f + rho2f));
        }

        // Local Mach number of the face states (Chen et al. 2020, Eq. 3.14)
        const scalar M = min
        (
            scalar(1),
            max
            (
                max(mag(U_pos), Ucut)/max(c_pos, small),
                max(mag(U_neg), Ucut)/max(c_neg, small)
            )
        );

        // Mach-number scaling function (Chen et al. 2020, Eq. 3.13)
        const scalar f = M*sqrt(4 + sqr(1 - sqr(M)))/(1 + sqr(M));

        // pf = (f - 1)*p_d with the acoustic part of the star pressures,
        // p_d = phiL*phiR*(qR - qL)/(phiR - phiL) (Chen et al. 2020,
        // Eqs. 3.6 and 3.12)
        pf = (f - 1)*rho_pos*(SL - qL)*rho_neg*(SR - qR)*(qR - qL)/denom;
    }


    // Assemble fluxes over the four regions  L | L* | R* | R

    if (SL >= 0)
    {
        // Supersonic from the left - use left physical flux
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
    else if (SStar >= 0)
    {
        // Left star region
        const scalar omegaL = 1.0/stabilise(SL - SStar, small);
        const scalar wL = (SL - qL)*omegaL;

        // Mixture quantities
        const scalar rhoStarL  = wL*rho_pos;
        const vector rhoUStarL =
            omegaL*((SL - qL)*rhoUL + (pStarL - p_pos)*n);

        // Contact-advected scalars (simple density-ratio scaling)
        const scalar alpha1Rho1StarL  = wL*alpha1Rho1_pos;
        const scalar alpha2Rho2StarL  = wL*alpha2Rho2_pos;
        const scalar alphaRhoE1StarL  = wL*alphaRhoE1_pos;
        const scalar alphaRhoE2StarL  = wL*alphaRhoE2_pos;

        // Mixture total energy star state (Bharate & Mandal 2025, Eq. 16)
        const scalar rhoEStarL =
            ((SL - qL)*rhoEL - p_pos*qL + pStarL*SStar)
          / stabilise(SL - SStar, small);

        const scalar qStarA = SStar*magSf;

        phi             = rhoStarL*qStarA;
        phiUp           = rhoUStarL*qStarA + (pStarL + pf)*Sf;
        phiAlpha1Rho1   = alpha1Rho1StarL*qStarA;
        phiAlpha2Rho2   = alpha2Rho2StarL*qStarA;

        // Non-conservative volume-fraction advection with the
        // HLLC-consistent face velocity (Johnsen & Colonius 2006)
        const scalar qHatA = wL*qStarA;

        phiAlpha1       = alpha1_pos*qHatA;
        phiv            = qHatA;
        alphaf          = alpha1_pos;

        // Per-phase internal energy: Toro form (purely advective)
        const scalar FL_alphaRhoE1 = alphaRhoE1_pos*qL*magSf;
        const scalar FL_alphaRhoE2 = alphaRhoE2_pos*qL*magSf;

        phiAlpha1RhoE1 =
            FL_alphaRhoE1 + SL*magSf*(alphaRhoE1StarL - alphaRhoE1_pos);
        phiAlpha2RhoE2 =
            FL_alphaRhoE2 + SL*magSf*(alphaRhoE2StarL - alphaRhoE2_pos);

        // Mixture total energy: Toro form.  Satisfies
        //
        //     phiRhoE(sigmaKappaf, p) = phiRhoE(0, Pi) + sigmaKappaf*phiAlpha1
        //
        // with Pi = p - sigmaKappaf*alpha1, as in HLL
        const scalar FL_rhoE = (rhoEL + p_pos)*qL*magSf;
        phiRhoE = FL_rhoE + SL*magSf*(rhoEStarL - rhoEL) + pf*qStarA;
    }
    else if (SR > 0)
    {
        // Right star region
        const scalar omegaR = 1.0/stabilise(SR - SStar, small);
        const scalar wR = (SR - qR)*omegaR;

        const scalar rhoStarR  = wR*rho_neg;
        const vector rhoUStarR =
            omegaR*((SR - qR)*rhoUR + (pStarR - p_neg)*n);

        const scalar alpha1Rho1StarR  = wR*alpha1Rho1_neg;
        const scalar alpha2Rho2StarR  = wR*alpha2Rho2_neg;
        const scalar alphaRhoE1StarR  = wR*alphaRhoE1_neg;
        const scalar alphaRhoE2StarR  = wR*alphaRhoE2_neg;

        // Mixture total energy star state (Bharate & Mandal 2025, Eq. 16)
        const scalar rhoEStarR =
            ((SR - qR)*rhoER - p_neg*qR + pStarR*SStar)
          / stabilise(SR - SStar, small);

        const scalar qStarA = SStar*magSf;

        phi             = rhoStarR*qStarA;
        phiUp           = rhoUStarR*qStarA + (pStarR + pf)*Sf;
        phiAlpha1Rho1   = alpha1Rho1StarR*qStarA;
        phiAlpha2Rho2   = alpha2Rho2StarR*qStarA;

        // Non-conservative volume-fraction advection
        const scalar qHatA = wR*qStarA;

        phiAlpha1       = alpha1_neg*qHatA;
        phiv            = qHatA;
        alphaf          = alpha1_neg;

        // Per-phase internal energy: Toro form
        const scalar FR_alphaRhoE1 = alphaRhoE1_neg*qR*magSf;
        const scalar FR_alphaRhoE2 = alphaRhoE2_neg*qR*magSf;

        phiAlpha1RhoE1 =
            FR_alphaRhoE1 + SR*magSf*(alphaRhoE1StarR - alphaRhoE1_neg);
        phiAlpha2RhoE2 =
            FR_alphaRhoE2 + SR*magSf*(alphaRhoE2StarR - alphaRhoE2_neg);

        // Mixture total energy: Toro form
        const scalar FR_rhoE = (rhoER + p_neg)*qR*magSf;
        phiRhoE = FR_rhoE + SR*magSf*(rhoEStarR - rhoER) + pf*qStarA;
    }
    else
    {
        // Supersonic from the right - use right physical flux
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

    // Maximum signal speed for the CFL condition
    amaxSf = max(mag(SL), mag(SR))*magSf;
}


// ************************************************************************* //
