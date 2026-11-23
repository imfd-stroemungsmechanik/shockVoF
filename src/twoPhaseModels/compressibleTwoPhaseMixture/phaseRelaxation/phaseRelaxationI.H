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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline void
Foam::phaseRelaxation::relaxCell
(
    // EoS parameters (consts)
    const scalar g1,
    const scalar g2,
    const scalar pInf1,
    const scalar pInf2,
    const scalar Cv1,
    const scalar Cv2,
    // Inputs (cell-local)
    const scalar rho,
    const scalar eMixture,
    const scalar alphaRho1,
    const scalar alphaRho2,
    // In/out
    scalar& pc,
    scalar& alpha1,
    scalar& alpha2,
    scalar& alphaRhoE1,
    scalar& alphaRhoE2,
    scalar& p,
    scalar& T
) const
{
    // Limit pc so that the Step B decode keeps p2 (pc > 0) or p1 = p2 + pc
    // (pc < 0) above the floor
    auto limitPc = [&](const scalar a1, const scalar a2)
    {
        const scalar Gamma = a1/(g1 - 1.0) + a2/(g2 - 1.0);
        const scalar Pi    = a1*g1*pInf1/(g1 - 1.0) + a2*g2*pInf2/(g2 - 1.0);
        const scalar eStar =
            max(rho*eMixture - Pi - Gamma*Foam::phaseFloors::p, scalar(0));

        pc = min
        (
            max(pc, -(g2 - 1.0)*eStar/max(a2, small)),
            (g1 - 1.0)*eStar/max(a1, small)
        );
    };

    limitPc(alpha1, 1.0 - alpha1);

    // Step A: volume-fraction relaxation to the equilibrium p1 - p2 = pc
    // with the pre-relaxation interface pressure pI_0 (Bharate & Mandal
    // 2026, Eqs. 37-40).  Skipped in a (nearly) single-phase cell, where
    // the relaxation is the identity but would decode a vanished phase's
    // pressure.
    const scalar alpha1_0 = alpha1;
    const scalar alpha2_0 = 1.0 - alpha1_0;

    if (min(alpha1_0, alpha2_0) > Foam::phaseFloors::relaxAlpha)
    {
        // Pre-relaxation phasic pressures from SG EoS
        const scalar p1_0 = max
        (
            (g1 - 1.0)*alphaRhoE1/alpha1_0 - g1*pInf1,
            -pInf1 + Foam::phaseFloors::p
        );
        const scalar p2_0 = max
        (
            (g2 - 1.0)*alphaRhoE2/alpha2_0 - g2*pInf2,
            -pInf2 + Foam::phaseFloors::p
        );

        // Pre-relaxation phasic densities
        const scalar rho1_0 = max(alphaRho1/alpha1_0, Foam::phaseFloors::rho);
        const scalar rho2_0 = max(alphaRho2/alpha2_0, Foam::phaseFloors::rho);

        // Acoustic impedances Z_k = rho_k*c_k
        const scalar Z1_0 = Foam::sqrt(rho1_0*g1*(p1_0 + pInf1));
        const scalar Z2_0 = Foam::sqrt(rho2_0*g2*(p2_0 + pInf2));

        // Interface pressure of the pre-relaxation state with the
        // impedance weight theta of p1 (Bharate & Mandal 2026, Eq. 3)
        const scalar theta = Z2_0/max(Z1_0 + Z2_0, small);
        const scalar pI_0 = theta*p1_0 + (1.0 - theta)*p2_0;

        // Coefficients C_k = gamma_k*pInf_k + (gamma_k - 1)*pI_0
        const scalar C1 = g1*pInf1 + (g1 - 1.0)*pI_0;
        const scalar C2 = g2*pInf2 + (g2 - 1.0)*pI_0;

        // Monic quadratic p*^2 + b*p* + d = 0 for the relaxed phase-2
        // pressure p* = p2*, with p1* = p* + pc (Bharate & Mandal 2026,
        // Eq. 39, b simplified with alpha1 + alpha2 = 1)
        const scalar b = alpha1_0*C2 + alpha2_0*C1
                       - alpha1_0*p1_0 - alpha2_0*p2_0 + pc;
        const scalar d = -(alpha1_0*(p1_0 - pc)*C2 + alpha2_0*p2_0*(C1 + pc));

        const scalar disc = max(b*b - 4.0*d, 0.0);
        const scalar sqrtDisc = Foam::sqrt(disc);

        // Root with the positive square root, the physically admissible
        // one (Bharate & Mandal 2026, Eq. 40), evaluated without
        // cancellation.
        // For pc = 0 and positive phase pressures, d < 0 and this root is
        // positive; for pc != 0 it can fall below the floor, which then
        // bounds it.
        const scalar pRelax = max
        (
            b > 0
          ? -2.0*d/(b + sqrtDisc)
          : 0.5*(sqrtDisc - b),
            Foam::phaseFloors::p
        );

        // Relaxed volume fraction alpha1* = alpha1*(p1 + C1)/(p1* + C1)
        // (Bharate & Mandal 2026, Eq. 38)
        const scalar denom = max(pRelax + pc + C1, small);
        scalar alpha1_new = alpha1_0*(p1_0 + C1)/denom;
        alpha1_new = max(min(alpha1_new, 1.0), 0.0);

        alpha1 = alpha1_new;
        alpha2 = 1.0 - alpha1_new;
    }

    // Step B: Kapila mixture-EoS decode of the phase-2 pressure from the
    // conservative mixture energy with p1 = p2 + pc (Bharate & Mandal
    // 2026, Eq. 41); the relaxed pressure of Step A is discarded.  The mixture
    // pressure is p = alpha1*p1 + alpha2*p2 = p2 + alpha1*pc
    limitPc(alpha1, alpha2);

    const scalar Gamma = alpha1/(g1 - 1.0) + alpha2/(g2 - 1.0);
    const scalar Pi    = alpha1*g1*pInf1/(g1 - 1.0)
                       + alpha2*g2*pInf2/(g2 - 1.0);

    const scalar p2 = max
    (
        (rho*eMixture - Pi - alpha1*pc/(g1 - 1.0))/max(Gamma, small),
        Foam::phaseFloors::p
    );
    const scalar p1 = p2 + pc;

    p = p2 + alpha1*pc;

    // Step C: reinitialize per-phase internal energies from p1, p2
    alphaRhoE1 = alpha1*(p1 + g1*pInf1)/(g1 - 1.0);
    alphaRhoE2 = alpha2*(p2 + g2*pInf2)/(g2 - 1.0);

    // Mixture temperature: the heat-capacity-weighted mean of the phase
    // temperatures, sum(alphaRho_k*Cv_k*T_k)/sum(alphaRho_k*Cv_k)
    const scalar thermalE = alpha1*(p1 + pInf1)/(g1 - 1.0)
        + alpha2*(p2 + pInf2)/(g2 - 1.0);

    const scalar rhoCvMix = max
    (
        alphaRho1*Cv1 + alphaRho2*Cv2,
        Foam::phaseFloors::rho
    );

    T = thermalE/rhoCvMix;
}


// ************************************************************************* //
