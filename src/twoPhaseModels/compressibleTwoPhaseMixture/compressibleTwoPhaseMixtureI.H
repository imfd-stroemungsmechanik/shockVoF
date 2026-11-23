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

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline void Foam::compressibleTwoPhaseMixture::primitivesToState
(
    const scalar p1,
    const scalar p2,
    const scalar alpha1,
    const scalar rho1,
    const scalar rho2,
    scalar& rho,
    scalar& c,
    scalar& alphaRho1,
    scalar& alphaRho2,
    scalar& alphaRhoE1,
    scalar& alphaRhoE2
) const
{
    // alpha_k * rho_k * e_k = alpha_k * (p + g_k*pInf_k)/(g_k - 1)
    const scalar g1    = sgEoS1_.gamma().value();
    const scalar g2    = sgEoS2_.gamma().value();
    const scalar pInf1 = sgEoS1_.pInf().value();
    const scalar pInf2 = sgEoS2_.pInf().value();

    const scalar alpha2 = 1.0 - alpha1;

    // Partial and mixture densities
    alphaRho1 = alpha1*rho1;
    alphaRho2 = alpha2*rho2;
    rho = alphaRho1 + alphaRho2;

    // Per-phase internal energy densities
    alphaRhoE1 = alpha1*(p1 + g1*pInf1)/(g1 - 1.0);
    alphaRhoE2 = alpha2*(p2 + g2*pInf2)/(g2 - 1.0);

    // Frozen mixture sound speed
    c = sqrt
    (
        max
        (
            (
                alpha1*g1*(p1 + pInf1)
              + alpha2*g2*(p2 + pInf2)
            )
          / rho,
            small
        )
    );
}


// ************************************************************************* //
