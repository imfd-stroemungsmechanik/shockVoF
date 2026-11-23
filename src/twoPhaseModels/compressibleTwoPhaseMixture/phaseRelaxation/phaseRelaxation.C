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

#include "phaseRelaxation.H"
#include "phaseModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseRelaxation::phaseRelaxation
(
    const stiffenedGasEoS& sgEoS1,
    const stiffenedGasEoS& sgEoS2
)
:
    sgEoS1_(sgEoS1),
    sgEoS2_(sgEoS2)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::phaseRelaxation::update
(
    phaseModel& phase1,
    phaseModel& phase2,
    const volScalarField& rho,
    const volScalarField& eMixture,
    volScalarField& pc,
    volScalarField& p,
    volScalarField& T
) const
{
    // EoS parameters as raw scalars (lifted out of the loop)
    const scalar g1    = sgEoS1_.gamma().value();
    const scalar g2    = sgEoS2_.gamma().value();
    const scalar pInf1 = sgEoS1_.pInf().value();
    const scalar pInf2 = sgEoS2_.pInf().value();
    const scalar Cv1   = sgEoS1_.Cv().value();
    const scalar Cv2   = sgEoS2_.Cv().value();

    // Field references
    volScalarField& alpha1     = phase1.alpha();
    volScalarField& alpha2     = phase2.alpha();
    volScalarField& alphaRhoE1 = phase1.alphaE();
    volScalarField& alphaRhoE2 = phase2.alphaE();

    // --- Internal field cell loop ---
    {
        const scalarField& rhoI       = rho.primitiveField();
        const scalarField& eMixtureI  = eMixture.primitiveField();
        const scalarField& alphaRho1I = phase1.alphaRho().primitiveField();
        const scalarField& alphaRho2I = phase2.alphaRho().primitiveField();

        scalarField& pcI         = pc.primitiveFieldRef();
        scalarField& alpha1I     = alpha1.primitiveFieldRef();
        scalarField& alpha2I     = alpha2.primitiveFieldRef();
        scalarField& alphaRhoE1I = alphaRhoE1.primitiveFieldRef();
        scalarField& alphaRhoE2I = alphaRhoE2.primitiveFieldRef();
        scalarField& pI          = p.primitiveFieldRef();
        scalarField& TI          = T.primitiveFieldRef();

        forAll(rhoI, celli)
        {
            relaxCell
            (
                g1,
                g2,
                pInf1,
                pInf2,
                Cv1,
                Cv2,
                rhoI[celli],
                eMixtureI[celli],
                alphaRho1I[celli],
                alphaRho2I[celli],
                pcI[celli],
                alpha1I[celli],
                alpha2I[celli],
                alphaRhoE1I[celli],
                alphaRhoE2I[celli],
                pI[celli],
                TI[celli]
            );
        }
    }

    // Boundary conditions.  Coupled patches are not relaxed in place:
    // correctBoundaryConditions() sets them from the relaxed cells
    alpha1.correctBoundaryConditions();
    alpha2.correctBoundaryConditions();
    alphaRhoE1.correctBoundaryConditions();
    alphaRhoE2.correctBoundaryConditions();
    p.correctBoundaryConditions();
    T.correctBoundaryConditions();
    pc.correctBoundaryConditions();
}


// ************************************************************************* //
