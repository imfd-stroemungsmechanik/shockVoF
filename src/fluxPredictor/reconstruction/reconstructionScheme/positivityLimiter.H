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

InNamespace
    Foam

Description
    Zhang-Shu scaling limiter for the reconstructed face states.

    A single factor theta per cell scales all six primitives toward the
    cell value so that p_k >= phaseFloors::p and rho_k >= phaseFloors::rho
    on every reconstructed face.  The density constraint of a phase with
    alpha_k <= phaseFloors::alpha is skipped and its face density clamped
    instead; alpha1 is not constrained.  Faces of non-coupled patches are not
    limited, only their phase pressures are floored.  The neighbor side of
    coupled patches is replaced by the remote owner side.

    References:
    \verbatim
        Zhang, X. & Shu, C.-W. (2010).
        On positivity-preserving high order discontinuous Galerkin
        schemes for compressible Euler equations on rectangular meshes.
        Journal of Computational Physics, 229(23), 8918-8934.
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef positivityLimiter_H
#define positivityLimiter_H

#include "volFields.H"
#include "surfaceFields.H"
#include "fvMesh.H"
#include "swapCoupledOwnToNei.H"
#include "phaseFloors.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Compute the Zhang-Shu scaling parameter for a single face value
//  with a one-sided lower bound (positivity)
//
//  Returns the maximum theta in [0, 1] such that
//      phi_C + theta * (phi_f - phi_C) >= phiMin
//
inline scalar zhangShuThetaLower
(
    const scalar phiFace,
    const scalar phiCell,
    const scalar phiMin
)
{
    if (phiFace >= phiMin)
    {
        return 1.0;
    }

    // phiFace < phiMin: solve phiCell + t*(phiFace - phiCell) = phiMin
    const scalar denom = phiCell - phiFace;

    return (denom > vSmall)
      ? max(scalar(0), (phiCell - phiMin) / denom)
      : 0.0;
}


//- Zhang-Shu theta of the phase densities; skips a phase with
//  alpha_k <= phaseFloors::alpha
inline scalar phaseDensityTheta
(
    const scalar alpha1Cell,
    const scalar rho1Face,
    const scalar rho1Cell,
    const scalar rho2Face,
    const scalar rho2Cell
)
{
    scalar theta = 1.0;

    // Phase 1 present: alpha1 > phaseFloors::alpha
    if (alpha1Cell > Foam::phaseFloors::alpha)
    {
        theta = min
        (
            theta,
            zhangShuThetaLower(rho1Face, rho1Cell, Foam::phaseFloors::rho)
        );
    }

    // Phase 2 present: alpha2 = 1 - alpha1 > phaseFloors::alpha
    if (alpha1Cell < 1.0 - Foam::phaseFloors::alpha)
    {
        theta = min
        (
            theta,
            zhangShuThetaLower(rho2Face, rho2Cell, Foam::phaseFloors::rho)
        );
    }

    return theta;
}


//- Zhang-Shu limiter for all primitives (see file description).
//  Min-accumulates theta into positivityLimiter if non-null
inline void positivityLimit
(
    surfaceScalarField& alpha1Own,
    surfaceScalarField& alpha1Nei,
    surfaceScalarField& p1Own,
    surfaceScalarField& p1Nei,
    surfaceScalarField& p2Own,
    surfaceScalarField& p2Nei,
    surfaceVectorField& UOwn,
    surfaceVectorField& UNei,
    surfaceScalarField& rho1Own,
    surfaceScalarField& rho1Nei,
    surfaceScalarField& rho2Own,
    surfaceScalarField& rho2Nei,
    const volScalarField& alpha1,
    const volScalarField& p1,
    const volScalarField& p2,
    const volVectorField& U,
    const volScalarField& rho1,
    const volScalarField& rho2,
    volScalarField* positivityLimiter = nullptr
)
{
    const fvMesh& mesh = alpha1.mesh();
    const labelUList& owner = mesh.owner();
    const labelUList& neighbor = mesh.neighbour();


    // Pass 1: per-cell theta = min over the reconstructed faces of the
    // (p1, p2, rho1, rho2) constraints

    scalarField theta(mesh.nCells(), 1.0);

    // --- Internal faces ---
    forAll(owner, facei)
    {
        const label o = owner[facei];
        const label n = neighbor[facei];

        // Owner-side constraints
        theta[o] = min
        (
            theta[o],
            min
            (
                min
                (
                    zhangShuThetaLower
                    (
                        p1Own[facei], p1[o], Foam::phaseFloors::p
                    ),
                    zhangShuThetaLower
                    (
                        p2Own[facei], p2[o], Foam::phaseFloors::p
                    )
                ),
                phaseDensityTheta
                (
                    alpha1[o],
                    rho1Own[facei], rho1[o],
                    rho2Own[facei], rho2[o]
                )
            )
        );

        // Neighbor-side constraints
        theta[n] = min
        (
            theta[n],
            min
            (
                min
                (
                    zhangShuThetaLower
                    (
                        p1Nei[facei], p1[n], Foam::phaseFloors::p
                    ),
                    zhangShuThetaLower
                    (
                        p2Nei[facei], p2[n], Foam::phaseFloors::p
                    )
                ),
                phaseDensityTheta
                (
                    alpha1[n],
                    rho1Nei[facei], rho1[n],
                    rho2Nei[facei], rho2[n]
                )
            )
        );
    }

    // Coupled patches, owner side: the local cell's reconstruction,
    // constrained like the internal faces.  A non-coupled patch carries the
    // boundary-condition state instead, which is not part of the cell's
    // reconstruction and so does not limit it
    forAll(alpha1.boundaryField(), patchi)
    {
        if (!alpha1.boundaryField()[patchi].coupled())
        {
            continue;
        }

        const labelUList& faceCells =
            alpha1.boundaryField()[patchi].patch().faceCells();

        const scalarField& bP1   = p1Own.boundaryField()[patchi];
        const scalarField& bP2   = p2Own.boundaryField()[patchi];
        const scalarField& bRho1 = rho1Own.boundaryField()[patchi];
        const scalarField& bRho2 = rho2Own.boundaryField()[patchi];

        forAll(faceCells, facei)
        {
            const label c = faceCells[facei];

            theta[c] = min
            (
                theta[c],
                min
                (
                    min
                    (
                        zhangShuThetaLower
                        (
                            bP1[facei], p1[c], Foam::phaseFloors::p
                        ),
                        zhangShuThetaLower
                        (
                            bP2[facei], p2[c], Foam::phaseFloors::p
                        )
                    ),
                    phaseDensityTheta
                    (
                        alpha1[c],
                        bRho1[facei], rho1[c],
                        bRho2[facei], rho2[c]
                    )
                )
            );
        }
    }


    // Diagnostics: min-accumulate the per-cell theta into the caller's
    // field (reset per time step in reconstructionScheme::reconstruct())
    if (positivityLimiter)
    {
        scalarField& pl = positivityLimiter->primitiveFieldRef();

        forAll(pl, celli)
        {
            pl[celli] = min(pl[celli], theta[celli]);
        }

        positivityLimiter->correctBoundaryConditions();
    }


    // Pass 2: apply the single per-cell theta to all primitives,
    // phi_f = phi_C + theta_C*(phi_f^rec - phi_C)

    // --- Internal faces ---
    forAll(owner, facei)
    {
        const label o = owner[facei];
        const label n = neighbor[facei];

        if (theta[o] < 1.0)
        {
            const scalar t = theta[o];

            alpha1Own[facei] = alpha1[o] + t*(alpha1Own[facei] - alpha1[o]);
            p1Own[facei]     = p1[o]     + t*(p1Own[facei]     - p1[o]);
            p2Own[facei]     = p2[o]     + t*(p2Own[facei]     - p2[o]);
            UOwn[facei]      = U[o]      + t*(UOwn[facei]      - U[o]);
            rho1Own[facei]   = rho1[o]   + t*(rho1Own[facei]   - rho1[o]);
            rho2Own[facei]   = rho2[o]   + t*(rho2Own[facei]   - rho2[o]);
        }

        if (theta[n] < 1.0)
        {
            const scalar t = theta[n];

            alpha1Nei[facei] = alpha1[n] + t*(alpha1Nei[facei] - alpha1[n]);
            p1Nei[facei]     = p1[n]     + t*(p1Nei[facei]     - p1[n]);
            p2Nei[facei]     = p2[n]     + t*(p2Nei[facei]     - p2[n]);
            UNei[facei]      = U[n]      + t*(UNei[facei]      - U[n]);
            rho1Nei[facei]   = rho1[n]   + t*(rho1Nei[facei]   - rho1[n]);
            rho2Nei[facei]   = rho2[n]   + t*(rho2Nei[facei]   - rho2[n]);
        }
    }

    // --- Boundary patches ---
    forAll(alpha1.boundaryField(), patchi)
    {
        scalarField& bAlpha1Own = alpha1Own.boundaryFieldRef()[patchi];
        scalarField& bP1Own     = p1Own.boundaryFieldRef()[patchi];
        scalarField& bP2Own     = p2Own.boundaryFieldRef()[patchi];
        vectorField& bUOwn      = UOwn.boundaryFieldRef()[patchi];
        scalarField& bRho1Own   = rho1Own.boundaryFieldRef()[patchi];
        scalarField& bRho2Own   = rho2Own.boundaryFieldRef()[patchi];

        if (alpha1.boundaryField()[patchi].coupled())
        {
            // Owner side: apply the local cell's theta.  The nei side is
            // left alone here: Pass 3 replaces it with the remote
            // processor's own side, limited by the remote cell's theta
            const labelUList& faceCells =
                alpha1.boundaryField()[patchi].patch().faceCells();

            forAll(faceCells, facei)
            {
                const label c = faceCells[facei];

                if (theta[c] < 1.0)
                {
                    const scalar t = theta[c];

                    bAlpha1Own[facei] =
                        alpha1[c] + t*(bAlpha1Own[facei] - alpha1[c]);
                    bP1Own[facei] =
                        p1[c]     + t*(bP1Own[facei]     - p1[c]);
                    bP2Own[facei] =
                        p2[c]     + t*(bP2Own[facei]     - p2[c]);
                    bUOwn[facei] =
                        U[c]      + t*(bUOwn[facei]      - U[c]);
                    bRho1Own[facei] =
                        rho1[c]   + t*(bRho1Own[facei]   - rho1[c]);
                    bRho2Own[facei] =
                        rho2[c]   + t*(bRho2Own[facei]   - rho2[c]);
                }
            }
        }
        else
        {
            // Non-coupled patches: both sides hold the boundary-condition
            // state, which is not scaled
            forAll(bP1Own, facei)
            {
                bP1Own[facei] = max(bP1Own[facei], Foam::phaseFloors::p);
                bP2Own[facei] = max(bP2Own[facei], Foam::phaseFloors::p);
            }

            scalarField& bAlpha1Nei = alpha1Nei.boundaryFieldRef()[patchi];
            scalarField& bP1Nei     = p1Nei.boundaryFieldRef()[patchi];
            scalarField& bP2Nei     = p2Nei.boundaryFieldRef()[patchi];
            vectorField& bUNei      = UNei.boundaryFieldRef()[patchi];
            scalarField& bRho1Nei   = rho1Nei.boundaryFieldRef()[patchi];
            scalarField& bRho2Nei   = rho2Nei.boundaryFieldRef()[patchi];

            bAlpha1Nei = bAlpha1Own;
            bP1Nei     = bP1Own;
            bP2Nei     = bP2Own;
            bUNei      = bUOwn;
            bRho1Nei   = bRho1Own;
            bRho2Nei   = bRho2Own;
        }
    }


    // Pass 3: replace the nei side on coupled patches with the limited
    // own side from the remote processor, guaranteeing bit-identical
    // L/R states regardless of floating-point ordering.
    swapCoupledOwnToNei(alpha1Own, alpha1Nei);
    swapCoupledOwnToNei(p1Own,     p1Nei);
    swapCoupledOwnToNei(p2Own,     p2Nei);
    swapCoupledOwnToNei(UOwn,      UNei);
    swapCoupledOwnToNei(rho1Own,   rho1Nei);
    swapCoupledOwnToNei(rho2Own,   rho2Nei);


    // Vanished-phase safety clamp: where a phase's density constraint was
    // skipped, its reconstructed face density may sit below the floor.
    // Pin every face density to rhoFloor so the Riemann/EoS decode stays
    // evaluable
    const dimensionedScalar rhoFloor
    (
        rho1Own.dimensions(),
        Foam::phaseFloors::rho
    );

    rho1Own == max(rho1Own, rhoFloor);
    rho1Nei == max(rho1Nei, rhoFloor);
    rho2Own == max(rho2Own, rhoFloor);
    rho2Nei == max(rho2Nei, rhoFloor);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
