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
    Single-pass flux-corrected-transport limiter (Zalesak 1979).

    Computes one blending factor theta_f in [0, 1] per face such that the
    partial densities, the phasic internal energies (including their
    pressure-work source) and the mixture total energy stay above a lower
    bound. Only lower bounds are enforced, so only Zalesak's R- is
    formed. No bound is imposed on alpha1; its fluxes are blended with
    the same theta.

    References:
    \verbatim
        Zalesak, S. T. (1979).
        Fully multidimensional flux-corrected transport algorithms for
        fluids.
        Journal of Computational Physics, 31(3), 335-362.
    \endverbatim

\*---------------------------------------------------------------------------*/

#include "shockVoF.H"
#include "fvcSurfaceIntegrate.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
namespace
{

//- Construct a dimensionless cell field initialized to 1, with default
//  (coupled-aware) boundaries so patchNeighbourField() returns the remote
//  value on processor patches.  Not registered
volScalarField unitField(const fvMesh& mesh, const word& name)
{
    return volScalarField
    (
        IOobject
        (
            name,
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedScalar(dimless, 1)
    );
}


//- Cell values of fvc::div(phi), without the processor halo exchange
//  of fvc::div: the limiter reads only cell values
tmp<scalarField> divCells(const fvMesh& mesh, const surfaceScalarField& phi)
{
    tmp<scalarField> tdiv(new scalarField(mesh.nCells(), scalar(0)));
    fvc::surfaceIntegrate(tdiv.ref(), phi);
    return tdiv;
}


//- Overflow-safe cell limiting fraction R- = clamp(Q-/P-, 0, 1) from the
//  antidiffusive sum Pm (<= 0) leaving each cell and the low-order updated
//  cell values U
void limiterFraction
(
    const scalarField& Pm,
    const scalarField& U,
    const scalar Umin,
    volScalarField& Rm
)
{
    scalarField& RmI = Rm.primitiveFieldRef();

    forAll(U, celli)
    {
        const scalar Qm = Umin - U[celli];   // headroom down (<= 0 in bounds)
        const scalar Pm_i = Pm[celli];       // <= 0

        // R- : limits the antidiffusion flowing OUT of the cell
        if (Pm_i >= -vSmall)     RmI[celli] = 1;   // no outgoing antidiffusion
        else if (Qm >= 0)        RmI[celli] = 0;   // at/below the lower bound
        else if (Qm <= Pm_i)     RmI[celli] = 1;   // room to give it all
        else                     RmI[celli] = Qm/Pm_i;         // 0 < r < 1
    }

    Rm.correctBoundaryConditions();
}


//- Zalesak FCT limiter for a lower bound.  Reduces theta (in place, taking
//  the minimum) so every cell average of U stays at or above Umin under the
//  blended update U = lowOrderU - dt/V*sum theta*A, dt the length of the
//  forward-Euler sub-stage
void zalesakLimit
(
    const fvMesh& mesh,
    surfaceScalarField& theta,
    const scalar dt,
    const scalarField& lowOrderU,
    const surfaceScalarField& A,
    const scalar Umin,
    const volScalarField* S = nullptr,
    const surfaceScalarField* Av = nullptr
)
{
    const labelUList& owner = mesh.owner();
    const labelUList& neighbor = mesh.neighbour();
    const scalarField& V = mesh.V();

    // Non-conservative source cell values, or null for a conservative field
    const scalarField* Sc = S ? &S->primitiveField() : nullptr;

    // Step 1: sum the maximum admissible (theta = 1) antidiffusive increment
    // out of each cell (P-)
    scalarField Pm(mesh.nCells(), scalar(0));

    forAll(owner, facei)
    {
        const label o = owner[facei];
        const label n = neighbor[facei];

        scalar aOwn = A[facei];
        scalar aNei = A[facei];
        if (Sc)
        {
            aOwn -= (*Sc)[o]*(*Av)[facei];
            aNei -= (*Sc)[n]*(*Av)[facei];
        }

        const scalar dOwn = -dt/V[o]*aOwn;
        if (dOwn < 0) Pm[o] += dOwn;

        const scalar dNei = dt/V[n]*aNei;
        if (dNei < 0) Pm[n] += dNei;
    }

    forAll(mesh.boundary(), patchi)
    {
        const fvsPatchScalarField& pA = A.boundaryField()[patchi];
        const labelUList& fc = mesh.boundary()[patchi].faceCells();

        forAll(fc, i)
        {
            const label o = fc[i];
            scalar aOwn = pA[i];
            if (Sc)
            {
                aOwn -= (*Sc)[o]*Av->boundaryField()[patchi][i];
            }
            const scalar dOwn = -dt/V[o]*aOwn;
            if (dOwn < 0) Pm[o] += dOwn;
        }
    }

    // Step 2: cell limiting fraction R-
    volScalarField Rm(unitField(mesh, "fct::Rm"));
    limiterFraction(Pm, lowOrderU, Umin, Rm);

    const scalarField& RmI = Rm.primitiveField();

    // Step 3: a face is limited by the R- of the cell its antidiffusion
    // leaves
    forAll(owner, facei)
    {
        const label o = owner[facei];
        const label n = neighbor[facei];

        scalar aOwn = A[facei];
        scalar aNei = A[facei];
        if (Sc)
        {
            aOwn -= (*Sc)[o]*(*Av)[facei];
            aNei -= (*Sc)[n]*(*Av)[facei];
        }

        scalar tfOwn = 1;
        if (aOwn > vSmall)       tfOwn = RmI[o];   // owner loses -> R-

        scalar tfNei = 1;
        if (aNei < -vSmall)      tfNei = RmI[n];   // neighbor loses -> R-

        theta[facei] = min(theta[facei], min(tfOwn, tfNei));
    }

    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& patch = mesh.boundary()[patchi];
        const fvsPatchScalarField& pA = A.boundaryField()[patchi];
        fvsPatchScalarField& pTheta = theta.boundaryFieldRef()[patchi];
        const labelUList& fc = patch.faceCells();

        if (patch.coupled())
        {
            const scalarField RmNei
            (
                Rm.boundaryField()[patchi].patchNeighbourField()
            );
            scalarField SNei;
            if (Sc)
            {
                SNei = S->boundaryField()[patchi].patchNeighbourField();
            }

            forAll(fc, i)
            {
                const label o = fc[i];

                scalar aOwn = pA[i];
                scalar aNei = pA[i];
                if (Sc)
                {
                    const scalar av = Av->boundaryField()[patchi][i];
                    aOwn -= (*Sc)[o]*av;
                    aNei -= SNei[i]*av;
                }

                scalar tfOwn = 1;
                if (aOwn > vSmall)       tfOwn = RmI[o];

                scalar tfNei = 1;
                if (aNei < -vSmall)      tfNei = RmNei[i];

                pTheta[i] = min(pTheta[i], min(tfOwn, tfNei));
            }
        }
        else
        {
            forAll(fc, i)
            {
                const label o = fc[i];

                scalar aOwn = pA[i];
                if (Sc)
                {
                    aOwn -= (*Sc)[o]*Av->boundaryField()[patchi][i];
                }

                scalar tfOwn = 1;
                if (aOwn > vSmall)       tfOwn = RmI[o];

                pTheta[i] = min(pTheta[i], tfOwn);
            }
        }
    }
}

} // End anonymous namespace
} // End namespace Foam


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField>
Foam::solvers::shockVoF::fluxLimiter(const dimensionedScalar& deltaT)
{
    tmp<surfaceScalarField> tTheta
    (
        surfaceScalarField::New
        (
            "fctTheta_f",
            mesh,
            dimensionedScalar(dimless, 1)
        )
    );
    surfaceScalarField& theta = tTheta.ref();

    // The bounds are enforced only when the flow is solved; with the flow
    // frozen only alpha1 is advected (see solveEquations)
    if (pimple.flow())
    {
        // Length of the forward-Euler sub-stage being limited
        const dimensionedScalar& dt = deltaT;

        // Relative admissibility tolerance
        const scalar tol = 1e-6;

        // Global maxima of the five constrained fields
        scalarList maxU
        ({
            max(alpha1Rho1_.primitiveField()),
            max(alpha2Rho2_.primitiveField()),
            max(alpha1RhoE1_.primitiveField()),
            max(alpha2RhoE2_.primitiveField()),
            max(rhoE_.primitiveField())
        });
        Pstream::listCombineGather(maxU, maxEqOp<scalar>());
        Pstream::listCombineScatter(maxU);

        const scalar maxAlpha1Rho1 = maxU[0];
        const scalar maxAlpha2Rho2 = maxU[1];
        const scalar maxAlpha1RhoE1 = maxU[2];
        const scalar maxAlpha2RhoE2 = maxU[3];
        const scalar maxRhoE = maxU[4];

        // Low-order volumetric-flux divergence and its antidiffusive part
        const scalarField lowOrderDivPhiv(divCells(mesh, flux_.phivLO()));
        const surfaceScalarField Av(flux_.phiv() - flux_.phivLO());

        // --- Partial densities: positivity, conservative ---
        {
            const scalarField lowOrderAlpha1Rho1
            (
                alpha1Rho1_.primitiveField()
              - dt.value()*divCells(mesh, flux_.phiAlpha1Rho1LO())
            );
            const surfaceScalarField A
            (
                flux_.phiAlpha1Rho1() - flux_.phiAlpha1Rho1LO()
            );
            const scalar Umin = -(small + tol*max(maxAlpha1Rho1, small));
            zalesakLimit
            (
                mesh, theta, dt.value(), lowOrderAlpha1Rho1, A, Umin
            );
        }
        {
            const scalarField lowOrderAlpha2Rho2
            (
                alpha2Rho2_.primitiveField()
              - dt.value()*divCells(mesh, flux_.phiAlpha2Rho2LO())
            );
            const surfaceScalarField A
            (
                flux_.phiAlpha2Rho2() - flux_.phiAlpha2Rho2LO()
            );
            const scalar Umin = -(small + tol*max(maxAlpha2Rho2, small));
            zalesakLimit
            (
                mesh, theta, dt.value(), lowOrderAlpha2Rho2, A, Umin
            );
        }

        // --- Phasic internal energies: positivity, non-conservative ---
        // d(alphaRhoEk) = -div(phiAlphaRhoEk) - alpha_k*p_k*div(phiv) with
        // the phase pressures p1 = p + alpha2*pc and p2 = p - alpha1*pc,
        // so S = -alpha_k*p_k
        const volScalarField alpha1P1(alpha1_*mixture_.p1());
        const volScalarField alpha2P2(alpha2_*mixture_.p2());
        {
            const scalarField lowOrderAlpha1RhoE1
            (
                alpha1RhoE1_.primitiveField()
              - dt.value()
               *(
                    divCells(mesh, flux_.phiAlpha1RhoE1LO())
                  + alpha1P1.primitiveField()*lowOrderDivPhiv
                )
            );
            const surfaceScalarField A
            (
                flux_.phiAlpha1RhoE1() - flux_.phiAlpha1RhoE1LO()
            );

            // Source multiplier; zalesakLimit() reads its coupled patch
            // values as the neighbor cells' values
            const volScalarField S1("fct::S1", -alpha1P1);
            const scalar Umin = -(small + tol*max(maxAlpha1RhoE1, small));

            zalesakLimit
            (
                mesh, theta, dt.value(), lowOrderAlpha1RhoE1, A, Umin,
                &S1, &Av
            );
        }
        {
            const scalarField lowOrderAlpha2RhoE2
            (
                alpha2RhoE2_.primitiveField()
              - dt.value()
               *(
                    divCells(mesh, flux_.phiAlpha2RhoE2LO())
                  + alpha2P2.primitiveField()*lowOrderDivPhiv
                )
            );
            const surfaceScalarField A
            (
                flux_.phiAlpha2RhoE2() - flux_.phiAlpha2RhoE2LO()
            );
            const volScalarField S2("fct::S2", -alpha2P2);
            const scalar Umin = -(small + tol*max(maxAlpha2RhoE2, small));

            zalesakLimit
            (
                mesh, theta, dt.value(), lowOrderAlpha2RhoE2, A, Umin,
                &S2, &Av
            );
        }

        // --- Mixture total energy: positivity, conservative flux ---
        {
            const scalarField lowOrderRhoE
            (
                rhoE_.primitiveField()
              - dt.value()*divCells(mesh, flux_.phiRhoELO())
            );
            const surfaceScalarField A(flux_.phiRhoE() - flux_.phiRhoELO());
            const scalar Umin = -(small + tol*max(maxRhoE, small));
            zalesakLimit
            (
                mesh, theta, dt.value(), lowOrderRhoE, A, Umin
            );
        }
    }

    return tTheta;
}


// ************************************************************************* //
