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

#include "fluxPredictor.H"
#include "fvMesh.H"
#include "surfaceFields.H"
#include "fvsPatchFields.H"
#include "IStringStream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fluxPredictor, 0);
}

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
namespace
{

//- Blend flux in place: flux = fluxLO + theta*(flux - fluxLO)
template<class Type>
void blendFlux
(
    SurfaceField<Type>& flux,
    const SurfaceField<Type>& fluxLO,
    const surfaceScalarField& theta
)
{
    Field<Type>& f = flux.primitiveFieldRef();
    const Field<Type>& fLO = fluxLO.primitiveField();
    const scalarField& t = theta.primitiveField();

    forAll(f, facei)
    {
        f[facei] = fLO[facei] + t[facei]*(f[facei] - fLO[facei]);
    }

    typename SurfaceField<Type>::Boundary& fBf = flux.boundaryFieldRef();

    forAll(fBf, patchi)
    {
        fvsPatchField<Type>& pf = fBf[patchi];
        const fvsPatchField<Type>& pfLO = fluxLO.boundaryField()[patchi];
        const fvsPatchScalarField& pt = theta.boundaryField()[patchi];

        forAll(pf, facei)
        {
            pf[facei] = pfLO[facei] + pt[facei]*(pf[facei] - pfLO[facei]);
        }
    }
}

} // End anonymous namespace
} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluxPredictor::fluxPredictor
(
    const volScalarField& p,
    const volVectorField& U,
    compressibleTwoPhaseMixture& mixture,
    const interfaceModel& interface
)
:
    mesh_(p.mesh()),
    mixture_(mixture),
    interface_(interface),
    alpha1_(mixture_.phase1().alpha()),
    U_(U),
    rho1_(mixture_.phase1().rho()),
    rho2_(mixture_.phase2().rho()),

    reconstruction_
    (
        reconstructionScheme::New
        (
            mesh_,
            mesh_.schemes().subDict("interpolationSchemes")
        )
    ),

    reconstructionLowOrder_
    (
        reconstructionScheme::New
        (
            mesh_,
            dictionary
            (
                IStringStream
                (
                    "reconstruction"
                    "{"
                    "    type                oneDim;"
                    "    reconstruct(alpha)  upwind;"
                    "    reconstruct(p)      upwind;"
                    "    reconstruct(U)      upwind;"
                    "    reconstruct(rho)    upwind;"
                    "    positivityLimiter   no;"
                    "}"
                )()
            )
        )
    ),

    fluxScheme_
    (
        fluxScheme::New
        (
            mesh_,
            mesh_.schemes().subDict("fluxScheme")
        )
    ),

    fluxSchemeLowOrder_
    (
        fluxScheme::New
        (
            "HLL",
            mesh_,
            mesh_.schemes().subDict("fluxScheme")
        )
    ),

    amaxSf_
    (
        IOobject
        (
            "amaxSf",
            mesh_.time().name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("amaxSf", dimVelocity*dimArea, Zero)
    ),

    phi_
    (
        IOobject
        (
            "phi",
            mesh_.time().name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("phi", dimMass/dimTime, Zero)
    ),
    phiUp_
    (
        IOobject
        (
            "phiUp",
            mesh_.time().name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("phiUp", dimMass*dimVelocity/dimTime, Zero)
    ),
    phiAlpha1Rho1_
    (
        IOobject
        (
            IOobject::groupName("phiAlpha1Rho1", alpha1_.group()),
            mesh_.time().name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("phiAlpha1Rho1", dimMass/dimTime, Zero)
    ),
    phiAlpha2Rho2_
    (
        IOobject
        (
            IOobject::groupName
            (
                "phiAlpha2Rho2",
                mixture_.phase2().alpha().group()
            ),
            mesh_.time().name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("phiAlpha2Rho2", dimMass/dimTime, Zero)
    ),
    phiAlpha1_
    (
        IOobject
        (
            "phiAlpha1",
            mesh_.time().name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("phiAlpha1", dimVolume/dimTime, Zero)
    ),
    phiAlpha1RhoE1_
    (
        IOobject
        (
            "phiAlpha1RhoE1",
            mesh_.time().name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("phiAlpha1RhoE1", dimEnergy/dimTime, Zero)
    ),
    phiAlpha2RhoE2_
    (
        IOobject
        (
            "phiAlpha2RhoE2",
            mesh_.time().name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("phiAlpha2RhoE2", dimEnergy/dimTime, Zero)
    ),
    phiRhoE_
    (
        IOobject
        (
            "phiRhoE",
            mesh_.time().name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("phiRhoE", dimEnergy/dimTime, Zero)
    ),
    phiv_
    (
        IOobject
        (
            "phiv",
            mesh_.time().name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("phiv", dimVolume/dimTime, Zero)
    ),
    alphaf_
    (
        IOobject
        (
            "alphaf",
            mesh_.time().name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphaf", dimless, Zero)
    ),

    phiLO_
    (
        IOobject
        (
            "phiLO",
            mesh_.time().name(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("phiLO", dimMass/dimTime, Zero)
    ),
    phiUpLO_
    (
        IOobject
        (
            "phiUpLO",
            mesh_.time().name(),
            mesh_
        ),
        mesh_,
        dimensionedVector("phiUpLO", dimMass*dimVelocity/dimTime, Zero)
    ),
    phiAlpha1Rho1LO_
    (
        IOobject
        (
            "phiAlpha1Rho1LO",
            mesh_.time().name(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("phiAlpha1Rho1LO", dimMass/dimTime, Zero)
    ),
    phiAlpha2Rho2LO_
    (
        IOobject
        (
            "phiAlpha2Rho2LO",
            mesh_.time().name(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("phiAlpha2Rho2LO", dimMass/dimTime, Zero)
    ),
    phiAlpha1LO_
    (
        IOobject
        (
            "phiAlpha1LO",
            mesh_.time().name(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("phiAlpha1LO", dimVolume/dimTime, Zero)
    ),
    phiAlpha1RhoE1LO_
    (
        IOobject
        (
            "phiAlpha1RhoE1LO",
            mesh_.time().name(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("phiAlpha1RhoE1LO", dimEnergy/dimTime, Zero)
    ),
    phiAlpha2RhoE2LO_
    (
        IOobject
        (
            "phiAlpha2RhoE2LO",
            mesh_.time().name(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("phiAlpha2RhoE2LO", dimEnergy/dimTime, Zero)
    ),
    phiRhoELO_
    (
        IOobject
        (
            "phiRhoELO",
            mesh_.time().name(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("phiRhoELO", dimEnergy/dimTime, Zero)
    ),
    phivLO_
    (
        IOobject
        (
            "phivLO",
            mesh_.time().name(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("phivLO", dimVolume/dimTime, Zero)
    ),
    alphafLO_
    (
        IOobject
        (
            "alphafLO",
            mesh_.time().name(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("alphafLO", dimless, Zero)
    ),
    curTimeIndex_(-1)
{
    if (debug)
    {
        fluxLimiterPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "fluxLimiter",
                    mesh_.time().name(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar(dimless, 1),
                "zeroGradient"
            )
        );
    }

    // Both flux sets from the start, so that the boundary conditions that
    // read the fluxes (waveTransmissive, inletOutlet, ...) can be evaluated
    // before the first time step.
    update();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fluxPredictor::update()
{
    // Phase pressures, reconstructed in place of the mixture pressure
    const volScalarField& p1 = mixture_.p1();
    const volScalarField& p2 = mixture_.p2();

    // Owner/neighbor face states, filled by each reconstruction in turn
    tmp<surfaceScalarField> alpha1_pos, alpha1_neg;
    tmp<surfaceScalarField> p1_pos, p1_neg;
    tmp<surfaceScalarField> p2_pos, p2_neg;
    tmp<surfaceVectorField> U_pos, U_neg;
    tmp<surfaceScalarField> rho1_pos, rho1_neg;
    tmp<surfaceScalarField> rho2_pos, rho2_neg;

    // Both flux sets raise the signal speed from zero
    amaxSf_ = dimensionedScalar(amaxSf_.dimensions(), 0);

    // --- High-order fluxes: face states from the runtime-selected
    //     reconstruction through the selected flux scheme.  Interface
    //     sharpening (rho-THINC), the alpha clamp and positivity limiting
    //     are all applied inside reconstruct(). ---
    reconstruction_->reconstruct
    (
        alpha1_, p1, p2, U_, rho1_, rho2_,
        interface_.normal(),
        alpha1_pos, alpha1_neg,
        p1_pos,     p1_neg,
        p2_pos,     p2_neg,
        U_pos,      U_neg,
        rho1_pos,   rho1_neg,
        rho2_pos,   rho2_neg
    );

    evaluateFluxes
    (
        *fluxScheme_,
        alpha1_pos(), alpha1_neg(),
        p1_pos(),     p1_neg(),
        p2_pos(),     p2_neg(),
        U_pos(),      U_neg(),
        rho1_pos(),   rho1_neg(),
        rho2_pos(),   rho2_neg(),
        phi_,
        phiUp_,
        phiAlpha1Rho1_,
        phiAlpha2Rho2_,
        phiAlpha1_,
        phiAlpha1RhoE1_,
        phiAlpha2RhoE2_,
        phiRhoE_,
        phiv_,
        alphaf_
    );

    // --- Low-order fluxes: the cell-center states, from the upwind
    //     reconstruction, through HLL (constructed by type name, not
    //     user-selectable).  Reconstructing into the same tmps releases
    //     the high-order states. ---
    reconstructionLowOrder_->reconstruct
    (
        alpha1_, p1, p2, U_, rho1_, rho2_,
        interface_.normal(),
        alpha1_pos, alpha1_neg,
        p1_pos,     p1_neg,
        p2_pos,     p2_neg,
        U_pos,      U_neg,
        rho1_pos,   rho1_neg,
        rho2_pos,   rho2_neg
    );

    evaluateFluxes
    (
        *fluxSchemeLowOrder_,
        alpha1_pos(), alpha1_neg(),
        p1_pos(),     p1_neg(),
        p2_pos(),     p2_neg(),
        U_pos(),      U_neg(),
        rho1_pos(),   rho1_neg(),
        rho2_pos(),   rho2_neg(),
        phiLO_,
        phiUpLO_,
        phiAlpha1Rho1LO_,
        phiAlpha2Rho2LO_,
        phiAlpha1LO_,
        phiAlpha1RhoE1LO_,
        phiAlpha2RhoE2LO_,
        phiRhoELO_,
        phivLO_,
        alphafLO_
    );
}


void Foam::fluxPredictor::evaluateFluxes
(
    const fluxScheme& scheme,
    const surfaceScalarField& alpha1_pos,
    const surfaceScalarField& alpha1_neg,
    const surfaceScalarField& p1_pos,
    const surfaceScalarField& p1_neg,
    const surfaceScalarField& p2_pos,
    const surfaceScalarField& p2_neg,
    const surfaceVectorField& U_pos,
    const surfaceVectorField& U_neg,
    const surfaceScalarField& rho1_pos,
    const surfaceScalarField& rho1_neg,
    const surfaceScalarField& rho2_pos,
    const surfaceScalarField& rho2_neg,
    surfaceScalarField& phi,
    surfaceVectorField& phiUp,
    surfaceScalarField& phiAlpha1Rho1,
    surfaceScalarField& phiAlpha2Rho2,
    surfaceScalarField& phiAlpha1,
    surfaceScalarField& phiAlpha1RhoE1,
    surfaceScalarField& phiAlpha2RhoE2,
    surfaceScalarField& phiRhoE,
    surfaceScalarField& phiv,
    surfaceScalarField& alphaf
)
{
    const labelUList& owner = mesh_.owner();

    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();

    // Face capillary coefficient sigma*kappa_f
    const surfaceScalarField& sigmaKappaf = interface_.sigmaKappaf();

    // Face mixture pressure alpha1*p1 + alpha2*p2
    const surfaceScalarField p_pos
    (
        "p_pos",
        p2_pos + alpha1_pos*(p1_pos - p2_pos)
    );
    const surfaceScalarField p_neg
    (
        "p_neg",
        p2_neg + alpha1_neg*(p1_neg - p2_neg)
    );

    // Signal speed of each Riemann problem
    scalar amaxSfi = 0;

    // --- Internal faces ---
    forAll(owner, facei)
    {
        // Left state
        scalar rhoL, cL, alpha1Rho1L, alpha2Rho2L,
               alphaRhoE1L, alphaRhoE2L;
        mixture_.primitivesToState
        (
            p1_pos[facei],
            p2_pos[facei],
            alpha1_pos[facei],
            rho1_pos[facei],
            rho2_pos[facei],
            rhoL,
            cL,
            alpha1Rho1L,
            alpha2Rho2L,
            alphaRhoE1L,
            alphaRhoE2L
        );

        // Right state
        scalar rhoR, cR, alpha1Rho1R, alpha2Rho2R,
               alphaRhoE1R, alphaRhoE2R;
        mixture_.primitivesToState
        (
            p1_neg[facei],
            p2_neg[facei],
            alpha1_neg[facei],
            rho1_neg[facei],
            rho2_neg[facei],
            rhoR,
            cR,
            alpha1Rho1R,
            alpha2Rho2R,
            alphaRhoE1R,
            alphaRhoE2R
        );

        scheme.evaluateFlux
        (
            phi[facei],
            phiUp[facei],
            phiAlpha1Rho1[facei],
            phiAlpha2Rho2[facei],
            phiAlpha1[facei],
            phiAlpha1RhoE1[facei],
            phiAlpha2RhoE2[facei],
            phiRhoE[facei],
            phiv[facei],
            amaxSfi,

            alphaf[facei],

            p_pos[facei],
            U_pos[facei],
            rhoL,
            cL,
            alpha1_pos[facei],
            alpha1Rho1L,
            alpha2Rho2L,
            alphaRhoE1L,
            alphaRhoE2L,

            p_neg[facei],
            U_neg[facei],
            rhoR,
            cR,
            alpha1_neg[facei],
            alpha1Rho1R,
            alpha2Rho2R,
            alphaRhoE1R,
            alphaRhoE2R,

            Sf[facei],
            magSf[facei],

            sigmaKappaf[facei]
        );

        amaxSf_[facei] = max(amaxSf_[facei], amaxSfi);
    }


    // --- Boundary faces (coupled AND non-coupled) ---
    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& patch = mesh_.boundary()[patchi];

        // Output flux patch references
        fvsPatchScalarField& pPhi = phi.boundaryFieldRef()[patchi];
        fvsPatchVectorField& pPhiUp = phiUp.boundaryFieldRef()[patchi];
        fvsPatchScalarField& pPhiA1R1 =
            phiAlpha1Rho1.boundaryFieldRef()[patchi];
        fvsPatchScalarField& pPhiA2R2 =
            phiAlpha2Rho2.boundaryFieldRef()[patchi];
        fvsPatchScalarField& pPhiAlpha1 = phiAlpha1.boundaryFieldRef()[patchi];
        fvsPatchScalarField& pPhiA1RE1 =
            phiAlpha1RhoE1.boundaryFieldRef()[patchi];
        fvsPatchScalarField& pPhiA2RE2 =
            phiAlpha2RhoE2.boundaryFieldRef()[patchi];
        fvsPatchScalarField& pPhiRhoE = phiRhoE.boundaryFieldRef()[patchi];
        fvsPatchScalarField& pPhiv = phiv.boundaryFieldRef()[patchi];
        fvsPatchScalarField& pAmaxSf = amaxSf_.boundaryFieldRef()[patchi];
        fvsPatchScalarField& pAlphaf = alphaf.boundaryFieldRef()[patchi];

        const fvsPatchVectorField& pSf    = Sf.boundaryField()[patchi];
        const fvsPatchScalarField& pMagSf = magSf.boundaryField()[patchi];

        const fvsPatchScalarField& pp_pos = p_pos.boundaryField()[patchi];
        const fvsPatchScalarField& pp1_pos = p1_pos.boundaryField()[patchi];
        const fvsPatchScalarField& pp2_pos = p2_pos.boundaryField()[patchi];
        const fvsPatchVectorField& pU_pos = U_pos.boundaryField()[patchi];
        const fvsPatchScalarField& pRho1_pos =
            rho1_pos.boundaryField()[patchi];
        const fvsPatchScalarField& pRho2_pos =
            rho2_pos.boundaryField()[patchi];
        const fvsPatchScalarField& pAlpha1_pos =
            alpha1_pos.boundaryField()[patchi];

        const fvsPatchScalarField& pp_neg = p_neg.boundaryField()[patchi];
        const fvsPatchScalarField& pp1_neg = p1_neg.boundaryField()[patchi];
        const fvsPatchScalarField& pp2_neg = p2_neg.boundaryField()[patchi];
        const fvsPatchVectorField& pU_neg = U_neg.boundaryField()[patchi];
        const fvsPatchScalarField& pRho1_neg =
            rho1_neg.boundaryField()[patchi];
        const fvsPatchScalarField& pRho2_neg =
            rho2_neg.boundaryField()[patchi];
        const fvsPatchScalarField& pAlpha1_neg =
            alpha1_neg.boundaryField()[patchi];
        const fvsPatchScalarField& pSigmaKappaf =
            sigmaKappaf.boundaryField()[patchi];

        forAll(patch, facei)
        {
            scalar rhoL, cL, alpha1Rho1L, alpha2Rho2L,
                   alphaRhoE1L, alphaRhoE2L;
            mixture_.primitivesToState
            (
                pp1_pos[facei],
                pp2_pos[facei],
                pAlpha1_pos[facei],
                pRho1_pos[facei],
                pRho2_pos[facei],
                rhoL,
                cL,
                alpha1Rho1L,
                alpha2Rho2L,
                alphaRhoE1L,
                alphaRhoE2L
            );

            scalar rhoR, cR, alpha1Rho1R, alpha2Rho2R,
                   alphaRhoE1R, alphaRhoE2R;
            mixture_.primitivesToState
            (
                pp1_neg[facei],
                pp2_neg[facei],
                pAlpha1_neg[facei],
                pRho1_neg[facei],
                pRho2_neg[facei],
                rhoR,
                cR,
                alpha1Rho1R,
                alpha2Rho2R,
                alphaRhoE1R,
                alphaRhoE2R
            );

            scheme.evaluateFlux
            (
                pPhi[facei],
                pPhiUp[facei],
                pPhiA1R1[facei],
                pPhiA2R2[facei],
                pPhiAlpha1[facei],
                pPhiA1RE1[facei],
                pPhiA2RE2[facei],
                pPhiRhoE[facei],
                pPhiv[facei],
                amaxSfi,

                pAlphaf[facei],

                pp_pos[facei],
                pU_pos[facei],
                rhoL,
                cL,
                pAlpha1_pos[facei],
                alpha1Rho1L,
                alpha2Rho2L,
                alphaRhoE1L,
                alphaRhoE2L,

                pp_neg[facei],
                pU_neg[facei],
                rhoR,
                cR,
                pAlpha1_neg[facei],
                alpha1Rho1R,
                alpha2Rho2R,
                alphaRhoE1R,
                alphaRhoE2R,

                pSf[facei],
                pMagSf[facei],

                pSigmaKappaf[facei]
            );

            pAmaxSf[facei] = max(pAmaxSf[facei], amaxSfi);
        }
    }
}


void Foam::fluxPredictor::blend(const surfaceScalarField& theta)
{
    // Blend of each transport flux toward its low-order value
    blendFlux(phi_, phiLO_, theta);
    blendFlux(phiUp_, phiUpLO_, theta);
    blendFlux(phiAlpha1Rho1_, phiAlpha1Rho1LO_, theta);
    blendFlux(phiAlpha2Rho2_, phiAlpha2Rho2LO_, theta);
    blendFlux(phiAlpha1_, phiAlpha1LO_, theta);
    blendFlux(phiAlpha1RhoE1_, phiAlpha1RhoE1LO_, theta);
    blendFlux(phiAlpha2RhoE2_, phiAlpha2RhoE2LO_, theta);
    blendFlux(phiRhoE_, phiRhoELO_, theta);
    blendFlux(phiv_, phivLO_, theta);
    blendFlux(alphaf_, alphafLO_, theta);

    // Write diagnostics field
    if (debug)
    {
        volScalarField& fluxLimiter = fluxLimiterPtr_();

        // Reset diagnostics field at the beginning of time step
        if (curTimeIndex_ != mesh_.time().timeIndex())
        {
            fluxLimiter == dimensionedScalar(dimless, 1);

            curTimeIndex_ = mesh_.time().timeIndex();
        }

        const labelUList& owner = mesh_.owner();
        const labelUList& neighbor = mesh_.neighbour();

        forAll(owner, facei)
        {
            fluxLimiter[owner[facei]] = min
            (
                fluxLimiter[owner[facei]],
                theta[facei]
            );
            fluxLimiter[neighbor[facei]] = min
            (
                fluxLimiter[neighbor[facei]],
                theta[facei]
            );
        }

        forAll(mesh_.boundary(), patchi)
        {
            const fvsPatchScalarField& pTheta = theta.boundaryField()[patchi];
            const labelUList& faceCells = mesh_.boundary()[patchi].faceCells();

            forAll(faceCells, i)
            {
                fluxLimiter[faceCells[i]] = min
                (
                    fluxLimiter[faceCells[i]],
                    pTheta[i]
                );
            }
        }

        fluxLimiter.correctBoundaryConditions();
    }
}


// ************************************************************************* //
