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

#include "phaseModel.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseModel::phaseModel
(
    const fvMesh& mesh,
    const word& phaseName,
    volScalarField& alpha,
    const volScalarField& p,
    const volScalarField& T
)
:
    IOdictionary
    (
        IOobject
        (
            IOobject::groupName("physicalProperties", phaseName),
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    name_(phaseName),
    eos_
    (
        phaseEquationOfState::New
        (
            phaseName,
            subDict("equationOfState")
        )
    ),
    transport_
    (
        new phaseTransport
        (
            phaseName,
            mesh,
            subDict("transport")
        )
    ),
    alpha_(alpha),

    alphaRho_
    (
        IOobject
        (
            IOobject::groupName("alphaRho", phaseName),
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("alphaRho", dimDensity, 0),
        "calculated"
    ),
    alphaE_
    (
        IOobject
        (
            IOobject::groupName("alphaE", phaseName),
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("alphaE", dimEnergy/dimVolume, 0),
        "calculated"
    ),

    rho_
    (
        IOobject
        (
            IOobject::groupName("rho", phaseName),
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("rho", dimDensity, 0),
        "calculated"
    ),

    c_
    (
        IOobject
        (
            IOobject::groupName("c", phaseName),
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("c", dimVelocity, 0),
        "calculated"
    ),

    rhoRead_(typeIOobject<volScalarField>(rho_).headerOk())
{
    if (rhoRead_)
    {
        Info<< "Reading phase density " << rho_.name() << endl;
    }
    else
    {
        Info<< "Initializing phase density " << rho_.name()
            << " from p and T" << endl;
    }

    // Seed the phase state at the mixture pressure;
    // compressibleTwoPhaseMixture::initialize() seeds it again at the phase
    // pressures once the capillary pressure is known
    initialize(p, T);

    // Initialize the temperature-dependent transport properties at the
    // mixture temperature (single-temperature closure)
    correctTransport(T);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseModel::~phaseModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::phaseModel::initialize
(
    const volScalarField& p,
    const volScalarField& T
)
{
    const phaseEquationOfState& eos = eos_();

    // Specific internal energy, used only to seed alphaE below
    volScalarField e
    (
        IOobject
        (
            IOobject::groupName("e", alpha_.group()),
            mesh_.time().name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedScalar("e", dimEnergy/dimMass, 0),
        "calculated"
    );

    // Internal field
    {
        const scalarField& pi = p.primitiveField();
        const scalarField& Ti = T.primitiveField();

        scalarField& rhoi = rho_.primitiveFieldRef();
        scalarField& ei   = e.primitiveFieldRef();
        scalarField& ci   = c_.primitiveFieldRef();

        forAll(pi, celli)
        {
            if (!rhoRead_)
            {
                rhoi[celli] = eos.rho(pi[celli], Ti[celli]);
            }

            ei[celli] = eos.e(pi[celli], rhoi[celli]);
            ci[celli] = eos.c(pi[celli], rhoi[celli]);
        }
    }

    // Boundary patches
    forAll(mesh_.boundary(), patchi)
    {
        const fvPatchScalarField& pPatch = p.boundaryField()[patchi];
        const fvPatchScalarField& TPatch = T.boundaryField()[patchi];

        fvPatchScalarField& rhoPatch = rho_.boundaryFieldRef()[patchi];
        fvPatchScalarField& ePatch   = e.boundaryFieldRef()[patchi];
        fvPatchScalarField& cPatch   = c_.boundaryFieldRef()[patchi];

        forAll(pPatch, facei)
        {
            if (!rhoRead_)
            {
                rhoPatch[facei] = eos.rho(pPatch[facei], TPatch[facei]);
            }

            ePatch[facei] = eos.e(pPatch[facei], rhoPatch[facei]);
            cPatch[facei] = eos.c(pPatch[facei], rhoPatch[facei]);
        }
    }

    // Update processor patches
    rho_.correctBoundaryConditions();
    e.correctBoundaryConditions();
    c_.correctBoundaryConditions();

    // Seed the RK-conserved fields from the EoS-consistent state
    alphaRho_ == alpha_*rho_;
    alphaE_   == alpha_*rho_*e;
}


void Foam::phaseModel::correct
(
    const volScalarField& p,
    const volScalarField& T
)
{
    const phaseEquationOfState& eos = eos_();

    // Internal field
    {
        const scalarField& pi        = p.primitiveField();
        const scalarField& Ti        = T.primitiveField();
        const scalarField& alphai    = alpha_.primitiveField();
        const scalarField& alphaRhoi = alphaRho_.primitiveField();

        scalarField& rhoi = rho_.primitiveFieldRef();
        scalarField& ci   = c_.primitiveFieldRef();

        forAll(pi, celli)
        {
            const scalar a_k  = alphai[celli];
            const scalar ar_k = alphaRhoi[celli];

            scalar rho_k = (a_k > Foam::phaseFloors::alpha)
                         ? ar_k/a_k
                         : eos.rho(pi[celli], Ti[celli]);
            rho_k = max(rho_k, Foam::phaseFloors::rho);

            const scalar c_k = eos.c(pi[celli], rho_k);

            rhoi[celli] = rho_k;
            ci[celli]   = c_k;
        }
    }

    // Update patches
    rho_.correctBoundaryConditions();
    c_.correctBoundaryConditions();
}


void Foam::phaseModel::correctBoundariesFromPrimitives
(
    const volScalarField& p,
    const volScalarField& T
)
{
    const phaseEquationOfState& eos = eos_();

    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& patch = mesh_.boundary()[patchi];

        if (patch.coupled())
        {
            continue;
        }

        const scalarField& pP     = p.boundaryField()[patchi];
        const scalarField& alphaP = alpha_.boundaryField()[patchi];

        scalarField& rhoP = rho_.boundaryFieldRef()[patchi];
        scalarField& cP   = c_.boundaryFieldRef()[patchi];
        scalarField& arP  = alphaRho_.boundaryFieldRef()[patchi];
        scalarField& aEP  = alphaE_.boundaryFieldRef()[patchi];

        if (T.boundaryField()[patchi].fixesValue())
        {
            // ---- Dirichlet T (inlet): reconstruct from (p, T) ----
            const scalarField& TP = T.boundaryField()[patchi];

            forAll(pP, facei)
            {
                const scalar rhof = eos.rho(pP[facei], TP[facei]);
                const scalar ef   = eos.e(pP[facei], rhof);
                const scalar cf   = eos.c(pP[facei], rhof);

                rhoP[facei] = rhof;
                cP[facei]   = cf;
                arP[facei]  = alphaP[facei]*rhof;
                aEP[facei]  = alphaP[facei]*rhof*ef;
            }
        }
        else
        {
            // ---- Non-Dirichlet T (outlet / zeroGradient / wave) ----
            const labelUList& faceCells =
                mesh_.boundary()[patchi].faceCells();

            forAll(pP, facei)
            {
                const scalar rhof = rho_.primitiveField()[faceCells[facei]];

                const scalar ef = eos.e(pP[facei], rhof);
                const scalar cf = eos.c(pP[facei], rhof);

                rhoP[facei] = rhof;
                cP[facei]   = cf;
                arP[facei]  = alphaP[facei]*rhof;
                aEP[facei]  = alphaP[facei]*rhof*ef;
            }
        }
    }
}


Foam::tmp<Foam::volScalarField>
Foam::phaseModel::nu() const
{
    return transport_().mu()/rho_;
}


Foam::tmp<Foam::scalarField>
Foam::phaseModel::nu(const label patchi) const
{
    return tmp<scalarField>
    (
        new scalarField
        (
            transport_().mu().boundaryField()[patchi]
          / rho_.boundaryField()[patchi]
        )
    );
}


void Foam::phaseModel::correctTransport(const volScalarField& T)
{
    transport_->correct(T);
}


// ************************************************************************* //
