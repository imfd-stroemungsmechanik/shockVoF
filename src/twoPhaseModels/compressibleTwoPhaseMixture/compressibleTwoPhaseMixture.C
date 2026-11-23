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

#include "compressibleTwoPhaseMixture.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::compressibleTwoPhaseMixture::updateMixtureTransportAndAcoustics()
{
    const volScalarField& alpha1 = twoPhaseMixture::alpha1();
    const volScalarField& alpha2 = twoPhaseMixture::alpha2();

    // Volumetric mixing rule for the viscosity.  The thermal conductivity
    // is mixed per phase with the turbulent contribution in
    // twoPhaseThermophysicalTransportModel::kappaEff(), so no mixture
    // value is formed here
    const dimensionedScalar rhoSmall("rhoSmall", dimDensity, small);

    mu_ = alpha1*phase1_.mu() + alpha2*phase2_.mu();
    nu_ = mu_/max(rho_, rhoSmall);

    // Frozen mixture sound speed
    const dimensionedScalar cSqrSmall("cSqrSmall", sqr(dimVelocity), small);

    const volScalarField rhoCSqrMix
    (
        alpha1*phase1_.rho()*sqr(phase1_.c())
      + alpha2*phase2_.rho()*sqr(phase2_.c())
    );

    c_ = sqrt(max(rhoCSqrMix/max(rho_, rhoSmall), cSqrSmall));

    psi_ = 1.0/sqr(c_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::compressibleTwoPhaseMixture::compressibleTwoPhaseMixture
(
    const fvMesh& mesh
)
:
    // twoPhaseMixture is the phaseProperties IOdictionary: it reads the
    // phase names and alpha1 and sets alpha2 = 1 - alpha1
    twoPhaseMixture(mesh),
    viscosity(),

    p_
    (
        IOobject
        (
            "p",
            mesh.time().name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    T_
    (
        IOobject
        (
            "T",
            mesh.time().name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    pc_
    (
        IOobject
        (
            "pc",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimPressure, 0),
        "zeroGradient"
    ),

    p1_
    (
        IOobject
        (
            IOobject::groupName("p", twoPhaseMixture::phase1Name()),
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        p_ + twoPhaseMixture::alpha2()*pc_
    ),

    p2_
    (
        IOobject
        (
            IOobject::groupName("p", twoPhaseMixture::phase2Name()),
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        p_ - twoPhaseMixture::alpha1()*pc_
    ),

    phase1_
    (
        mesh,
        twoPhaseMixture::phase1Name(),
        twoPhaseMixture::alpha1(),
        p_,
        T_
    ),
    phase2_
    (
        mesh,
        twoPhaseMixture::phase2Name(),
        twoPhaseMixture::alpha2(),
        p_,
        T_
    ),

    sgEoS1_(refCast<const stiffenedGasEoS>(phase1_.eos())),
    sgEoS2_(refCast<const stiffenedGasEoS>(phase2_.eos())),

    rho_
    (
        IOobject
        (
            "rho",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("rho", dimDensity, 0),
        "calculated"
    ),

    mu_
    (
        IOobject
        (
            "mu",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("mu", dimDynamicViscosity, 0),
        "calculated"
    ),

    nu_
    (
        IOobject
        (
            "nu",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("nu", dimKinematicViscosity, 0),
        "calculated"
    ),
    c_
    (
        IOobject
        (
            "c",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("c", dimVelocity, 0),
        "calculated"
    ),
    psi_
    (
        IOobject
        (
            "psi",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("psi", dimless/sqr(dimVelocity), 0),
        "calculated"
    ),

    relaxation_(sgEoS1_, sgEoS2_)
{
    correctRho();

    updateMixtureTransportAndAcoustics();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::compressibleTwoPhaseMixture::~compressibleTwoPhaseMixture()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::compressibleTwoPhaseMixture::correctRho()
{
    rho_ = max
    (
        phase1_.alphaRho() + phase2_.alphaRho(),
        dimensionedScalar(dimDensity, phaseFloors::rho)
    );
    rho_.correctBoundaryConditions();
}


void Foam::compressibleTwoPhaseMixture::correctPhasePressures()
{
    const volScalarField& alpha1 = twoPhaseMixture::alpha1();
    volScalarField& alpha2 = twoPhaseMixture::alpha2();

    // alpha2 = 1 - alpha1 on the non-coupled patches: alpha2 is calculated
    // there, so evaluating its boundary conditions does not follow alpha1's
    forAll(alpha2.boundaryField(), patchi)
    {
        if (!alpha2.boundaryField()[patchi].coupled())
        {
            alpha2.boundaryFieldRef()[patchi] ==
                1.0 - alpha1.boundaryField()[patchi];
        }
    }

    p1_ == p_ + alpha2*pc_;
    p2_ == p_ - alpha1*pc_;
}


void Foam::compressibleTwoPhaseMixture::initialize
(
    const volScalarField& sigmaKappa
)
{
    // Capillary pressure limited so that p1 and p2 stay above the floor
    const dimensionedScalar aSmall(dimless, small);
    const volScalarField pMargin
    (
        max(p_ - dimensionedScalar(dimPressure, phaseFloors::p), 0*p_)
    );

    pc_ = min
    (
        max(sigmaKappa, -pMargin/max(twoPhaseMixture::alpha2(), aSmall)),
        pMargin/max(twoPhaseMixture::alpha1(), aSmall)
    );
    pc_.correctBoundaryConditions();

    correctPhasePressures();

    phase1_.initialize(p1_, T_);
    phase2_.initialize(p2_, T_);

    correctRho();

    updateMixtureTransportAndAcoustics();
}


void Foam::compressibleTwoPhaseMixture::correct
(
    const volScalarField& eMixture,
    const volScalarField& sigmaKappa
)
{
    // Step 1: pressure relaxation (see phaseRelaxation.H) at the mixture
    // density of correctRho(), which limits the capillary pressure
    pc_ = sigmaKappa;

    relaxation_.update
    (
        phase1_,
        phase2_,
        rho_,
        eMixture,
        pc_,
        p_,
        T_
    );

    // Step 2: phase pressures of the relaxed state
    correctPhasePressures();

    // Step 3: conservative decode of per-phase rho and c
    phase1_.correct(p1_, T_);
    phase2_.correct(p2_, T_);

    // Step 4: boundary state from the relaxed primitives
    correctBoundariesFromPrimitives();

    // Refresh the temperature-dependent phase transport
    phase1_.correctTransport(T_);
    phase2_.correctTransport(T_);

    // Step 5: mixing-rule transport and frozen mixture acoustics
    updateMixtureTransportAndAcoustics();
}


Foam::tmp<Foam::volScalarField>
Foam::compressibleTwoPhaseMixture::nu() const
{
    return nu_;
}


Foam::tmp<Foam::scalarField>
Foam::compressibleTwoPhaseMixture::nu(const label patchi) const
{
    return nu_.boundaryField()[patchi];
}


void Foam::compressibleTwoPhaseMixture::correctBoundariesFromPrimitives()
{
    // Re-impose per-phase boundary state from primitive variables
    phase1_.correctBoundariesFromPrimitives(p1_, T_);
    phase2_.correctBoundariesFromPrimitives(p2_, T_);

    // Mixture density of the re-imposed partial densities, floored as in
    // correctRho(). The coupled patches keep the values exchanged by
    // correctRho()
    forAll(rho_.boundaryField(), patchi)
    {
        if (!rho_.boundaryField()[patchi].coupled())
        {
            rho_.boundaryFieldRef()[patchi] ==
                max
                (
                    phase1_.alphaRho().boundaryField()[patchi]
                  + phase2_.alphaRho().boundaryField()[patchi],
                    phaseFloors::rho
                );
        }
    }
}


// ************************************************************************* //
