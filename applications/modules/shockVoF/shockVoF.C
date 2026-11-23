/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2026 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2026 Martin Heinrich
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

#include "shockVoF.H"
#include "fvcSurfaceIntegrate.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(shockVoF, 0);
    addToRunTimeSelectionTable(solver, shockVoF, fvMesh);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::solvers::shockVoF::nRungeKuttaStages() const
{
    const label nStages =
        pimple.dict().lookupOrDefault<label>("nStages", 2);

    if (nStages < 2)
    {
        FatalIOErrorInFunction(pimple.dict())
            << "nStages must be at least 2, got " << nStages
            << exit(FatalIOError);
    }

    return nStages;
}


void Foam::solvers::shockVoF::correctCoNum()
{
    // Courant numbers of the full time step, as in other OpenFOAM solvers
    const scalar deltaT = runTime.deltaTValue();

    // Total mesh volume for both mean Courant numbers
    const scalar sumV = gSum(mesh.V().primitiveField());

    // Courant number from the maximum signal speed.  With the flow frozen
    // no wave travels and only alpha1 is advected, so the convective
    // Courant number of the volumetric flux is used instead
    const scalarField sumSpeedSf
    (
        pimple.flow()
      ? fvc::surfaceSum(flux_.amaxSf())().primitiveField()
      : fvc::surfaceSum(mag(flux_.phiv()))().primitiveField()
    );

    CoNum_ =
        0.5*gMax(sumSpeedSf/mesh.V().primitiveField())*deltaT;

    const scalar meanCoNum =
        0.5
       *(gSum(sumSpeedSf)/sumV)
       *deltaT;

    Info<< "Courant Number mean: " << meanCoNum
        << " max: " << CoNum_ << endl;

    const scalarField sumPhi
    (
        interface_.nearInterface()().primitiveField()
       *fvc::surfaceSum(mag(flux_.phiv()))().primitiveField()
    );

    alphaCoNum_ =
        0.5*gMax(sumPhi/mesh.V().primitiveField())*deltaT;

    const scalar meanAlphaCoNum =
        0.5
       *(gSum(sumPhi)/sumV)
       *deltaT;

    Info<< "Interface Courant Number mean: " << meanAlphaCoNum
        << " max: " << alphaCoNum_ << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::shockVoF::shockVoF(fvMesh& mesh)
:
    basicFluidSolver(mesh),

    U_
    (
        IOobject
        (
            "U",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    mixture_(mesh),

    p_(mixture_.p()),
    rho_(mixture_.rho()),

    alpha1_(mixture_.phase1().alpha()),
    alpha2_(mixture_.phase2().alpha()),

    alpha1Rho1_(mixture_.phase1().alphaRho()),
    alpha2Rho2_(mixture_.phase2().alphaRho()),

    alpha1RhoE1_(mixture_.phase1().alphaE()),
    alpha2RhoE2_(mixture_.phase2().alphaE()),

    rhoU_
    (
        IOobject
        (
            "rhoU",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho_*U_
    ),

    rhoE_
    (
        IOobject
        (
            "rhoE",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha1RhoE1_ + alpha2RhoE2_ + 0.5*rho_*magSqr(U_)
    ),

    interface_(alpha1_, mixture_),

    flux_(p_, U_, mixture_, interface_),

    momentumTransport_
    (
        rho_,
        U_,
        flux_.phi(),
        flux_.phiv(),
        flux_.phiAlpha1Rho1(),
        flux_.phiAlpha2Rho2(),
        mixture_
    ),

    thermophysicalTransport_(mixture_, momentumTransport_),

    alphaCoNum_(0)
{
    // Seed the phase states at the phase pressures of the initial
    // capillary pressure, and the conserved mixture fields from them
    interface_.correct();
    mixture_.initialize(interface_.sigmaKappa());
    rhoU_ == rho_*U_;
    rhoE_ == alpha1RhoE1_ + alpha2RhoE2_ + 0.5*rho_*magSqr(U_);

    // Create old-time step copies of the conserved fields
    alpha1_.oldTime();
    alpha1Rho1_.oldTime();
    alpha2Rho2_.oldTime();
    rhoU_.oldTime();
    alpha1RhoE1_.oldTime();
    alpha2RhoE2_.oldTime();
    rhoE_.oldTime();

    // Evaluate the boundary conditions of U before the first time step
    U_.correctBoundaryConditions();

    // Re-evaluate the inviscid fluxes for the initial time-step /
    // Courant-number estimate.
    flux_.update();

    if (transient())
    {
        correctCoNum();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::shockVoF::~shockVoF()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::solvers::shockVoF::maxDeltaT() const
{
    // Signal-speed Courant limit (and fvModels/maxDeltaT) from the base
    scalar deltaT = basicFluidSolver::maxDeltaT();

    // Additional interface Courant limit
    const scalar maxAlphaCo =
        runTime.controlDict().lookupOrDefault<scalar>("maxAlphaCo", vGreat);

    if (maxAlphaCo < vGreat && alphaCoNum_ > small)
    {
        deltaT = min(deltaT, maxAlphaCo/alphaCoNum_*runTime.deltaTValue());
    }

    return deltaT;
}


void Foam::solvers::shockVoF::preSolve()
{
    const bool meshChanged = mesh_.update();

    // Update primitive fields and fluxes
    if (meshChanged)
    {
        updatePrimitives();
        interface_.correct();

        flux_.update();
    }

    // Signal-speed and interface Courant numbers for time-step control
    correctCoNum();
}


void Foam::solvers::shockVoF::moveMesh()
{}


void Foam::solvers::shockVoF::motionCorrector()
{}


void Foam::solvers::shockVoF::prePredictor()
{
    if (!pimple.flow())
    {
        momentumPredictor();
    }
}


void Foam::solvers::shockVoF::momentumTransportPredictor()
{
    momentumTransport_.predict();
}


void Foam::solvers::shockVoF::thermophysicalTransportPredictor()
{
    // foamRun calls this whenever thermophysics is solved, but only the
    // energy equation reads the thermophysical transport, and it is not
    // advanced with the flow frozen
    if (pimple.flow())
    {
        thermophysicalTransport_.predict();
    }
}


void Foam::solvers::shockVoF::thermophysicalPredictor()
{
    // The energy equations are advanced inside momentumPredictor()
}


void Foam::solvers::shockVoF::pressureCorrector()
{
    // The pressure is recovered by relaxation inside momentumPredictor()
}


void Foam::solvers::shockVoF::momentumTransportCorrector()
{
    momentumTransport_.correct();
}


void Foam::solvers::shockVoF::thermophysicalTransportCorrector()
{
    // Guarded on pimple.flow() for the reason given in
    // thermophysicalTransportPredictor()
    if (pimple.flow())
    {
        thermophysicalTransport_.correct();
    }
}


void Foam::solvers::shockVoF::postSolve()
{}


// ************************************************************************* //
