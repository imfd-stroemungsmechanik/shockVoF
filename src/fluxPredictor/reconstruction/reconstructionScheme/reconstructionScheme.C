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

#include "reconstructionScheme.H"
#include "positivityLimiter.H"
#include "swapCoupledOwnToNei.H"
#include "Switch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(reconstructionScheme, 0);
    defineRunTimeSelectionTable(reconstructionScheme, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reconstructionScheme::reconstructionScheme
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    positivityLimiter_
    (
        dict.lookupOrDefault<Switch>("positivityLimiter", true)
    ),
    rhoThincPtr_(rhoThinc::New(mesh, dict)),
    curTimeIndex_(-1)
{
    Info<< "    positivityLimiter: " << Switch(positivityLimiter_) << endl;

    if (debug && positivityLimiter_)
    {
        positivityLimiterPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "positivityLimiter",
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
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::reconstructionScheme::~reconstructionScheme()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::reconstructionScheme::reconstruct
(
    const volScalarField& alpha1,
    const volScalarField& p1,
    const volScalarField& p2,
    const volVectorField& U,
    const volScalarField& rho1,
    const volScalarField& rho2,
    const interfaceNormal& normal,
    tmp<surfaceScalarField>& alpha1Own,
    tmp<surfaceScalarField>& alpha1Nei,
    tmp<surfaceScalarField>& p1Own,
    tmp<surfaceScalarField>& p1Nei,
    tmp<surfaceScalarField>& p2Own,
    tmp<surfaceScalarField>& p2Nei,
    tmp<surfaceVectorField>& UOwn,
    tmp<surfaceVectorField>& UNei,
    tmp<surfaceScalarField>& rho1Own,
    tmp<surfaceScalarField>& rho1Nei,
    tmp<surfaceScalarField>& rho2Own,
    tmp<surfaceScalarField>& rho2Nei
) const
{
    // Scheme-specific bulk reconstruction of every primitive
    reconstructStates
    (
        alpha1, p1, p2, U, rho1, rho2,
        alpha1Own, alpha1Nei,
        p1Own,     p1Nei,
        p2Own,     p2Nei,
        UOwn,      UNei,
        rho1Own,   rho1Nei,
        rho2Own,   rho2Nei
    );

    if (rhoThincPtr_.valid())
    {
        // Interface sharpening: overwrite alpha (tanh profile) and the
        // phasic densities (cell value, as in Garrick et al. 2017b) in the
        // interface band.
        rhoThincPtr_->correct
        (
            alpha1, normal.gradC(), rho1, rho2,
            alpha1Own.ref(), alpha1Nei.ref(),
            rho1Own.ref(), rho1Nei.ref(),
            rho2Own.ref(), rho2Nei.ref()
        );
    }

    // Clamp alpha
    alpha1Own.ref() = max(min(alpha1Own(), 1.0), 0.0);
    alpha1Nei.ref() = max(min(alpha1Nei(), 1.0), 0.0);

    // Zhang-Shu positivity limiter, incl. the coupled-patch swap
    if (positivityLimiter_)
    {
        // Diagnostics: reset the theta field at the start of a new time
        // step, then min-accumulate over this step's reconstruct() calls
        volScalarField* positivityLimiterDiag = nullptr;
        if (debug)
        {
            volScalarField& pl = positivityLimiterPtr_();

            if (curTimeIndex_ != mesh_.time().timeIndex())
            {
                pl == dimensionedScalar(dimless, 1);

                curTimeIndex_ = mesh_.time().timeIndex();
            }

            positivityLimiterDiag = &pl;
        }

        positivityLimit
        (
            alpha1Own.ref(),
            alpha1Nei.ref(),
            p1Own.ref(),
            p1Nei.ref(),
            p2Own.ref(),
            p2Nei.ref(),
            UOwn.ref(),
            UNei.ref(),
            rho1Own.ref(),
            rho1Nei.ref(),
            rho2Own.ref(),
            rho2Nei.ref(),
            alpha1,
            p1,
            p2,
            U,
            rho1,
            rho2,
            positivityLimiterDiag
        );
    }
    else
    {
        // Each side of a coupled face reconstructs and sharpens only its
        // own (owner) state, so the neighbor state is taken from the other
        // side, all six primitives.
        swapCoupledOwnToNei(alpha1Own(), alpha1Nei.ref());
        swapCoupledOwnToNei(p1Own(),     p1Nei.ref());
        swapCoupledOwnToNei(p2Own(),     p2Nei.ref());
        swapCoupledOwnToNei(UOwn(),      UNei.ref());
        swapCoupledOwnToNei(rho1Own(),   rho1Nei.ref());
        swapCoupledOwnToNei(rho2Own(),   rho2Nei.ref());
    }
}


// ************************************************************************* //
