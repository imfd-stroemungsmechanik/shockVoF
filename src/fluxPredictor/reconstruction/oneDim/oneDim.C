/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2025 OpenFOAM Foundation
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

#include "oneDim.H"
#include "fvMeshStitcher.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(oneDim, 0);
    addToRunTimeSelectionTable(reconstructionScheme, oneDim, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::oneDim::oneDim(const fvMesh& mesh, const dictionary& dict)
:
    reconstructionScheme(mesh, dict),
    dict_(dict)
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::oneDim::updateDirections() const
{
    if (mesh_.topoChanging() || mesh_.stitcher().stitches())
    {
        pos_.clear();
        neg_.clear();
    }

    if (!pos_.valid())
    {
        pos_ = surfaceScalarField::New
        (
            "pos",
            mesh_,
            dimensionedScalar(dimless, 1)
        );

        neg_ = surfaceScalarField::New
        (
            "neg",
            mesh_,
            dimensionedScalar(dimless, -1)
        );
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::oneDim::reconstructStates
(
    const volScalarField& alpha1,
    const volScalarField& p1,
    const volScalarField& p2,
    const volVectorField& U,
    const volScalarField& rho1,
    const volScalarField& rho2,
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
    updateDirections();

    // Limited face reconstruction of every primitive (scheme per role,
    // selected under reconstruction as reconstruct(alpha|p|U|rho)).  The
    // phase pressures share the reconstruct(p) scheme, the phase densities
    // the reconstruct(rho) scheme.
    alpha1Own = interpolate(alpha1, pos_(), "alpha");
    alpha1Nei = interpolate(alpha1, neg_(), "alpha");
    p1Own     = interpolate(p1,     pos_(), "p");
    p1Nei     = interpolate(p1,     neg_(), "p");
    p2Own     = interpolate(p2,     pos_(), "p");
    p2Nei     = interpolate(p2,     neg_(), "p");
    UOwn      = interpolate(U,      pos_(), "U");
    UNei      = interpolate(U,      neg_(), "U");
    rho1Own   = interpolate(rho1,   pos_(), "rho");
    rho1Nei   = interpolate(rho1,   neg_(), "rho");
    rho2Own   = interpolate(rho2,   pos_(), "rho");
    rho2Nei   = interpolate(rho2,   neg_(), "rho");
}


// ************************************************************************* //
