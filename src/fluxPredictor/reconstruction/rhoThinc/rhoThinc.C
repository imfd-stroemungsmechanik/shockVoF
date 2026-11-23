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

#include "rhoThinc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rhoThinc, 0);
    defineRunTimeSelectionTable(rhoThinc, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rhoThinc::rhoThinc(const fvMesh& mesh, const dictionary& dict)
:
    mesh_(mesh),
    beta_(dict.lookupOrDefault<scalar>("beta", 1.0)),
    delta_(dict.lookupOrDefault<scalar>("delta", 1e-4))
{
    if (beta_ <= 0 || beta_ > 4)
    {
        FatalIOErrorInFunction(dict)
            << "rhoThinc beta must lie in (0, 4], got " << beta_
            << exit(FatalIOError);
    }

    if (delta_ <= 0 || delta_ >= 0.5)
    {
        FatalIOErrorInFunction(dict)
            << "rhoThinc delta must lie in (0, 0.5), got " << delta_
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::rhoThinc::~rhoThinc()
{}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::rhoThinc::buildBoundaryNeighborData
(
    const volScalarField& alpha1,
    labelList& bndOffset,
    vectorField& bndDelta,
    scalarField& bndVal
) const
{
    const vectorField& cc = mesh_.C();
    const label nCells = mesh_.nCells();

    PtrList<scalarField> coupledNei(mesh_.boundary().size());

    forAll(alpha1.boundaryField(), patchi)
    {
        if (alpha1.boundaryField()[patchi].coupled())
        {
            coupledNei.set
            (
                patchi,
                new scalarField
                (
                    alpha1.boundaryField()[patchi].patchNeighbourField()
                )
            );
        }
    }

    labelList bndCount(nCells, 0);
    forAll(alpha1.boundaryField(), patchi)
    {
        const labelUList& fcs =
            alpha1.boundaryField()[patchi].patch().faceCells();
        forAll(fcs, fi)
        {
            ++bndCount[fcs[fi]];
        }
    }

    bndOffset.setSize(nCells + 1);
    bndOffset[0] = 0;
    forAll(bndCount, cI)
    {
        bndOffset[cI + 1] = bndOffset[cI] + bndCount[cI];
    }
    const label nBndEntries = bndOffset[nCells];

    bndDelta.setSize(nBndEntries);
    bndDelta = Zero;
    bndVal.setSize(nBndEntries);
    bndVal = 0;

    bndCount = 0;
    forAll(alpha1.boundaryField(), patchi)
    {
        const fvPatch& pp = mesh_.boundary()[patchi];
        const fvPatchScalarField& pphi = alpha1.boundaryField()[patchi];
        const labelUList& fcs = pp.faceCells();

        if (pphi.coupled())
        {
            const vectorField pdelta(pp.delta());
            const scalarField& nei = coupledNei[patchi];

            forAll(fcs, fi)
            {
                const label cI = fcs[fi];
                const label k  = bndOffset[cI] + bndCount[cI]++;
                bndDelta[k] = pdelta[fi];
                bndVal[k]   = nei[fi];
            }
        }
        else
        {
            const vectorField& Cf = mesh_.Cf().boundaryField()[patchi];

            forAll(fcs, fi)
            {
                const label cI = fcs[fi];
                const label k  = bndOffset[cI] + bndCount[cI]++;
                bndDelta[k] = Cf[fi] - cc[cI];
                bndVal[k]   = pphi[fi];
            }
        }
    }
}


// ************************************************************************* //
