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

#include "THINC.H"
#include "addToRunTimeSelectionTable.H"

#include <cmath>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(THINC, 0);
    addToRunTimeSelectionTable(rhoThinc, THINC, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::THINC::THINC(const fvMesh& mesh, const dictionary& dict)
:
    rhoThinc(mesh, dict)
{
    Info<< "    beta             : " << beta_ << nl
        << "    delta            : " << delta_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::THINC::~THINC()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::THINC::computeInterface
(
    const volScalarField& alpha1,
    const volVectorField& gradC,
    boolList& admissible,
    vectorField& nUnit
) const
{
    const vectorField& cc = mesh_.C();

    const scalarField& alphaI = alpha1.primitiveField();
    const vectorField& gradI  = gradC.primitiveField();

    // Boundary-neighbor data for the monotonicity gate (G3)
    labelList bndOffset;
    vectorField bndDelta;
    scalarField bndVal;
    buildBoundaryNeighborData(alpha1, bndOffset, bndDelta, bndVal);

    forAll(alphaI, cI)
    {
        const scalar alphaC = alphaI[cI];

        // G1: interface band
        if (alphaC <= delta_ || alphaC >= 1.0 - delta_)
        {
            continue;
        }

        // G2: well-defined interface normal
        const scalar magG = mag(gradI[cI]);
        if (magG < small)
        {
            continue;
        }

        // Unit normal, up-gradient
        const vector n = gradI[cI]/magG;

        // G3: monotonicity along n: some forward/backward neighbor pair
        // must bracket alphaC, else alpha is a local extremum along n and
        // the monotone tanh profile is inappropriate
        {
            scalar maxFwd = -great;
            scalar minFwd =  great;
            scalar maxBwd = -great;
            scalar minBwd =  great;
            bool hasFwd = false;
            bool hasBwd = false;

            const labelList& cellNb = mesh_.cellCells()[cI];
            forAll(cellNb, k)
            {
                const label nbI = cellNb[k];
                const scalar proj = (n & (cc[nbI] - cc[cI]));
                if (proj > small)
                {
                    maxFwd = max(maxFwd, alphaI[nbI]);
                    minFwd = min(minFwd, alphaI[nbI]);
                    hasFwd = true;
                }
                else if (proj < -small)
                {
                    maxBwd = max(maxBwd, alphaI[nbI]);
                    minBwd = min(minBwd, alphaI[nbI]);
                    hasBwd = true;
                }
            }

            for (label k = bndOffset[cI]; k < bndOffset[cI + 1]; ++k)
            {
                const scalar proj = (n & bndDelta[k]);
                if (proj > small)
                {
                    maxFwd = max(maxFwd, bndVal[k]);
                    minFwd = min(minFwd, bndVal[k]);
                    hasFwd = true;
                }
                else if (proj < -small)
                {
                    maxBwd = max(maxBwd, bndVal[k]);
                    minBwd = min(minBwd, bndVal[k]);
                    hasBwd = true;
                }
            }

            if (!hasFwd || !hasBwd)
            {
                continue;
            }

            const bool monoInc =
                (maxFwd > alphaC && minBwd < alphaC);
            const bool monoDec =
                (minFwd < alphaC && maxBwd > alphaC);

            if (!monoInc && !monoDec)
            {
                continue;
            }
        }

        nUnit[cI]      = n;
        admissible[cI] = true;
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::THINC::correct
(
    const volScalarField& alpha1,
    const volVectorField& gradC,
    const volScalarField& rho1,
    const volScalarField& rho2,
    surfaceScalarField& alpha1Own,
    surfaceScalarField& alpha1Nei,
    surfaceScalarField& rho1Own,
    surfaceScalarField& rho1Nei,
    surfaceScalarField& rho2Own,
    surfaceScalarField& rho2Nei
) const
{
    const label nCells = mesh_.nCells();

    // Interface detection and per-cell normal
    boolList    admissible(nCells, false);
    vectorField nUnit(nCells, Zero);

    computeInterface(alpha1, gradC, admissible, nUnit);

    const labelUList& owner     = mesh_.owner();
    const labelUList& neighbor  = mesh_.neighbour();
    const scalarField& alphaI   = alpha1.primitiveField();
    const scalarField& rho1I    = rho1.primitiveField();
    const scalarField& rho2I    = rho2.primitiveField();


    // Overwrite the owner/neighbor face states of a reconstructed field
    // in the admissible band with value(cell, face); non-band faces keep
    // their bulk reconstruction.  On coupled patches only the owner side
    // is written: the caller makes the neighbor side the remote owner
    // state, for bit-exact L/R (see reconstructionScheme::reconstruct()).
    auto overwrite =
        [&](surfaceScalarField& fOwn, surfaceScalarField& fNei, auto value)
    {
        scalarField& ownI = fOwn.primitiveFieldRef();
        scalarField& neiI = fNei.primitiveFieldRef();

        forAll(owner, faceI)
        {
            const label o = owner[faceI];
            const label n = neighbor[faceI];

            if (admissible[o])
            {
                ownI[faceI] = value(o, faceI);
            }
            if (admissible[n])
            {
                neiI[faceI] = value(n, faceI);
            }
        }

        forAll(mesh_.boundary(), patchi)
        {
            if (!mesh_.boundary()[patchi].coupled())
            {
                continue;
            }

            const polyPatch& ppoly = mesh_.boundary()[patchi].poly();
            const labelUList& fcs  = mesh_.boundary()[patchi].faceCells();
            fvsPatchScalarField& pOwn = fOwn.boundaryFieldRef()[patchi];

            forAll(fcs, fi)
            {
                const label c = fcs[fi];
                if (admissible[c])
                {
                    pOwn[fi] = value(c, ppoly.start() + fi);
                }
            }
        }
    };


    // Face volume fraction from the one-dimensional profile along the
    // outward face normal (Garrick et al. 2017b).  The face sits at
    // X = 1 of the local coordinate.
    const vectorField& fAreas = mesh_.faceAreas();
    const labelUList& fOwner  = mesh_.faceOwner();

    auto alphaThinc = [&](const label cI, const label fI) -> scalar
    {
        const scalar alphaC = alphaI[cI];

        // Outward unit face normal of cell cI at face fI.  faceAreas()
        // is indexed by global face label and points owner -> neighbor,
        // so it needs flipping when cI is the neighbor.
        vector dHat = fAreas[fI]/max(mag(fAreas[fI]), small);

        if (fOwner[fI] != cI)
        {
            dHat = -dHat;
        }

        const scalar nd = (nUnit[cI] & dHat);

        // Eq. 38: steepness weighted by the interface orientation, so a
        // direction tangential to the interface is left unsharpened
        const scalar betaF = beta_*mag(nd) + 0.01;

        // Garrick's sigma_i: the sign of the volume-fraction gradient
        // along dHat
        const scalar sgm = (nd >= 0) ? 1.0 : -1.0;

        // Eq. 40 substituted into Eq. 42, written with expm1 for accuracy
        // at small betaF:
        //   ratio   = (1 - exp(-u*a))/(exp(u*(1 - a)) - 1),  u = 2*sgm*betaF
        //   alpha_f = 0.5*(1 + tanh(sgm*betaF + 0.5*ln(ratio)))
        const scalar u = 2.0*sgm*betaF;

        const scalar num = -std::expm1(-u*alphaC);
        const scalar den =  std::expm1(u*(1.0 - alphaC));

        // Both share the sign of u, so the ratio is positive; guard the
        // degenerate ends anyway and fall back to the cell average there
        if (mag(num) < small || mag(den) < small)
        {
            return alphaC;
        }

        const scalar ratio = num/den;

        if (ratio <= 0)
        {
            return alphaC;
        }

        const scalar t = sgm*betaF + 0.5*Foam::log(ratio);

        return max(scalar(0), min(scalar(1), 0.5*(1.0 + Foam::tanh(t))));
    };

    // alpha : the tanh profile in the band
    overwrite(alpha1Own, alpha1Nei, alphaThinc);

    // rho1, rho2 : cell-center value in the band (Garrick et al. 2017b)
    overwrite
    (
        rho1Own, rho1Nei,
        [&](const label cI, const label) -> scalar { return rho1I[cI]; }
    );
    overwrite
    (
        rho2Own, rho2Nei,
        [&](const label cI, const label) -> scalar { return rho2I[cI]; }
    );
}


// ************************************************************************* //
