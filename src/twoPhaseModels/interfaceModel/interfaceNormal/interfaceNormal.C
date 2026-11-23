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

#include "interfaceNormal.H"
#include "gaussGrad.H"
#include "volPointInterpolation.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::interfaceNormal::color() const
{
    if (nSmooth_ == 0 && psiExp_ == 1)
    {
        return tmp<volScalarField>(alpha1_);
    }

    // The copy keeps alpha1's patch types, so coupled patches stay coupled
    tmp<volScalarField> tc
    (
        new volScalarField
        (
            IOobject
            (
                "interfaceNormal:C",
                mesh_.time().name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            alpha1_
        )
    );
    volScalarField& c = tc.ref();

    const labelUList& owner = mesh_.owner();
    const labelUList& neighbor = mesh_.neighbour();

    scalarField sum;
    labelList count;

    for (label pass = 0; pass < nSmooth_; ++pass)
    {
        // Average over the cell and its face neighbors
        const scalarField& prev = c.primitiveField();

        sum = prev;
        count.setSize(prev.size());
        count = 1;

        forAll(owner, facei)
        {
            sum[owner[facei]] += prev[neighbor[facei]];
            ++count[owner[facei]];

            sum[neighbor[facei]] += prev[owner[facei]];
            ++count[neighbor[facei]];
        }

        // alpha1's halo is current, so the first pass needs no exchange
        forAll(c.boundaryField(), patchi)
        {
            const fvPatchScalarField& pc = c.boundaryField()[patchi];

            if (pc.coupled())
            {
                const scalarField pnf(pc.patchNeighbourField());
                const labelUList& faceCells = pc.patch().faceCells();

                forAll(faceCells, fi)
                {
                    sum[faceCells[fi]] += pnf[fi];
                    ++count[faceCells[fi]];
                }
            }
        }

        scalarField& cur = c.primitiveFieldRef();

        forAll(cur, celli)
        {
            cur[celli] = sum[celli]/count[celli];
        }

        c.correctBoundaryConditions();
    }

    // The color function maps values pointwise, so it is applied to the
    // patch values too; re-evaluating fixed-value patches would undo it
    if (psiExp_ != 1)
    {
        const scalar e = psiExp_;

        auto C = [e](const scalar alpha)
        {
            const scalar a = min(max(alpha, scalar(0)), scalar(1));
            const scalar aE = pow(a, e);

            return aE/(aE + pow(1.0 - a, e));
        };

        scalarField& ci = c.primitiveFieldRef();

        forAll(ci, celli)
        {
            ci[celli] = C(ci[celli]);
        }

        volScalarField::Boundary& cBf = c.boundaryFieldRef();

        forAll(cBf, patchi)
        {
            forAll(cBf[patchi], facei)
            {
                cBf[patchi][facei] = C(cBf[patchi][facei]);
            }
        }
    }

    return tc;
}


Foam::tmp<Foam::volVectorField>
Foam::interfaceNormal::pointGrad(const volScalarField& c) const
{
    // Cell to point interpolation brings in the diagonal neighbors,
    // which a compact face interpolation does not see
    const volPointInterpolation& vpi = volPointInterpolation::New(mesh_);
    const pointScalarField cP(vpi.interpolate(c));

    // Face values as the average of the face's vertex values
    surfaceScalarField cf
    (
        IOobject
        (
            "interfaceNormal:cf",
            mesh_.time().name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedScalar("0", c.dimensions(), 0)
    );

    const faceList& faces = mesh_.faces();

    scalarField& cfi = cf.primitiveFieldRef();

    forAll(cfi, facei)
    {
        const face& f = faces[facei];

        scalar s = 0;

        forAll(f, fp)
        {
            s += cP[f[fp]];
        }

        cfi[facei] = s/f.size();
    }

    surfaceScalarField::Boundary& cfbf = cf.boundaryFieldRef();

    forAll(cfbf, patchi)
    {
        fvsPatchScalarField& pcf = cfbf[patchi];
        const label start = mesh_.boundary()[patchi].start();

        forAll(pcf, i)
        {
            const face& f = faces[start + i];

            scalar s = 0;

            forAll(f, fp)
            {
                s += cP[f[fp]];
            }

            pcf[i] = s/f.size();
        }
    }

    // gradf exchanges the halo of the gradient
    return fv::gaussGrad<scalar>::gradf(cf, "interfaceNormal:gradC");
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceNormal::interfaceNormal
(
    const fvMesh& mesh,
    const volScalarField& alpha1,
    const dictionary& dict
)
:
    mesh_(mesh),
    alpha1_(alpha1),
    psiExp_(dict.lookupOrDefault<scalar>("psiExp", 0.1)),
    nSmooth_(dict.lookupOrDefault<label>("nSmooth", 0))
{
    if (psiExp_ <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "interfaceNormal psiExp must be positive, got " << psiExp_
            << exit(FatalIOError);
    }

    if (nSmooth_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "interfaceNormal nSmooth must be non-negative, got "
            << nSmooth_ << exit(FatalIOError);
    }

    Info<< "Interface normal" << nl
        << "    psiExp           : " << psiExp_ << nl
        << "    nSmooth          : " << nSmooth_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interfaceNormal::~interfaceNormal()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::interfaceNormal::correct()
{
    gradCPtr_.clear();
}


const Foam::volVectorField& Foam::interfaceNormal::gradC() const
{
    if (!gradCPtr_.valid())
    {
        const tmp<volScalarField> tc(color());

        gradCPtr_.reset(pointGrad(tc()).ptr());
    }

    return gradCPtr_();
}


// ************************************************************************* //
