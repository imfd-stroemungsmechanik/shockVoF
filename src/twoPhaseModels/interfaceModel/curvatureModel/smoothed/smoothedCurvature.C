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

#include "smoothedCurvature.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcDiv.H"
#include "linear.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace curvatureModels
{
    defineTypeNameAndDebug(smoothed, 0);
    addToRunTimeSelectionTable(curvatureModel, smoothed, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::curvatureModels::smoothed::pointCellAverage(const volScalarField& c) const
{
    // The composite stencil is the full vertex neighborhood (9-point in 2D,
    // 27-point in 3D)
    const volPointInterpolation& vpi = volPointInterpolation::New(mesh_);
    const pointScalarField cP(vpi.interpolate(c));

    tmp<volScalarField> tavg
    (
        new volScalarField("curvatureModel::box(" + c.name() + ")", c)
    );
    volScalarField& avg = tavg.ref();

    // Point -> cell: mean over each cell's vertices
    const labelListList& cellPoints = mesh_.cellPoints();

    scalarField& avgi = avg.primitiveFieldRef();
    forAll(avgi, celli)
    {
        const labelList& pts = cellPoints[celli];

        scalar sum = 0;
        forAll(pts, i)
        {
            sum += cP[pts[i]];
        }

        avgi[celli] = sum/pts.size();
    }

    // Non-coupled patch values feed the boundary-point values of the next
    // cell-to-point pass; setting them to the adjacent cell value (zero
    // gradient) lets the filter act up to the boundary.
    volScalarField::Boundary& avgBf = avg.boundaryFieldRef();
    forAll(avgBf, patchi)
    {
        if (!avgBf[patchi].coupled())
        {
            avgBf[patchi] == avgBf[patchi].patchInternalField();
        }
    }

    return tavg;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::curvatureModels::smoothed::smoothed
(
    const fvMesh& mesh,
    const volScalarField& alpha1,
    const interfaceNormal& normal,
    const dictionary& dict
)
:
    curvatureModel(mesh, alpha1, normal, dict),
    deltaN_("deltaN", 1e-8/pow(average(mesh.V()), 1.0/3.0)),
    alphaTol_(dict.lookupOrDefault<scalar>("alphaTol", 1e-12)),
    nFilter_(dict.lookupOrDefault<label>("nFilter", 10)),
    wSmall_("wSmall", dimless, small)
{
    Info<< "    nFilter          : " << nFilter_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::curvatureModels::smoothed::~smoothed()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::curvatureModels::smoothed::correct()
{
    const surfaceVectorField& Sf = mesh_.Sf();

    // Bounded volume fraction and a*(1 - a), for the interface band and
    // the filter weights
    const volScalarField a(min(max(alpha1_, scalar(0)), scalar(1)));
    const volScalarField a1a(a*(1.0 - a));

    // Color gradient, interpolated linearly to the faces
    const surfaceVectorField gradCf(linearInterpolate(normal_.gradC()));

    // Face unit interface normal (deltaN-stabilized) and its area flux
    const surfaceVectorField nHatfv(gradCf/(mag(gradCf) + deltaN_));
    const surfaceScalarField nHatf(nHatfv & Sf);

    // Interface band: the cells in which the interface is resolved (1) or
    // not (0)
    const volScalarField band
    (
        "band",
        pos0(a1a - dimensionedScalar(dimless, alphaTol_))
    );

    // Garrick's interfacial filter (Garrick et al. 2017a):
    // box average weighted by w = (a(1-a))^2
    volScalarField kappaBand("kappaBand", band*(-fvc::div(nHatf)));

    const volScalarField w(band*sqr(a1a));
    const volScalarField wBox(pointCellAverage(w));

    for (label i = 0; i < nFilter_; ++i)
    {
        kappaBand == band*pointCellAverage(w*kappaBand)/max(wBox, wSmall_);
    }

    // Evaluate the zeroGradient boundary so that sigmaKappaf is well
    // defined on boundary faces; the masked expression above carries only
    // calculated patch values.
    kappa_ == kappaBand;
    kappa_.correctBoundaryConditions();
}


// ************************************************************************* //
