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

#include "interfaceModel.H"
#include "linear.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceModel::interfaceModel
(
    const volScalarField& alpha1,
    const dictionary& phaseProperties
)
:
    alpha1_(alpha1),
    sigma_
    (
        "sigma",
        dimForce/dimLength,
        phaseProperties.subDict("surfaceTension").lookup<scalar>("sigma")
    ),
    normal_
    (
        alpha1.mesh(),
        alpha1,
        phaseProperties.subOrEmptyDict("interfaceNormal")
    ),
    curvature_
    (
        curvatureModel::New
        (
            alpha1.mesh(),
            alpha1,
            normal_,
            phaseProperties.subDict("surfaceTension")
        )
    ),
    sigmaKappaf_
    (
        IOobject
        (
            "sigmaKappaf",
            alpha1.time().name(),
            alpha1.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha1.mesh(),
        dimensionedScalar("sigmaKappaf", dimPressure, 0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interfaceModel::~interfaceModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::interfaceModel::nearInterface() const
{
    return pos0(alpha1_ - 0.01)*pos0(0.99 - alpha1_);
}


void Foam::interfaceModel::correct()
{
    normal_.correct();
    curvature_->correct();
    sigmaKappaf_ == sigma_*linearInterpolate(curvature_->kappa());
}


// ************************************************************************* //
