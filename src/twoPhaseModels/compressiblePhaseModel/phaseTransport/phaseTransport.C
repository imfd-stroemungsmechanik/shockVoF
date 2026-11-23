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

#include "phaseTransport.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseTransport::phaseTransport
(
    const word& phaseName,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    muModel_
    (
        Function1<scalar>::New
        (
            "mu",
            dimTemperature,
            dimDynamicViscosity,
            dict
        )
    ),
    kappaModel_
    (
        Function1<scalar>::New
        (
            "kappa",
            dimTemperature,
            dimPower/(dimLength*dimTemperature),
            dict
        )
    ),

    mu_
    (
        IOobject
        (
            IOobject::groupName("mu", phaseName),
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("mu", dimDynamicViscosity, 0),
        "calculated"
    ),
    kappa_
    (
        IOobject
        (
            IOobject::groupName("kappa", phaseName),
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("kappa", dimPower/(dimLength*dimTemperature), 0),
        "calculated"
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::phaseTransport::correct(const volScalarField& T)
{
    // Internal field
    mu_.primitiveFieldRef() = muModel_->value(T.primitiveField());
    kappa_.primitiveFieldRef() = kappaModel_->value(T.primitiveField());

    // Boundary field
    volScalarField::Boundary& muBf = mu_.boundaryFieldRef();
    volScalarField::Boundary& kappaBf = kappa_.boundaryFieldRef();
    const volScalarField::Boundary& TBf = T.boundaryField();

    forAll(muBf, patchi)
    {
        muBf[patchi] = muModel_->value(TBf[patchi]);
        kappaBf[patchi] = kappaModel_->value(TBf[patchi]);
    }
}


// ************************************************************************* //
