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

#include "fluxScheme.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fluxScheme> Foam::fluxScheme::New
(
    const fvMesh& mesh,
    const dictionary& dict
)
{
    const word modelType(dict.lookup("type"));

    Info<< "Selecting high-order flux scheme " << modelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction(dict)
            << "Unknown " << fluxScheme::typeName << " type "
            << modelType << nl << nl
            << "Valid fluxScheme types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return
        autoPtr<fluxScheme>
        (
            cstrIter()(mesh, dict)
        );
}


Foam::autoPtr<Foam::fluxScheme> Foam::fluxScheme::New
(
    const word& type,
    const fvMesh& mesh,
    const dictionary& dict
)
{
    Info<< "Selecting low-order flux scheme " << type << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->find(type);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction(dict)
            << "Unknown " << fluxScheme::typeName << " type "
            << type << nl << nl
            << "Valid fluxScheme types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return
        autoPtr<fluxScheme>
        (
            cstrIter()(mesh, dict)
        );
}


// ************************************************************************* //
