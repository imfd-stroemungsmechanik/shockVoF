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

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::rhoThinc> Foam::rhoThinc::New
(
    const fvMesh& mesh,
    const dictionary& dict
)
{
    // No rhoTHINC sub-dictionary is the same as type none
    word modelType("none");

    if (dict.found("rhoTHINC"))
    {
        modelType = word(dict.subDict("rhoTHINC").lookup("type"));
    }

    if (modelType == "none")
    {
        Info<< "Selecting rho-THINC interface sharpening none" << endl;

        return autoPtr<rhoThinc>(nullptr);
    }

    const dictionary& thincDict = dict.subDict("rhoTHINC");

    Info<< "Selecting rho-THINC interface sharpening " << modelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        // none is handled above and never registered, but it is valid
        wordList validTypes(dictionaryConstructorTablePtr_->sortedToc());
        validTypes.setSize(validTypes.size() + 1);
        validTypes[validTypes.size() - 1] = word("none");

        FatalIOErrorInFunction(thincDict)
            << "Unknown " << rhoThinc::typeName << " type "
            << modelType << nl << nl
            << "Valid rhoThinc types are:" << nl
            << validTypes
            << exit(FatalIOError);
    }

    return
        autoPtr<rhoThinc>
        (
            cstrIter()(mesh, thincDict)
        );
}


// ************************************************************************* //
