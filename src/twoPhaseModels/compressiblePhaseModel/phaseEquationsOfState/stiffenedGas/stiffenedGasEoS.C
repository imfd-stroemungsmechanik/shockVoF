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

#include "stiffenedGasEoS.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(stiffenedGasEoS, 0);

    addToRunTimeSelectionTable
    (
        phaseEquationOfState,
        stiffenedGasEoS,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::stiffenedGasEoS::stiffenedGasEoS
(
    const word& phaseName,
    const dictionary& dict,
    const scalar pInf
)
:
    phaseEquationOfState(),
    gamma_
    (
        "gamma",
        dimless,
        dict.lookup<scalar>("gamma")
    ),
    pInf_
    (
        "pInf",
        dimPressure,
        pInf
    ),
    Cv_
    (
        "Cv",
        dimEnergy/dimMass/dimTemperature,
        dict.lookup<scalar>("Cv")
    ),
    R_
    (
        "R",
        Cv_*(gamma_ - 1.0)
    )
{
    if (gamma_.value() <= 1)
    {
        FatalIOErrorInFunction(dict)
            << "gamma must be greater than 1 for phase "
            << phaseName << ", got " << gamma_.value()
            << exit(FatalIOError);
    }

    if (pInf_.value() < 0)
    {
        FatalIOErrorInFunction(dict)
            << "pInf must be non-negative for phase "
            << phaseName << ", got " << pInf_.value()
            << exit(FatalIOError);
    }

    if (Cv_.value() <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "Cv must be positive for phase "
            << phaseName << ", got " << Cv_.value()
            << exit(FatalIOError);
    }
}


Foam::stiffenedGasEoS::stiffenedGasEoS
(
    const word& phaseName,
    const dictionary& dict
)
:
    stiffenedGasEoS(phaseName, dict, dict.lookup<scalar>("pInf"))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::stiffenedGasEoS::rho
(
    const scalar p,
    const scalar T
) const
{
    return (p + pInf_.value())/(R_.value()*T);
}


Foam::scalar Foam::stiffenedGasEoS::e
(
    const scalar p,
    const scalar rho
) const
{
    const scalar g    = gamma_.value();
    const scalar pInf = pInf_.value();
    return (p + g*pInf)/(rho*(g - 1.0));
}


Foam::scalar Foam::stiffenedGasEoS::T
(
    const scalar p,
    const scalar rho
) const
{
    return (p + pInf_.value())/(R_.value()*rho);
}


Foam::scalar Foam::stiffenedGasEoS::c
(
    const scalar p,
    const scalar rho
) const
{
    return Foam::sqrt(gamma_.value()*(p + pInf_.value())/rho);
}


// ************************************************************************* //
