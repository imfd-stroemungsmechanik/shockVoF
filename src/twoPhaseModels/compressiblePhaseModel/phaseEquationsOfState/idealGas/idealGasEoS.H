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

Class
    Foam::idealGasEoS

Description
    Ideal-gas (calorically perfect gas) equation of state: the stiffened-gas
    equation of state with the stiffness pressure pInf = 0,

        p   = (gamma - 1)*rho*e
        rho = p/(R*T),   with R = (gamma - 1)*Cv
        c   = sqrt(gamma*p/rho)
        T   = p/(R*rho)
        e   = p/(rho*(gamma - 1)) = Cv*T

    It is a stiffenedGasEoS, so components that require the stiffened gas
    accept it as well.

    Read from the equationOfState sub-dictionary of
    physicalProperties.<phase>, e.g., for air:
    \verbatim
    equationOfState
    {
        type    idealGas;
        gamma   1.4;
        Cv      717.5;
    }
    \endverbatim

    Cv = R/(gamma - 1) is the isochoric heat capacity [J/(kg K)].

SourceFiles
    idealGasEoS.C

See also
    Foam::stiffenedGasEoS
    Foam::phaseEquationOfState

\*---------------------------------------------------------------------------*/

#ifndef idealGasEoS_H
#define idealGasEoS_H

#include "stiffenedGasEoS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class idealGasEoS Declaration
\*---------------------------------------------------------------------------*/

class idealGasEoS
:
    public stiffenedGasEoS
{
public:

    //- Runtime type information
    TypeName("idealGas");


    // Constructors

        //- Construct from phase name and sub-dictionary
        idealGasEoS(const word& phaseName, const dictionary& dict);


    //- Destructor
    virtual ~idealGasEoS() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
