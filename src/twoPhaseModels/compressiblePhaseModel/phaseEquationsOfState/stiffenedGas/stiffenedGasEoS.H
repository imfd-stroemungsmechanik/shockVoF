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
    Foam::stiffenedGasEoS

Description
    Stiffened-gas equation of state,

        p   = (gamma - 1)*rho*e - gamma*pInf
        rho = (p + pInf)/(R*T),   with R = (gamma - 1)*Cv
        c   = sqrt(gamma*(p + pInf)/rho)
        T   = (p + pInf)/(R*rho)
        e   = (p + gamma*pInf)/(rho*(gamma - 1))

    Setting pInf = 0 recovers the ideal-gas law, which is also available
    as type idealGas (see idealGasEoS).

    Read from the equationOfState sub-dictionary of
    physicalProperties.<phase>, e.g., for water:
    \verbatim
    equationOfState
    {
        type    stiffenedGas;
        gamma   4.4;
        pInf    6e8;
        Cv      602.1;
    }
    \endverbatim

SourceFiles
    stiffenedGasEoS.C

See also
    Foam::phaseEquationOfState
    Foam::idealGasEoS

\*---------------------------------------------------------------------------*/

#ifndef stiffenedGasEoS_H
#define stiffenedGasEoS_H

#include "phaseEquationOfState.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class stiffenedGasEoS Declaration
\*---------------------------------------------------------------------------*/

class stiffenedGasEoS
:
    public phaseEquationOfState
{
    // Private Data

        // Equation of state parameters

            //- Heat capacity ratio [-]
            dimensionedScalar gamma_;

            //- Stiffness pressure [Pa]
            dimensionedScalar pInf_;

            //- Isochoric heat capacity [J/(kg K)]
            dimensionedScalar Cv_;

            //- Gas constant R = (gamma - 1)*Cv [J/(kg K)]
            dimensionedScalar R_;


protected:

    // Protected Constructors

        //- Construct from phase name and sub-dictionary with the given
        //  stiffness pressure pInf
        stiffenedGasEoS
        (
            const word& phaseName,
            const dictionary& dict,
            const scalar pInf
        );


public:

    //- Runtime type information
    TypeName("stiffenedGas");


    // Constructors

        //- Construct from phase name and sub-dictionary
        stiffenedGasEoS(const word& phaseName, const dictionary& dict);


    //- Destructor
    virtual ~stiffenedGasEoS() = default;


    // Member Functions

        // EoS-specific parameter access (dimensioned, for field
        // expressions in EoS-aware components)

            //- Heat capacity ratio [-]
            const dimensionedScalar& gamma() const
            {
                return gamma_;
            }

            //- Stiffness pressure [Pa]
            const dimensionedScalar& pInf() const
            {
                return pInf_;
            }

            //- Gas constant R = (gamma - 1)*Cv [J/(kg K)]
            const dimensionedScalar& R() const
            {
                return R_;
            }


        // Point-wise thermodynamic evaluations (virtual overrides)

            //- Density from pressure and temperature
            virtual scalar rho(const scalar p, const scalar T) const;

            //- Specific internal energy from pressure and density
            virtual scalar e(const scalar p, const scalar rho) const;

            //- Temperature from pressure and density
            virtual scalar T(const scalar p, const scalar rho) const;

            //- Speed of sound from pressure and density
            virtual scalar c(const scalar p, const scalar rho) const;


        // Thermodynamic property access

            //- Isobaric heat capacity Cp = gamma*Cv [J/(kg K)]
            virtual dimensionedScalar Cp() const
            {
                return gamma_*Cv_;
            }

            //- Isochoric heat capacity [J/(kg K)]
            virtual const dimensionedScalar& Cv() const
            {
                return Cv_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
