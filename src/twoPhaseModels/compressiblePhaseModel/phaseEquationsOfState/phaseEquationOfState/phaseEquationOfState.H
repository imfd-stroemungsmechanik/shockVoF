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
    Foam::phaseEquationOfState

Description
    Abstract base class for per-phase equations of state: pointwise
    rho(p, T), e(p, rho), T(p, rho) and c(p, rho), plus the constant Cp and
    Cv for mixing rules.  EoS-specific parameters stay on the concrete class
    and are reached by EoS-aware consumers via a downcast.  The solver
    requires both phase equations of state to be stiffenedGasEoS types
    (stiffenedGas, or idealGas, its special case with pInf = 0):
    compressibleTwoPhaseMixture binds them with refCast, so any other type is
    rejected at construction.

SourceFiles
    phaseEquationOfState.C
    phaseEquationOfStateNew.C

\*---------------------------------------------------------------------------*/

#ifndef phaseEquationOfState_H
#define phaseEquationOfState_H

#include "dictionary.H"
#include "dimensionedScalar.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class phaseEquationOfState Declaration
\*---------------------------------------------------------------------------*/

class phaseEquationOfState
{

public:

    //- Runtime type information
    TypeName("phaseEquationOfState");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            phaseEquationOfState,
            dictionary,
            (
                const word& phaseName,
                const dictionary& dict
            ),
            (phaseName, dict)
        );


    // Constructors

        //- Construct null
        phaseEquationOfState() = default;

        //- Disallow default bitwise copy construction
        phaseEquationOfState(const phaseEquationOfState&) = delete;


    // Selectors

        //- Select and construct a phaseEquationOfState
        static autoPtr<phaseEquationOfState> New
        (
            const word& phaseName,
            const dictionary& dict
        );


    //- Destructor
    virtual ~phaseEquationOfState() = default;


    // Member Functions

        // Point-wise thermodynamic evaluations

            //- Density from pressure and temperature [kg/m^3]
            virtual scalar rho
            (
                const scalar p,
                const scalar T
            ) const = 0;

            //- Specific internal energy from pressure and density [J/kg]
            virtual scalar e
            (
                const scalar p,
                const scalar rho
            ) const = 0;

            //- Temperature from pressure and density [K]
            virtual scalar T
            (
                const scalar p,
                const scalar rho
            ) const = 0;

            //- Speed of sound from pressure and density [m/s]
            virtual scalar c
            (
                const scalar p,
                const scalar rho
            ) const = 0;


        // Thermodynamic properties

            //- Isobaric heat capacity [J/(kg K)]
            virtual dimensionedScalar Cp() const = 0;

            //- Isochoric heat capacity [J/(kg K)]
            virtual const dimensionedScalar& Cv() const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const phaseEquationOfState&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
