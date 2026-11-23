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
    Foam::phaseTransport

Description
    Per-phase transport-property closure: the dynamic viscosity and
    thermal conductivity as Function1s of temperature.

    Both properties are read from the transport sub-dictionary of
    constant/physicalProperties.<phase> as OpenFOAM Function1 entries, so
    any Function1 type is accepted (constant, table, polynomial, ...):

    \verbatim
    transport
    {
        mu      1e-3;
        kappa   table
        (
            (273 0.56)
            (373 0.68)
        );
    }
    \endverbatim

    A bare scalar (mu 1e-3;) is read by Function1 as a constant.

    correct(T) evaluates mu(T) and kappa(T) into cached fields at the
    supplied temperature, which is the shared mixture temperature
    (single-temperature closure; see compressibleTwoPhaseMixture::correct).

SourceFiles
    phaseTransport.C

\*---------------------------------------------------------------------------*/

#ifndef phaseTransport_H
#define phaseTransport_H

#include "volFields.H"
#include "dictionary.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class phaseTransport Declaration
\*---------------------------------------------------------------------------*/

class phaseTransport
{
    // Private Data

        //- Dynamic viscosity as a function of temperature [Pa s]
        autoPtr<Function1<scalar>> muModel_;

        //- Thermal conductivity as a function of temperature [W/(m K)]
        autoPtr<Function1<scalar>> kappaModel_;

        //- Cached dynamic viscosity field [Pa s]
        volScalarField mu_;

        //- Cached thermal conductivity field [W/(m K)]
        volScalarField kappa_;


public:

    // Constructors

        //- Construct from phase name, mesh and transport sub-dictionary
        phaseTransport
        (
            const word& phaseName,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        phaseTransport(const phaseTransport&) = delete;


    //- Destructor
    ~phaseTransport() = default;


    // Member Functions

        //- Dynamic viscosity field [Pa s]
        const volScalarField& mu() const
        {
            return mu_;
        }

        //- Thermal conductivity field [W/(m K)]
        const volScalarField& kappa() const
        {
            return kappa_;
        }

        //- Evaluate mu(T) and kappa(T) into the cached fields at the
        //  supplied temperature
        void correct(const volScalarField& T);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const phaseTransport&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
