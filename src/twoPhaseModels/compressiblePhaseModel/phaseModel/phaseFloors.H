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

Namespace
    Foam::phaseFloors

Description
    Floors and thresholds shared by the phase models, the pressure
    relaxation and the reconstruction.

\*---------------------------------------------------------------------------*/

#ifndef phaseFloors_H
#define phaseFloors_H

#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseFloors
{
    //- Volume-fraction floor
    constexpr scalar alpha = 1.0e-8;

    //- Volume-fraction threshold below which a phase is treated as
    //  vanished by the pressure relaxation
    constexpr scalar relaxAlpha = 1.0e-6;

    //- Minimum admissible per-phase density
    constexpr scalar rho = 1.0e-6;

    //- Minimum admissible pressure
    constexpr scalar p = 1.0e-6;

} // End namespace phaseFloors
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
