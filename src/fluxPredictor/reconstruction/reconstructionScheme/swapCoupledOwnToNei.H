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

InNamespace
    Foam

Description
    Replace the neighbor side of every coupled patch with the owner side
    of the cell across the coupling, so both sides see identical L/R states.

\*---------------------------------------------------------------------------*/

#ifndef swapCoupledOwnToNei_H
#define swapCoupledOwnToNei_H

#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Overwrite the coupled-patch values of fNei with the values of fOwn on the
//  other side of the coupling
template<class Type>
inline void swapCoupledOwnToNei
(
    const SurfaceField<Type>& fOwn,
    SurfaceField<Type>& fNei
)
{
    const tmp<typename SurfaceField<Type>::Boundary> tOwnNbr
    (
        fOwn.boundaryField().boundaryNeighbourField()
    );
    const typename SurfaceField<Type>::Boundary& ownNbr = tOwnNbr();

    typename SurfaceField<Type>::Boundary& neiBf = fNei.boundaryFieldRef();

    forAll(neiBf, patchi)
    {
        if (fOwn.boundaryField()[patchi].coupled())
        {
            neiBf[patchi] == ownNbr[patchi];
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
