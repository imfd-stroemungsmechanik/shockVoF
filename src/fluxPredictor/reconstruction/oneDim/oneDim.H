/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2025 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2026 Martin Heinrich
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
    Foam::oneDim

Description
    One-dimensional (rhoCentralFoam/shockFluid style) face reconstruction.

    Each primitive is reconstructed to its owner/neighbor face state with
    an independently selectable, direction-biased limited interpolation
    scheme.  The bias is imposed by interpolating against the +1 (owner)
    and -1 (neighbor) face-direction fields, so any flux-limited surface
    interpolation scheme may be used per role:

    \verbatim
    interpolationSchemes
    {
        reconstruction
        {
            type                oneDim;

            reconstruct(alpha)  vanLeer;
            reconstruct(p)      vanLeer;
            reconstruct(U)      vanLeerV;
            reconstruct(rho)    vanLeer;

            positivityLimiter   yes;

            rhoTHINC
            {
                type            THINC;
                beta            1.0;
            }
        }
    }
    \endverbatim

    All four reconstruct(...) entries are required; reconstruct(p) applies
    to both phase pressures and reconstruct(rho) to both phase densities.

    Interface sharpening (rho-THINC) and positivity limiting are provided
    by the base class; oneDim implements only the bulk reconstruction.

SourceFiles
    oneDim.C

See also
    Foam::reconstructionScheme
    Foam::rhoThinc

\*---------------------------------------------------------------------------*/

#ifndef oneDim_H
#define oneDim_H

#include "reconstructionScheme.H"
#include "surfaceInterpolationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class oneDim Declaration
\*---------------------------------------------------------------------------*/

class oneDim
:
    public reconstructionScheme
{
    // Private Data

        //- The reconstruction sub-dictionary, holding the per-role
        //  reconstruct(alpha|p|U|rho) scheme specifications
        const dictionary dict_;

        //- Cached +1/-1 face-direction fields
        mutable tmp<surfaceScalarField> pos_;
        mutable tmp<surfaceScalarField> neg_;


    // Private Member Functions

        //- (Re)build the cached +1/-1 face-direction fields, clearing them
        //  first if the mesh topology has changed
        void updateDirections() const;

        //- Direction-biased limited interpolation of field vf, using the
        //  scheme named reconstruct(<key>) from the reconstruction dict
        template<class Type>
        inline tmp<SurfaceField<Type>> interpolate
        (
            const VolField<Type>& vf,
            const surfaceScalarField& dir,
            const word& key
        ) const
        {
            tmp<SurfaceField<Type>> tsf
            (
                surfaceInterpolationScheme<Type>::New
                (
                    mesh_,
                    dir,
                    dict_.lookup("reconstruct(" + key + ')')
                )().interpolate(vf)
            );

            tsf.ref().rename(vf.name() + '_' + dir.name());

            return tsf;
        }


protected:

    // Protected Member Functions

        //- Reconstruct the bulk owner/neighbor face states of every
        //  primitive with the per-role limited schemes
        virtual void reconstructStates
        (
            const volScalarField& alpha1,
            const volScalarField& p1,
            const volScalarField& p2,
            const volVectorField& U,
            const volScalarField& rho1,
            const volScalarField& rho2,
            tmp<surfaceScalarField>& alpha1Own,
            tmp<surfaceScalarField>& alpha1Nei,
            tmp<surfaceScalarField>& p1Own,
            tmp<surfaceScalarField>& p1Nei,
            tmp<surfaceScalarField>& p2Own,
            tmp<surfaceScalarField>& p2Nei,
            tmp<surfaceVectorField>& UOwn,
            tmp<surfaceVectorField>& UNei,
            tmp<surfaceScalarField>& rho1Own,
            tmp<surfaceScalarField>& rho1Nei,
            tmp<surfaceScalarField>& rho2Own,
            tmp<surfaceScalarField>& rho2Nei
        ) const;


public:

    //- Runtime type information
    TypeName("oneDim");


    // Constructors

        //- Construct from mesh and the reconstruction sub-dictionary
        oneDim(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~oneDim() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
