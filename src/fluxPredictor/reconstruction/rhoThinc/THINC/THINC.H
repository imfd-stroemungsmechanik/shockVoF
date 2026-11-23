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
    Foam::THINC

Description
    THINC face reconstruction of Garrick et al. (2017b).

    The face volume fraction follows the 1D tanh profile along the outward
    face normal, with the profile offset in closed form and the
    steepness weighted by the interface orientation,

        beta_f = beta*|n & dHat| + 0.01.

    \verbatim
    rhoTHINC
    {
        type            THINC;
        beta            1.0;
        delta           1e-4;
    }
    \endverbatim

    References:
    \verbatim
        Garrick, D. P., Hagen, W. A. & Regele, J. D. (2017b).
        An interface capturing scheme for modeling atomization in
        compressible flows.
        Journal of Computational Physics, 344, 260-280.
    \endverbatim

SourceFiles
    THINC.C

See also
    Foam::rhoThinc

\*---------------------------------------------------------------------------*/

#ifndef THINC_H
#define THINC_H

#include "rhoThinc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class THINC Declaration
\*---------------------------------------------------------------------------*/

class THINC
:
    public rhoThinc
{
    // Private Member Functions

        //- Detect the interface band (gates G1-G3) and cache the unit
        //  normal
        void computeInterface
        (
            const volScalarField& alpha1,
            const volVectorField& gradC,
            boolList& admissible,
            vectorField& nUnit
        ) const;


public:

    //- Runtime type information
    TypeName("THINC");


    // Constructors

        //- Construct from mesh and the rhoTHINC sub-dictionary
        THINC(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~THINC();


    // Member Functions

        //- Sharpen the interface in place (Garrick THINC method)
        virtual void correct
        (
            const volScalarField& alpha1,
            const volVectorField& gradC,
            const volScalarField& rho1,
            const volScalarField& rho2,
            surfaceScalarField& alpha1Own,
            surfaceScalarField& alpha1Nei,
            surfaceScalarField& rho1Own,
            surfaceScalarField& rho1Nei,
            surfaceScalarField& rho2Own,
            surfaceScalarField& rho2Nei
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
