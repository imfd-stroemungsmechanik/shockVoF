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
    Foam::rhoThinc

Description
    Abstract base class for rho-THINC interface sharpening (Garrick et al.
    2017b).

    correct() overwrites the bulk-reconstructed face states in admissible
    interface cells: alpha1 with a THINC profile along the interface normal
    (the color gradient of interfaceNormal), rho1 and rho2 with their cell
    values.  A cell is admissible if it passes the three gates
      G1: delta < alpha1 < 1 - delta,
      G2: |grad(C)| > small,
      G3: alpha1 is monotone along the normal across its neighbors.

    Selected by the rhoTHINC sub-dictionary of the reconstruction
    dictionary.

    \verbatim
    rhoTHINC
    {
        type    THINC;
        beta    1.0;
        delta   1e-4;
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
    rhoThinc.C
    rhoThincNew.C

See also
    Foam::THINC
    Foam::fluxPredictor

\*---------------------------------------------------------------------------*/

#ifndef rhoThinc_H
#define rhoThinc_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class rhoThinc Declaration
\*---------------------------------------------------------------------------*/

class rhoThinc
{
    // Private Member Functions

        //- Disallow default bitwise copy construction
        rhoThinc(const rhoThinc&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const rhoThinc&) = delete;


protected:

    // Protected Data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- THINC steepness parameter beta (dimensionless)
        scalar beta_;

        //- Band threshold: sharpen only where delta < alpha1 < 1 - delta
        scalar delta_;


    // Protected Member Functions

        //- Build per-cell CSR-packed boundary-neighbor data for the
        //  monotonicity gate (G3).  Coupled-patch neighbor values are
        //  prefetched (one halo exchange per coupled patch)
        void buildBoundaryNeighborData
        (
            const volScalarField& alpha1,
            labelList& bndOffset,
            vectorField& bndDelta,
            scalarField& bndVal
        ) const;


public:

    //- Runtime type information
    TypeName("rhoThinc");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            rhoThinc,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& dict
            ),
            (mesh, dict)
        );


    // Constructors

        //- Construct from mesh and the rhoTHINC sub-dictionary
        rhoThinc(const fvMesh& mesh, const dictionary& dict);


    // Selector

        //- Select from the reconstruction sub-dictionary.  Returns
        //  nullptr if there is no rhoTHINC sub-dictionary or its type is
        //  none (no sharpening)
        static autoPtr<rhoThinc> New
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~rhoThinc();


    // Member Functions

        //- Overwrite the alpha1, rho1 and rho2 face states in admissible
        //  cells.  On coupled patches only the owner side is written
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
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
