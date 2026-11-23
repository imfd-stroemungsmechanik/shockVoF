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
    Foam::interfaceNormal

Description
    Cell volume-fraction gradient, whose direction is the interface
    normal.  It is the one normal of the solver, used by the smoothed
    curvature model and by rho-THINC.

    The color function (Shukla et al. 2010; Garrick et al. 2017a) of the
    volume fraction a after nSmooth Jacobi passes is

        C = a^e/(a^e + (1 - a)^e),

    with the exponent e = psiExp and C = a for e = 1.  Its gradient is the
    Gauss gradient of face values averaged from the vertex values.  On a
    uniform Cartesian mesh this is the central difference weighted 1:2:1
    over the transverse neighbors, the wide 9-point (2D) or 27-point (3D)
    stencil of Garrick et al. (2017a).  It is hard-coded, independent of
    the gradSchemes of fvSchemes.

    The gradient is computed on the first request after correct(), which
    interfaceModel calls whenever alpha1 has changed.  The request involves
    communication, so it has to be made on all processors.

    \verbatim
    interfaceNormal
    {
        psiExp      0.1;
        nSmooth     0;
    }
    \endverbatim

    References:
    \verbatim
        Shukla, R. K., Pantano, C. & Freund, J. B. (2010).
        An interface capturing method for the simulation of multi-phase
        compressible flows.
        Journal of Computational Physics, 229(19), 7411-7439.

        Garrick, D. P., Owkes, M. & Regele, J. D. (2017a).
        A finite-volume HLLC-based scheme for compressible interfacial
        flows with surface tension.
        Journal of Computational Physics, 339, 46-67.
    \endverbatim

SourceFiles
    interfaceNormal.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceNormal_H
#define interfaceNormal_H

#include "volFields.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class interfaceNormal Declaration
\*---------------------------------------------------------------------------*/

class interfaceNormal
{
    // Private Data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Reference to the phase-1 volume fraction
        const volScalarField& alpha1_;

        //- Color-function exponent
        const scalar psiExp_;

        //- Jacobi passes on alpha1
        const label nSmooth_;

        //- Color gradient, empty while out of date
        mutable autoPtr<volVectorField> gradCPtr_;


    // Private Member Functions

        //- Color function
        tmp<volScalarField> color() const;

        //- Gauss gradient from face values averaged from the vertex values
        tmp<volVectorField> pointGrad(const volScalarField& c) const;


public:

    // Constructors

        //- Construct from mesh, phase-1 fraction and the interfaceNormal
        //  sub-dictionary
        interfaceNormal
        (
            const fvMesh& mesh,
            const volScalarField& alpha1,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        interfaceNormal(const interfaceNormal&) = delete;


    //- Destructor
    ~interfaceNormal();


    // Member Functions

        //- Mark the gradient out of date
        void correct();

        //- Color gradient
        const volVectorField& gradC() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const interfaceNormal&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
