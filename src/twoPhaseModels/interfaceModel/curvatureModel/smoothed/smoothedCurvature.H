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
    Foam::curvatureModels::smoothed

Description
    Smoothed-color-function curvature.

    The interface normal is the gradient of the smoothed color function
    (Shukla et al. 2010) provided by interfaceNormal, with its wide,
    point-based stencil.  Interpolated linearly to the faces, independent of
    fvSchemes, its divergence is compact second order, and the curvature is
    then smoothed across the diffuse band using a vertex-neighbor box
    average (the interfacial filter of Garrick et al. 2017a).

    kappa is evaluated in the interface band only, the cells with
    alpha1*(1 - alpha1) >= alphaTol, and is zero elsewhere.  The capillary
    force left unbalanced at the band edge scales with alphaTol.

    Dictionary entry:
    \verbatim
    surfaceTension
    {
        sigma           0.072;

        curvatureModel  smoothed;
        nFilter         10;
        alphaTol        1e-12;
    }
    \endverbatim

    The color function is set in phaseProperties/interfaceNormal.

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
    smoothedCurvature.C

\*---------------------------------------------------------------------------*/

#ifndef smoothedCurvature_H
#define smoothedCurvature_H

#include "curvatureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace curvatureModels
{

/*---------------------------------------------------------------------------*\
                          Class smoothed Declaration
\*---------------------------------------------------------------------------*/

class smoothed
:
    public curvatureModel
{
    // Private Data

        //- Stabilization of the interface unit normal against a vanishing
        //  volume-fraction gradient [1/m]
        dimensionedScalar deltaN_;

        //- Threshold on alpha1*(1 - alpha1) marking the cells in which the
        //  interface is resolved well enough for -div(nHat) to be a
        //  curvature
        scalar alphaTol_;

        //- Number of interfacial-filter iterations
        label nFilter_;

        //- Floor in the interfacial-filter denominator
        dimensionedScalar wSmall_;


    // Private Member Functions

        //- Vertex-neighbor box average via cell -> point -> cell
        //  interpolation (Garrick 2017a)
        tmp<volScalarField> pointCellAverage(const volScalarField& c) const;


public:

    //- Runtime type information
    TypeName("smoothed");


    // Constructors

        //- Construct from mesh, phase-1 fraction, interface normal and
        //  the surfaceTension sub-dictionary
        smoothed
        (
            const fvMesh& mesh,
            const volScalarField& alpha1,
            const interfaceNormal& normal,
            const dictionary& dict
        );


    //- Destructor
    virtual ~smoothed();


    // Member Functions

        //- Curvature from the smoothed color function, masked to the
        //  interface band
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace curvatureModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
