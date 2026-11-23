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
    Foam::reconstructionScheme

Description
    Abstract base class for the reconstruction of the primitive variables
    (alpha1, p1, p2, U, rho1, rho2) to owner/neighbor face states.

    The phase pressures p1 and p2 are reconstructed instead of the mixture
    pressure.  At the static capillary balance p1 - p2 = sigma*kappa both
    are uniform, so that the limiters do not act on the capillary pressure
    jump; the face mixture pressure alpha1_f*p1_f + alpha2_f*p2_f carries
    it through the face volume fraction.

    reconstruct() calls the scheme-specific reconstructStates() and then
    applies, in order:

      1. rho-THINC interface sharpening, if selected,
      2. clamping of the face volume fraction to [0, 1],
      3. the Zhang-Shu positivity limiter, if enabled,
      4. replacement of the neighbor state on coupled patches by the
         remote owner state.

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

    If rhoTHINC is absent or its type is none, no sharpening is applied.

SourceFiles
    reconstructionScheme.C
    reconstructionSchemeNew.C

See also
    Foam::oneDim
    Foam::rhoThinc
    Foam::fluxPredictor

\*---------------------------------------------------------------------------*/

#ifndef reconstructionScheme_H
#define reconstructionScheme_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "rhoThinc.H"
#include "interfaceNormal.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class reconstructionScheme Declaration
\*---------------------------------------------------------------------------*/

class reconstructionScheme
{
    // Private Member Functions

        //- Disallow default bitwise copy construction
        reconstructionScheme(const reconstructionScheme&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const reconstructionScheme&) = delete;


protected:

    // Protected Data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Apply the Zhang-Shu positivity-preserving limiter to the
        //  reconstructed face states
        const bool positivityLimiter_;

        //- rho-THINC interface sharpening (null if no rhoTHINC dict,
        //  or if its type is none)
        autoPtr<rhoThinc> rhoThincPtr_;


    // Positivity-limiter diagnostics

        //- Cell-centered positivity-limiter (Zhang-Shu theta) field for
        //  diagnostics, allocated iff the debug switch is set
        mutable autoPtr<volScalarField> positivityLimiterPtr_;

        //- Current time index (for resetting the diagnostics field)
        mutable label curTimeIndex_;


    // Protected Member Functions

        //- Reconstruct the bulk owner/neighbor face states of every
        //  primitive.  The only scheme-specific step; implemented by the
        //  derived class
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
        ) const = 0;


public:

    //- Runtime type information
    TypeName("reconstructionScheme");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            reconstructionScheme,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& dict
            ),
            (mesh, dict)
        );


    // Constructors

        //- Construct from mesh and the reconstruction sub-dictionary
        reconstructionScheme(const fvMesh& mesh, const dictionary& dict);


    // Selector

        //- Select from the interpolationSchemes dictionary, which requires
        //  a "reconstruction" sub-dictionary with a "type" entry
        static autoPtr<reconstructionScheme> New
        (
            const fvMesh& mesh,
            const dictionary& interpolationSchemes
        );


    //- Destructor
    virtual ~reconstructionScheme();


    // Member Functions

        //- Reconstruct the face states of all primitives and apply the
        //  steps listed in the class description
        void reconstruct
        (
            const volScalarField& alpha1,
            const volScalarField& p1,
            const volScalarField& p2,
            const volVectorField& U,
            const volScalarField& rho1,
            const volScalarField& rho2,
            const interfaceNormal& normal,
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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
