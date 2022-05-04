/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::temperatureCoupledBase

Description
    Common functions used in temperature coupled boundaries.

    The thermal conductivity \c kappa is obtained from the region fluidThermo
    or solidThermo as appropriate.

    For turbulent fluid regions \c kappa is obtained from the
    thermophysicalTransportModel \c kappaEff.

    For solid regions kappa may be obtained from the anisotropic alpha field
    by specifying the optional \c alphaAni field name.

    \par Keywords provided by this class:
      \table
        Property | Description                     | Required | Default value
        alphaAni | Name of the non-isotropic alpha | no  |
      \endtable

SourceFiles
    temperatureCoupledBase.C

\*---------------------------------------------------------------------------*/

#ifndef NeutronFluxCoupledBase_H
#define NeutronFluxCoupledBase_H

#include "scalarField.H"
#include "NamedEnum.H"
#include "fvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class NeutronFluxCoupledBase Declaration
\*---------------------------------------------------------------------------*/

class NeutronFluxCoupledBase
{
    // Private data

        //- Underlying patch
        const fvPatch& patch_;

   //- Name of D
        const word myDFieldName_;

  //- Name of D
        const word nbrDFieldName_;


public:

    // Constructors

        //- Construct from patch
        NeutronFluxCoupledBase
        (
            const fvPatch& patch
        );

        //- Construct from patch and dictionary
        NeutronFluxCoupledBase
        (
            const fvPatch& patch,
            const dictionary& dict
        );

         //- Construct from patch and NeutronFluxCoupledBase
        NeutronFluxCoupledBase
        (
            const fvPatch& patch,
            const NeutronFluxCoupledBase& base
        );


    // Member Functions

        //- Given patch Neutron Flux calculate corresponding D field
        tmp<scalarField> myD(const fvPatchScalarField& Fluxp) const;

        tmp<scalarField> nbrD(const fvPatchScalarField& Fluxp) const;

        //- Write
        void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
