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

\*---------------------------------------------------------------------------*/

#include "NeutronFluxCoupledBase.H"
#include "fluidThermo.H"
#include "solidThermo.H"
#include "thermophysicalTransportModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::NeutronFluxCoupledBase::NeutronFluxCoupledBase
(
    const fvPatch& patch
)
:
    patch_(patch),
    myDFieldName_("undefined-myDFieldName"),
    nbrDFieldName_("undefined-nbrDFieldName")
{}


Foam::NeutronFluxCoupledBase::NeutronFluxCoupledBase
(
    const fvPatch& patch,
    const dictionary& dict
)
:
    patch_(patch),
    myDFieldName_(dict.lookup("myDiffFieldName")),
    nbrDFieldName_(dict.lookup("nbrDiffFieldName"))
{}


Foam::NeutronFluxCoupledBase::NeutronFluxCoupledBase
(
    const fvPatch& patch,
    const NeutronFluxCoupledBase& base

)
:
    patch_(patch),
    myDFieldName_(base.myDFieldName_),
    nbrDFieldName_(base.nbrDFieldName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::NeutronFluxCoupledBase::Diff
(
    const fvPatchScalarField& Fluxp
) const
{

    const fvMesh& mesh = patch_.boundaryMesh().mesh();
    const label patchi = patch_.index();

    if (mesh.foundObject<volScalarField>(myDFieldName_))
    {
        return patch_.lookupPatchField<volScalarField, scalar>
        (
          myDFieldName_
        );
    }
    else if (mesh.foundObject<volScalarField>(nbrDFieldName_))
    {
        return patch_.lookupPatchField<volScalarField, scalar>
        (
          nbrDFieldName_
        );
    }
    else
    {
        FatalErrorInFunction
              << "Did not find field " << myDFieldName_
              << "or field " << nbrDFieldName_
              << " on mesh " << mesh.name() << " patch " << patch_.name()
              << exit(FatalError);
    }

}




void Foam::NeutronFluxCoupledBase::write(Ostream& os) const
{}


// ************************************************************************* //
