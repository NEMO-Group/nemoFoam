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
    DFieldName_("undefined-DFieldName")
{}


Foam::NeutronFluxCoupledBase::NeutronFluxCoupledBase
(
    const fvPatch& patch,
    const dictionary& dict
)
:
    patch_(patch),
    DFieldName_(dict.lookup("DiffFieldName"))
{}


Foam::NeutronFluxCoupledBase::NeutronFluxCoupledBase
(
    const fvPatch& patch,
    const NeutronFluxCoupledBase& base

)
:
    patch_(patch),
    DFieldName_(base.DFieldName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::NeutronFluxCoupledBase::D
(
    const fvPatchScalarField& Fluxp
) const
{
    const fvMesh& mesh = patch_.boundaryMesh().mesh();
    const label patchi = patch_.index();
//    const objectRegistry& db();
//   const volScalarField& D(patchi);//= mesh.lookupObject<volScalarField>("D");
//   const word Dname = dict.Dname;
//   const word& phase(Tp.internalField().group())
   const word Dname = DFieldName_;
   const volScalarField& D = mesh.lookupObject<volScalarField>(Dname);
   const fvPatchScalarField& D_p = D.boundaryField()[patchi];

    return D_p;//mesh.lookupObject<volScalarField>("D2_2").value(patchi);//(patchi);

 /*   const word& phase(Tp.internalField().group());

    const word fluidThermoName
    (
        IOobject::groupName(basicThermo::dictName, phase)
    );

    if (mesh.foundObject<fluidThermo>(fluidThermoName))
    {
        static word ttmName
        (
            IOobject::groupName
            (
                thermophysicalTransportModel::typeName,
                phase
            )
        );

        if (mesh.foundObject<thermophysicalTransportModel>(ttmName))
        {
            const thermophysicalTransportModel& ttm =
                mesh.lookupObject<thermophysicalTransportModel>(ttmName);

            return ttm.kappaEff(patchi);
        }
        else
        {
            const fluidThermo& thermo =
                mesh.lookupObject<fluidThermo>(fluidThermoName);

            return thermo.kappa(patchi);
        }
    }
    else if (mesh.foundObject<solidThermo>(basicThermo::dictName))
    {
        const solidThermo& thermo =
            mesh.lookupObject<solidThermo>(basicThermo::dictName);

        if (!thermo.isotropic())
        {
            const symmTensorField kappa(thermo.KappaLocal(patchi));
            const vectorField n(patch_.nf());

            return n & kappa & n;
        }
        else
        {
            return thermo.kappa(patchi);
        }
    }
*/
/*    else
    {
        FatalErrorInFunction
            << "Cannot find a fluidThermo or solidThermo instance"
            << exit(FatalError);

        return scalarField::null();
    }
    */
}


void Foam::NeutronFluxCoupledBase::write(Ostream& os) const
{}


// ************************************************************************* //
