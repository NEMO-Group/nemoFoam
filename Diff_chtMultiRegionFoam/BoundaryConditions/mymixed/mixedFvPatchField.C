/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "mymixedFvPatchField.H"
#include "fvc.H"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::mymixedFvPatchField<Type>::mymixedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF),
    J_in_(p.size()),
    D_coeff_(p.size()),
    Grad_Phi_(p.size()),
    n_bdr_(p.size())




{}


template<class Type>
Foam::mymixedFvPatchField<Type>::mymixedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF, dict, false),
    J_in_("J_in", dict, p.size()),
    D_coeff_("D_coeff", dict, p.size()),
    Grad_Phi_("Grad_Phi", dict, p.size()),
    n_bdr_("n_bdr", dict, p.size())
{
    evaluate();
}


template<class Type>
Foam::mymixedFvPatchField<Type>::mymixedFvPatchField
(
    const mymixedFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper,
    const bool mappingRequired
)
:
    fvPatchField<Type>(ptf, p, iF, mapper, mappingRequired),
    J_in_(mapper(ptf.J_in_)),
    n_bdr_(mapper(ptf.n_bdr_))

{}


template<class Type>
Foam::mymixedFvPatchField<Type>::mymixedFvPatchField
(
    const mymixedFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf, iF),
    J_in_(ptf.J_in_),
    D_coeff_(ptf.D_coeff_),
    Grad_Phi_(ptf.Grad_Phi_),
    n_bdr_(ptf.n_bdr_)

{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::mymixedFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fvPatchField<Type>::autoMap(m);
    m(J_in_, J_in_);
    m(n_bdr_, n_bdr_);
}


template<class Type>
void Foam::mymixedFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fvPatchField<Type>::rmap(ptf, addr);

    const mymixedFvPatchField<Type>& mptf =
        refCast<const mymixedFvPatchField<Type>>(ptf);

    J_in_.rmap(mptf.J_in, addr);
    n_bdr_.rmap(mptf.n_bdr_, addr);
}


template<class Type>
void Foam::mymixedFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    const fvMesh& mesh = patch_.boundaryMesh().mesh();
    const label patchi = patch_.index();
    const word& Dname = D_coeff_;
    const volScalarField& DFld = mesh.lookupObject<volScalarField>(Dname);
    const fvPatchScalarField& Dbndr = DFld.boundaryField()[patchi];

    const volScalarField& FluxFld = mesh.lookupObject<volScalarField>(Grad_Phi_);
    const fvPatchScalarField& Fluxbndr = FluxFld.boundaryField()[patchi];

    const volScalarField& grad = fvc::grad(Fluxbndr);
    n_bdr_ = mesh.Sf;


    Field<Type>::operator=
    (
        4 * J_in_
      -
        2 * Dbndr *
        (
            (grad * this->patch().deltaCoeffs()) & n_bdr_
        )
    );

    fvPatchField<Type>::evaluate();
}




template<class Type>
void Foam::mymixedFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    writeEntry(os, "J_in", J_in_);
    writeEntry(os, "D_coeff", D_coeff_);
    writeEntry(os, "Grad_Phi", Grad_Phi_);
    writeEntry(os, "n_bdr", n_bdr_);
}


// ************************************************************************* //
