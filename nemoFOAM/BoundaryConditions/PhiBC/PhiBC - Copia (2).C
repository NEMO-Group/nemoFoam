/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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

#include "PhiBC.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
/*
Foam::scalar Foam::PhiBC::t() const
{
    return db().time().timeOutputValue();
}
*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PhiBC::
PhiBC
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    beta_(p.size(), Zero),
    diffCoeffName_("DiffCoeffName"),
    J0_(p.size(), Zero),
    Direction_("Undefined")

{
    refValue() = 0;
    refGrad() = 0;
    valueFraction() = 1.0;
}


Foam::PhiBC::
PhiBC
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    beta_("beta", dict, p.size()),
    diffCoeffName_(dict.lookupOrDefault<word>("DiffCoeffName", "wordDefault")),
    J0_("Current", dict, p.size()),
    Direction_(dict.lookupOrDefault<word>("CurrentDirection", "Not defined"))


{
    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0;
        valueFraction() = 1;
    }

    /*
    // Initialise with the value entry if evaluation is not possible
    fvPatchScalarField::operator=
    (
        scalarField("value", dict, p.size())
    );
    refValue() = *this;
    */
}


Foam::PhiBC::
PhiBC
(
    const PhiBC& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    beta_(mapper(ptf.beta_)),
    diffCoeffName_(ptf.diffCoeffName_),
    J0_(mapper(ptf.J0_)),
    Direction_(ptf.Direction_)

{}


Foam::PhiBC::
PhiBC
(
    const PhiBC& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    beta_(ptf.beta_),
    diffCoeffName_(ptf.diffCoeffName_),
    J0_(ptf.J0_),
    Direction_(ptf.Direction_)

{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::PhiBC::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
    m(beta_, beta_);
    m(J0_, J0_);

}


void Foam::PhiBC::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const PhiBC& tiptf =
        refCast<const PhiBC>(ptf);

    beta_.rmap(tiptf.beta_, addr);
    J0_.rmap(tiptf.J0_, addr);

}


void Foam::PhiBC::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchScalarField& diffCoeffp =
        patch().lookupPatchField<volScalarField, scalar>(diffCoeffName_);

    const scalarField gamma_ = beta_;


    refGrad() = 0;
    refValue() = 2*J0_;


    if (Direction_ == "Positive")
    {
      valueFraction() = (gamma_)/(gamma_ - diffCoeffp*patch().deltaCoeffs());
    }
    else if (Direction_ == "Negative")
    {
      valueFraction() = (gamma_)/(gamma_ + diffCoeffp*patch().deltaCoeffs());
    }

    mixedFvPatchScalarField::updateCoeffs();
}



void Foam::PhiBC::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);


    writeEntry(os, "beta", beta_);
    writeEntry(os, "diffCoeffName", diffCoeffName_);
    writeEntry(os, "Current", J0_);
    writeEntry(os, "refValue", refValue());
    writeEntry(os, "refGradient", refGrad());
    writeEntry(os, "valueFraction", valueFraction());
//    writeEntry(os, "value", value);

}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        PhiBC
    );
}

// ************************************************************************* //
