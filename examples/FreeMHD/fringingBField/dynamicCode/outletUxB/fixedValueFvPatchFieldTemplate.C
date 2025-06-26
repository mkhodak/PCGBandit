/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
    Copyright (C) YEAR AUTHOR, AFFILIATION
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

#include "fixedValueFvPatchFieldTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
#include "PatchFunction1.H"

//{{{ begin codeInclude
#line 21 "/scratch/gpfs/bw0594/FOAM_RUN/lmxCh13.2_ClosedChannel_mesh1mm/system/liquid/codeDict.outletUxB"
#include "scalar.H"
        #include "fvCFD.H"
        #include "mappedPatchBase.H"
        #include "fvPatchFieldMapper.H"
        #include "volFields.H"
//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// dynamicCode:
// SHA1 = d724b23ade7e72ed88804038169381cbd93eab02
//
// unique function name that can be checked if the correct library version
// has been loaded
extern "C" void outletUxB_d724b23ade7e72ed88804038169381cbd93eab02(bool load)
{
    if (load)
    {
        // Code that can be explicitly executed after loading
    }
    else
    {
        // Code that can be explicitly executed before unloading
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makeRemovablePatchTypeField
(
    fvPatchScalarField,
    outletUxBFixedValueFvPatchScalarField
);

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
outletUxBFixedValueFvPatchScalarField::
outletUxBFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    parent_bctype(p, iF)
{
    if (false)
    {
        printMessage("Construct outletUxB : patch/DimensionedField");
    }
}


Foam::
outletUxBFixedValueFvPatchScalarField::
outletUxBFixedValueFvPatchScalarField
(
    const outletUxBFixedValueFvPatchScalarField& rhs,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    parent_bctype(rhs, p, iF, mapper)
{
    if (false)
    {
        printMessage("Construct outletUxB : patch/DimensionedField/mapper");
    }
}


Foam::
outletUxBFixedValueFvPatchScalarField::
outletUxBFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    parent_bctype(p, iF, dict)
{
    if (false)
    {
        printMessage("Construct outletUxB : patch/dictionary");
    }
}


Foam::
outletUxBFixedValueFvPatchScalarField::
outletUxBFixedValueFvPatchScalarField
(
    const outletUxBFixedValueFvPatchScalarField& rhs
)
:
    parent_bctype(rhs),
    dictionaryContent(rhs)
{
    if (false)
    {
        printMessage("Copy construct outletUxB");
    }
}


Foam::
outletUxBFixedValueFvPatchScalarField::
outletUxBFixedValueFvPatchScalarField
(
    const outletUxBFixedValueFvPatchScalarField& rhs,
    const DimensionedField<scalar, volMesh>& iF
)
:
    parent_bctype(rhs, iF)
{
    if (false)
    {
        printMessage("Construct outletUxB : copy/DimensionedField");
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::
outletUxBFixedValueFvPatchScalarField::
~outletUxBFixedValueFvPatchScalarField()
{
    if (false)
    {
        printMessage("Destroy outletUxB");
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::
outletUxBFixedValueFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        printMessage("updateCoeffs outletUxB");
    }

//{{{ begin code
    #line 35 "/scratch/gpfs/bw0594/FOAM_RUN/lmxCh13.2_ClosedChannel_mesh1mm/system/liquid/codeDict.outletUxB"
//get this patch's mesh
    const fvMesh& mesh 	= patch().boundaryMesh().mesh();
    const label id 		= patch().index();
   
   //get mesh variable
    const volVectorField& Bmesh 	= mesh.lookupObject<volVectorField>("B");
    const volVectorField& Umesh 	= mesh.lookupObject<volVectorField>("U");
    const volScalarField& Phimesh 	= mesh.lookupObject<volScalarField>("potE");
    
    //get patch internal varialbe and delta
    const vectorField B      = Bmesh.boundaryField().boundaryInternalField()[id];
    const vectorField U      = Umesh.boundaryField().boundaryInternalField()[id];
    const scalarField phi    = Phimesh.boundaryField().boundaryInternalField()[id];
    const scalarField delta  = patch().deltaCoeffs();

    
    scalarField result(patch().size(), Zero);
    vectorField nf = patch().nf();
    
    forAll(result, i)
    {
		result[i] = phi[i]+((U[i]^B[i])&nf[i]/delta[i]);

	}
	
    
    operator==(result);
//}}} end code

    this->parent_bctype::updateCoeffs();
}


// ************************************************************************* //

