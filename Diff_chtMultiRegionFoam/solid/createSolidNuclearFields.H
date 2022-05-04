//************************** Open and Read the nuclearPropertiesDict******************************//

#include "readSolidNuclearProperties.H"

//************************** Read the nuclear Fields from latest time ******************************//

// Initialise fields

PtrList<volScalarField> fluxTotSolid(solidRegions.size());
PtrList<volScalarField> fluxSolid(nEnergyGroupsInit*solidRegions.size());
PtrList<volScalarField> ExternalNeutronSource(nEnergyGroupsInit*solidRegions.size());

PtrList<volScalarField> precSolid(nPrecGroupsInit*solidRegions.size());

PtrList<volScalarField> decSolid(nDecGroupsInit*solidRegions.size());

forAll(solidRegions,jj)
{

     // Define the total flux (it will be printed by the code) and read N-neutron flux fields from Initial fields in case/0 or latestTime

	Info<< "\n Reading "<< nEnergyGroupsInit << " neutron flux fields for region: "<< solidRegions[jj].name() << "\n" << endl;

	fluxTotSolid.set
    	(
        		jj,
        		new volScalarField
        		(
            		IOobject
            		(
									"fluxTot",
									runTime.timeName(),
									solidRegions[jj],
									IOobject::NO_READ,
									IOobject::AUTO_WRITE
								),
								solidRegions[jj],
								dimensionedScalar("", dimless/dimArea/dimTime, 0)
							)
			);

	for (label i=0; i<nEnergyGroupsInit; i++)
	{
		label i_jjInd = nEnergyGroupsInit*jj + i;

		fluxSolid.set
		(
			i_jjInd,
			new volScalarField
			(
				IOobject
				(
					"flux" + name(i+1),
					runTime.timeName(),
					solidRegions[jj],
					IOobject::MUST_READ,
					IOobject::AUTO_WRITE
				),
			solidRegions[jj]
			)
		);

		fluxTotSolid[jj] += fluxSolid[i_jjInd]; // initialize fluxTot with initial flux fields

		ExternalNeutronSource.set
		(
			i_jjInd,
			new volScalarField
			(
				IOobject
				(
					"ExternalNeutronSource" + name(i+1),
					runTime.timeName(),
					solidRegions[jj],
					IOobject::MUST_READ,
					IOobject::AUTO_WRITE
				),
			solidRegions[jj]
			)
		);

	}

	// Read R - delayed neutron precursors fields

	Info<< "\n Reading" << nPrecGroupsInit << " neutron precursor fields for region: "<< solidRegions[jj].name() << ".\n" << endl;

	for (label i=0; i<nPrecGroupsInit; i++)
	{
		label i_jjInd = nPrecGroupsInit*jj + i;

		precSolid.set
		(
			i_jjInd,
			new volScalarField
			(
				IOobject
				(
					"prec" + name(i+1),
					runTime.timeName(),
					solidRegions[jj],
					IOobject::READ_IF_PRESENT,
					IOobject::AUTO_WRITE
				),
			solidRegions[jj]
			)
		);
	}



// Read M - decay heat precursors fields

	Info<< "\n Reading" << nDecGroupsInit << " decay heat precursor fields for region: "<< solidRegions[jj].name() << ".\n" << endl;

	for (label i=0; i<nDecGroupsInit; i++)
	{
		label i_jjInd = nDecGroupsInit*jj + i;

		decSolid.set
		(
			i_jjInd,
			new volScalarField
			(
				IOobject
				(
					"dec" + name(i+1),
					runTime.timeName(),
					solidRegions[jj],
					IOobject::READ_IF_PRESENT,
					IOobject::AUTO_WRITE
				),
			solidRegions[jj]
			)
		);
	}

}
// *********************************** Create neutron transport fields  ************************************************************************* //

// Temperature-dependent properties are declared as fields, which populate corresponding pointer lists. However these are no written to files.
// Initialisation to reference values

PtrList<volScalarField> D(nEnergyGroupsInit*nRegions);
PtrList<volScalarField> Stot(nEnergyGroupsInit*nRegions);
PtrList<volScalarField> Sa(nEnergyGroupsInit*nRegions);
PtrList<volScalarField> Sf(nEnergyGroupsInit*nRegions);

PtrList<volScalarField> Ss(nEnergyGroupsInit*nEnergyGroupsInit*nRegions);

PtrList<volScalarField> Sink(nEnergyGroupsInit*nRegions); // Disappearing of neutron operator due to Scattering out from the group + absorption (==fission + capture), to be used later  <------ VERY IMPORTANT
PtrList<volScalarField> neutronSource(nEnergyGroupsInit*nRegions);



//With this loop the pointers defined in the previous lines are going to point the nuclear data contained in the scalar lists defined in "readNuclearFields.H"
forAll(solidRegions,jj) // jj count on regions
{

	Info<< "\n Creating neutron transport properties fields for region: " <<solidRegions[jj].name() << endl;

	for(label i=0; i<nEnergyGroupsInit; i++)  // i counts on groups
	{
		label i_jjInd = nEnergyGroupsInit*jj+i;
		word iName = name(i + 1);
		word jjName = name(jj+1);
		word i_jjName = jjName + "_" + iName;  //RegionIndex_GroupIndex

		D.set
		(
			i_jjInd,
			new volScalarField
			(
				IOobject
				(
					"D" + i_jjName,
					runTime.timeName(),
					solidRegions[jj],
					IOobject::NO_READ,
					IOobject::NO_WRITE
				),
				solidRegions[jj],
				D_ref[i_jjInd]
			)
		);

		Stot.set
		(
			i_jjInd,
			new volScalarField
			(
				IOobject
				(
					"Stot" + i_jjName,
					runTime.timeName(),
					solidRegions[jj],
					IOobject::NO_READ,
					IOobject::NO_WRITE
				),
				solidRegions[jj],
				Stot_ref[i_jjInd]
			)
		);

		Sa.set
		(
			i_jjInd,
			new volScalarField
			(
				IOobject
				(
					"Sa" + i_jjName,
					runTime.timeName(),
					solidRegions[jj],
					IOobject::NO_READ,
					IOobject::NO_WRITE
				),
				solidRegions[jj],
				Sa_ref[i_jjInd]
			)
		);

		Sf.set
		(
			i_jjInd,
			new volScalarField
			(
				IOobject
				(
					"Sf" + i_jjName,
					runTime.timeName(),
					solidRegions[jj],
					IOobject::NO_READ,
					IOobject::NO_WRITE
				),
				solidRegions[jj],
				Sf_ref[i_jjInd]
			)
		);

		for(label ll=0 ; ll<nEnergyGroupsInit ; ll++)
		{

			label RegionIndex4Scattering = nEnergyGroupsInit*nEnergyGroupsInit*jj;
			label EnergyIndex4Scattering = nEnergyGroupsInit * i;
			label Index = RegionIndex4Scattering + EnergyIndex4Scattering + ll;
			word jj_i_llName = name(jj + 1) + "_" + name(i+1) + "TO" + name(ll+1); // Order is: Region_StartEnergyTOArrivalEnergy


			Ss.set
			(
				Index,
				new volScalarField
				(
					IOobject
					(
						"Ss" + jj_i_llName,
						runTime.timeName(),
						solidRegions[jj],
						IOobject::NO_READ,
						IOobject::NO_WRITE
					),
					solidRegions[jj],
					Ss_ref[Index]
				)
			);
		} // end of scattering loop

		Sink.set
		(
			i_jjInd,
			new volScalarField
			(
				IOobject
				(
					"Sink" + i_jjName,
					runTime.timeName(),
					solidRegions[jj],
					IOobject::NO_READ,
					IOobject::NO_WRITE
				),
				solidRegions[jj],
				dimensionedScalar("", dimless/dimLength, 0)
			)
		);

		neutronSource.set
		(
			i_jjInd,
			new volScalarField
			(
				IOobject
				(
					"neutronSource" + i_jjName,
					runTime.timeName(),
					solidRegions[jj],
					IOobject::NO_READ,
					IOobject::NO_WRITE
				),
			solidRegions[jj],
			dimensionedScalar("", dimless/dimVolume/dimTime, 0)
			)
		);

	} //end of for loop on energy groups
} // end of loop on Regions
// ****************************************************Create other fields**************************************************************************** //


Switch Debug(runTime.controlDict().lookup("Debug"));

PtrList<volScalarField> unitField(nRegions);
PtrList<volScalarField> zeroField(nRegions);

PtrList<volScalarField> RfTot(nRegions);
PtrList<volScalarField> RiTot(nRegions);
PtrList<volScalarField> Rf(nEnergyGroupsInit*nRegions);  // Fission rate
PtrList<volScalarField> Ri(nEnergyGroupsInit*nRegions);  // Interaction rate = total interactions - fission interaction
PtrList<volScalarField> NuRfTot(nRegions);
PtrList<volScalarField> EfRfTot(nRegions);


PtrList<volScalarField> qTot(nRegions); // Total volumetric power = q_fission + q_kerma [W/m3]
PtrList<volScalarField> QTot(nRegions); // Total power = Q_fission + Q_kerma [W]

PtrList<volScalarField> q_Decay(nRegions*nDecGroupsInit); // Volumetric power due to decay per decay heat group [W/m3]

PtrList<volScalarField> q_PromptFission(nRegions*nEnergyGroupsInit); // Volumetric power due to prompt fission per energy group [W/m3]
PtrList<volScalarField> q_Kerma(nRegions*nEnergyGroupsInit); // Volumetric power due to kerma per energy group [W/m3]

forAll(solidRegions,jj)
{

// Definition of a unit dimless field to help declaration of other fields

	unitField.set
	(
		jj,
		new volScalarField
		(
			IOobject
			(
				"unitField" + name(jj+1),
				runTime.timeName(),
				solidRegions[jj],
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			solidRegions[jj],
			dimensionedScalar("", dimless, 1)
		)
	);

	zeroField.set
	(
		jj,
		new volScalarField
		(
			IOobject
			(
				"zeroField" + name(jj+1),
				runTime.timeName(),
				solidRegions[jj],
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
		solidRegions[jj],
		dimensionedScalar("", dimless, 0)
		)
	);

	RfTot.set // fission rate per unit volume and per unit time
	(
		jj,
		new volScalarField
		(
			IOobject
			(
				"RfTot",
				runTime.timeName(),
				solidRegions[jj],
				IOobject::READ_IF_PRESENT,
				IOobject::AUTO_WRITE
			),
		solidRegions[jj],
		dimensionedScalar("", dimless/dimVolume/dimTime, 0)
		)
	);

	RiTot.set // Interaction rate per unit volume and per unit time
	(
		jj,
		new volScalarField
		(
			IOobject
			(
				"RiTot",
				runTime.timeName(),
				solidRegions[jj],
				IOobject::READ_IF_PRESENT,
				IOobject::AUTO_WRITE
			),
		solidRegions[jj],
		dimensionedScalar("", dimless/dimVolume/dimTime, 0)
		)
	);

	for (label i=0; i<nEnergyGroupsInit; i++)
	{
		label i_jjInd = nEnergyGroupsInit*jj + i;

		Rf.set
		(
			i_jjInd,
			new volScalarField
			(
				IOobject
				(
					"Rf" + name(i+1),
					runTime.timeName(),
					solidRegions[jj],
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
			solidRegions[jj],
			dimensionedScalar("", dimless/dimVolume/dimTime, 0)
			)
		);

		Ri.set
		(
			i_jjInd,
			new volScalarField
			(
				IOobject
				(
					"Ri" + name(i+1),
					runTime.timeName(),
					solidRegions[jj],
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
			solidRegions[jj],
			dimensionedScalar("", dimless/dimVolume/dimTime, 0)
			)
		);
	}

	NuRfTot.set // fission rate per unit volume and per unit time
	(
		jj,
		new volScalarField
		(
			IOobject
			(
				"NuRfTot",
				runTime.timeName(),
				solidRegions[jj],
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
		solidRegions[jj],
		dimensionedScalar("", dimless/dimVolume/dimTime, 0)
		)
	);

	EfRfTot.set // fission rate per unit volume and per unit time
	(
		jj,
		new volScalarField
		(
			IOobject
			(
				"EfRfTot",
				runTime.timeName(),
				solidRegions[jj],
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
		solidRegions[jj],
		dimensionedScalar("", dimEnergy/dimVolume/dimTime, 0)
		)
	);

	qTot.set
	(
		jj,
		new volScalarField
		(
			IOobject
			(
				"qTot",
				runTime.timeName(),
				solidRegions[jj],
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
		solidRegions[jj],
		dimensionedScalar("", dimPower/dimVolume, 0)
		)
	);

	QTot.set
	(
		jj,
		new volScalarField
		(
			IOobject
			(
				"QTot",
				runTime.timeName(),
				solidRegions[jj],
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
		solidRegions[jj],
		dimensionedScalar("", dimPower, 0)
		)
	);

	for (label i=0; i<nDecGroupsInit; i++)
	{
		label i_jjInd = nDecGroupsInit*jj + i;

		q_Decay.set
		(
			i_jjInd,
			new volScalarField
			(
				IOobject
				(
					"q_Decay" + name(i+1),
					runTime.timeName(),
					solidRegions[jj],
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
			solidRegions[jj],
			dimensionedScalar("", dimPower/dimVolume, 0)
			)
		);
	}

	for (label i=0; i<nEnergyGroupsInit; i++)
	{
		label i_jjInd = nEnergyGroupsInit*jj + i;

		q_PromptFission.set
		(
			i_jjInd,
			new volScalarField
			(
				IOobject
				(
					"q_PromptFission" + name(i+1),
					runTime.timeName(),
					solidRegions[jj],
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
			solidRegions[jj],
			dimensionedScalar("", dimPower/dimVolume, 0)
			)
		);

		q_Kerma.set
		(
			i_jjInd,
			new volScalarField
			(
				IOobject
				(
					"q_Kerma" + name(i+1),
					runTime.timeName(),
					solidRegions[jj],
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
			solidRegions[jj],
			dimensionedScalar("", dimPower/dimVolume, 0)
			)
		);
	}
}

/*
    volScalarField Qprompt = Q;
    volScalarField Qdecay = Q;

    volScalarField qPrompt = q;
    volScalarField qDecay = q;


// Read values related to reactor control
    Info<< "\nReading reactor control parameters" << endl;

    Switch tempCoupling(runTime.controlDict().lookup("tempCoupling"));

    dimensionedScalar KeffInit("KeffInit", dimless, runTime.controlDict().lookupOrDefault<scalar>("Keff", 1));


    dimensionedScalar Qnominal("Qnominal", dimPower, runTime.controlDict().lookup("nominalPower"));



    volScalarField Keff
    (
        IOobject
        (
            "Keff",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        KeffInit
    );

    volScalarField reactivity
    (
        IOobject
        (
            "reactivity",
            runTime.timeName(),
            mesh
        ),
        (Keff - 1)/Keff
    );

    Info<< "    Nominal reactor power  : " << Qnominal.value()/1E+06 << " MW" << endl;
    Info<< "    Keff                   : " << max(Keff).value() << endl;
    Info<< "    Reactivity             : " << max(reactivity).value()*1e5 << " pcm" << endl;
*/
// Initialise fields
//    #include "updateCrossSections.H"
//    #include "updateFissionRate.H"
//   #include "updatePowerSource.H"
//    #include "updateNeutronSource.H"
