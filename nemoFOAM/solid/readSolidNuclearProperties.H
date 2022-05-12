// Open and define the Dictionary for Nuclear Properties (For Fluxeqns, NeutronsPrecursors and HeatPrecursors)

    IOdictionary nuclearPropertiesDictInit	//Define the reference dictionary to read the number of neutron groups, precursor groups and decay precursors;
    (
        IOobject
        (
            "nuclearProperties",
            runTime.constant(),
            solidRegions[0],
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::AUTO_WRITE
        )
    );

    dictionary neutronTransportDictInit(nuclearPropertiesDictInit.subDict("neutronTransport")); // define reference subDict inside reference Dict file

    const label nEnergyGroupsInit(readInt(neutronTransportDictInit.lookup("energyGroups")));
    const label nRegions = solidRegions.size();

// Open and define the delayed neutron precursors SubDict. Spoiler: Same number of precursor groups for all regions.

    dictionary neutronPrecursorsDictInit(nuclearPropertiesDictInit.subDict("delayedNeutronPrecursors"));

    const label nPrecGroupsInit(readInt(neutronPrecursorsDictInit.lookup("groups")));

// Open and definition of decay heat precursors SubDict. Spoiler: Same number of decay heat precursors groups for all regions.

    dictionary decayHeatPrecursorsDictInit(nuclearPropertiesDictInit.subDict("decayHeatPrecursors"));

    const label nDecGroupsInit(readInt(decayHeatPrecursorsDictInit.lookup("groups")));

// * * * * * * * * * * * * * * * * * * * * * * * * DECLARATION * * * * * * * * * * * * * * * * * * * * * * * * * * *

// We need nRegions*nEnergyGroups data, for what regards fluxEqns, despite of Scattering that is much more complicated. Let's see after what is done;

    PtrList<dimensionedScalar> DopplerControl(nRegions);
    PtrList<dimensionedScalar> rhoRefXS(nRegions);
    PtrList<dimensionedScalar> TRefXS(nRegions);


    PtrList<dimensionedScalar> D_ref(nEnergyGroupsInit*nRegions);
    PtrList<dimensionedScalar> alphaD_ref(nEnergyGroupsInit*nRegions);

    PtrList<dimensionedScalar> Stot_ref(nEnergyGroupsInit*nRegions);
    PtrList<dimensionedScalar> alphaStot_ref(nEnergyGroupsInit*nRegions);

    PtrList<dimensionedScalar> Sa_ref(nEnergyGroupsInit*nRegions);
    PtrList<dimensionedScalar> alphaSa_ref(nEnergyGroupsInit*nRegions);

    PtrList<dimensionedScalar> Sf_ref(nEnergyGroupsInit*nRegions);
    PtrList<dimensionedScalar> alphaSf_ref(nEnergyGroupsInit*nRegions);

    PtrList<dimensionedScalar> Ss_ref(nEnergyGroupsInit*nRegions*nEnergyGroupsInit);
    PtrList<dimensionedScalar> alphaSs_ref(nEnergyGroupsInit*nRegions*nEnergyGroupsInit);



    PtrList<dimensionedScalar> invVel(nEnergyGroupsInit*nRegions);
    PtrList<dimensionedScalar> Nu(nEnergyGroupsInit*nRegions);
    PtrList<dimensionedScalar> Ef(nEnergyGroupsInit*nRegions);
    PtrList<dimensionedScalar> Kerma(nEnergyGroupsInit*nRegions);
    PtrList<dimensionedScalar> chiPrompt(nEnergyGroupsInit*nRegions);
    PtrList<dimensionedScalar> chiDelayed(nEnergyGroupsInit*nRegions);


// Now we need nPrecGroups*nRegions for what regards precEqns, despite of precBetaTot that is a number per each region.

    PtrList<dimensionedScalar> precLambda(nPrecGroupsInit*nRegions);
    PtrList<dimensionedScalar> precBeta(nPrecGroupsInit*nRegions);
    PtrList<dimensionedScalar> precBetaTot(nRegions);

// Now we need nDecGroupsInit*nRegions for what regards decEqns, despite of decBetaTot that is a single number per each region.

    PtrList<dimensionedScalar> decLambda(nDecGroupsInit*nRegions);
    PtrList<dimensionedScalar> decBeta(nDecGroupsInit*nRegions);
    PtrList<dimensionedScalar> decBetaTot(nRegions);



// *******************************************************************************************************************************************************//

forAll(solidRegions,jj)
{
	word jjName = name(jj + 1);


 // ****************************** FLUXEQNS part ************************************************************************//

Info<< "\n\n*** Opening nuclear properties dict for region: " <<solidRegions[jj].name() << nl << endl;


    IOdictionary nuclearPropertiesDict	//Define the dictionary
    (
		IOobject
      (
        "nuclearProperties",
        runTime.constant(),
        solidRegions[jj],
	      IOobject::MUST_READ_IF_MODIFIED,
	      IOobject::NO_WRITE
      )
    );

 // Read neutron transport constants
    Info<< "\n***** Reading neutron transport constants" << endl;

    dictionary neutronTransportDict(nuclearPropertiesDict.subDict("neutronTransport")); // define the subDict inside the big Dict file

    label nEnergyGroups(readInt(neutronTransportDict.lookup("energyGroups")));

 // Other parameters
    scalar rhoRefXS_(readScalar(neutronTransportDict.lookup("rhoRef")));
    scalar TRefXS_(readScalar(neutronTransportDict.lookup("TRef")));

  // Diffusion parameters
    List<scalar> D_ref_(neutronTransportDict.lookup("DRef"));
    List<scalar> alphaD_ref_(neutronTransportDict.lookup("alphaDRef"));

  //Total cross sections
    List<scalar> Stot_ref_(neutronTransportDict.lookup("StotRef"));
    List<scalar> alphaStot_ref_(neutronTransportDict.lookup("alphaStotRef"));

  //Absorption cross sections
    List<scalar> Sa_ref_(neutronTransportDict.lookup("SaRef"));
    List<scalar> alphaSa_ref_(neutronTransportDict.lookup("alphaSaRef"));

  //Fission cross sections
    List<scalar> Sf_ref_(neutronTransportDict.lookup("SfRef"));
    List<scalar> alphaSf_ref_(neutronTransportDict.lookup("alphaSfRef"));

  //Scattering cross sections. These are written in "matrices" form in dictionary because they exploit the fact that i-Energy neutron can become (i->N)-Energy neutron, but the solver read them as VECTOR!!!!
    List<scalar> Ss_ref_(neutronTransportDict.lookup("SsRef"));
    List<scalar> alphaSs_ref_(neutronTransportDict.lookup("alphaSsRef"));

  //Inverse of neutron mean velocity
    List<scalar> invVel_(neutronTransportDict.lookup("inverseVelocity"));

  //Average neutrons produced by fission
    List<scalar> Nu_(neutronTransportDict.lookup("Nu"));

  //Average energy released per fission
    List<scalar> Ef_(neutronTransportDict.lookup("Ef"));

  //Average Kerma (?)
    List<scalar> Kerma_(neutronTransportDict.lookup("Kerma"));

  //Prompt neutron spectrum
    List<scalar> chiPrompt_(neutronTransportDict.lookup("promptYeld"));

  //Delayed neutron spectrum
    List<scalar> chiDelayed_(neutronTransportDict.lookup("delayedYeld"));





  // Enable Doppler effect switch for cross sections update
    Switch DopplerControl_(neutronTransportDict.lookup("useDoppler"));

    if (!DopplerControl_)
    {
        DopplerControl.set(jj, new dimensionedScalar("DopplerControl" + jjName, dimless, 0));
    }
    else
    {
        DopplerControl.set(jj, new dimensionedScalar("DopplerControl" + jjName, dimless, 1));
    }

    rhoRefXS.set(jj, new dimensionedScalar("rhoRefXS" + jjName, dimDensity, rhoRefXS_));
    TRefXS.set(jj, new dimensionedScalar("TRefXS" + jjName, dimTemperature, TRefXS_));


  // Control on input data size, made only on diffusion parameters
    if( D_ref_.size() == nEnergyGroupsInit )
    {
        Info<< "    Correctly found " << nEnergyGroups << " energy groups." << endl;
    }
    else
    {
        FatalError
            << "In " << neutronTransportDict.name() << " of region: "<< solidRegions[jj].name() << ", invalid number of neutron transport constants. Check energy groups value. " << exit(FatalError);
    }



    // Populate pointer lists with dimensioned values

    //Remember the sintax used here: X_:list of scalar defined previously. X : list of pointer (up to now empty)

    forAll(D_ref_, i)  // Remember the Foam sintax: forAll(list, j) is a shortcut to the classical for(int i=0; i<(list).size(); i++). j is automatically a label (=integer).
    {
        label i_jjInd = nEnergyGroupsInit*jj + i;  // jj is Region counter, i is energy counter
        word i_jjName = name(jj + 1) + "_" + name(i + 1);

        D_ref.set(i_jjInd, new dimensionedScalar("D" + i_jjName + "_ref", dimLength, D_ref_[i]));
        Stot_ref.set(i_jjInd, new dimensionedScalar("Stot" + i_jjName + "_ref", dimless/dimLength, Stot_ref_[i]));
        Sa_ref.set(i_jjInd, new dimensionedScalar("Sa" + i_jjName + "_ref", dimless/dimLength, Sa_ref_[i]));
        Sf_ref.set(i_jjInd, new dimensionedScalar("Sf" + i_jjName + "_ref", dimless/dimLength, Sf_ref_[i]));

        alphaD_ref.set(i_jjInd, new dimensionedScalar("alphaD" + i_jjName + "_ref", dimLength, alphaD_ref_[i]));
        alphaStot_ref.set(i_jjInd, new dimensionedScalar("alphaStot" + i_jjName + "_ref", dimless/dimLength, alphaStot_ref_[i]));
        alphaSa_ref.set(i_jjInd, new dimensionedScalar("alphaSa" + i_jjName + "_ref", dimless/dimLength, alphaSa_ref_[i]));
        alphaSf_ref.set(i_jjInd, new dimensionedScalar("alphaSf" + i_jjName + "_ref", dimless/dimLength, alphaSf_ref_[i]));


        forAll(D_ref_, ll)
        {
            label RegionIndex4Scattering = nEnergyGroupsInit*nEnergyGroupsInit*jj;
            label EnergyIndex4Scattering = nEnergyGroupsInit * i;
            label Index = RegionIndex4Scattering + EnergyIndex4Scattering + ll;
            label ReadingSsIndex = EnergyIndex4Scattering + ll;
            word jj_i_llName = name(jj + 1) + "_" + name(i+1) + "TO" + name(ll+1); // Order is: Region_StartEnergyTOArrivalEnergy

            Ss_ref.set(Index, new dimensionedScalar("Ss" + jj_i_llName + "_ref", dimless/dimLength, Ss_ref_[ReadingSsIndex]));	//necessary trick to use a vector representation for Ss
            alphaSs_ref.set(Index, new dimensionedScalar("alphaSs" + jj_i_llName + "_ref", dimless/dimLength, alphaSs_ref_[ReadingSsIndex]));
        }

        Nu.set(i_jjInd, new dimensionedScalar("Nu" + i_jjName, dimless, Nu_[i]));
        Ef.set(i_jjInd, new dimensionedScalar("Ef" + name(i + 1), dimEnergy, Ef_[i]));
        Kerma.set(i_jjInd, new dimensionedScalar("Kerma" + name(i + 1), dimEnergy/dimLength, Kerma_[i]));
        invVel.set(i_jjInd, new dimensionedScalar("invVel" + i_jjName, dimTime/dimLength, invVel_[i]));
        chiPrompt.set(i_jjInd, new dimensionedScalar("chiPrompt" + i_jjName, dimless, chiPrompt_[i]));
        chiDelayed.set(i_jjInd, new dimensionedScalar("chiDelayed" + name(i + 1), dimless, chiDelayed_[i]));
    }


 // ****************************** precEQNS part ********************************************************************************************************//


Info<< "\nReading delayed neutron precursors constants for region: " << solidRegions[jj].name() << nl << endl;

    dictionary neutronPrecursorsDict(nuclearPropertiesDict.subDict("delayedNeutronPrecursors"));

    label nPrecGroups(readInt(neutronPrecursorsDict.lookup("groups")));
    List<scalar> precLambda_(neutronPrecursorsDict.lookup("lambda"));;
    List<scalar> precBeta_(neutronPrecursorsDict.lookup("beta"));
    scalar precBetaTot_=0;


    // Check input data size
    if( precLambda_.size() == nPrecGroups && precBeta_.size() == nPrecGroups && nPrecGroups == nPrecGroupsInit )
    {
        Info<< "\n   In region: " << solidRegions[jj].name() <<" found " << nPrecGroups << " groups of delayed neutron precursors." << endl;
    }
    else
    {
        FatalError
            << "in " << neutronPrecursorsDict.name() << " of region: "<< solidRegions[jj].name() << " , invalid number of delayed neutron precursors constants." << exit(FatalError);
    }

    // Assign dimensioned values to pointer list
    forAll(precLambda_, i)
    {
	label i_jjInd = nPrecGroupsInit*jj + i;
	word i_jjName = name(jj + 1) + "_" + name(i + 1);  //RegionIndex_precIndex


        precLambda.set(i_jjInd, new dimensionedScalar("precLambda" + i_jjName, dimless/dimTime, precLambda_[i]));
        precBeta.set(i_jjInd, new dimensionedScalar("precBeta" + i_jjName, dimless, precBeta_[i]));
        precBetaTot_ += precBeta[i].value();
    }

    precBetaTot.set(jj, new dimensionedScalar("precBetaTot" + name(jj+1), dimless, precBetaTot_ ));

 // ****************************** decEQNS part ********************************************************************************************************//

// Read decay heat precursors constants

	Info<< "\nReading decay heat precursors constants for region: " << solidRegions[jj].name() << nl << endl;

	dictionary decayHeatPrecursorsDict(nuclearPropertiesDict.subDict("decayHeatPrecursors"));

	label nDecGroups(readInt(decayHeatPrecursorsDict.lookup("groups")));

	// Scalar lists to lookup dictionary
	List<scalar> decLambda_(decayHeatPrecursorsDict.lookup("lambda"));
	List<scalar> decBeta_(decayHeatPrecursorsDict.lookup("beta"));
	scalar decBetaTot_ = 0;


    // Check input data size
	if( decLambda_.size() == nDecGroups && decBeta_.size() == nDecGroups && nDecGroups == nDecGroupsInit )
	{
		Info<< "\n   In region: " << solidRegions[jj].name() <<" found " << nDecGroups << " groups of delayed neutron precursors." << endl;
	}
	else
	{
		FatalError
		  << "in " << decayHeatPrecursorsDict.name() << " of region: "<< solidRegions[jj].name() << ", invalid number of decay heat precursors constants." << exit(FatalError);
	}

	// Assign dimensioned values to pointer list
	forAll(decLambda_, i)
	{
		label i_jjInd = nDecGroupsInit*jj + i;
		word i_jjName = name(jj + 1) + "_" + name(i + 1);  //RegionIndex_precIndex


		decLambda.set(i_jjInd, new dimensionedScalar("decLambda" + i_jjName, dimless/dimTime, decLambda_[i]));
		decBeta.set(i_jjInd, new dimensionedScalar("decBeta" + i_jjName, dimless, decBeta_[i]));
		decBetaTot_ += decBeta[i].value();
	}

	decBetaTot.set(jj, new dimensionedScalar("decBetaTot" + name(jj+1), dimless, decBetaTot_ ));
}




Info<< "\n End of readSolidNuclearProperties" << endl;


Info<< "Read Nuclear field correct ended." << endl;
