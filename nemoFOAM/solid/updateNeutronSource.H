{
if (Debug)  {Info << "\nupdateNeutronSource: START" << endl;}
  for(label ii = 0 ; ii<nEnergyGroupsInit; ii++)
   {
	label posEqn = ii+i*nEnergyGroupsInit;

     // Initialize the SourceTerms as zero fields
        volScalarField tScattering_NS = zeroField[i]*dimensionedScalar("", dimless/dimVolume/dimTime, 0); //These are the neutrons that enter into the i-th group from more energetic groups due to scattering
		if (Debug)  {Info << "\n *** tScattering["<<posEqn<<"] = zeroField["<<i<<"]" << endl;}
		if (Debug)  {Info << "\n *** tScattering["<<posEqn<<"] = "<< tScattering_NS << endl;}

	volScalarField tFission_NS = zeroField[i]*dimensionedScalar("", dimless/dimVolume/dimTime, 0); //These are the neutrons that enter into the i-th group due to fission produced by neutrons not in group i-th
		if (Debug)  {Info << "\n *** tFission["<<posEqn<<"] = zeroField["<<i<<"]" << endl;}

	volScalarField tDelayed_NS = zeroField[i]*dimensionedScalar("", dimless/dimVolume/dimTime, 0);  //These are the delayed neutrons that enter into the i-th group due to emission of precursors
		if (Debug)  {Info << "\n *** tDelayed["<<posEqn<<"] = zeroField["<<i<<"]" << endl;}

     // Evaluate neutron sources based on current fluxes and cross sections
	for(label iii = 0; iii<nEnergyGroupsInit; iii++)
	{
		label Pos = iii + i*nEnergyGroupsInit;
			if (Debug)  {Info << "\n *** Pos = "<< Pos << endl;}

		if (ii!=iii)
		{
			if (Debug)  {Info << "\n *** i = "<< i << "; ii = "<< ii << "; iii = " << iii << ". Region is: "<<solidRegions[i].name() << endl;}
			if (Debug)  {Info << "\n *** posEqn = " << posEqn << ";  Pos = " << Pos << endl;}

			tFission_NS += chiPrompt[posEqn]*(1 - precBetaTot[i])*Nu[Pos]*Sf[Pos]*fluxSolid[Pos]; // Should be OK
				if (Debug)  {Info << "\n *** tFission["<<posEqn<<"] += chiPrompt["<<posEqn<<"]*(1 - precBetaTot["<<i<<"])*Nu["<<Pos<<"]*Sf["<<Pos<<"]*fluxSolid["<<Pos<<"] " <<endl;}

	       		tScattering_NS += Ss[i*nEnergyGroupsInit*nEnergyGroupsInit+ii+nEnergyGroupsInit*iii]*fluxSolid[Pos];  // Should be OK
				if (Debug)  {Info << "\n *** tScattering["<<posEqn<<"] += Ss["<<i*nEnergyGroupsInit*nEnergyGroupsInit+ii+nEnergyGroupsInit*iii<<"]*fluxSolid["<<Pos<<"];" <<endl;}
				if (Debug)  {Info << "\n *** tScattering["<<posEqn<<"] = "<< tScattering_NS << endl;}
		}
	}
	if (Debug)  {Info <<"\n ** Fission source & Scattering Source updated."<< endl;}


	for(label iii = 0; iii<nPrecGroupsInit; iii++)
	{
		label Pos = iii + i*nPrecGroupsInit;
			if (Debug)  {Info << "\n *** Pos_prec = "<< Pos << endl;}

		tDelayed_NS += chiDelayed[posEqn]*precLambda[Pos]*precSolid[Pos];
			if (Debug)  {Info << "\n *** tDelayed_NS["<< posEqn <<"] += chiDelayed["<<posEqn<<"]*precLambda["<<Pos<<"]*precSolid["<<Pos<<"] "<< endl;}
	}
	if (Debug)  {Info <<"\n ** Delayed source updated."<< endl;}

    // Repopulate effectively neutronSource
	neutronSource[posEqn] = tScattering_NS + tFission_NS + tDelayed_NS;

//	if(posEqn == 0)   {Info << "\n * neutronSource["<<posEqn<<"] = "<< neutronSource[posEqn] << endl;}
		if (Debug)  {Info << "\n * neutronSource["<<posEqn<<"] correctly updated" <<endl;}
		if (Debug)  {Info << "\n * neutronSource["<<posEqn<<"] = "<< neutronSource[posEqn] << endl;}
    }

if (Debug)  {Info << "\n updateNeutronSource: END" << endl;}
}
