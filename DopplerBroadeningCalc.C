/// @file 
/// Contains a method that generates a plot of energy resolution as a function 
///	of angle due to doppler broadening.
/// @author Karl Smith

#include "TROOT.h"
#include "TClass.h"
#include "TPad.h"
#include "TAxis.h"

#if !defined(__CINT__)
#include "DopplerBroadening.cpp"
#endif

/// @brief Plots the broadening due to Doppler effects including energy shift, 
///	opening angle, and spread in beta values. Also, plots the individual
///	components.
/// @param[in] energyMeV The energy of the emitted gamma-ray in MeV.
/// @param[in] beta The fraction of the speed of light of incoming beam.
/// @param[in] dThetaDeg The angular coverage of the detector in degrees.
/// @param[in] resolutionConst The constant term in the 1/sqrt(e) energy 
///	resolution function. This should be in units of sqrt(MeV).
///             
/// @param[in] dBeta The width of the distribution of the beta values.
/// @return Nothing.
void DopplerBroadeningCalc(const float& energyMeV, const float& beta, const float& dThetaDeg = 0, const float resolutionConst = 1, const float& dBeta = 0) {
#ifdef __CINT__
	if (!TClass::GetDict("DopplerBroadening")) {
		gROOT->ProcessLine(".L DopplerBroadening.cpp+");
   }
#endif

	DopplerBroadening broadening(energyMeV, beta, dThetaDeg, resolutionConst, dBeta);

	TF1 *totalBroadening = broadening.GetTotalBroadening();
	totalBroadening -> SetMinimum(0);
	totalBroadening -> SetLineColor(kBlack);
	totalBroadening ->GetXaxis()->SetTitle("Angle [#circ]");
	totalBroadening ->GetYaxis()->SetTitle("Resolution [dE/E]");
	totalBroadening -> Draw("");

	//Broadening due to change in gamma energy.
	TF1 * energyBroadening = broadening.GetEnergyBroadening();
	energyBroadening -> SetMinimum(0);
	energyBroadening -> SetLineColor(kBlue);
	energyBroadening -> Draw("SAME");

	//Broadening due to solid angle coverage.
	broadening.GetSolidAngleBroadening() -> Draw("SAME");

	//Broadening due to beta value distribution.
	TF1 *betaBroadening = broadening.GetBetaBroadening();
	betaBroadening -> SetLineColor(kGreen + 1);
	betaBroadening -> Draw("SAME");

	gPad->BuildLegend(0.5,0.33,0.88,0.12);
}
