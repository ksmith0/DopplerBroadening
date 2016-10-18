/// @file 
/// Contains declarations for DopplerBroadening class.
/// @author Karl Smith

#include "DopplerBroadening.hpp"
#include "TAxis.h"

/// @brief Default constructor.
/// @param[in] energyMeV The energy of the emitted gamma-ray in MeV.
/// @param[in] beta The fraction of the speed of light of incoming beam.
/// @param[in] dThetaDeg The angular coverage of the detector in degrees.
/// @param[in] resolutionConst The constant term in the 1/sqrt(e) resolution 
///             term.
/// @param[in] dBeta The change in beta. 
DopplerBroadening::DopplerBroadening(const float& energyMeV, const float& beta, const float& dThetaDeg /*= 0*/, const float& resolutionConst /*= 1*/, const float& dBeta /*= 0*/) : 
	energyMeV_(energyMeV), beta_(beta), dThetaDeg_(dThetaDeg), resolutionConst_(resolutionConst), dBeta_(dBeta)
{	
	energyBroadening_ = new TF1("energyBroadening", DopplerBroadening::EnergyBroadening, 0, 180, 5);
	solidAngleBroadening_ = new TF1("solidAngleBroadening" , DopplerBroadening::SolidAngleBroadening , 0, 180, 5);
	betaBroadening_ = new TF1("betaBroadening" , DopplerBroadening::BetaBroadening , 0, 180, 5);
	totalBroadening_ = new TF1("totalBroadening" , DopplerBroadening::TotalBroadening , 0, 180, 5);

	energyBroadening_->SetTitle("Energy Broadening");
	solidAngleBroadening_->SetTitle("Solid Angle Broadening");
	betaBroadening_->SetTitle("Beta Broadening");
	totalBroadening_->SetTitle("Total Broadening");

	UpdateParameters();
}

/// @brief Updates the parameters for each TF1 when called.
void DopplerBroadening::UpdateParameters() {
	energyBroadening_->SetParameters(energyMeV_, beta_, dThetaDeg_ * TMath::DegToRad(), resolutionConst_, dBeta_);
	solidAngleBroadening_->SetParameters(energyMeV_, beta_, dThetaDeg_ * TMath::DegToRad(), resolutionConst_, dBeta_);
	betaBroadening_->SetParameters(energyMeV_, beta_, dThetaDeg_ * TMath::DegToRad(), resolutionConst_, dBeta_);
	totalBroadening_->SetParameters(energyMeV_, beta_, dThetaDeg_ * TMath::DegToRad(), resolutionConst_, dBeta_);
}

