/// @file 
/// Contains header information for DopplerBroadening class.
/// @author Karl Smith

#ifndef DOPPLERBROADENING_HPP
#define DOPPLERBROADENING_HPP

#include "TMath.h"
#include "TF1.h"

/// @brief Produces a number of TF1 objects describing the change in energy
///	resolution due to energy shift, opening angle and change in beta value.
/// @author Karl Smith
class DopplerBroadening {
	public:
		/// @brief Default Constructor
		/// @param[in] energyMeV The energy of the emitted gamma-ray in MeV.
		/// @param[in] beta The fraction of the speed of light of incoming beam.
		/// @param[in] dThetaDeg The angular coverage of the detector in degrees.
		/// @param[in] resolutionConst The constant term in the 1/sqrt(e) resolution 
		///             term.
		/// @param[in] dBeta The change in beta. 
		DopplerBroadening(const float& energyMeV, const float& beta, const float& dThetaDeg = 0, const float& resolutionConst = 1, const float& dBeta = 0);

		/// @brief Returns a TF1 describing the change in energy resolution due 
		///	to gamma energy doppler shift.
		/// @return A TF1 object with the change in energy resolution due to 
		///	energy shift as a function of polar angle \f$\theta\f$ in degrees.
		TF1* GetEnergyBroadening() { return energyBroadening_; };
		
		/// @brief Returns a TF1 describing the change in energy resolution due 
		///	to the solid angle coverage of the detector.
		/// @return A TF1 object with the change in energy resolution due to 
		///	opening angle of the detector as a function of polar angle 
		///	\f$\theta\f$ in degrees.
		TF1* GetSolidAngleBroadening() { return solidAngleBroadening_; };

		/// @brief Returns a TF1 describing the change in energy resolution due
		///	to the distribution of beta values.
		/// @return A TF1 object with the change in energy resolution due to 
		///	beta distribution width as a function of polar angle \f$\theta\f$
		///	in degrees.
		TF1* GetBetaBroadening() { return betaBroadening_; };

		/// @brief Returns a TF1 describing the change in energy resolution due
		///	to gamma energy doppler shift, the solid angle coverage of the 
		///	detector and the distribution of beta values.
		/// @return A TF1 object with the change in energy resolution due to 
		///	energy shift, opening angle of the detector, and beta distribution 
		///	width as a function of polar angle \f$\theta\f$ in degrees.
		TF1* GetTotalBroadening() { return totalBroadening_; };

	private:
		/// @brief Updates the parameters for each TF1 when called.
		void UpdateParameters();

		float energyMeV_; ///< Energy of the gamma-ray in MeV.
		float beta_; ///< Beta value, ratio of velocity to speed of light.
		float dThetaDeg_; ///< Opening angle of the detector.
		float resolutionConst_; ///< Constant term in \f$ \frac{1}{\sqrt{E}} \f$ of resolution function.
		float dBeta_; ///< Width of beta distribtuion.

		TF1 *energyBroadening_; ///< Pointer to energy broadening TF1.
		TF1 *solidAngleBroadening_;///< Pointer to solid angle broadening TF1.
		TF1 *betaBroadening_; ///< Pointer to beta distribution broadening TF1.
		TF1 *totalBroadening_; ///< Pointer to total broadening TF1.

		/// @brief Computes the Doppler shift for the given angle and parameters.
		/// @param[in] ang The angle value in degrees, index corresponds to power - 1.
		/// @param[in] par The parameters of the function.
		/// @return The shifted energy of the gamma-ray do to Doppler effect in MeV.
		/// The parameters are:
		/// 	0: Emitted gamma ray energy.
		///	1: The beta value of the frame of reference.
		static Double_t DopplerShift(Double_t *ang, Double_t *par) {
			return ( par[0] * ( 1 - pow(par[1], 2) ) / 
					( 1 - par[1] * cos( ang[0] * TMath::DegToRad() ) ) );
		};

		/// @brief Computes the Doppler  broadening in resolution for the given angle 
		///         and parameters.
		/// @param[in] ang The angle value in degrees, index corresponds to power - 1.
		/// @param[in] par The parameters of the function.
		/// @return The energy resolution expected at the given angle.
		/// The parameters are:
		/// 	0: Emitted gamma ray energy.
		///	1: The beta value of the frame of reference.
		///	3: The constant in \f$\frac{1}{\sqrt{E}}\f$ term that defines 
		///		detector resolution in units of \f$\sqrt{MeV}\f$.
		static Double_t EnergyBroadening(Double_t *ang, Double_t *par) {
			return ( par[3] / sqrt( DopplerShift(ang, par) ) );
		};

		/// @brief Computes the resolution broadening due to opening angle of the 
		///         detector.
		/// @param[in] ang The angle value in degrees, index corresponds to power - 1.
		/// @param[in] par The parameters of the function.
		/// @return The broadening for given parameters and angle.
		/// The parameters are:
		/// 	0: Emitted gamma ray energy.
		///	1: The beta value of the frame of reference.
		///	2: The opening angle of the detector in radians.
		static Double_t SolidAngleBroadening(Double_t *ang, Double_t *par) {
			Double_t angleRad = ang[0] * TMath::DegToRad();
			return ( par[2] * par[1] * sin( angleRad ) / 
					( 1 - par[1] * cos ( angleRad ) ) );
		};

		/// @brief Returns the total broadening due to energy shift, opening angle and 
		///         spread in beta values.
		/// @param[in] ang The angle value in degrees, index corresponds to power - 1.
		/// @param[in] par The parameters of the function.
		/// @return The total broadening expected.
		/// The parameters are:
		///	1: The beta value of the frame of reference.
		///	4: The width of the beta distribution.
		static Double_t BetaBroadening(Double_t *ang, Double_t *par) {
			Double_t angleRad = ang[0] * TMath::DegToRad();
			return ( par[4] * abs( cos(angleRad) - par[1]) / 
					( 1 - par[1] * cos (angleRad) ) /
					( 1 - pow(par[1], 2) ) );

		};

		/// @brief Returns the total broadening due to energy shift, opening angle and 
		///         spread in beta values.
		/// @param[in] ang The angle value in degrees, index corresponds to power - 1.
		/// @param[in] par The parameters of the function.
		/// @return The total broadening expected.
		static Double_t TotalBroadening(Double_t *ang, Double_t *par) {
			return ( sqrt( pow(EnergyBroadening(ang, par),2) + 
						pow(SolidAngleBroadening(ang, par), 2 ) + 
						pow(BetaBroadening(ang, par), 2 ) ) );
		};


};

#endif //DOPPLERBROADENING_HPP

