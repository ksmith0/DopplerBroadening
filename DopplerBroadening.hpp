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
///
/// Doppler broadening is caused by doppler shift of an emitted \f$\gamma\f$-ray
/// due to a moving frame of reference. This shift is determined by the emitted
/// gamma energy, \f$E_\gamma\f$, the speed of the frame of reference written
/// with respect to the speed of light, \f$\beta\f$, and the polar angle of the
/// detector relative to the direction of the reference frame, \f$\theta\f$. The
/// detected energy, \f$E'\f$, is then written as 
/// \f[
///	E' = E_\gamma \frac{\sqrt{1-\beta^2}}{1-\beta cos\theta}
/// \f]
/// 
/// The resolution of a detector is dependent on the energy of the 
/// \f$\gamma\f$-ray detected. This intrinsic resolution can typically be 
/// written as
/// \f[
/// 	\frac{\delta E_{int}}{E'}(E') = \frac{const.}{\sqrt{E'}}
/// \f]
/// As the doppler shift causes the \f$\gamma\f$-ray energy to change it adds a 
/// broadening to the resolution. The total broadening can then be computed from
/// the intrinsic resolution plus the variance of the doppler shift equation 
/// computed via the sum of the partial derivatives squared. For simplification
/// one should divide the broadening by the detected energy giving the following
/// \f[
///	\frac{\delta E'}{E'} = \sqrt{
///		\left(\frac{\delta E_{int}}{E'}(E')\right)^2 +
///		\left(\frac{\partial E'}{\partial E_\gamma}\frac{\delta E_\gamma}{E'}\right)^2 +
///		\left(\frac{\partial E'}{\partial \theta}\frac{\delta \theta}{E'}\right)^2 +
///		\left(\frac{\partial E'}{\partial \beta}\frac{\delta \beta}{E'}\right)^2
///	}
/// \f]
/// Each contribution can be computed separately. The first term from the 
/// intrinsic resolution has been explained above. Computing the second term, 
/// derivative of the doppler shift equation with respect to the energy of the
/// \f$\gamma\f$-ray we find
/// \f[
///	\frac{\partial E'}{\partial E_\gamma}\frac{\delta E_\gamma}{E'} = 
///		\frac{\sqrt{1-\beta^2}}{1-\beta cos\theta}\frac{\delta E_\gamma}{E'} = 
///		\frac{\delta E_\gamma}{E_\gamma}
/// \f]
/// which is simply the uncertainty in the \f$\gamma\f$-ray emission.
///
/// The third term, the derivative of the doppler shift with respect to the 
/// polar angle can be computed as follows
/// \f[
///	\frac{\partial E'}{\partial \theta}\frac{\delta \theta}{E'} = 
///		\frac{E_\gamma \beta \sqrt{1-\beta^2} sin\theta}
///			{(1-\beta cos\theta)^2} \frac{\delta \theta}{E'} = 
///		\frac{\beta sin\theta}{(1-\beta cos\theta)}\delta \theta
/// \f]
/// 
/// Finally, the fourth term, the derivative of the doppler shift with respect
/// to the beta value can be determined as
/// \f[
///	\frac{\partial E'}{\partial \beta}\frac{\delta \beta}{E'} = 
///		\frac{E_\gamma |cos\theta - \beta|}
///			{\sqrt{1-\beta^2}(1-\beta cos\theta)^2}\frac{\delta \beta}{E'} =
///		\frac{|cos \theta - \beta|}{(1-\beta^2)(1-\beta cos\theta)}\delta\beta
/// \f]
///		
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
		/// The parameters are:
		///	* par[0] - Emitted gamma ray energy, \f$ E_\gamma \f$.
		///	* par[1] - The beta value, \f$\beta\f$, of the frame of reference.
		/// @return The shifted energy of the gamma-ray do to Doppler effect in MeV.
		/// \f[
		///	E' = E_\gamma \frac{1-\beta^2}{1-\beta cos\theta}
		/// \f]
		static Double_t DopplerShift(Double_t *ang, Double_t *par) {
			return ( par[0] * ( 1 - pow(par[1], 2) ) / 
					( 1 - par[1] * cos( ang[0] * TMath::DegToRad() ) ) );
		};

		/// @brief Computes the Doppler  broadening in resolution for the given angle 
		///         and parameters.
		/// @param[in] ang The angle value in degrees, index corresponds to power - 1.
		/// @param[in] par The parameters of the function.
		/// The parameters are:
		///	* par[0] - Emitted gamma ray energy, \f$ E_\gamma \f$.
		///	* par[1] - The beta value, \f$\beta\f$, of the frame of reference.
		///	* par[3] - The constant in \f$\frac{const.}{\sqrt{E}}\f$ term that defines 
		///		detector resolution in units of \f$\sqrt{MeV}\f$.
		/// @return The energy resolution expected at the given angle from the doppler
		///	shifted energy \f$E'\f$.
		/// \f[
		///	\frac{\delta E_{int}}{E'}(E') = \frac{const.}{\sqrt{E'}}
		/// \f]
		static Double_t EnergyBroadening(Double_t *ang, Double_t *par) {
			return ( par[3] / sqrt( DopplerShift(ang, par) ) );
		};

		/// @brief Computes the resolution broadening due to opening angle of the 
		///         detector.
		/// @param[in] ang The angle value in degrees, index corresponds to power - 1.
		/// @param[in] par The parameters of the function.
		/// The parameters are:
		///	* par[0] - Emitted gamma ray energy, \f$ E_\gamma \f$.
		///	* par[1] - The beta value, \f$\beta\f$, of the frame of reference.
		///	* par[2] - The opening angle of the detector, \f$\delta \theta\f$, in radians.
		/// @return The broadening for given parameters and angle.
		/// \f[
		///	\frac{\partial E'}{\partial \theta}\frac{\delta\theta}{E'} = 
		///		\frac{\beta sin\theta}{1-\beta cos\theta}\delta \theta
		/// \f]
		static Double_t SolidAngleBroadening(Double_t *ang, Double_t *par) {
			Double_t angleRad = ang[0] * TMath::DegToRad();
			return ( par[2] * par[1] * sin( angleRad ) / 
					( 1 - par[1] * cos ( angleRad ) ) );
		};

		/// @brief Returns the total broadening due to energy shift, opening angle and 
		///         spread in beta values.
		/// @param[in] ang The angle value in degrees, index corresponds to power - 1.
		/// @param[in] par The parameters of the function.
		/// The parameters are:
		///	* par[1] - The beta value, \f$\beta\f$, of the frame of reference.
		///	* par[4] - The width of the beta distribution, \f$\delta\beta\f$.
		/// @return The total broadening expected.
		/// \f[
		///	\frac{\partial E'}{\partial \beta}\frac{\delta\beta}{E'} =
		///		\frac{\delta\beta |cos\theta - \beta|}{(1-\beta^2)(1 - \beta cos\theta)}
		/// \f]
		static Double_t BetaBroadening(Double_t *ang, Double_t *par) {
			Double_t angleRad = ang[0] * TMath::DegToRad();
			return ( par[4] * abs( cos(angleRad) - par[1]) / 
					(( 1 - pow(par[1], 2) ) * ( 1 - par[1] * cos (angleRad) ));

		};

		/// @brief Returns the total broadening due to energy shift, opening angle and 
		///         spread in beta values.
		/// @param[in] ang The angle value in degrees, index corresponds to power - 1.
		/// @param[in] par The parameters of the function.
		/// The parameters are:
		///	* par[0] - Emitted gamma ray energy, \f$ E_\gamma \f$.
		///	* par[1] - The beta value, \f$\beta\f$, of the frame of reference.
		///	* par[2] - The opening angle of the detector, \f$\delta \theta\f$, in radians.
		///	* par[3] - The constant in \f$\frac{const.}{\sqrt{E}}\f$ term that defines 
		///		detector resolution in units of \f$\sqrt{MeV}\f$.
		///	* par[4] - The width of the beta distribution, \f$\delta\beta\f$.
		/// @return The total broadening expected.
		/// \f[
		///	\frac{\delta E'}{E'} = \sqrt{
		///		\left(\frac{\delta E_{int}}{E'}(E')\right)^2 +
		///		\left(\frac{\partial E'}{\partial \theta}\frac{\delta \theta}{E'}\right)^2 +
		///		\left(\frac{\partial E'}{\partial \beta}\frac{\delta \beta}{E'}\right)^2
		///	}
		/// \f]
		/// \note The following term has been omitted:
		/// \f[ 
		///	\left(\frac{\partial E'}{\partial E_\gamma}\frac{\delta E_\gamma}{E'}\right)^2
		/// \f]
		static Double_t TotalBroadening(Double_t *ang, Double_t *par) {
			return ( sqrt( pow(EnergyBroadening(ang, par),2) + 
						pow(SolidAngleBroadening(ang, par), 2 ) + 
						pow(BetaBroadening(ang, par), 2 ) ) );
		};


};

#endif //DOPPLERBROADENING_HPP

