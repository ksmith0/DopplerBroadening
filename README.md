Doppler Broadening
===
The project provide the DopplerBroadening class as well as a simple method to create plots. This method can be envoked from root by the following:
```
.x DopplerBroadeningCalc.C(energyMeV, beta, dThetaDeg, resolutionConst, dBeta)
```
where the arguments are as follows:
 * energyMeV The energy of the emitted gamma-ray in MeV.
 * beta The fraction of the speed of light of incoming beam.
 * dThetaDeg The angular coverage of the detector in degrees.
 * resolutionConst The constant term in the 1/sqrt(e) resolution term.
 * dBeta The change in beta. 

Further documentation is available at https://ksmith0.github.io/DopplerBroadening.
