/**
 * \class CRTDetSim
 *
 * \ingroup crt
 *
 * \brief Provides CRTData from simulations
 *
 * Converts IDEs from largeant (or whichever producer) to 
 * CRTData. This is meant to mimic the physical detector as much as
 * possible.
 *
 *
 * \author $Author: Kevin Wierman<kevin.wierman@pnnl.gov> 
 *
 * \version $Revision: 1.0 $
 *
 * \date $Date: 2016/12/12 $
 *
 * Contact: kevin.wierman@pnnl.gov
 *
 * Created on: Tuesday, December 13, 2016
 *
**/

#ifndef CRTDetSim_HH_
#define CRTDetSim_HH_

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "lardataalg/DetectorInfo/ElecClock.h"
#include "CLHEP/Random/RandomEngine.h"
#include "fhiclcpp/ParameterSet.h"

#include <string>

namespace crt{
  class CRTDetSim :  public art:: EDProducer{
    /// Name of the producer of the IDEs
    std::string fG4ModuleLabel;
    double fGlobalT0Offset;  //!< Time delay fit: Gaussian normalization
    double fTDelayNorm;  //!< Time delay fit: Gaussian normalization
    double fTDelayShift;  //!< Time delay fit: Gaussian x shift
    double fTDelaySigma;  //!< Time delay fit: Gaussian width
    double fTDelayOffset;  //!< Time delay fit: Gaussian baseline offset
    double fTDelayRMSGausNorm;  //!< Time delay RMS fit: Gaussian normalization
    double fTDelayRMSGausShift;  //!< Time delay RMS fit: Gaussian x shift
    double fTDelayRMSGausSigma;  //!< Time delay RMS fit: Gaussian width
    double fTDelayRMSExpNorm;  //!< Time delay RMS fit: Exponential normalization
    double fTDelayRMSExpShift;  //!< Time delay RMS fit: Exponential x shift
    double fTDelayRMSExpScale;  //!< Time delay RMS fit: Exponential scale
    double fNpeScaleNorm;  //!< Npe vs. distance: 1/r^2 scale
    double fNpeScaleShift;  //!< Npe vs. distance: 1/r^2 x shift
    double fQ0;  // Average energy deposited for mips, for charge scaling [GeV]
    double fQPed;  // ADC offset for the single-peak peak mean [ADC]
    double fQSlope;  // Slope in mean ADC / Npe [ADC]
    double fQRMS;  // ADC single-pe spectrum width [ADC]
    double fQThreshold;  //!< ADC charge threshold [ADC]
    double fTResInterpolator;  // Interpolator time resolution [ns]
    double fPropDelay;  // Delay in pulse arrival time [ns/m]
    double fPropDelayError;  // Delay in pulse arrival time, uncertainty [ns/m]
    bool fUseEdep;  //!< Use the true G4 energy deposited, assume mip if false.
    bool fModelTransAtten;  //!< simplify electronics response
    bool fModelLongAtten;  //!< simplify electronics response
    double fAbsLenEff;  // Effective abs. length for transverse Npe scaling [cm]
    double fStripCoincidenceWindow;
    double fTaggerPlaneCoincidenceWindow;
    double fSipmTimeResponse;
    double fCRTClockFreq;
    bool fverbose;
    bool fSumThresh;

    /**
     * Get the channel trigger time relative to the start of the MC event.
     *
     * @param engine The random number generator engine
     * @param clock The clock to count ticks on
     * @param t0 The starting time (which delay is added to)
     * @param npe Number of observed photoelectrons
     * @param r Distance between the energy deposit and strip readout end [mm]
     * @return The channel trigger time [ns]
     */
    double getChannelTriggerTicks(CLHEP::HepRandomEngine* engine,
                                  detinfo::ElecClock& clock,
                                  float t0, float npeMean, float r);

  public:

    /// Default ctor
    CRTDetSim(const fhicl::ParameterSet&);

    /// Default dtor
    ~CRTDetSim();

    /// art::EDProducer::produce implementation
    virtual void produce (art::Event&);

    /// Set up the Configuration Parameters
    void reconfigure(fhicl::ParameterSet const & p) ;

  };
}


#endif  //CRTDetSim_HH_
