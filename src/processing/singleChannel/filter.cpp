#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <memory>
#define RTSEIS_LOGGING 1
#include "rtseis/log.h"
//#include "rtseis/postProcessing/singleChannel/lowpass.hpp"
#include "rtseis/utilities/filterRepresentations/ba.hpp"
#include "rtseis/utilities/filterRepresentations/sos.hpp"
#include "rtseis/utilities/filterRepresentations/fir.hpp"
#include "rtseis/utilities/filterImplementations/firFilter.hpp"
#include "rtseis/utilities/filterImplementations/iirFilter.hpp"
#include "rtseis/utilities/filterImplementations/sosFilter.hpp"
#include "rtseis/utilities/filterImplementations/iiriirFilter.hpp"

//namespace FR = RTSeis::Utilities::FilterRepresentations;
//using namespace RTSeis::PostProcessing::SingleChannel;

enum FilterJob
{
    FIR,
    BA,
    SOS,
    NONE
};

class FilterParameters
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Default constructor.
     */
    FilterParameters(void);
    /*!
     * @brief Initializes a digital FIR filter.
     * @param[in] fir  The FIR filter to set.
     * @param[in] lzeroPhase  If true then the filter will be applied in both
     *                        the forwards and backwards direction to the data.
     *                        If false then the filter will be applied causally.
     */ 
    FilterParameters(const RTSeis::Utilities::FilterRepresentations::FIR &fir,
                     const bool lzeroPhase = false);
    /*! 
     * @brief Initializes a digital IIR filter held in second order sections.
     * @param[in] sos   The IIR filter held as second-order-sections to set.
     * @param[in] lzeroPhase  If true then the filter will be applied in both
     *                        the forwards and backwards direction to the data.
     *                        If false then the filter will be applied causally.
     */ 
    FilterParameters(const RTSeis::Utilities::FilterRepresentations::SOS &sos,
                     const bool lzeroPhase = false);
    /*! 
     * @brief Initializes a digital IIR filter has as numerator and denominator
     *        coefficients.
     * @param[in] ba  The IIR filter held as to set.
     * @param[in] lzeroPhase  If true then the filter will be applied in both
     *                        the forwards and backwards direction to the data.
     *                        If false then the filter will be applied causally.
     */ 
    FilterParameters(const RTSeis::Utilities::FilterRepresentations::BA &ba,
                     const bool lzeroPhase = false);
    /*! @} */

    /*! @name Zero Phase Filtering
     * @{
     */
    /*!
     * @brief Enables or disables zero-phase filtering.
     * @param[in] lzeroPhase  If true then the filter will be applied in both
     *                        the forwards and backwards direction to the data.
     *                        If false then the filter will be applied causally.
     */
    void setZeroPhase(const bool lzeroPhase);
    /*!
     * @brief Determines how the filter is to be applied.
     * @result True indicates that the filter will be applied both the forward
     *         and reverse directions to the signal.
     */
    bool isZeroPhase(void) const;
    /*! @} */
    /*!
     * @brief Determines if the filter parameters are valid.
     * @result True indicates that the filter parameters are valid.
     */
    bool isValid(void) const;
private:
    class FilterParametersImpl;
    std::unique_ptr<FilterParametersImpl> pImpl;
};

class Filter
{
public:
    Filter(void);
    Filter(const FilterParameters &parms);
private:
    class FilterImpl;
    std::unique_ptr<FilterImpl> pImpl;
};

class FilterParameters::FilterParametersImpl
{
public:
    /// Default constructor
    FilterParametersImpl(void)
    {
        return;
    }
    /// Initialize FIR filter
    FilterParametersImpl(const RTSeis::Utilities::FilterRepresentations::FIR &firIn,
                         const bool lzeroPhaseIn) :
        fir(firIn),
        job(FilterJob::FIR),
        lzeroPhase(lzeroPhaseIn),
        lhaveFilter(true)
    {
         return;
    }
    /// Initialize IIR BA filter 
    FilterParametersImpl(const RTSeis::Utilities::FilterRepresentations::BA &baIn,
                         const bool lzeroPhaseIn) :
        ba(baIn),
        job(FilterJob::BA),
        lzeroPhase(lzeroPhaseIn),
        lhaveFilter(true) 
    {
        return;
    }
    /// Initialize IIR SOS filter
    FilterParametersImpl(const RTSeis::Utilities::FilterRepresentations::SOS &sosIn,
                         const bool lzeroPhaseIn) :
        sos(sosIn),
        job(FilterJob::SOS),
        lzeroPhase(lzeroPhaseIn),
        lhaveFilter(true)
    {
        return;
    }
    /// Copy operator
    FilterParametersImpl& operator=(const FilterParametersImpl &parms)
    {
        if (&parms == this){return *this;}
        ba  = parms.ba;
        fir = parms.fir;
        sos = parms.sos;
        job = parms.job;
        lhaveFilter = parms.lhaveFilter;
        return *this;
    }
public:
    class RTSeis::Utilities::FilterRepresentations::BA ba;
    class RTSeis::Utilities::FilterRepresentations::FIR fir;
    class RTSeis::Utilities::FilterRepresentations::SOS sos;
    FilterJob job = FilterJob::NONE;
    bool lzeroPhase = false;
    bool lhaveFilter = false;
};

FilterParameters::FilterParameters(void) :
    pImpl(new FilterParametersImpl())
{

}

FilterParameters::FilterParameters(
    const RTSeis::Utilities::FilterRepresentations::FIR &fir,
    const bool lzeroPhase) :
    pImpl(new FilterParametersImpl(fir, lzeroPhase))
{
    return;
}

void FilterParameters::setZeroPhase(const bool lzeroPhase)
{
    pImpl->lzeroPhase = lzeroPhase;
    return;
}

bool FilterParameters::isZeroPhase(void) const
{
    return pImpl->lzeroPhase;
}
//============================================================================//


