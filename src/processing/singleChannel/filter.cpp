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

enum FilterImplementation
{
    FIR,    // FIR filtering
    BA,     // IIR filtering
    SOS,    // SOS IIR filtering
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
     * @brief Copy constructor.
     * @param[in] parms  The filter parameters class from which to intiialize.
     */
    FilterParameters(const FilterParameters &parms);
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

    /*!
     * @brief Copy operator.
     * @param[in] parms  The filter parameters to copy.
     * @result A deep copy of the filter parameters.
     */
     FilterParameters& operator=(const FilterParameters &parms);

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Default destructor.
     */
    ~FilterParameters(void);
    /*!
     * @brief Resets the filter parameters.
     */
    void clear(void);
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
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Default constructor.  Note, this class can not yet be used to
     *        filter data until the parameters have been set with
     *        \c setParameters().
     */
    Filter(void);
    /*!
     * @brief Constructs the filter class from the parameters.
     * @param[in] parms  The filter parameters from which to construct
     *                   the class.
     * @throws std::invalid_argument if parms is invalid.
     */
    Filter(const FilterParameters &parms);
    /*!
     * @brief Copy constructor.
     * @param[in] filter  Filter to class from which to construct this class.
     */ 
    Filter(const Filter &filter);
    /*! @} */
   
    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy assignment operator.
     * @param[in] filter  The filter class to copy.
     * @result A deep copy of the filter.
     */
    Filter &operator=(const Filter &filter);
    /*! @} */

    /*! @name Destructors.
     * @{
     */
    /*!
     * @brief Default destructor.
     */
    ~Filter(void);
    /*!
     * @brief Resets the filters and clears the parameters.
     */
    void clear(void);
    /*! @} */

    /*!
     * @brief Sets the filtering parameters.
     * @param[in] parms  The filter parameters to set.
     * @throws std::invalid_argument if parms is invalid.
     */
    void setParameters(const FilterParameters &parms);
    /*!
     * @brief Applies the filter to the signal.
     * 
     * @throws std::invalid_argument if the class is not yet set.
     */
    void apply(const std::vector<double> &x, std::vector<double> &y);
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
        job(FilterImplementation::FIR),
        lzeroPhase(lzeroPhaseIn),
        lhaveFilter(true)
    {
         return;
    }
    /// Initialize IIR BA filter 
    FilterParametersImpl(const RTSeis::Utilities::FilterRepresentations::BA &baIn,
                         const bool lzeroPhaseIn) :
        ba(baIn),
        job(FilterImplementation::BA),
        lzeroPhase(lzeroPhaseIn),
        lhaveFilter(true) 
    {
        return;
    }
    /// Initialize IIR SOS filter
    FilterParametersImpl(const RTSeis::Utilities::FilterRepresentations::SOS &sosIn,
                         const bool lzeroPhaseIn) :
        sos(sosIn),
        job(FilterImplementation::SOS),
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
    FilterImplementation getImplementation(void) const
    {
        return job;
    }
public:
    class RTSeis::Utilities::FilterRepresentations::BA ba;
    class RTSeis::Utilities::FilterRepresentations::FIR fir;
    class RTSeis::Utilities::FilterRepresentations::SOS sos;
    FilterImplementation job = FilterImplementation::NONE;
    bool lzeroPhase = false;
    bool lhaveFilter = false;
};

class Filter::FilterImpl
{
public:
    FilterImpl(void){return;}
    FilterParameters parms;
    RTSeis::Utilities::FilterImplementations::FIRFilter fir;
    RTSeis::Utilities::FilterImplementations::IIRFilter iir;
    RTSeis::Utilities::FilterImplementations::IIRIIRFilter iiriir;
    RTSeis::Utilities::FilterImplementations::SOSFilter sos;
    bool linit = false;
};

FilterParameters::FilterParameters(void) :
    pImpl(new FilterParametersImpl())
{
    return;
}

FilterParameters::FilterParameters(
    const RTSeis::Utilities::FilterRepresentations::FIR &fir,
    const bool lzeroPhase) :
    pImpl(new FilterParametersImpl(fir, lzeroPhase))
{
    return;
}

FilterParameters::FilterParameters(const FilterParameters &parms)
{
    *this = parms;
}

FilterParameters& FilterParameters::operator=(const FilterParameters &parms)
{
    if (&parms == this){return *this;}
    if (pImpl){clear();}
    pImpl = std::unique_ptr<FilterParametersImpl> (new FilterParametersImpl());
    pImpl->ba = parms.pImpl->ba;
    pImpl->fir = parms.pImpl->fir;
    pImpl->sos = parms.pImpl->sos;
    pImpl->job = parms.pImpl->job;
    pImpl->lzeroPhase = parms.pImpl->lzeroPhase;
    pImpl->lhaveFilter = parms.pImpl->lhaveFilter;
    return *this; 
}

FilterParameters::~FilterParameters(void)
{
    clear();
    return;
}

void FilterParameters::clear(void)
{
    if (pImpl)
    {
        pImpl->ba.clear();
        pImpl->fir.clear();
        pImpl->sos.clear();
        pImpl->job = FilterImplementation::NONE;
        pImpl->lzeroPhase = false;
        pImpl->lhaveFilter = false;
    }
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

bool FilterParameters::isValid(void) const
{
    if (!pImpl->lhaveFilter){return false;}
    // N.B. I don't need to scrutinize the filter coefficients, e.g.,
    // as[3*i] = 0 because, for the filters to have non-zero length, then
    // the calling applicaiton must have passed the internal checks in the
    // respective classes.
    if (pImpl->job == FilterImplementation::SOS)
    {
        if (pImpl->sos.getNumberOfSections() < 1){return false;}
    }
    else if (pImpl->job == FilterImplementation::FIR)
    {
        if (pImpl->fir.getNumberOfFilterTaps() < 1){return false;}
    }
    else if (pImpl->job == FilterImplementation::BA)
    {
        if (pImpl->ba.getNumberOfNumeratorCoefficients() < 1 ||
            pImpl->ba.getNumberOfDenominatorCoefficients() < 1)
        {
            return false;
        }
    }
    else
    {
        return false;
    }
    return true;
}

//============================================================================//


Filter::Filter(void) :
    pImpl(new FilterImpl())
{
    return;
}

Filter& Filter::operator=(const Filter &filter)
{
    if (&filter == this){return *this;}
    pImpl->fir    = filter.pImpl->fir;
    pImpl->iir    = filter.pImpl->iir;
    pImpl->iiriir = filter.pImpl->iiriir;
    pImpl->sos    = filter.pImpl->sos;
    pImpl->linit  = filter.pImpl->linit;
    return *this;
}

Filter::~Filter(void)
{
    clear();
}

void Filter::clear(void)
{
    if (pImpl)
    {
        pImpl->parms.clear();
        pImpl->fir.clear();
        pImpl->iir.clear();
        pImpl->iiriir.clear();
        pImpl->sos.clear();
        pImpl->linit = false;
    }
    return;
}

void Filter::setParameters(const FilterParameters &parms)
{
    clear();
    if (!parms.isValid())
    {
        throw std::invalid_argument("Parameters are not valid");
    }
    // Initialize the filter
    pImpl->parms = parms;
    FilterImplementation job;// = pImpl->parms->pImpl.getImplementation();
    if (pImpl->parms.isZeroPhase() && job == FilterImplementation::BA)
    {
        
    } 
    return;
}

void Filter::apply(const std::vector<double> &x, std::vector<double> &y)
{
    // Nothing to do
    if (x.size() == 0)
    {
        y.resize(0);
        return;
    }
    if (!pImpl->linit)
    {
        y.resize(0);
        throw std::invalid_argument("Filter was never set");
    }
    // Zero-phase filtering is its own animal
    return;
}

