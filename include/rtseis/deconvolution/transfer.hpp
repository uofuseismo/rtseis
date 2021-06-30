#ifndef RTSEIS_DECONVOLUTION_TRANSFER_HPP
#define RTSEIS_DECONVOLUTION_TRANSFER_HPP 1
#include <memory>
using namespace RTSeis::Deconvolution
{
// Forward declaration
class InstrumentResponse;
/// @class Transfer transfer.hpp "rtseis/utilities/deconvolution/transfer.hpp"
/// @brief Removes and/or adds an instrument response.
/// @note This is a frequency domain based algorithm.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
template<class T = double>
class Transfer
{
public:
    /// @name Constructors
    /// @{
    /// @brief Constructor.
    Transfer()
    /// @brief Copy constructor.
    Transfer(const Transfer &transfer);
    /// @brief Move constructor.
    Transfer(Transfer &&transfer) noexcept;
    /// @}

    /// @name Destructors
    /// @{
    /// @brief Destructor
    ~Transfer();
    /// @brief Releases memory and resets the class.
    void clear() noexcept;
    /// @}

    /// @brief Initializes the transfer class.  This will remove the current
    ///        response from the data.
    /// @param[in] npts             The number of samples in the signal.
    /// @param[in] samplingRate     The sampling rate in Hz.
    /// @param[in] currentResponse  The instrument response that will be
    ///                             deconvolved from the data.
    void initialize(int npts,
                    double samplingRate,
                    const InstrumentResponse &currentResponse);
    /// @brief Initializes the transfer class. This will remove the current
    ///        response from the data then convolve in the target response.
    /// @param[in] npts             The number of samples in the signal.
    /// @param[in] samplingRate     The sampling rate in Hz.
    /// @param[in] currentResponse  The instrument response that will be 
    ///                             deconvolved from the data.
    /// @param[in] targetResponse   The instrument response to be convolved
    ///                             into the data.
    void initialize(int npts,
                    double samplingRate,
                    const InstrumentResponse &currentResponse,
                    const InstrumentResponse &targetResponse);
    /// @brief Determines if the class is initialized.
    /// @param[in] npts   The number of points in the signal.
    void apply(const int npts, const double x[], double y[]);
private:
    class TransferImpl;
    std::unique_ptr<TransferImpl> pImpl;
};
}
#endif
