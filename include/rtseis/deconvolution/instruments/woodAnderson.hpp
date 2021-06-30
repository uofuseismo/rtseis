#ifndef RTSEIS_DECONVOLUTION_INSTRUMENTS_WOODANDERSON_HPP
#define RTSEIS_DECONVOLUTION_INSTRUMENTS_WOODANDERSON_HPP
#include <memory>
namespace RTSeis::FilterRepresentations
{
class BA;
}
namespace RTSeis::Deconvolution::Instruments
{
/// @brief Defines a Wood Anderson instrument response.
///        This assumes a constant of constant of 2800 and 
class WoodAnderson
{
public:
    /// @result The transfer function for the Wood Anderson instrument. 
    [[nodiscard]]
    static RTSeis::FilterRepresentations::BA getTransferFunction() noexcept;
};
}
#endif
