#ifndef PRIVATE_PAD_HPP
#define PRIVATE_PAD_HPP
namespace
{
/// @brief Padding utility.
/// @param[in] n              The number of samples.  This must be positive.
/// @param[in] precisionSize  The size of the precision: ex: sizeof(float) = 4.
/// @param[in] alignment      The bit alignment.  This should be a power of 2.
/// @result The padded length of a row of a row major matrix so
///         that the next row begins on the desired alignment.
[[maybe_unused]]
int padLength(const int n,
                     const size_t precisionSize = sizeof(double),
                     const int alignment=64)
{
    auto size = static_cast<int> (precisionSize);
    int padLength = 0;
    auto xmod = (n*size)%alignment;
    if (xmod != 0){padLength = (alignment - xmod)/size;}
    auto nptsPadded = n + padLength;
    return nptsPadded;
}
/// @brief Padding utility for double matrices.
/// @param[in] n   The number of samples.  
/// @param[in] alignment  The bit alignment.  This should be a 
///                       multiple of 8.
/// @result The padded length of a row of a row major matrix so
///         that the next row begins on the desired alignment. 
[[maybe_unused]]
int padLength64f(const int n, const int alignment=64)
{
    return padLength(n, sizeof(double), alignment);
}
/// @brief Padding utility for float matrices.
/// @param[in] n   The number of samples.
/// @param[in] alignment  The bit alignment.  This should be a
///                       multiple of 8.
/// @result The padded length of a row of a row major matrix so
///         that the next row begins on the desired alignment.
[[maybe_unused]]
int padLength32f(const int n, const int alignment=64)
{
    return padLength(n, sizeof(float), alignment);
}
}
#endif
