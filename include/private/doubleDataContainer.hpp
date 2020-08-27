#ifndef RTSEIS_PRIVATE_VECTOR_HPP
#define RTSEIS_PRIVATE_VECTOR_HPP 1
#include <memory>
#include <ipps.h>

namespace RTSeis
{
namespace Private
{

/*!
 * @brief Container for the filtered data.
 */
class DoubleDataContainer
{
public:
    /*! @name Constructors
     * @{
     */
    /*!
     * @brief Default constructor.
     */
    DoubleDataContainer(void)
    {
        return;
    }
    /*!
     * @brief Copy constructor.
     * @param[in] data  Data container from which to initialize this class.
     * @note The data pointer, if set, cannot be copied.
     */
    DoubleDataContainer(const DoubleDataContainer &data)
    {
        *this = data;
    }
    /*!
     * @brief Move constructor.
     * @param[in,out] data  On input the data container whose memory is to be
     *                      moved to this class.
     * @param[in,out] data  On exit this is no longer a valid data class.
     */
    DoubleDataContainer(DoubleDataContainer &&data)
    {
        *this = std::move(data);
    }

    /*! @name Operators
     * @{
     */
    /*!
     * @brief Copy operator.
     * @param[in] data  Data container to copy.
     * @result A deep copy of the input data container.  Note, the data pointer
     *         cannot be copied.  Thus, the input data retains ownership.
     */
    DoubleDataContainer& operator=(const DoubleDataContainer &data)
    {
        if (&data == this){return *this;}
        maxx_ = data.maxx_;
        nx_ = data.nx_;
        maxy_ = data.maxy_;
        ny_ = data.ny_;
        if (maxx_ > 0)
        {
            x_ = ippsMalloc_64f(maxx_);
            ippsCopy_64f(data.x_, x_, maxx_);
        }
        if (maxy_ > 0)
        {
            y_ = ippsMalloc_64f(maxy_);
            ippsCopy_64f(data.y_, y_, maxy_);
        }
        //xptr_ = std::move(data.xptr_;
        return *this;
    }
    /*!
     * @brief Move operator.
     * @param[in,out] data  Data container to move.  On exit, data's behavior
     *                      will be undefined.
     * @result A deep copy of the input data container.
     */
    DoubleDataContainer& operator=(DoubleDataContainer &&data)
    {
        maxx_ = data.maxx_;
        nx_ = data.nx_;
        maxy_ = data.maxy_;
        ny_ = data.ny_;
        x_ = std::move(data.x_);
        y_ = std::move(data.y_);
        if (data.xptr_){xptr_ = std::move(data.xptr_);}
        // Clean up
        data.x_ = 0;
        data.y_ = 0;
        data.xptr_.release(); // = nullptr;
        data.maxx_ = 0;
        data.nx_ = 0;
        data.maxy_ = 0;
        data.ny_ = 0;
    }
    /*! @} */

    /*! @name Destructors
     * @{
     */
    /*!
     * @brief Default destructor.
     */
    ~DoubleDataContainer(void)
    {
        clear();
    }
    /*!
     * @brief Resets the data container and releases memory.
     */
    void clear(void)
    {
        if (x_){ippsFree(x_);}
        if (y_){ippsFree(y_);}
        x_ = nullptr;
        y_ = nullptr;
        xptr_.release(); // = nullptr;
        maxx_ = 0;
        nx_ = 0;
        maxy_ = 0;
        ny_ = 0;
        return;
    }
    /*! @} */


    /*!
     * @brief Gets the number of samples in the input signal.
     * @result The number of samples in the input signal.
     */
    int getNumberOfInputSamples(void) const noexcept
    {
        return nx_;
    }

    /*!
     * @brief Gets the number of samples in the output (filtered) signal.
     * @result The number of samples in the output signal.
     */
    int getNumberOfOutputSamples(void) const noexcept
    {
        return ny_;
    }

    /*!
     * @brief Sets a pointer to the input data.
     * @param[in] nx  The number of samples in x.
     * @param[in] x   The data pointer to set.
     */
    void setInputDataPointer(const int nx, std::unique_ptr<const double> x) noexcept
    {
        xptr_.release(); // = nullptr;
        // Set pointer if possible
        if (nx > 0 && x)
        {
            xptr_ = std::move(x);
            //xptr_ = x;
            nx_ = nx;
        }
        else
        {
            if (nx_ > 0){RTSEIS_ERRMSG("%s", "x is NULL");}
            nx_ = 0;
        }
        return;
    }
    /*!
     * @brief Releases the input data pointer.
     * @note This will also set the number of input samples to 0.
     */ 
    void releaseInputDataPointer(void)
    {
        nx_ = 0;
        xptr_.release(); // = nullptr;
        return;
    }
    /*!
     * @brief Sets the input data.
     * @param[in] nx  The number of samples in x.
     * @param[in] x   The signal to filter.  This is an array of dimension [nx]. 
     */
    void setInputData(const int nx, const double x[]) noexcept
    {
        xptr_.release();// = nullptr; // Invalidate old pointer
        // Resize x
        if (nx > maxx_)
        {
            if (x_){ippsFree(x_);}
            maxx_ = nx;
            x_ = ippsMalloc_64f(maxx_);
        }
        // Copy if possible
        if (nx > 0 && x)
        {
            nx_ = nx;
            ippsCopy_64f(x, x_, nx_);
        }
        else
        {
            if (nx_ > 0){RTSEIS_ERRMSG("%s", "x is NULL");}
            nx_ = 0;
        }
        return;
    }
    /*!
     * @brief Gets a pointer to the input data.
     * @result If this is a nullptr then there is no input data.
     * @result Otherwise, this is a pointer to the data to filter whose
     *         length can be determined by \c getNumberOfInputSamples().
     * @note If this class goes out of scope then this pointer will be unsafe
     *       for use.
     */
    const double *getInputDataPointer(void) noexcept
    {
        const double *xptr = nullptr;
        if (nx_ > 0)
        {
            // "get" means I still own this pointer
            xptr = xptr_.get(); // if xptr_ is null then xptr will be null
        }
        else
        {
            if (x_ != nullptr){xptr = x_;}
        }
        return xptr;
    }


    /*!
     * @brief Sets the output (filtered) data.
     * @param[in] ny  The number of samples in y.
     * @param[in] y   The filtered signal to set.  This is an array
     *                of dimension [ny]. 
     */
    void setOutputData(const int ny, const double y[]) noexcept
    {
        // Resize y
        if (ny > maxy_)
        {
            if (y_){ippsFree(y_);}
            maxy_ = ny;
            y_ = ippsMalloc_64f(maxy_);
        }
        // Copy if possible
        if (ny > 0 && y)
        {
            ny_ = ny;
            ippsCopy_64f(y, y_, ny_);
        }
        else
        {
            if (ny_ > 0){RTSEIS_ERRMSG("%s", "y is NULL");}
            ny_ = 0;
        }
        return;
    }
    /*!
     * @brief Gets a pointer to the output data.
     * @result If this is a nullptr then there are no output data points.
     * @result Otherwise, this is a pointer to the filtered signal whose
     *         length can be determined by \c getNumberOfOutputSamples().
     * @note If this class goes out of scope then this pointer will be
     *       unsafe for use.
     */
    const double *getOutputDataPointer(void) const noexcept
    {
        const double *yptr = nullptr;
        if (ny_ > 0){yptr = y_;}
        return yptr;
    }
    
private:
    /// A pointer to the input data
    std::unique_ptr<const double> xptr_ = nullptr;
    /// The input data
    double *x_ = nullptr;
    /// The output data
    double *y_ = nullptr;
    /// Max space reserved to x
    int maxx_ = 0;
    /// Number of samples in x
    int nx_ = 0;
    /// Max space reserved to y
    int maxy_ = 0;
    /// Number of samples in y
    int ny_ = 0;
}; 


};
};

#endif
