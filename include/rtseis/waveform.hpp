#ifndef RTSEIS_WAVEFORM_H
#define RTSEIS_WAVEFORM_H 1
#include <stdlib.h>
#include "rtseis/config.h"
#include "rtseis/modules.h"
#include <string>

enum rtseisCommand_enum
{
    RTSEIS_COMMAND_NONE = 0,
    RTSEIS_COMMAND_DEMEAN = 1, 
    RTSEIS_COMMAND_DETREND = 2
};

namespace RTSeis
{
namespace Data
{

class Waveform : 
    public RTSeis::Modules::Demean, 
           RTSeis::Modules::Detrend
{
    public:
        Waveform(void);
        int setNetworkName(const std::string network);
        int setStationName(const std::string station);
        int setChannelName(const std::string channel);
        int setLocationCode(const std::string location);
        int getLocationCode(std::string &location);

        int setData(const int nx, const double x[])
        {
            if (x_ != nullptr){free(x_); x_ = nullptr;}
            xptr_ = xptr_;
            lxptr_ = false;
            lhaveX_ = false;
            int ierr = setNumberOfSamples_(nx);
            if (ierr != 0)
            {
                //RTSEIS_ERRMSG("%s", "No points to set on x");
                return -1;
            }
            if (x == nullptr)
            {
                //RTSEIS_ERRMSG("%s", "x is null");
                return -1;
            }
            // Set the data
            nx_ = nx;
            x_ = static_cast<double *> (calloc(nx_, sizeof(double))); 
            xptr_ = x_;
            lhaveX_ = true;
            return 0;
        }

        int detrend(void)
        {
            int ierr;
            if (nx_ == 0){return 0;}
            if (!lhaveX_)
            {
                //RTSEIS_ERRMSG("%s", "No input data");
                return -1;
            }
            resizeY_(nx_);
            job_ = RTSEIS_COMMAND_DETREND;
            ierr = compute_();
            return ierr;
        }
        int demean(void)
        {
            int ierr;
            if (nx_ == 0){return 0;}
            if (!lhaveX_)
            {
                return -1;
            }
            if (nx_ < 2)
            {
                return -1;
            }
            resizeY_(ny_); 
            job_ = RTSEIS_COMMAND_DEMEAN;
            ierr = compute_();
            return ierr;
        }
        /*!
         * @brief Sets the parameters for the demeaning.
         * @result 0 indicates success.
         * @ingroup rtseis_data_waveform
         */
        int setDemeanParameters(void)
        {
            int ierr = RTSeis::Modules::Demean::setParameters(precision_);
            return ierr;
        }
        /*!
         * @brief Sets the parameters for the detrending.
         * @result 0 indicates success.
         * @ingroup rtseis_data_waveform
         */
        int setDetrendParameters(const DetrendParameters &parameters)
        {
            int ierr = RTSeis::Modules::Detrend::setParameters(parameters);
            return ierr;
        }
    private:
        int compute_(void)
        {
            int ierr = 0;
            if (job_ == RTSEIS_COMMAND_NONE)
            {
                return 0;
            }
            else if (job_ == RTSEIS_COMMAND_DEMEAN)
            {
                ierr = RTSeis::Modules::Detrend::detrend(nx_, xptr_, y_);
            }
            else if (job_ == RTSEIS_COMMAND_DETREND)
            {
                ierr = RTSeis::Modules::Demean::demean(nx_, xptr_, y_);
            }
            else
            {
                //RTSEIS_ERRMSG("%s", "Unknown job");
                return -1;
            }
            return ierr;
        }

        int resizeY_(const int n)
        {
            if (n > ny_)
            {
                if (y_ != nullptr)
                {
                    free(y_);
                    y_ = nullptr;
                }
                ny_ = n;
            }
            if (y_ == nullptr)
            {
                y_ = static_cast<double *> (calloc(ny_, sizeof(double)));
            }
            return 0;
        }
        double *xptr_ = nullptr;
        double *x_ = nullptr;
        double *y_ = nullptr;
        bool lxptr_ = false;
        bool lhaveX_ = false;
        int setNumberOfSamples_(const int nx)
        {
            nx_ = 0;
            if (nx < 1)
            {
                //RTSEIS_ERRMSG("%s", "Number of samples must be positive");
                return -1;
            }
            nx_ = nx;
            return 0;
        }
        /*!< Network code. */
        std::string network_;
        /*!< Station code. */
        std::string station_;
        /*!< Channel code. */
        std::string channel_;
        /*!< Location code. */
        std::string location_;
        /*!< Sampling period. */
        double dt_ = 0;
        /*!< Number of samples in input signal. */
        int nx_ = 0; 
        /*!< Number of samples in output signals. */
        int ny_ = 0;
        /*!< Precision. */
        enum rtseisPrecision_enum precision_ = RTSEIS_DOUBLE;
        /*!< Command to apply. */
        enum rtseisCommand_enum job_ = RTSEIS_COMMAND_NONE;
};

}; /* Data */
}; /* RTSeis */


#endif
