#ifndef RTSEIS_WAVEFORM_H
#define RTSEIS_WAVEFORM_H 1
#include <stdlib.h>
#include "rtseis/config.h"
#include "rtseis/modules.hpp"
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

class Waveform;
class Command;

class Command
{
    enum class FilterType
    {
        /*!< No action. */
        NONE = 0,
        /*!< Remove the mean from the data. */
        DEMEAN = 1,
        /*!< Remove the trend fro the data. */
        DETREND = 2
    };
    public:
        /*!
         * @brief Default constructor.
         */
        Command(void)
        {
            clear();
        }
        /*!
         * @brief Initializes a demean command from the given parameters.
         * @param[in] parms  The demean parameters.
         */
        Command(const Modules::DemeanParameters &parms)
        {
            clear();
            int ierr = demean_.setParameters(parms);
            if (ierr != 0)
            {
                //RTSEIS_ERRMSG("%s", "Failed to set demean parameters");
                filterType_ = FilterType::NONE;
                demean_.clear();
                return;
            }
            filterType_ = FilterType::DEMEAN;
            return;
        }
        /*!
         * @brief Initializes a detrend command from the given parameters.
         * @param[in] parms  The detrend parameters.
         */
        Command(const Modules::DetrendParameters &parms)
        {
            clear();
            int ierr = detrend_.setParameters(parms);
            if (ierr != 0)
            {
                //RTSEIS_ERRMSG("%s", "Failed to set detrend parameters");
                filterType_ = FilterType::NONE;
                detrend_.clear();
                return;
            }
            filterType_ = FilterType::DETREND;
            return;
        }
        /*!
         * @brief Deletes the current command and sets the filter to none.
         */
        void clear(void)
        {
            if (filterType_ != FilterType::NONE)
            {
                if (filterType_ == FilterType::DEMEAN)
                {
                    demean_.clear();
                }
                else if (filterType_ == FilterType::DETREND)
                {
                    detrend_.clear();
                }
            }
            filterType_ = FilterType::NONE;
        }
        //int setDetrendCommand(const DetrendParameters );
    private:
        /*!< Module for removing mean from data. */
        Modules::Demean demean_;
        /*!< Module for removing trend from data. */
        Modules::Detrend detrend_;
        /*!< Defines the filter type to apply. */
        FilterType filterType_ = FilterType::NONE;
    friend class Waveform;
};

class Waveform : 
    public RTSeis::Modules::Demean//, RTSeis::Modules::Detrend
{
    public:
        Waveform(void);
        int setNetworkName(const std::string network);
        int setStationName(const std::string station);
        int setChannelName(const std::string channel);
        int setLocationCode(const std::string location);
        int getLocationCode(std::string &location);

        /*!
         * @brief Sets the data to process.
         * @param[in] nx   Number of points in time series.
         * @param[in] x    Time series to process.
         * @result 0 indicates success.
         * @ingroup rtseis_data_waveform
         */
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
        /*!
         * @brief Applies the demean command to the data.
         * @result 0 indicates succes.
         * @ingroup rtseis_data_waveform
         */
        int demean(void)
        {
            int ierr;
            Modules::DemeanParameters parms;
            command_ = Command(parms);
            if (nx_ == 0){return 0;} 
            if (!lhaveX_)
            {
                //if (!lhaveX_){RTSEIS_ERRMSG("%s", "No data");}
                return -1; 
            }
            resizeY_(ny_);
            ierr = command_.demean_.demean(nx_, x_, y_);
            if (ierr != 0)
            {
                //RTSEIS_ERRMSG("%s", "Failed to apply demean");
            }
            return ierr;
        }
        /*!
         * @brief Applies the detrend command to the data.
         * @result 0 indicates succes.
         * @ingroup rtseis_data_waveform
         */
        int detrend(void)
        {
            int ierr;
            Modules::DetrendParameters parms;
            command_ = Command(parms);
            if (nx_ == 0){return 0;}
            if (!lhaveX_ || nx_ < 2)
            {
                //if (!lhaveX_){RTSEIS_ERRMSG("%s", "No data");}
                //if (nx_ < 2){RTSEIS_ERRMSG("%s", "At least 2 points needed");}
                return -1;
            }
            resizeY_(ny_);
            ierr = command_.detrend_.detrend(nx_, x_, y_);
            if (ierr != 0)
            {
                //RTSEIS_ERRMSG("%s", "Failed to apply detrend");
            }
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
        int setDetrendParameters(const Modules::DetrendParameters &parameters)
        {
            detrend_.setParameters(parameters);
            return 0; //ierr;
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
                ierr = detrend_.detrend(nx_, xptr_, y_);
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
        /*!< Classes. */
        Modules::Detrend detrend_;
        /*!< Sampling period. */
        double dt_ = 0;
        /*!< Number of samples in input signal. */
        int nx_ = 0; 
        /*!< Number of samples in output signals. */
        int ny_ = 0;
        /*!< Command to apply. */
        Data::Command command_;
        /*!< Precision. */
        enum rtseisPrecision_enum precision_ = RTSEIS_DOUBLE;
        /*!< Command to apply. */
        enum rtseisCommand_enum job_ = RTSEIS_COMMAND_NONE;
};

}; /* Data */
}; /* RTSeis */


#endif
