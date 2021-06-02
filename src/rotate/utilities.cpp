#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <stdexcept>
#include "rtseis/rotate/utilities.hpp"

/// Convert north/east to radial/transverse
template<typename T>
void RTSeis::Rotate::northEastToRadialTransverse(
    const int nSamples,
    const T backAzimuth,
    const T north[],
    const T east[],
    T *radialIn[],
    T *transverseIn[])
{
    // Checks
    if (nSamples < 1){return;}
    T *radial = *radialIn;
    T *transverse = *transverseIn;
    if (north == nullptr || east == nullptr || 
        radial == nullptr || transverse == nullptr)
    {
        if (north == nullptr){throw std::invalid_argument("north is NULL");}
        if (east == nullptr){throw std::invalid_argument("east is NULL");}
        if (radial == nullptr){throw std::invalid_argument("radial is NULL");}
        throw std::invalid_argument("transverse is NULL"); 
    }
    // Use a counterclockwise rotation matrix.  Flip back-azimuth to
    // theta and make it increase negatively.  Comparing with
    //  {T} = [ cos(az) sin(az)]{E}
    //  {R} = [-sin(az) cos(az)]{N}
    // with 
    //  {T} = [ cos(180-az) sin(180-az)]{E}
    //  {R} = [-sin(180-az) cos(180-az)]{N}
    // becomes 
    //  {T} = [-cos(baz)  sin(baz)]{E}
    //  {R} = [-sin(baz) -cos(baz)]{N}
    // Reordering:
    //  R = -N cos(baz) - E sin(baz) 
    //  T =  N sin(baz) - E cos(baz)
    // where the -az means az increases clockwise.  This is exactly the equation
    // proposed by (4.19-4.20) of Haskov - Routine Data Processing in Earthquake
    // Seismology.  Stein and Wysession Eqn 44 on pg 58 results in a
    // polarity flip on the transvserse which corresponds to a right-handed 
    // system which I don't want because stations usually use a left-handed
    // system.
    T cb = std::cos(backAzimuth);
    T sb = std::sin(backAzimuth);
    T ncb =-cb;
    T nsb =-sb;
    #pragma omp simd
    for (int i=0; i<nSamples; ++i)
    {
        radial[i]     = ncb*north[i] + nsb*east[i];
        transverse[i] =  sb*north[i] + ncb*east[i]; 
         
    }
}

/// Convert radial/transverse to north/east
template<typename T>
void RTSeis::Rotate::radialTransverseToNorthEast(
    const int nSamples,
    const T backAzimuth,
    const T radial[],
    const T transverse[],
    T *northIn[],
    T *eastIn[])
{
    // Checks
    if (nSamples < 1){return;}
    T *north = *northIn;
    T *east = *eastIn;
    if (north == nullptr || east == nullptr ||
        radial == nullptr || transverse == nullptr)
    {
        if (north == nullptr){throw std::invalid_argument("north is NULL");}
        if (east == nullptr){throw std::invalid_argument("east is NULL");}
        if (radial == nullptr){throw std::invalid_argument("radial is NULL");}
        throw std::invalid_argument("transverse is NULL");
    }
    // Compute Inverse of pg 95 Haskov 4.19 - 4.20
    T cb = std::cos(backAzimuth);
    T sb = std::sin(backAzimuth);
    T ncb =-cb;
    T nsb =-sb;
    #pragma omp simd
    for (int i=0; i<nSamples; ++i)
    {
        north[i] = ncb*radial[i] +  sb*transverse[i];
        east[i]  = nsb*radial[i] + ncb*transverse[i];
    }
}

/// Convert vertical/north/east to longitudinal/radial/transverse
template<typename T>
void RTSeis::Rotate::verticalNorthEastToLongitudinalRadialTransverse(
    const int nSamples,
    const T backAzimuth,
    const T incidenceAngle,
    const T vertical[],
    const T north[],
    const T east[],
    T *longitudinalIn[],
    T *radialIn[],
    T *transverseIn[]
    )
{
    // Checks
    if (nSamples < 1){return;}
    T *longitudinal = *longitudinalIn;
    T *radial = *radialIn;
    T *transverse = *transverseIn;
    if (vertical == nullptr || north == nullptr || east == nullptr ||
        longitudinal == nullptr || radial == nullptr || transverse == nullptr)
    {
        if (vertical == nullptr)
        {
            throw std::invalid_argument("vertical is NULL");
        }
        if (north == nullptr){throw std::invalid_argument("north is NULL");}
        if (east == nullptr){throw std::invalid_argument("east is NULL");}
        if (longitudinal == nullptr)
        {
            throw std::invalid_argument("longitudinal is NULL");
        }
        if (radial == nullptr){throw std::invalid_argument("radial is NULL");}
        throw std::invalid_argument("transverse is NULL");
    }
    T ci = std::cos(incidenceAngle);
    T si = std::sin(incidenceAngle);
    T cb = std::cos(backAzimuth);
    T sb = std::sin(backAzimuth);
    T sisb = si*sb;
    T cisb = ci*sb;
    T sicb = si*cb;
    T cicb = ci*cb; 
    #pragma omp simd
    for (int i=0; i<nSamples; ++i)
    {
        longitudinal[i] = ci*vertical[i] - sisb*east[i] - sicb*north[i];
        radial[i]       =-si*vertical[i] - cisb*east[i] - cicb*north[i];
        transverse[i]   =                -   cb*east[i] +   sb*north[i];
    }
}

/// Convert longitudinal/radial/transverse to vertical/north/east
template<typename T>
void RTSeis::Rotate::longitudinalRadialTransverseToVerticalNorthEast(
    const int nSamples,
    const T backAzimuth,
    const T incidenceAngle,
    const T longitudinal[],
    const T radial[],
    const T transverse[],
    T *verticalIn[],
    T *northIn[],
    T *eastIn[]
    )
{
    // Checks
    if (nSamples < 1){return;}
    T *vertical = *verticalIn;
    T *north = *northIn;
    T *east = *eastIn;
    if (vertical == nullptr || north == nullptr || east == nullptr ||
        longitudinal == nullptr || radial == nullptr || transverse == nullptr)
    {
        if (vertical == nullptr)
        {
            throw std::invalid_argument("vertical is NULL");
        }
        if (north == nullptr){throw std::invalid_argument("north is NULL");}
        if (east == nullptr){throw std::invalid_argument("east is NULL");}
        if (longitudinal == nullptr)
        {
            throw std::invalid_argument("longitudinal is NULL");
        }
        if (radial == nullptr){throw std::invalid_argument("radial is NULL");}
        throw std::invalid_argument("transverse is NULL");
    }
    T ci = std::cos(incidenceAngle);
    T si = std::sin(incidenceAngle);
    T cb = std::cos(backAzimuth);
    T sb = std::sin(backAzimuth);
    T sisb = si*sb;
    T cisb = ci*sb;
    T sicb = si*cb;
    T cicb = ci*cb;
    #pragma omp simd
    for (int i=0; i<nSamples; ++i)
    {
        vertical[i] =   ci*longitudinal[i] - si*radial[i];
        east[i]     =-sisb*longitudinal[i] - cisb*radial[i] - cb*transverse[i];
        north[i]    =-sicb*longitudinal[i] - cicb*radial[i] + sb*transverse[i];
    }
}

///--------------------------------------------------------------------------///
///                    Template function instantiation                       ///
///--------------------------------------------------------------------------///
template
void RTSeis::Rotate::northEastToRadialTransverse<double>(
    const int nSamples,
    const double backAzimuth,
    const double north[],
    const double east[],
    double *radial[],
    double *transverse[]);
template
void RTSeis::Rotate::northEastToRadialTransverse<float>(
    const int nSamples,
    const float backAzimuth,
    const float north[],
    const float east[],
    float *radial[],
    float *transverse[]);

template
void RTSeis::Rotate::radialTransverseToNorthEast<double>(
    const int nSamples,
    const double backAzimuth,
    const double radial[],
    const double transverse[],
    double *north[],
    double *east[]);
template
void RTSeis::Rotate::radialTransverseToNorthEast<float>(
    const int nSamples,
    const float backAzimuth,
    const float radial[],
    const float transverse[],
    float *north[],
    float *east[]);

template
void RTSeis::Rotate::verticalNorthEastToLongitudinalRadialTransverse<double>(
    const int nSamples,
    const double backAzimuth,
    const double incidenceAngle,
    const double vertical[],
    const double north[],
    const double east[],
    double *longitudinalIn[],
    double *radialIn[],
    double *transverseIn[]
    );
template
void RTSeis::Rotate::verticalNorthEastToLongitudinalRadialTransverse<float>(
    const int nSamples,
    const float backAzimuth,
    const float incidenceAngle,
    const float vertical[],
    const float north[],
    const float east[],
    float *longitudinalIn[],
    float *radialIn[],
    float *transverseIn[]
    );

template
void RTSeis::Rotate::longitudinalRadialTransverseToVerticalNorthEast<double>(
    const int nSamples,
    const double backAzimuth,
    const double incidenceAngle,
    const double longitudinal[],
    const double radial[],
    const double transverse[],
    double *verticalIn[],
    double *northIn[],
    double *eastIn[]
    );
template
void RTSeis::Rotate::longitudinalRadialTransverseToVerticalNorthEast<float>(
    const int nSamples,
    const float backAzimuth,
    const float incidenceAngle,
    const float longitudinal[],
    const float radial[],
    const float transverse[],
    float *verticalIn[],
    float *northIn[],
    float *eastIn[]
    );
