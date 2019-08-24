#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <stdexcept>
#include "rtseis/private/throw.hpp"
#include "rtseis/utilities/rotate/utilities.hpp"

/// Convert north/east to radial/transverse
template<typename T>
void RTSeis::Utilities::Rotate::northEastToRadialTransverse(
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
        if (north == nullptr){RTSEIS_THROW_IA("%s", "north is NULL");}
        if (east == nullptr){RTSEIS_THROW_IA("%s", "east is NULL");}
        if (radial == nullptr){RTSEIS_THROW_IA("%s", "radial is NULL");}
        RTSEIS_THROW_IA("%s", "transverse is NULL"); 
    }
    // Stein and Wysession Eqn 44.
    T theta = 3*M_PI_2 - backAzimuth;
    T ca = std::cos(theta);
    T sa = std::sin(theta);
    #pragma omp simd
    for (int i=0; i<nSamples; ++i)
    {
        radial[i]     = ca*east[i] + sa*north[i];
        transverse[i] =-sa*east[i] + ca*north[i];
    }
}

/// Convert radial/transverse to north/east
template<typename T>
void RTSeis::Utilities::Rotate::radialTransverseToNorthEast(
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
        if (north == nullptr){RTSEIS_THROW_IA("%s", "north is NULL");}
        if (east == nullptr){RTSEIS_THROW_IA("%s", "east is NULL");}
        if (radial == nullptr){RTSEIS_THROW_IA("%s", "radial is NULL");}
        RTSEIS_THROW_IA("%s", "transverse is NULL");
    }
    // Stein and Wysession Eqn 44.
    T theta = 3*M_PI_2 - backAzimuth;
    T ca = std::cos(theta);
    T sa = std::sin(theta);
    #pragma omp simd
    for (int i=0; i<nSamples; ++i)
    {
        north[i] = ca*radial[i] - sa*transverse[i];
        east[i]  = sa*radial[i] + ca*transverse[i];
    }
}

/// Convert vertical/north/east to longitudinal/radial/transverse
template<typename T>
void RTSeis::Utilities::Rotate::verticalNorthEastToLongitudinalRadialTransverse(
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
        if (vertical == nullptr){RTSEIS_THROW_IA("%s", "vertical is NULL");}
        if (north == nullptr){RTSEIS_THROW_IA("%s", "north is NULL");}
        if (east == nullptr){RTSEIS_THROW_IA("%s", "east is NULL");}
        if (longitudinal == nullptr)
        {
            RTSEIS_THROW_IA("%s", "longitudinal is NULL");
        }
        if (radial == nullptr){RTSEIS_THROW_IA("%s", "radial is NULL");}
        RTSEIS_THROW_IA("%s", "transverse is NULL");
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
        radial[i]       = si*vertical[i] - cisb*east[i] - cicb*north[i];
        transverse[i]   =                -   cb*east[i] +   sb*north[i];
    }
}

/// Convert longitudinal/radial/transverse to vertical/north/east
template<typename T>
void RTSeis::Utilities::Rotate::longitudinalRadialTransverseToVerticalNorthEast(
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
        if (vertical == nullptr){RTSEIS_THROW_IA("%s", "vertical is NULL");}
        if (north == nullptr){RTSEIS_THROW_IA("%s", "north is NULL");}
        if (east == nullptr){RTSEIS_THROW_IA("%s", "east is NULL");}
        if (longitudinal == nullptr)
        {
            RTSEIS_THROW_IA("%s", "longitudinal is NULL");
        }
        if (radial == nullptr){RTSEIS_THROW_IA("%s", "radial is NULL");}
        RTSEIS_THROW_IA("%s", "transverse is NULL");
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
        vertical[i] =   ci*longitudinal[i] - sisb*radial[i] - sicb*transverse[i];
        east[i]     =   si*longitudinal[i] - cisb*radial[i] - cicb*transverse[i];
        north[i]    =-sicb*longitudinal[i] - sb*radial[i] - cb*transverse[i];
    }
}

/// Template function instantiation
template
void RTSeis::Utilities::Rotate::northEastToRadialTransverse<double>(
    const int nSamples,
    const double backAzimuth,
    const double north[],
    const double east[],
    double *radial[],
    double *transverse[]);
template
void RTSeis::Utilities::Rotate::northEastToRadialTransverse<float>(
    const int nSamples,
    const float backAzimuth,
    const float north[],
    const float east[],
    float *radial[],
    float *transverse[]);

template
void RTSeis::Utilities::Rotate::verticalNorthEastToLongitudinalRadialTransverse<double>(
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
void RTSeis::Utilities::Rotate::verticalNorthEastToLongitudinalRadialTransverse<float>(
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
