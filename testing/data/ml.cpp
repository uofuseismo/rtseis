#include <iostream>
#include <algorithm>
#include <cassert>
#include <limits>
#include <numeric>
#include <cmath>
#include <vector>
#include <fstream>

namespace
{

template<class T>
std::vector<T> cosineTaper(const std::vector<T> &x,
                           const double pct = 5)
{
    constexpr T zero = 0;
    std::vector<T> xt(x.size(), zero);
    if (x.empty()){return xt;}
    auto nx = static_cast<int> (x.size());
    auto nPct = static_cast<int> (std::round((nx*pct)/100)) + 1;
    auto m = std::max(2, std::min(nx, nPct));
    if (nx < 3){return xt;} // All zeros
    std::copy(x.begin(), x.end(), xt.begin());
    int mp12 = m/2; // Taper first (m+1)/2 points
    for (int i = 0; i < mp12; ++i)
    {
        auto scale = static_cast<T> (std::sin(i*M_PI/(m - 1)));
        //std::cout << i << " " << scale << std::endl;
        xt[i] = x[i]*scale;
        xt[nx - mp12 + i] = x[nx - mp12 + i]*scale;
    }
    return xt;
}

void getWoodAndersonConstants(double samplingPeriod, const bool isVelocity,
                              double *h, double *f, double *g)
{
    auto samplingRate = static_cast<int> (std::round(1/samplingPeriod)); 
    // defaults : "vsp with delta-t = 0.01 sec" = short period velocity
    // h0   = 0.781; f0  = 1.29; g0 = 2963.0;
    double h0 =  0.78262E+00;
    double f0 =  0.12886E+01;
    double g0 =  0.29694E+04; // 7 Hz
    if (isVelocity)
    {
        if (samplingRate == 100)
        {
            // Velocity : "vsp with delta-t = 0.01 sec" = short period velocity
            //h0   = 0.781; f0  = 1.29; g0 = 2963.0; // from paper, like 10 Hz
            //h0 =  0.78067E+00; f0 =  0.12884E+01; g0 =  0.29635E+04; // 10 Hz
            h0 =  0.78262E+00; f0 =  0.12886E+01; g0 =  0.29694E+04; // 7 Hz
        }
        else if (samplingRate == 200) // NEW added 09/14/2007 aww 
        {
            //h0 =  0.79245E+00; f0 =  0.12691E+01; g0 =  0.28861E+04; // 10 Hz
            h0 =  0.79251E+00; f0 =  0.12692E+01; g0 =  0.28865E+04; // 7 Hz
        }
        else if (samplingRate == 80) // NEW added 09/14/2007 aww
        {
            //h0 =  0.77412E+00; f0 =  0.12976E+01; g0 =  0.29996E+04; // 10 Hz
            h0 =  0.77721E+00; f0 =  0.12980E+01; g0 =  0.30093E+04; // 7 Hz
        }
        else if (samplingRate == 40)
        {
            //h0   = 0.743; f0  = 1.342; g0 = 3186.0; //  from Hiroo 08/18/06, like 7 HZ value
            //h0 =  0.72958E+00; f0 =  0.13398E+01; g0 =  0.31401E+04; // 10 Hz
            h0 =  0.74306E+00; f0 =  0.13422E+01; g0 =  0.31859E+04; // 7 Hz
        }
        else if (samplingRate == 50) // NEW added 09/14/2007 aww
        {
            //h0 =  0.74985E+00; f0 =  0.13238E+01; g0 =  0.30929E+04; // 10 Hz
            h0 =  0.75820E+00; f0 =  0.13251E+01; g0 =  0.31203E+04; // 7 Hz
        }
        else if (samplingRate == 20)
        {
            // Velocity : "vbb with delta-t = 0.05 sec" = broad-band velocity
            //h0   = 0.568; f0  = 1.39; g0 = 3110.0; // from paper, like 10 Hz
            //h0 =  0.57034E+00; f0 =  0.13924E+01; g0 =  0.31175E+04; // 10 Hz
            h0 =  0.63382E+00; f0 =  0.14094E+01; g0 =  0.33662E+04; // 7 Hz
        }
        else if (samplingRate == 250)
        {
            h0 =  0.79384E+00; f0 =  0.12656E+01; g0 =  0.28695E+04; // VEL fmax=7 250sps WA2800
        }
        else if (samplingRate == 500)
        {
            h0 =  0.79702E+00; f0 =  0.12578E+01; g0 =  0.28349E+04; // VEL fmax=7 500sps WA2800
        }
        else
        {
            std::cerr << "Unhandled sampling rate: " << samplingRate
                      << std::endl;
        }
    }
    else // accelerometer values
    {
        if (samplingRate == 100)
        {
            // Acceleration : "lg with delta-t = 0.01 sec" = "low gain" acceleration
            //h0   = 0.781; f0  = 1.29; g0 = 2963.0; // from paper, like 10 Hz?
            //h0 =  0.78401E+00; f0 =  0.12873E+01; g0 =  0.29683E+04; // 10 Hz
            h0 =  0.78427E+00; f0 =  0.12876E+01; g0 =  0.29700E+04; // 7 Hz
        }
        else if (samplingRate == 200)
        {
            //BUG ? g0 in line below appears to be 2080 value for 10Hz not the 2800 value 09/14/2007
            //h0 =  0.793; f0  = 1.27; g0 = 2144.0; // from Hiroo 9/1/04
            //h0 =  0.79245E+00; f0 =  0.12691E+01; g0 =  0.28861E+04; // 10 Hz
            h0 =  0.79251E+00; f0 =  0.12692E+01; g0 =  0.28865E+04; // 7 Hz
        }
        else if (samplingRate == 40)
        {
            //h0   = 0.754; f0  = 1.336; g0 = 3191.0; // from Hiroo 08/18/06, like 7 HZ 
            //h0 =  0.75220E+00; f0 =  0.13331E+01; g0 =  0.31774E+04; // 10 Hz
            h0 =  0.75422E+00; f0 =  0.13357E+01; g0 =  0.31913E+04; // 7 Hz
        }
        else if (samplingRate == 50) // NEW added 09/14/2007 aww
        {
            //h0 =  0.76394E+00; f0 =  0.13194E+01; g0 =  0.31152E+04; // 10 Hz
            h0 =  0.76517E+00; f0 =  0.13210E+01; g0 =  0.31235E+04; // 7 Hz
        }
        else if (samplingRate == 80)
        {
            //h0   = 0.781;  f0  = 1.29; g0 =  2963.0; // from paper
            //h0   = 0.774;  f0  = 1.30; g0 =  2999.0; // from RT code
            //h0 =  0.77937E+00; f0 =  0.12958E+01; g0 =  0.30073E+04; // 10 Hz
            h0 =  0.77980E+00; f0 =  0.12963E+01; g0 =  0.30101E+04; // 7 Hz
        }
        else if (samplingRate == 20) // NEW added 09/14/2007 aww
        {
            //h0 =  0.67044E+00; f0 =  0.13636E+01; g0 =  0.32957E+04; // 10 Hz
            h0 =  0.68405E+00; f0 =  0.13829E+01; g0 =  0.34011E+04; // 7 Hz
        }
        else if (samplingRate == 250) // NEW added 09/14/2007 aww
        {
            h0 =  0.79409E+00; f0 =  0.12654E+01; g0 =  0.28695E+04; // ACC fmax=7 250sps WA2800
        }
        else if (samplingRate == 500) // NEW added 09/14/2007 aww
        {
            h0 =  0.79711E+00; f0 =  0.12578E+01; g0 =  0.28349E+04; // ACC fmax=7 500sps WA2800
        }
        else
        {
            std::cerr << "Unhandled sampling rate: " << samplingRate
                      << std::endl;
        }
    }
    *h = h0;
    *f = f0;
    *g = g0;
}

std::vector<double> readSeismogram(const std::string &fileName)
{
    std::ifstream infl(fileName, std::ios::in);
    std::vector<double> x;
    if (infl.is_open())
    {
        std::string line;
        x.reserve(14000);
        while (std::getline(infl, line))
        {
            double t, xi;
            sscanf(line.c_str(), "%lf, %lf\n", &t, &xi);
            x.push_back(xi);
        }
        infl.close();
    }
    return x;
}

/// Removes the mean from a dataset
template<class T>
std::vector<T> removeMean(const std::vector<T> &x)
{
    const T zero = 0;
    std::vector<T> y(x.size());
    if (x.empty()){return y;}
    auto xsum = std::accumulate(x.begin(), x.end(), zero);
    auto xmean = xsum/static_cast<T> (x.size());
    const T *__restrict__ xPtr = x.data();
    T *__restrict__ yPtr = y.data();
    for (int i = 0; i < static_cast<int> (x.size()); ++i)
    {
        yPtr[i] = xPtr[i] - xmean;
    } 
    return y;
}

/// Converts all elements in x from meters to centimeters.
template<class T>
std::vector<T> metersToCentimeters(const std::vector<T> &x)
{
    std::vector<T> y(x.size());
    // y = x*100
    const T scalar = static_cast<T> (100);
    std::transform(x.begin(), x.end(), y.begin(),
                   [&scalar](auto &c)
                   {
                      return scalar*c;
                   });
    return y;
}

}

/// Applies the Wood Anderson gain correction in-place
template<class T>
void correctWoodAndersonGain(std::vector<T> &x)
{
    // y = y*(2080/2800)
    const T fac2080_2800 = static_cast<T> (2080/2800.);
    std::transform(x.begin(), x.end(), x.begin(),
                   [&fac2080_2800](auto &c)
                   {
                       return c*fac2080_2800;
                   });
}

/// @brief Computes the Wood Anderson time domain filter for a
///        signal proportional to acceleration.
template<class T>
std::vector<T> accelerationFilter(const std::vector<T> &x,
                                  const double gain,
                                  const double dt,
                                  const bool applyFactor = false)
{
    constexpr T zero = 0;
    std::vector<T> y(x.size(), zero);
    constexpr T qHighPass = static_cast<T> (0.998);
    constexpr bool isVelocity = false;
    double h0, f0, g0;
    getWoodAndersonConstants(dt, isVelocity, &h0, &f0, &g0);    
    constexpr auto b0 = 2/(1 + qHighPass);
    double wdt = (2*M_PI)*(f0*dt);
    double c1 = 1 + h0*wdt;
    double c2 = 1 + 2*(h0*wdt) + wdt*wdt;
    auto c1x2 = static_cast<T> (2*c1);
    auto c2Inv = static_cast<T> (1/c2);
    auto gwa_dt2 = static_cast<T> (g0*(dt*dt));
    // This will remove the gain from acceleration and apply a highpass filter.
    auto divisor = static_cast<T> (1./(b0*gain));
    // Compute the acceleration and apply the Wood Anderson filter in one shot.
    // This is the composite of the two loops in Jiggle.
    y[0] = 0;
    y[1] = 0;
    T acc_i_1 = (x[1] - x[0])*divisor; // a[1] (note qHighPass*a[0] = 0)
    T acc_i = 0;
    for (int i = 2; i < static_cast<int> (y.size()); ++i)
    {
        // Compute the acceleration.  The is Eqn 10 from Kanamori et al., 1999
        acc_i = (x[i] - x[i-1])*divisor + qHighPass*acc_i_1;
        // Apply Wood-Anderson filter.  Eqn 3 assumes an input velocity signal
        // that is differentiated on the fly (v[k] - v[k-1])/dt.  That means,
        // we have to multiply by dt.  Additionally, the division by the gain
        // was also performed so all that's left to do is multiply by optimized
        // gain times the dt.  This is what gwa_dt2 does.
        y[i] = (acc_i*gwa_dt2 + c1x2*y[i-1] - y[i-2])*c2Inv;
        acc_i_1 = acc_i;
    } 
    constexpr bool debug = true;//false;
    if (debug)
    {
        // Compute the acceleration.  The is Eqn 10 from Kanamori et al., 1999.
        std::vector<T> acc(x.size(), zero);
        for (int i = 1; i < static_cast<int> (x.size()); ++i)
        {
            acc[i] = (x[i] - x[i-1])*divisor + qHighPass*acc[i-1];
        }   
        // Apply Wood-Anderson filter.  This is Eqn 3 from Kanamori et al., 1999
        // but reformed for an input acceleration signal.
        std::vector<T> yref(y.size(), zero);
        for (int i = 2; i < static_cast<int> (y.size()); ++i)
        {
            yref[i] = (acc[i]*gwa_dt2 + c1x2*yref[i-1] - yref[i-2])*c2Inv;
        }
        // Ensure my filter implementation matches the Jiggle variant exactly.
        T errmax = 0;
        for (int i = 0; i < static_cast<int> (y.size()); ++i)
        {
            errmax = std::max(errmax, std::abs(yref[i] - y[i]));
        }
        assert(errmax < std::numeric_limits<T>::epsilon()*10);
        //std::cout << "Max error: " << errmax << std::endl;
    }
    // At one point the gain on a Wood Anderson filter was revised from 2800
    // to 2080.  This section will remove the 2800 and rescale by 2080.
    if (applyFactor){correctWoodAndersonGain(y);}
    return y;
}

/// @brief Computes the Wood Anderson time domain filter for a signal
///        proportional to velocity.
template<class T>
std::vector<T> velocityFilter(const std::vector<T> &x,
                              const double gain,
                              const double dt,
                              const bool applyFactor = false)
{
    constexpr T zero = 0;
    std::vector<T> y(x.size(), 0);
    constexpr bool isVelocity = true;
    double h0, f0, g0;
    getWoodAndersonConstants(dt, isVelocity, &h0, &f0, &g0);
    double wdt = (2*M_PI)*(f0*dt);
    double c1 = 1 + h0*wdt;
    double c2 = 1 + 2*(h0*wdt) + wdt*wdt;
    auto scalar = static_cast<T> ((g0/gain)*dt);
    auto c1x2 = static_cast<T> (2*c1);
    auto c2Inv = static_cast<T> (1/c2);
    // This is Equation 3 of Kanamori et al., 1999. 
    // Note - Jiggle's implementation is wrong and doesn't correctly account for
    // the zero initial conditions. 
    y[0] = ( (x[0] - 0)*scalar + (c1x2*0 - 0) )*c2Inv;
    y[1] = ( (x[1] - x[0])*scalar + (c1x2*y[0] - 0) )*c2Inv;
    for (int i = 2; i < static_cast<int> (x.size()) - 3; ++i)
    {
        y[i] = ( (x[i] - x[i-1])*scalar + (c1x2*y[i-1] - y[i-2]) )*c2Inv;
    }
    // At one point the gain on a Wood Anderson filter was revised from 2800
    // to 2080.  This section will remove the 2800 and rescale by 2080.
    if (applyFactor){correctWoodAndersonGain(y);}
    return y;
} 

int main()
{
    double dt = 0.01;
    double t0 = 1624835706.4;
    double tp = 1624835745.64;
    auto idx = static_cast<int> ((tp - t0)/dt);
    double gain = 1274800764.8712056;
    double accGain = 321168.428435;
    double ampObs = 0.0137;
    double ampNextCycle = -0.0134;
    auto velTrace = readSeismogram("UU.SPU.HHN.01.txt");
    auto velDemeaned = metersToCentimeters(removeMean(velTrace));
    //auto velDemeaned = cosineTaper(removeMean(velTrace), 5);
    auto y = velocityFilter(velDemeaned, gain, dt);
    //y = metersToCentimeters(y);
    std::cout << velTrace.size() << std::endl;
    std::ofstream ofl("spu_wa.txt", std::ios::out);
    for (int i = 0; i < static_cast<int> (y.size()); ++i)
    {
        ofl << i*dt << " " << y[i] << std::endl;
    } 
    ofl.close();
    std::cout <<-(-0.009155 - 0.008697) << std::endl;
    std::cout << "Amp from Jiggle: " << (ampObs - ampNextCycle) << " " << (ampObs - ampNextCycle)/2 << std::endl;
    std::cout << "observed amplitude: " << ampObs << " amp/2: " << ampObs/2 << std::endl; 
    std::cout << y[idx-1] << " " << y[idx] << " " << y[idx+1] << std::endl;
    std::cout << y[idx+11-1] << " " << y[idx+11] << " " << y[idx+11+1] << std::endl;
    std::cout << "Recovered amplitude: " << y[idx] - y[idx+11] << " " << (y[idx] - y[idx+11])/2 << std::endl;
    std::cout << std::endl;

    t0 = 1617656190.35; 
    tp = 1617656234.91;
    ampObs = 0.4739763736724854;
    idx = static_cast<int> ((tp - t0)/dt);
    auto accTrace = readSeismogram("UU.TMU.ENN.01.txt");
    auto accDemeaned = metersToCentimeters(cosineTaper(removeMean(accTrace), 5));
    y = accelerationFilter(accDemeaned, accGain, dt);
    ofl.open("tmu_wa.txt", std::ios::out);
    for (int i = 0; i < static_cast<int> (y.size()); ++i)
    {
        ofl << i*dt << " " << y[i] << std::endl;
    }
    std::cout << (0.16 - -0.311) << std::endl;
    std::cout << "observed acc amp: " << ampObs << " amp/2 " << ampObs/2 << std::endl;
    std::cout << y[idx-1] << " " << y[idx] << " " << y[idx+1] << std::endl;
    std::cout << y[idx+15] << " " << y[idx+15+1] << " " << y[idx+15+2] << std::endl;
    std::cout << "Recovered acc amplitude: " << y[idx+15+1] - y[idx] << " " << (y[idx+15+1] - y[idx])/2 << std::endl; 
std::cout << y[idx-15] << " "  << y[idx-14] << " " << y[idx-13] << " " << y[idx-12] << std::endl;
    ofl.close();
    
    return EXIT_SUCCESS;    
}
