#include <cstdio>
#include <cstdlib>
#include <string>
#include "rtseis/version.hpp"

using namespace RTSeis;

int Version::getMajor() noexcept
{
    return RTSEIS_MAJOR;
}

int Version::getMinor() noexcept
{
    return RTSEIS_MINOR;
}

int Version::getPatch() noexcept
{
    return RTSEIS_PATCH;
}

bool Version::isAtLeast(const int major, const int minor,
                        const int patch) noexcept
{
    if (RTSEIS_MAJOR < major){return false;}
    if (RTSEIS_MAJOR > major){return true;}
    if (RTSEIS_MINOR < minor){return false;}
    if (RTSEIS_MINOR > minor){return true;}
    if (RTSEIS_PATCH < patch){return false;}
    return true;
}

std::string Version::getVersion() noexcept
{
    std::string version(RTSEIS_VERSION);
    return version;
}
