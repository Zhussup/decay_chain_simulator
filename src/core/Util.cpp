#include "Util.hpp"
#include <cstdio>

namespace decay {

std::string humanTime(double seconds) {
    char buf[64];
    if (seconds < 1.0) {
        std::snprintf(buf, sizeof buf, "%.2g s", seconds);
        return buf;
    }
    if (seconds < 60.0)      { std::snprintf(buf, sizeof buf, "%.1f s",   seconds); return buf; }
    if (seconds < 3600.0)    { std::snprintf(buf, sizeof buf, "%.1f min", seconds / 60.0); return buf; }
    if (seconds < 86400.0)   { std::snprintf(buf, sizeof buf, "%.1f h",   seconds / 3600.0); return buf; }
    if (seconds < 31557600.0){ std::snprintf(buf, sizeof buf, "%.1f d",   seconds / 86400.0); return buf; }

    const double years = seconds / 31557600.0;
    if (years < 1e3)  { std::snprintf(buf, sizeof buf, "%.2f yr", years); return buf; }
    if (years < 1e6)  { std::snprintf(buf, sizeof buf, "%.2f ky", years / 1e3); return buf; }
    if (years < 1e9)  { std::snprintf(buf, sizeof buf, "%.2f My", years / 1e6); return buf; }
    if (years < 1e15) { std::snprintf(buf, sizeof buf, "%.2f Gy", years / 1e9); return buf; }
    std::snprintf(buf, sizeof buf, "%.1e yr", years);
    return buf;
}

// счёт частиц/распадов: 42; 12.4k; 2.5M
std::string humanCount(double count) {
    char buf[32];
    if (count < 1e4) { std::snprintf(buf, sizeof buf, "%.0f", count); return buf; }
    if (count < 1e6) { std::snprintf(buf, sizeof buf, "%.1fk", count / 1e3); return buf; }
    if (count < 1e9) { std::snprintf(buf, sizeof buf, "%.1fM", count / 1e6); return buf; }
    std::snprintf(buf, sizeof buf, "%.1fB", count / 1e9);
    return buf;
}

} // namespace decay