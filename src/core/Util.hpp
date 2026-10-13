// Человекочитаемое время (и его же производная «время/с» для акселератора).
#pragma once
#include <string>

namespace decay {

// "153 s", "26.8 d", "30.1 ky", "4.47 Gy"
std::string humanTime(double seconds);

// particle/decay counts: 42; 12.4k; 2.5M
std::string humanCount(double count);

// accelerator rate: model seconds per real second
inline std::string humanRatePerSec(double simSecondsPerRealSecond) {
    return humanTime(simSecondsPerRealSecond) + "/s";
}

} // namespace decay