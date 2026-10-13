// Панель управления (ImGui): цепочка, число атомов, пауза и — главное — акселератор.
#pragma once
#include "core/Chains.hpp"
#include "core/DecaySim.hpp"
#include <cmath>
#include <cstdint>

namespace ui {

struct Panel {
    int    chainIdx = 0;
    float  atomsLog = 5.0f;    // атомов = 10^atomsLog
    bool   paused    = false;
    float  accelExp  = 16.0f;  // порядок множителя времени: время × 10^accelExp
    bool   wantReset = false;
    uint64_t lastEvents = 0;

    uint32_t atomsTarget() const {
        return (uint32_t)std::lround(std::pow(10.0, atomsLog));
    }
};

// palette — rgb по 3 float на нуклид текущей цепочки
void drawPanel(Panel& p, const decay::DecaySim& sim, const float* palette, int palCount,
               double fps, uint64_t frameEvents);

} // namespace ui