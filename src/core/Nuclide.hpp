// Модель нуклида и цепочки распада. Данные цепочек — в Chains.cpp.
#pragma once
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace decay {

struct Branch {
    int   to   = 0;    // индекс дочернего нуклида в той же цепочке
    float prob = 1.f;  // вероятность ветки (сумма по нуклиду = 1)
};

struct Nuclide {
    std::string name = "?";      // "Cs-137"
    int         Z        = 0;
    int         A        = 0;
    double      halfLife = std::numeric_limits<double>::infinity(); // секунд; inf = стабилен
    const char* mode     = "";  // преобладающий канал: "α", "β⁻", "γ", "" у стабильных
    std::vector<Branch> branches;

    bool   isStable() const { return !std::isfinite(halfLife); }
    double lambda()   const { return isStable() ? 0.0 : 0.69314718055994530942 / halfLife; }
};

struct Chain {
    const char*          id = "?";   // короткое имя для UI: "U-238"
    const char*          title = "?";
    std::vector<Nuclide> nucs;       // nucs[0] — исходное ядро; ветвление — через branches
    int                  start = 0;  // с какого нуклида начинаем симуляцию

    int find(const char* name) const {
        for (int i = 0; i < (int)nucs.size(); ++i)
            if (nucs[i].name == name) return i;
        return -1;
    }
};

} // namespace decay