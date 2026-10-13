// Стохастическое моделирование: 1 атом = 1 частица, событийно-ориентированный МК.
// У каждого атома своё расписание: следующий распад в момент tNext ~ Exp(λ).
// Атом независимо продвигается от t до конца окна (схема Гиллеспи на атом),
// поэтому жёсткость цепочки (разброс λ на ~30 порядков) схеме не страшна —
// ровно того, чего требует акселератор времени ×N.
// Состояние атома: индекс нуклида + состояние PCG — 1:1 переносится на GPU.
#pragma once
#include "Nuclide.hpp"
#include "Rng.hpp"
#include <cstdint>
#include <vector>

namespace decay {

class DecaySim {
public:
    void reset(const Chain& chain, uint32_t atoms, uint64_t seed);
    void advanceTo(double simT); // обработать все события атомов до simT

    const Chain* chain() const { return chain_; }
    uint32_t atoms() const { return atoms_; }
    double simTime() const { return t_; }
    uint64_t events() const { return events_; }

    const std::vector<int64_t>& counts() const { return counts_; }
    const std::vector<uint16_t>& nuc() const { return nuc_; } // для рендера, id на атом

private:
    const Chain*    chain_ = nullptr;
    uint32_t        atoms_ = 0;
    double          t_ = 0.0;
    uint64_t        events_ = 0;
    std::vector<uint16_t> nuc_;
    std::vector<double>   tNext_;   // время следующего распада атома (+inf у стабильных)
    std::vector<uint64_t> rng_;     // PCG state per atom
    std::vector<int64_t>  counts_;
};

} // namespace decay