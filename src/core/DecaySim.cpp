#include "DecaySim.hpp"

namespace decay {

void DecaySim::reset(const Chain& chain, uint32_t atoms, uint64_t seed) {
    chain_ = &chain;
    atoms_ = atoms;
    t_ = 0.0;
    events_ = 0;
    nuc_.assign(atoms, (uint16_t)chain.start); // все атомы стартуют в начале цепочки
    tNext_.assign(atoms, std::numeric_limits<double>::infinity());
    rng_.resize(atoms);
    counts_.assign(chain.nucs.size(), 0);
    counts_[(size_t)chain.start] = atoms;

    uint64_t x = seed;
    const double lam0 = chain.nucs[(size_t)chain.start].lambda();
    for (uint32_t i = 0; i < atoms; ++i) {
        rng_[i] = rng::splitmix64(x);
        uint64_t& st = rng_[i];
        tNext_[i] = std::isfinite(lam0) ? rng::Pcg32::expo(st, lam0)
                                        : std::numeric_limits<double>::infinity();
    }
}

void DecaySim::advanceTo(double simT) {
    const double end = simT;
    if (end <= t_)
        return;

    for (uint32_t i = 0; i < atoms_; ++i) {
        uint16_t cur = nuc_[i];
        double tEv = tNext_[i];
        uint64_t st = rng_[i]; // держим состояние атома в регистре до конца его окна

        while (tEv <= end) {
            const Nuclide& n = chain_->nucs[cur];
            // выбор ветки по кумулятивным вероятностям; последний вариант — по умолчанию
            double r = rng::Pcg32::u01(st);
            double acc = 0.0;
            int to = n.branches.back().to;
            for (auto& b : n.branches) {
                acc += b.prob;
                if (r < acc) { to = b.to; break; }
            }
            counts_[(size_t)cur]--;
            counts_[(size_t)to]++;
            events_++;

            if (chain_->nucs[to].isStable()) {
                cur = (uint16_t)to;
                tEv = std::numeric_limits<double>::infinity();
                break;
            }
            cur = (uint16_t)to;
            tEv += rng::Pcg32::expo(st, chain_->nucs[to].lambda());
        }

        nuc_[i] = cur;
        tNext_[i] = tEv;
        rng_[i] = st;
    }
    t_ = end;
}

} // namespace decay