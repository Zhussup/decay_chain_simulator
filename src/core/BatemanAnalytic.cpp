#include "BatemanAnalytic.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace decay {

BatemanAnalytic::BatemanAnalytic(const Chain& chain, double n0)
    : chain_(&chain), source_(chain.start), n0_(n0) {
    pathsTo_.assign(chain.nucs.size(), {});
    std::vector<NodeRef> stack;
    stack.push_back(NodeRef{source_, 1.0});
    dfs(chain, source_, stack, pathsTo_);
}

void BatemanAnalytic::dfs(const Chain& c, int cur, std::vector<NodeRef>& stack,
                          std::vector<std::vector<Path>>& out) const {
    out[cur].push_back(Path{stack}); // путь из источника в текущий нуклид
    if (c.nucs[cur].isStable())
        return;
    for (auto& b : c.nucs[cur].branches) {
        stack.push_back(NodeRef{b.to, double(b.prob)});
        dfs(c, b.to, stack, out);
        stack.pop_back();
    }
}

void BatemanAnalytic::pathTerms(const Path& p, double t,
                                std::vector<std::pair<double,int>>& terms) const {
    const size_t L = p.node.size();
    // лог произведения λ·w по переходам (в каждом члене суммы)
    double lnProd = 0.0;
    for (size_t r = 1; r < L; ++r) {
        double lamParent = chain_->nucs[p.node[r - 1].idx].lambda();
        double w         = p.node[r].w;
        lnProd += std::log(lamParent) + std::log(w);
    }
    for (size_t k = 0; k < L; ++k) {
        // член Бейтмана: (Π λ·w) · n0 · e^{−λk t} / Π_{j≠k}(λj−λk) — именно деление
        double lamK = chain_->nucs[p.node[k].idx].lambda();
        int    sign = 1;
        double lnAbs = lnProd + std::log(n0_) - lamK * t;
        for (size_t j = 0; j < L; ++j) {
            if (j == k) continue;
            double lamJ = chain_->nucs[p.node[j].idx].lambda();
            double d = lamJ - lamK;
            sign *= (d > 0.0) ? 1 : -1;
            lnAbs -= std::log(std::abs(d));
        }
        terms.push_back({lnAbs, sign});
    }
}

double BatemanAnalytic::N(int idx, double t) const {
    const Nuclide& nuc = chain_->nucs[idx];
    if (nuc.isStable()) {
        // сохранение массы: всё, что ещё не в нестабильных членах, уже в стабильных
        double sum = n0_;
        for (int i = 0; i < (int)chain_->nucs.size(); ++i)
            if (!chain_->nucs[i].isStable())
                sum -= N(i, t); // нестабильные не вызывают стабильные — рекурсии нет
        return sum;
    }
    std::vector<std::pair<double,int>> terms;
    for (auto& p : pathsTo_[idx])
        pathTerms(p, t, terms);

    double m = -std::numeric_limits<double>::infinity();
    for (auto& [lnAbs, sign] : terms)
        m = std::max(m, lnAbs);
    if (m == -std::numeric_limits<double>::infinity())
        return 0.0; // все члены ушли в ноль

    long double acc = 0.0L;
    for (auto& [lnAbs, sign] : terms)
        acc += (long double)sign * std::expl((long double)(lnAbs - m));

    if (acc <= 0.0L)
        return 0.0; // исчезающе малая величина или потеря точности
    return std::exp((double)(m + std::log((double)acc)));
}

double BatemanAnalytic::condition(int idx, double t) const {
    if (chain_->nucs[idx].isStable())
        return 1.0; // считается через сохранение — своего подавления нет
    std::vector<std::pair<double,int>> terms;
    for (auto& p : pathsTo_[idx])
        pathTerms(p, t, terms);
    double m = -std::numeric_limits<double>::infinity();
    for (auto& [lnAbs, sign] : terms)
        m = std::max(m, lnAbs);
    if (m == -std::numeric_limits<double>::infinity())
        return 1.0;
    long double acc = 0.0L;
    for (auto& [lnAbs, sign] : terms)
        acc += (long double)sign * std::expl((long double)(lnAbs - m));
    // доминирующий член масштабирует сумму в |acc| раз меньше себя:
    // подавление = |член_max| / |Σ| = 1/|acc|
    const double la = std::fabs((double)acc);
    if (la <= 0.0)
        return std::numeric_limits<double>::infinity();
    return 1.0 / la;
}

} // namespace decay