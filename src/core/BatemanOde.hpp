// Численная реализация кинетики цепочки: обратный (неявный) Эйлер.
// Неотрицательность и сохранение массы при ЛЮБОМ шаге, устойчива при любом
// разбросе λ. Это важно: "жёсткая" цепочка (U-238 против долей секунды
// у дочерних) убивает явный Эйлер, а этот шаг — нет.
//
// Правило шага: dt ~ dtFactor/λ ограничивается только "активными" нуклидами:
//  - население заметно (выше relThreshold·n0);
//  - и оно НЕ в квазистационарном режиме (население ≈ приток/λ): равновесные
//    члены при dt ≫ 1/λ сама схема выводит точно на равновесие, поэтому не
//    затыкать шаг в их микро-масштабе — иначе счёт шагов уходит в бесконечность.
#pragma once
#include "Nuclide.hpp"
#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

namespace decay {

class BatemanOde {
public:
    explicit BatemanOde(const Chain& chain, double relThreshold = 1e-6,
                        double dtFactor = 0.005, double maxSteps = 2e7)
        : chain_(&chain), rel_(relThreshold), dtF_(dtFactor), maxSteps_(maxSteps) {
        const int n = (int)chain.nucs.size();
        in_.assign(n, {});
        out_.resize(n, 0.0);
        for (int i = 0; i < n; ++i) {
            if (chain.nucs[i].isStable()) continue;
            out_[i] = chain.nucs[i].lambda();
            for (auto& b : chain.nucs[i].branches)
                in_[b.to].push_back({i, chain.nucs[i].lambda() * double(b.prob)});
        }
        n_.assign((size_t)n, 0.0);
    }

    void reset(int source, double n0) {
        n_.assign(chain_->nucs.size(), 0.0);
        n_[(size_t)source] = n0;
        source_ = source;
        n0_ = n0;
        t_ = 0.0;
        steps_ = 0.0;
    }

    void advanceTo(double tTarget) {
        const size_t n = chain_->nucs.size();
        while (t_ < tTarget && steps_ < maxSteps_) {
            double lmax = 0.0;
            for (size_t i = 0; i < n; ++i) {
                if (n_[i] <= rel_ * n0_) continue;
                const double lam = chain_->nucs[i].lambda();
                if (lam <= 0.0) continue;
                double flux = 0.0;
                for (auto& [from, rate] : in_[i]) flux += rate * n_[from];
                const double qs = flux / lam;
                if (std::fabs(n_[i] - qs) <= 0.01 * std::max(qs, rel_ * n0_))
                    continue; // квазистационарный: его масштаб времени не важен
                lmax = std::max(lmax, lam);
            }
            double dt = (lmax > 0.0) ? dtF_ / lmax : (tTarget - t_);
            if (dt > tTarget - t_) dt = tTarget - t_;
            if (dt <= 0.0) break;

            // Обратный Эйлер по треугольной системе (родители раньше детей).
            // x_i' = (x_i + dt·Σ_{k<i} rate_{k→i}·x_k') / (1 + dt·λ_out_i):
            // консервативность точная, положительность безусловная.
            for (size_t i = 0; i < n; ++i) {
                double num = n_[i];
                for (auto& [from, rate] : in_[i])
                    num += dt * rate * n_[from]; // n_[from] уже обновлён — неявно по притоку
                n_[i] = num / (1.0 + dt * out_[i]);
            }
            t_ += dt;
            steps_ += 1.0;
        }
    }

    const std::vector<double>& N() const { return n_; }
    double time() const { return t_; }
    double steps() const { return steps_; }

private:
    const Chain* chain_;
    double rel_, dtF_, maxSteps_;
    int source_ = 0;
    double n0_ = 0.0;
    double t_ = 0.0, steps_ = 0.0;
    std::vector<double> n_;
    std::vector<std::vector<std::pair<int, double>>> in_;
    std::vector<double> out_;
};

} // namespace decay