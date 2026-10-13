// Проверка физики ядра (консоль): аналитика Бейтмана в лог-домене ↔
// схема Патанкара ↔ событийный МК. Запуск: ./build/bateman_test
// Правила сравнения:
//  - МК — статистически точный судья (σ = sqrt(E), допуск 4.5σ);
//  - виды с ожиданием < 50 атомов не проверяем (шум МК);
//  - аналитику сверяем только там, где её условие (подавление в сумме) хорошее;
//    иначе она ненадёжна, и этот вид проверяется парой МК ↔ ОП;
//  - Патанкар: сохранение массы точно; точность 2% (источник/стабильный),
//    5% (промежуточные) — эталон на общих основаниях.
#include "core/Chains.hpp"
#include "core/BatemanAnalytic.hpp"
#include "core/BatemanOde.hpp"
#include "core/DecaySim.hpp"
#include "core/Util.hpp"
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

using namespace decay;

static int gPass = 0, gFail = 0, gSkip = 0;

static void verdict(bool ok, const std::string& what) {
    std::printf("  [%s] %s\n", ok ? "PASS" : "FAIL", what.c_str());
    if (ok) ++gPass; else ++gFail;
}
static void skipRow(const std::string& what) {
    std::printf("  [SKIP] %s\n", what.c_str());
    ++gSkip;
}

static Chain makeBranchedChain() {
    // X0 →(0.7) X1 → X3; X0 →(0.3) X2 → X3; X3 → X4 → X5(стаб.)
    Chain c;
    c.id = "synthB";
    c.title = "синтетическая ветвящаяся";
    c.start = 0;
    const double lam[4] = {2.0, 1.5, 0.8, 0.5};
    c.nucs.resize(5);
    for (int i = 0; i < 4; ++i) {
        Nuclide& n = c.nucs[i];
        n.name = "X" + std::to_string(i);
        n.A = i;
        n.halfLife = 0.69314718055994530942 / lam[i];
        n.mode = "synth";
    }
    c.nucs[4].name = "X5";
    c.nucs[4].A = 5;
    c.nucs[4].halfLife = std::numeric_limits<double>::infinity();
    c.nucs[0].branches.assign({Branch{1, 0.7f}, Branch{2, 0.3f}});
    c.nucs[1].branches.assign({Branch{3, 1.f}});
    c.nucs[2].branches.assign({Branch{3, 1.f}});
    c.nucs[3].branches.assign({Branch{4, 1.f}});
    return c;
}

static void compareAt(const Chain& c, double t, double nAtoms, uint64_t seed, bool withOde) {
    BatemanAnalytic an(c, (double)nAtoms);

    DecaySim mc;
    mc.reset(c, (uint32_t)nAtoms, seed);
    mc.advanceTo(t);

    BatemanOde ode(c);
    ode.reset(c.start, (double)nAtoms);
    ode.advanceTo(t);

    std::printf("  %s, t = %s\n", c.id, humanTime(t).c_str());

    // консервация общей массы у всех трёх
    double sumAn = 0.0, sumOde = 0.0;
    long long sumMc = 0;
    for (int i = 0; i < (int)c.nucs.size(); ++i) {
        sumAn  += an.N(i, t);
        sumOde += withOde ? ode.N()[i] : 0.0;
        sumMc  += mc.counts()[(size_t)i];
    }
    verdict(std::fabs(sumMc - nAtoms) < 0.5, std::string(c.id) + ": МК сохранил число атомов");
    if (withOde)
        // структурную утечку схема не даёт; остаток — только накопление
        // округлений double за ~10^5-10^6 шагов (~1e-9 относит.)
        verdict(std::fabs(sumOde - nAtoms) <= 1e-7 * nAtoms,
                std::string(c.id) + ": ОП сохранил массу (ост. " + std::to_string(sumOde - nAtoms) +
                ", шагов " + std::to_string((long long)ode.steps()) + ")");

    const int n = (int)c.nucs.size();
    for (int i = 0; i < n; ++i) {
        const std::string nn   = c.nucs[i].name;
        const bool       stabl = c.nucs[i].isStable();
        const double     expect = an.N(i, t);
        const double     mcv    = double(mc.counts()[i]);
        const double     cond   = an.condition(i, t);

        // стабильный член считается вычитанием суммы нестабильных: испорченные
        // (отрицательные, несходящиеся) нестабильные тянут за собой стабильный
        bool anOK = an.valid(i, t);
        if (stabl && anOK) {
            double sun = 0.0;
            bool   bad = false;
            for (int j = 0; j < n; ++j)
                if (!c.nucs[j].isStable()) {
                    const double nv = an.N(j, t);
                    if (nv < 0.0) bad = true;
                    sun += nv;
                }
            if (bad || sun > nAtoms * (1.0 + 1e-6))
                anOK = false;
        }

        if (!anOK) {
            // аналитика плохо обусловлена: сверяем МК с ОП, если там что-то видно
            if (withOde && ode.N()[i] >= 100.0) {
                double sig = std::sqrt(ode.N()[i]);
                double dev = std::fabs(mcv - ode.N()[i]) / sig;
                verdict(dev <= 4.5, nn + ": МК " + std::to_string(mcv) + " vs ОП " +
                                        std::to_string(ode.N()[i]) + " (" + std::to_string(dev) + "s), усл.аналитики " + std::to_string(cond));
            } else {
                skipRow(nn + ": условие аналитики " + std::to_string(cond) +
                        " — надёжной сверки нет");
            }
            continue;
        }
        if (expect < 50.0) {
            skipRow(nn + ": ожидание " + std::to_string(expect) + " < 50 (шум МК)");
            continue;
        }
        // σ для МК-счёта: у нестабильного — Пуассона; у стабильного — биномиальный
        // (счёт = сумма независимых индикаторов «атом дошёл до конца»)
        const double sig = stabl
            ? std::sqrt(std::max(nAtoms * (expect / nAtoms) * (1.0 - expect / nAtoms), 1.0))
            : std::sqrt(std::max(expect, 1.0));
        const double dev = std::fabs(mcv - expect) / sig;
        char row[200];
        std::snprintf(row, sizeof row, "%s: МК %.5g vs ожид %.5g (%.1fσ), ОП %.5g, усл %.0e",
                      nn.c_str(), mcv, expect, dev, withOde ? ode.N()[i] : 0.0, cond);
        verdict(dev <= 4.5, row);

        if (withOde && expect >= 1e-3 * nAtoms) {
            const double rel  = std::fabs(ode.N()[i] - expect) / expect;
            const double tol  = (i == c.start || c.nucs[i].isStable()) ? 0.02 : 0.05;
            std::snprintf(row, sizeof row, "%s: ОП против аналитики, rel %.2e", nn.c_str(), rel);
            verdict(rel <= tol, row);
        }
    }
}

int main() {
    std::printf("=== Синтетика: однолинейные цепочки ===\n");
    {
        Chain c1 = makeLineChain({1.0, 0.3, 0.05});
        for (double t : {0.5, 3.0, 30.0})
            compareAt(c1, t, 2e5, 11, true);

        Chain c2 = makeLineChain({1e8, 1.0, 1e-9}); // жёсткая: τ от 1e-8 до 1e9 с
        for (double t : {1e-3, 10.0, 1e3})
            compareAt(c2, t, 2e5, 22, false); // Патанкар на таком разбросе шагов не прогонишь — только МК ↔ аналитика
    }

    std::printf("=== Синтетика: ветвление ===\n");
    {
        Chain cb = makeBranchedChain();
        for (double t : {0.5, 2.0})
            compareAt(cb, t, 2e5, 33, true);
    }

    std::printf("=== Реальные цепочки ===\n");
    const double YR = 365.25 * 86400.0;
    {
        const auto& lib = library();
        auto byId = [&](const char* id) -> const Chain* {
            for (auto& c : lib)
                if (std::string(c.id) == id) return &c;
            return nullptr;
        };

        const Chain* po = byId("Po-210");
        for (double t : {30.0 * 86400.0, 276.752 * 86400.0})
            compareAt(*po, t, 5e5, 41, true);

        // старт прямо с изомера Ba-137m — короткая чистая цепочка из 2 звеньев
        {
            Chain cc = *byId("Cs-137");
            int iso = cc.find("Ba-137m");
            cc.start = iso;
            compareAt(cc, 153.1, 5e5, 51, true);
            compareAt(cc, 306.2, 5e5, 52, true);
        }

        const Chain* cs = byId("Cs-137");
        for (double t : {86400.0, 3.0 * YR})
            compareAt(*cs, t, 2e6, 61, true);

        const Chain* u238 = byId("U-238");
        for (double t : {1e9 * YR, 4.4683e9 * YR})
            compareAt(*u238, t, 5e5, 71, true);

        const Chain* th232 = byId("Th-232");
        for (double t : {1e9 * YR, 4e9 * YR})
            compareAt(*th232, t, 5e5, 81, true);
    }

    std::printf("ИТОГО: pass %d, fail %d, skip %d\n", gPass, gFail, gSkip);
    return gFail > 0 ? 1 : 0;
}