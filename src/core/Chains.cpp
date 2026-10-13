// Справочные данные: цепочки Po-210, Cs-137, Th-232, U-238.
// Полураспады — по данным IAEA/BNL; малозначимые ветки (<0.1%) опущены
// (β-ветки Po-218 и Bi-214, изомер Pa-234) — так проще, на картину не влияет.
#include "Chains.hpp"
#include <algorithm>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace decay {
namespace units {
    inline constexpr double s(double v)   { return v; }
    inline constexpr double us(double v)  { return v * 1e-6; }
    inline constexpr double min(double v) { return v * 60.0; }
    inline constexpr double h(double v)   { return v * 3600.0; }
    inline constexpr double d(double v)   { return v * 86400.0; }
    inline constexpr double yr(double v)  { return v * 365.25 * 86400.0; }
} // namespace units

namespace {

using NamedBranch = std::pair<const char*, float>; // имя дочернего, вероятность

struct Nl {
    const char* name;
    int         Z;
    int         A;
    double      T;     // период полураспада, с
    const char* mode;
    std::vector<NamedBranch> b; // пусто у стабильных
};

Chain build(const char* id, const char* title, const std::vector<Nl>& raw) {
    Chain c;
    c.id = id;
    c.title = title;
    c.start = 0;
    c.nucs.resize(raw.size());
    for (size_t i = 0; i < raw.size(); ++i) {
        Nuclide& n = c.nucs[i];
        n.name = raw[i].name;
        n.Z = raw[i].Z;
        n.A = raw[i].A;
        n.halfLife = raw[i].T;
        n.mode = raw[i].mode;
    }
    // имена -> индексы, затем вероятности веток
    for (size_t i = 0; i < raw.size(); ++i) {
        for (auto& [to, p] : raw[i].b) {
            int r = c.find(to);
            if (r < 0) throw std::runtime_error("unknown nuclide: " + std::string(to));
            c.nucs[i].branches.push_back({r, p});
        }
    }
    // валидация: у нестабильных должны быть ветки, сумма вероятностей = 1
    for (const Nuclide& n : c.nucs) {
        if (n.isStable()) continue;
        if (n.branches.empty())
            throw std::runtime_error("nuclide has no decay branch: " + n.name);
        double sum = 0.0;
        for (auto& b : n.branches) sum += b.prob;
        if (std::abs(sum - 1.0) > 1e-3)
            throw std::runtime_error("branch probs sum != 1: " + n.name + " (" + std::to_string(sum) + ")");
    }
    // аналитика Бейтмана требует попарно различных λ (иначе вырождение частных производных)
    for (size_t i = 0; i < c.nucs.size(); ++i)
        for (size_t j = i + 1; j < c.nucs.size(); ++j) {
            double li = c.nucs[i].lambda(), lj = c.nucs[j].lambda();
            if (li > 0 && lj > 0 && std::abs(li - lj) < 1e-9 * std::max(li, lj))
                throw std::runtime_error("degenerate lambdas: " + c.nucs[i].name + " vs " + c.nucs[j].name);
        }
    // решатели (обратный Эйлер по подстановке, МК) опираются на топологический
    // порядок: ветка распада уменьшает A, значит должна идти только вперёд по индексам
    for (size_t i = 0; i < c.nucs.size(); ++i)
        for (auto& b : c.nucs[i].branches)
            if (b.to <= (int)i)
                throw std::runtime_error("branch not forward (" + std::to_string(b.to) + " <= " +
                                         std::to_string(i) + "): " + c.nucs[i].name);
    return c;
}

} // namespace

Chain makeLineChain(const std::vector<double>& lambdas) {
    Chain c;
    c.id = "synth";
    c.title = "синтетическая линия";
    c.start = 0;
    c.nucs.resize(lambdas.size() + 1);
    char buf[16];
    for (size_t i = 0; i <= lambdas.size(); ++i) {
        std::snprintf(buf, sizeof buf, "X%zu", i);
        Nuclide& n = c.nucs[i];
        n.name = buf;
        n.Z = 0;
        n.A = int(i);
        if (i < lambdas.size()) {
            n.halfLife = 0.69314718055994530942 / lambdas[i];
            n.mode = "synth";
        } else {
            n.halfLife = std::numeric_limits<double>::infinity();
            n.mode = "";
        }
        if (i + 1 < c.nucs.size()) n.branches.push_back({int(i + 1), 1.f});
    }
    return c;
}

const std::vector<Chain>& library() {
    static const std::vector<Chain> LIB = {
        build("Po-210", "Po-210 to Pb-206 (1 step)", {
            {"Po-210", 84, 210, units::d(138.376), "a", {{"Pb-206", 1.f}}},
            {"Pb-206", 82, 206, units::s(std::numeric_limits<double>::infinity()), "", {}},
        }),
        build("Cs-137", "Cs-137 to Ba-137 (branching)", {
            {"Cs-137",  55, 137, units::yr(30.08), "b", {{"Ba-137m", 0.946f}, {"Ba-137", 0.054f}}},
            {"Ba-137m", 56, 137, units::s(153.1), "isomer", {{"Ba-137", 1.f}}},
            {"Ba-137",  56, 137, units::s(std::numeric_limits<double>::infinity()), "", {}},
        }),
        build("Th-232", "Th-232 to Pb-208 (12, branching)", {
            {"Th-232", 90, 232, units::yr(1.405e10), "a",    {{"Ra-228", 1.f}}},
            {"Ra-228", 88, 228, units::yr(5.75),     "b",   {{"Ac-228", 1.f}}},
            {"Ac-228", 89, 228, units::h(6.15),      "b",   {{"Th-228", 1.f}}},
            {"Th-228", 90, 228, units::yr(1.9116),   "a",    {{"Ra-224", 1.f}}},
            {"Ra-224", 88, 224, units::d(3.6319),    "a",    {{"Rn-220", 1.f}}},
            {"Rn-220", 86, 220, units::s(55.6),      "a",    {{"Po-216", 1.f}}},
            {"Po-216", 84, 216, units::s(0.145),     "a",    {{"Pb-212", 1.f}}},
            {"Pb-212", 82, 212, units::h(10.64),     "b",   {{"Bi-212", 1.f}}},
            {"Bi-212", 83, 212, units::min(60.55),   "b/a", {{"Po-212", 0.6406f}, {"Tl-208", 0.3594f}}},
            {"Po-212", 84, 212, units::us(0.299),    "a",    {{"Pb-208", 1.f}}},
            {"Tl-208", 81, 208, units::min(3.053),   "b",   {{"Pb-208", 1.f}}},
            {"Pb-208", 82, 208, units::s(std::numeric_limits<double>::infinity()), "", {}},
        }),
        build("U-238", "U-238 to Pb-206 (15 nuclides)", {
            {"U-238",   92, 238, units::yr(4.4683e9), "a",  {{"Th-234", 1.f}}},
            {"Th-234",  90, 234, units::d(24.10),     "b", {{"Pa-234m", 1.f}}},
            {"Pa-234m", 91, 234, units::min(1.159),   "b", {{"U-234", 1.f}}},
            {"U-234",   92, 234, units::yr(2.455e5),  "a",  {{"Th-230", 1.f}}},
            {"Th-230",  90, 230, units::yr(7.54e4),   "a",  {{"Ra-226", 1.f}}},
            {"Ra-226",  88, 226, units::yr(1600.0),   "a",  {{"Rn-222", 1.f}}},
            {"Rn-222",  86, 222, units::d(3.8235),    "a",  {{"Po-218", 1.f}}},
            {"Po-218",  84, 218, units::min(3.098),   "a",  {{"Pb-214", 1.f}}},
            {"Pb-214",  82, 214, units::min(26.8),    "b", {{"Bi-214", 1.f}}},
            {"Bi-214",  83, 214, units::min(19.9),    "b", {{"Po-214", 1.f}}},
            {"Po-214",  84, 214, units::us(164.3),    "a",  {{"Pb-210", 1.f}}},
            {"Pb-210",  82, 210, units::yr(22.20),    "b", {{"Bi-210", 1.f}}},
            {"Bi-210",  83, 210, units::d(5.012),     "b", {{"Po-210", 1.f}}},
            {"Po-210",  84, 210, units::d(138.376),   "a",  {{"Pb-206", 1.f}}},
            {"Pb-206",  82, 206, units::s(std::numeric_limits<double>::infinity()), "", {}},
        }),
    };
    return LIB;
}

} // namespace decay