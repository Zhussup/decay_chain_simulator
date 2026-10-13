// Точное решение цепочки Бейтмана в лог-домене.
// Стойко к разбросу λ (жёсткие цепочки), ветвление учитывается суммой по путям
// (цепочки распада — DAG: массовое число только убывает).
// Стабильные члены отдельно не решаются: N(стаб) = n0 − Σ N(нестаб) — точно.
//
// Обусловленность: сумма Бейтмана — знакопеременная; при малых t·λ доминирующий
// член может на много порядков превышать результат. condition(idx) даёт оценку
// |доминирующий член| / |результат|; если она велика, double начал терять цифры
// и надо сверяться с ODE/МК (это ровно та «жёсткость», о которой была оговорка).
#pragma once
#include "Nuclide.hpp"

namespace decay {

class BatemanAnalytic {
public:
    explicit BatemanAnalytic(const Chain& chain, double n0);

    // ожидаемое число атомов нуклида idx в момент t
    double N(int idx, double t) const;

    // оценка подавления в знакопеременной сумме; 1.0 — идеал, >>1 — риск потери точности
    double condition(int idx, double t) const;

    // выше этого порога результат считается ненадёжным
    static constexpr double kMaxCondition = 1e12;

    bool valid(int idx, double t) const { return condition(idx, t) < kMaxCondition; }

private:
    struct NodeRef { int idx; double w; };   // вход в нуклид с весом ветки (у источника w=1)
    struct Path   { std::vector<NodeRef> node; };

    void dfs(const Chain& c, int cur, std::vector<NodeRef>& stack,
             std::vector<std::vector<Path>>& out) const;

    // точное решение для одного пути (без n0): сумма (знак, ln|член|)
    void pathTerms(const Path& p, double t, std::vector<std::pair<double,int>>& terms) const;

    const Chain* chain_ = nullptr;
    int    source_ = 0;
    double n0_     = 1.0;
    std::vector<std::vector<Path>> pathsTo_; // по целевому нуклиду
};

} // namespace decay