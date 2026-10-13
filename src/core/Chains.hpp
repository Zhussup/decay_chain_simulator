// Библиотека цепочек с данными, забитыми вручную (без ENDF-парсинга).
#pragma once
#include "Nuclide.hpp"

namespace decay {

// Готовый набор цепочек; порядок индексов фиксирован.
const std::vector<Chain>& library();

// Вспомогательный конструктор для тестов: однолинейная цепочка по λ, последний член — стабильный.
Chain makeLineChain(const std::vector<double>& lambdas /* 1/с */);

} // namespace decay