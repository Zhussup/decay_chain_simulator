// Рендер частиц: инстансинг точек (gl_Points), цвет — из палитры по индексу нуклида.
// Позиция атома статична, каждый кадр обновляется только массив индексов (2 байта/атом).
#pragma once
#include <cstdint>

namespace decay {

class Renderer {
public:
    bool init();
    void setPalette(const float* rgb, int nucCount);    // по 3 float на нуклид
    void setPositions(const float* xy, uint32_t atoms); // по 2 float на атом
    void updateNucs(const std::uint16_t* ids, uint32_t atoms);
    void setViewport(int fbW, int fbH);
    void beginFrame() const; // очистка + состояние blend
    void draw(uint32_t atoms) const;

private:
    static constexpr int kMaxNuc = 24;

    unsigned int prog_ = 0;
    unsigned int vao_ = 0, vboPos_ = 0, vboNuc_ = 0;
    float pal_[kMaxNuc * 3] = {};
    int   palCount_ = 0;
    int   fbW_ = 0, fbH_ = 0;
    int   locAspect_ = -1, locSize_ = -1, locPal_ = -1, locAlpha_ = -1;
};

} // namespace decay