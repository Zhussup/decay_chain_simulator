#include "Renderer.hpp"
#include <glad/gl.h>
#include <algorithm>
#include <cstdio>

namespace {

const char* kVertexSrc = R"(#version 330 core
layout (location = 0) in vec2  aPos;
layout (location = 1) in uint  aNuc;

uniform float uAspect;
uniform vec3  uPal[24];
uniform float uSize;

out vec3 vCol;

void main() {
    gl_Position = vec4(aPos.x * uAspect, aPos.y, 0.0, 1.0);
    gl_PointSize = uSize;
    vCol = uPal[int(min(aNuc, uint(23)))];
})";

const char* kFragmentSrc = R"(#version 330 core
in vec3 vCol;
out vec4 o;

uniform float uAlpha;

void main() {
    vec2 d = gl_PointCoord - 0.5;
    float r2 = dot(d, d);
    if (r2 > 0.25) discard;
    o = vec4(vCol, uAlpha * (1.0 - 4.0 * r2));
})";

unsigned int compile(unsigned int type, const char* src, char err[512]) {
    unsigned int sh = glCreateShader(type);
    glShaderSource(sh, 1, &src, nullptr);
    glCompileShader(sh);
    int ok = 0;
    glGetShaderiv(sh, GL_COMPILE_STATUS, &ok);
    if (!ok) {
        glGetShaderInfoLog(sh, 512, nullptr, err);
        glDeleteShader(sh);
        return 0;
    }
    return sh;
}

} // namespace

namespace decay {

bool Renderer::init() {
    char err[512] = {};

    unsigned int vs = compile(GL_VERTEX_SHADER, kVertexSrc, err);
    unsigned int fs = vs ? compile(GL_FRAGMENT_SHADER, kFragmentSrc, err) : 0;
    if (!fs) glDeleteShader(vs);
    if (!fs) {
        std::fprintf(stderr, "shader: %s\n", err);
        return false;
    }

    prog_ = glCreateProgram();
    glAttachShader(prog_, vs);
    glAttachShader(prog_, fs);
    glLinkProgram(prog_);
    glDeleteShader(vs);
    glDeleteShader(fs);
    int ok = 0;
    glGetProgramiv(prog_, GL_LINK_STATUS, &ok);
    if (!ok) {
        glGetProgramInfoLog(prog_, 512, nullptr, err);
        std::fprintf(stderr, "program: %s\n", err);
        return false;
    }
    locAspect_ = glGetUniformLocation(prog_, "uAspect");
    locSize_   = glGetUniformLocation(prog_, "uSize");
    locPal_    = glGetUniformLocation(prog_, "uPal");
    locAlpha_  = glGetUniformLocation(prog_, "uAlpha");

    glGenVertexArrays(1, &vao_);
    glGenBuffers(1, &vboPos_);
    glGenBuffers(1, &vboNuc_);
    glBindVertexArray(vao_);
    glBindBuffer(GL_ARRAY_BUFFER, vboPos_);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), nullptr);
    glEnableVertexAttribArray(0);
    glVertexAttribDivisor(0, 1);
    glBindBuffer(GL_ARRAY_BUFFER, vboNuc_);
    glVertexAttribIPointer(1, 1, GL_UNSIGNED_SHORT, 2, nullptr);
    glEnableVertexAttribArray(1);
    glVertexAttribDivisor(1, 1);
    glBindVertexArray(0);
    palCount_ = 0;
    return true;
}

void Renderer::setPalette(const float* rgb, int nucCount) {
    palCount_ = std::min(nucCount, kMaxNuc);
    if (palCount_ > 0)
        std::copy(rgb, rgb + palCount_ * 3, pal_);
}

void Renderer::setPositions(const float* xy, uint32_t atoms) {
    glBindBuffer(GL_ARRAY_BUFFER, vboPos_);
    glBufferData(GL_ARRAY_BUFFER, (GLsizeiptr)atoms * 2 * sizeof(float), xy, GL_STATIC_DRAW);
}

void Renderer::updateNucs(const std::uint16_t* ids, uint32_t atoms) {
    glBindBuffer(GL_ARRAY_BUFFER, vboNuc_);
    glBufferData(GL_ARRAY_BUFFER, (GLsizeiptr)atoms * sizeof(uint16_t), ids, GL_STREAM_DRAW);
}

void Renderer::setViewport(int fbW, int fbH) {
    fbW_ = fbW;
    fbH_ = fbH;
    glViewport(0, 0, fbW, fbH);
}

void Renderer::beginFrame() const {
    glClearColor(0.055f, 0.062f, 0.078f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

void Renderer::draw(uint32_t atoms) const {
    glUseProgram(prog_);
    const float aspect = fbH_ > 0 ? float(fbW_) / float(fbH_) : 1.0f;
    glUniform1f(locAspect_, aspect);
    glUniform1f(locSize_, std::clamp(float(fbH_) / 300.0f, 2.0f, 6.0f));
    glUniform1f(locAlpha_, 0.9f);
    glUniform3fv(locPal_, kMaxNuc, pal_);
    glBindVertexArray(vao_);
    glDrawArraysInstanced(GL_POINTS, 0, 1, (GLint)atoms);
    glBindVertexArray(0);
}

} // namespace decay