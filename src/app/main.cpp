// Окно, главный цикл, сшивка: цепочка -> DecaySim -> Renderer, панель ImGui.
// Запуск: ./build/decay_sim [--seed N]
#include <glad/gl.h>
#include <GLFW/glfw3.h>

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include "core/Chains.hpp"
#include "core/DecaySim.hpp"
#include "core/Util.hpp"
#include "gfx/Renderer.hpp"
#include "ui/Ui.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <vector>

namespace {

void glfwErrorCb(int, const char* msg) {
    std::fprintf(stderr, "GLFW: %s\n", msg);
}

void hsvToRgb(float hDeg, float s, float v, float* out) {
    const float h = std::fmod(hDeg, 360.0f) / 60.0f; // [0..6)
    const int   i = int(h);
    const float f = h - float(i);
    const float p = v * (1.0f - s);
    const float q = v * (1.0f - f * s);
    const float t = v * (1.0f - (1.0f - f) * s);
    float r, g, b;
    switch (i % 6) {
        case 0: r = v; g = t; b = p; break;
        case 1: r = q; g = v; b = p; break;
        case 2: r = p; g = v; b = t; break;
        case 3: r = p; g = q; b = v; break;
        case 4: r = t; g = p; b = v; break;
        default: r = v; g = p; b = q; break;
    }
    out[0] = r; out[1] = g; out[2] = b;
}

} // namespace

int main(int argc, char** argv) {
    uint32_t seed = 42;
    for (int i = 1; i + 1 < argc; ++i)
        if (std::strcmp(argv[i], "--seed") == 0)
            seed = (uint32_t)std::strtoul(argv[i + 1], nullptr, 0);

    glfwSetErrorCallback(glfwErrorCb);
    if (!glfwInit()) return 1;
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_SAMPLES, 4);
    GLFWwindow* win = glfwCreateWindow(1280, 800, "Decay Chain Simulator", nullptr, nullptr);
    if (!win) {
        glfwTerminate();
        std::fprintf(stderr, "нет окна\n");
        return 1;
    }
    glfwMakeContextCurrent(win);
    glfwSwapInterval(1);
    int glVer = gladLoadGL((GLADloadfunc)glfwGetProcAddress);
    if (glVer < 3300) {
        std::fprintf(stderr, "нужен OpenGL 3.3 (получено %d)\n", glVer);
        glfwTerminate();
        return 1;
    }

    const auto& lib = decay::library();
    ui::Panel panel;
    for (int i = 0; i < (int)lib.size(); ++i)
        if (std::string(lib[i].id) == "U-238") panel.chainIdx = i; // стартовая цепочка

    decay::DecaySim sim;
    decay::Renderer renderer;
    if (!renderer.init()) {
        glfwTerminate();
        return 1;
    }

    std::vector<float> positions; // vec2 на атом
    std::vector<float> palette;   // vec3 на нуклид
    uint32_t  atoms = 0;
    int       builtChain = -1;
    double    simT = 0.0;

    const auto rebuildSim = [&]() {
        const decay::Chain& c = lib[(size_t)panel.chainIdx];
        atoms = panel.atomsTarget();
        sim.reset(c, atoms, seed);
        simT = 0.0;
        panel.lastEvents = 0;

        positions.resize((size_t)atoms * 2);
        {
            uint64_t x = (uint64_t)seed ^ 0x9e3779b97f4a7c15ull ^ uint64_t(panel.chainIdx);
            for (uint32_t i = 0; i < atoms; ++i) {
                const double u = decay::rng::Pcg32::u01(x);
                const double v = decay::rng::Pcg32::u01(x);
                positions[2 * (size_t)i + 0] = float(-1.0 + 2.0 * u);
                positions[2 * (size_t)i + 1] = float(-1.0 + 2.0 * v);
            }
        }
        const int nucs = (int)c.nucs.size();
        palette.assign((size_t)nucs * 3, 1.0f);
        for (int i = 0; i < nucs; ++i) {
            if (c.nucs[(size_t)i].isStable()) {
                palette[3 * (size_t)i + 0] = 0.80f;
                palette[3 * (size_t)i + 1] = 0.85f;
                palette[3 * (size_t)i + 2] = 0.90f;
            } else {
                hsvToRgb(360.0f * std::fmod(float(i) * 0.618033989f, 1.0f), 0.62f, 0.95f,
                         &palette[3 * (size_t)i]);
            }
        }
        renderer.setPalette(palette.data(), nucs);
        renderer.setPositions(positions.data(), atoms);
        builtChain = panel.chainIdx;
    };

    // ImGui
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(win, true);
    ImGui_ImplOpenGL3_Init("#version 330 core");
    {
        int fbW = 0, fbH = 0;
        glfwGetFramebufferSize(win, &fbW, &fbH);
        const float dpi = std::clamp(float(fbH) / 900.0f, 1.0f, 2.5f);
        static const ImWchar ranges[] = {
            0x0020, 0x007E,  // латиница
            0x00D7, 0x00D7,  // ×
            0x00F7, 0x00F7,  // ÷
            0x0391, 0x03C9,  // греческие (α, β, γ)
            0x0400, 0x045F,  // кириллица + ё
            0x2116, 0x2116,  // №
            0
        };
        const char* fontPaths[] = {
            "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
            "/usr/share/fonts/dejavu/DejaVuSans.ttf",
        };
        bool fontLoaded = false;
        for (const char* f : fontPaths)
            if (io.Fonts->AddFontFromFileTTF(f, 15.0f * dpi, nullptr, ranges)) {
                fontLoaded = true;
                break;
            }
        if (!fontLoaded) io.Fonts->AddFontDefault(); // без кириллицы, но работает
    }

    rebuildSim();

    double lastT = glfwGetTime();
    int    fbW = 0, fbH = 0;
    while (!glfwWindowShouldClose(win)) {
        glfwPollEvents();

        int nw = 0, nh = 0;
        glfwGetFramebufferSize(win, &nw, &nh);
        if (nw != fbW || nh != fbH) {
            fbW = nw;
            fbH = nh;
            renderer.setViewport(fbW, fbH);
        }

        // ==== время: настенный ход * акселератор ====
        const double now = glfwGetTime();
        double dtReal = now - lastT;
        lastT = now;
        dtReal = std::min(dtReal, 0.1); // защита от скачка после лага

        if (panel.wantReset || panel.chainIdx != builtChain ||
            panel.atomsTarget() != atoms) {
            rebuildSim();
            panel.wantReset = false;
        }
        if (!panel.paused)
            simT += dtReal * std::pow(10.0, double(panel.accelExp));
        sim.advanceTo(simT);

        const uint64_t evNow = sim.events();
        const uint64_t frameEvents = evNow - panel.lastEvents;
        panel.lastEvents = evNow;

        // ==== ImGui ====
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
        if (!io.WantTextInput) {
            if (ImGui::IsKeyReleased(ImGuiKey_Space)) panel.paused = !panel.paused;
            if (ImGui::IsKeyReleased(ImGuiKey_RightArrow) || ImGui::IsKeyReleased(ImGuiKey_Equal))
                panel.accelExp = std::min(18.0f, panel.accelExp + 1.0f);
            if (ImGui::IsKeyReleased(ImGuiKey_LeftArrow) || ImGui::IsKeyReleased(ImGuiKey_Minus))
                panel.accelExp = std::max(0.0f, panel.accelExp - 1.0f);
            if (ImGui::IsKeyReleased(ImGuiKey_R)) panel.wantReset = true;
        }
        ui::drawPanel(panel, sim, palette.data(), (int)palette.size() / 3, (double)io.Framerate,
                      frameEvents);

        // ==== рендер ====
        if (frameEvents > 0)
            renderer.updateNucs(sim.nuc().data(), atoms);
        renderer.beginFrame();
        renderer.draw(atoms);
        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(win);
    }

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(win);
    glfwTerminate();
    return 0;
}