#include "Ui.hpp"
#include "core/Util.hpp"
#include "imgui.h"
#include <cstdio>
#include <string>
#include <vector>

namespace ui {

void drawPanel(Panel& p, const decay::DecaySim& sim, const float* palette, int palCount,
               double fps, uint64_t frameEvents) {
    const auto& lib = decay::library();

    ImGuiIO& io = ImGui::GetIO();
    const float w = 370.0f, hgt = io.DisplaySize.y - 24.0f;
    ImGui::SetNextWindowPos(ImVec2(io.DisplaySize.x - w - 12.0f, 12.0f), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(w, hgt), ImGuiCond_FirstUseEver);
    ImGui::Begin("Control", nullptr, ImGuiWindowFlags_NoCollapse);
    ImGui::TextDisabled("Space = pause, left/right accel order, R or Reset button = restart");
    ImGui::Separator();

    // chain list (titles collected once)
    static std::vector<std::string> titleStorage;
    static std::vector<const char*> titlePtrs;
    if (titlePtrs.empty()) {
        for (const auto& c : lib) {
            titleStorage.push_back(c.title);
            titlePtrs.push_back(titleStorage.back().c_str());
        }
    }
    if (ImGui::Combo("Chain", &p.chainIdx, titlePtrs.data(), (int)titlePtrs.size()))
        p.wantReset = true;

    ImGui::SliderFloat("Particles", &p.atomsLog, 2.0f, 6.0f, "10^%.1f");
    if (ImGui::IsItemDeactivatedAfterEdit())
        p.wantReset = true;
    ImGui::TextDisabled("in view: %s particles", decay::humanCount(p.atomsTarget()).c_str());

    // ==== ACCELERATOR ====
    ImGui::SeparatorText("Time accelerator");
    ImGui::SliderFloat("order##acc", &p.accelExp, 0.0f, 18.0f, "x10^%.0f");
    if (ImGui::Button("/10")) p.accelExp = std::max(0.0f, p.accelExp - 1.0f);
    ImGui::SameLine();
    if (ImGui::Button("x10")) p.accelExp = std::min(18.0f, p.accelExp + 1.0f);
    ImGui::SameLine();
    if (ImGui::Button(p.paused ? "Resume" : "Pause")) p.paused = !p.paused;

    const double factor = std::pow(10.0, double(p.accelExp));
    ImGui::Text("rate: %s", decay::humanRatePerSec(factor).c_str());
    ImGui::Text("sim time: %s%s", decay::humanTime(sim.simTime()).c_str(),
                p.paused ? "  (paused)" : "");
    ImGui::Text("decays: %s total, %s this frame",
                decay::humanCount(double(sim.events())).c_str(),
                decay::humanCount(double(frameEvents)).c_str());

    // restart from scratch: fresh atoms scattered anew
    ImGui::SeparatorText("State");
    ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.62f, 0.20f, 0.18f, 1.0f));
    ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.78f, 0.28f, 0.24f, 1.0f));
    ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(0.88f, 0.34f, 0.28f, 1.0f));
    if (ImGui::Button("Reset"))
        p.wantReset = true;
    if (ImGui::IsItemHovered())
        ImGui::SetTooltip("rebuild the simulation: fresh particles, t = 0");
    ImGui::PopStyleColor(3);
    ImGui::SameLine();
    ImGui::TextDisabled("same seed and chain, back to t = 0");

    ImGui::SeparatorText("Nuclides");
    if (ImGui::BeginTable("counts", 4, ImGuiTableFlags_SizingFixedFit)) {
        ImGui::TableSetupColumn("clr",  ImGuiTableColumnFlags_WidthFixed, 14.0f);
        ImGui::TableSetupColumn("name", ImGuiTableColumnFlags_WidthFixed, 62.0f);
        ImGui::TableSetupColumn("N",    ImGuiTableColumnFlags_WidthFixed, 74.0f);
        ImGui::TableSetupColumn("half-life", ImGuiTableColumnFlags_WidthFixed, 96.0f);

        for (int i = 0; i < (int)sim.chain()->nucs.size() && i < palCount; ++i) {
            const decay::Nuclide& n = sim.chain()->nucs[(size_t)i];
            ImGui::PushID(i);
            ImGui::TableNextRow();
            // color chip: same palette the particles use
            ImGui::TableNextColumn();
            const ImVec2 pmin = ImGui::GetCursorScreenPos();
            ImU32 col = ImGui::GetColorU32(ImVec4(palette[3 * i + 0], palette[3 * i + 1],
                                                  palette[3 * i + 2], 1.0f));
            ImGui::GetWindowDrawList()->AddRectFilled(
                pmin, ImVec2(pmin.x + 11.0f, pmin.y + 11.0f), col, 3.0f);

            ImGui::TableNextColumn();
            ImGui::TextUnformatted(n.name.c_str());
            ImGui::TableNextColumn();
            char cnt[24];
            std::snprintf(cnt, sizeof cnt, "%lld", (long long)sim.counts()[(size_t)i]);
            ImGui::TextUnformatted(cnt);
            ImGui::TableNextColumn();
            if (n.isStable()) {
                ImGui::TextUnformatted("stable");
            } else {
                char ttl[48];
                std::snprintf(ttl, sizeof ttl, "%s %s", decay::humanTime(n.halfLife).c_str(), n.mode);
                ImGui::TextUnformatted(ttl);
            }
            ImGui::PopID();
            ImGui::Dummy(ImVec2(0.0f, 1.0f));
        }
        ImGui::EndTable();
    }

    ImGui::Text("%.0f FPS", fps);
    ImGui::End();
}

} // namespace ui