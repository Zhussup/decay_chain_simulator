# Decay Chain Simulator

A college project: radioactive decay physics (the Bateman equations) with an
interactive OpenGL visualization. The key idea: **one atom = one particle**.
Probabilistic Monte-Carlo simulation and the analytic Bateman solution are
computed for the same initial conditions and overlaid — the curves matching is
physical validation of both models. The centerpiece feature is the **time
accelerator ×10^N** (a dial up to ×10^18): the event-driven MC is built so that
chain stiffness (U-238 lives for 4.5 billion years while its daughters decay
in fractions of a second) slows down neither the simulation nor the visuals.

## Running

```bash
cmake -S . -B build
cmake --build build -j4
./build/decay_sim            # graphics: particles + ImGui panel
./build/bateman_test         # console physics check (exit 0 = everything agrees)
```

Nuclide data (half-lives and branching) — IAEA/BNL; negligible branches (<0.1%)
are omitted. Dependencies: system glfw3, OpenGL; glad/ImGui/glm are vendored
in `external/`.

## Controls

| Input                  | Action                                    |
|------------------------|-------------------------------------------|
| Space                  | pause / resume                            |
| Left / Right arrows    | accelerator order ÷10 / ×10               |
| `=` / `-`              | same (keyboard)                           |
| R                      | reset simulation                          |
| Reset button (panel)   | reset simulation                          |
| Chain combo            | Po-210, Cs-137, Th-232, U-238             |
| Particles slider       | 10^2 … 10^6 particles                     |
| Accelerator slider     | ×10^0 … ×10^18                            |

## Structure

```
src/core/               physics (no graphics)
  Nuclide.hpp             nuclide/branch/chain model (+find())
  Chains.hpp/.cpp         chain library: Po-210, Cs-137, Th-232, U-238;
                          makeLineChain() for synthetic tests; validation
                          (branches only forward = topological order,
                          branch probabilities sum to 1, distinct lambdas)
  BatemanAnalytic.hpp/.cpp  analytic Bateman solution in log domain
  BatemanOde.hpp          backward (implicit) Euler over the triangular system
  DecaySim.hpp/.cpp       event-driven MC: one atom = one particle, tNext ~ Exp(lambda)
  Rng.hpp                 PCG32 (XSH-RR 64/32) per atom + splitmix64 seeding
  Util.hpp/.cpp           human-readable times/rates/counts
src/gfx/Renderer.hpp/.cpp  OpenGL 3.3: instanced points, per-nuclide palette
src/ui/Ui.hpp/.cpp          ImGui panel (accelerator, nuclide table, reset)
src/app/main.cpp            window, main loop, wiring
tests/test_bateman.cpp      three-referee cross-check (analytic <-> ODE <-> MC)
external/                   glad2 (generated), imgui 1.92.1, glm (unused so far)
```

## How it works

**Bateman analytics.** For every nuclide we enumerate all decay paths from the
source (a chain is a DAG: the mass number only decreases, so branching is
explicit). Each sum term is `(prod lambda*w) * n0 * exp(-lambda_k t) /
prod_{j!=k}(lambda_j - lambda_k)`; the sum is evaluated in the log domain
(signed log-sum-exp with a `long double` accumulator) so that huge mutually
cancelling terms do not destroy each other in double precision. Stable members
follow from mass conservation. A condition estimate
(|largest term| / |result|) with a 1e12 cutoff decides where the analytics can
be trusted: there it acts as the referee; elsewhere the MC does.

**The ODE reference.** Backward Euler: the chain is already in topological
order, so the system is triangular and solves by forward substitution — mass
conservation is exact, positivity is unconditional, any lambda spread is fine.
The step `dt = 0.005/lambda_max` is limited only by "active" nuclides
(population above 1e-6*n0 and not in quasi-steady state) — otherwise
secular-equilibrium members (Ba-137m) stall the step at their micro-scale.
U-238 over 4.5 Gyr takes ~156 steps.

**The MC.** Every atom carries its own schedule: its next decay happens at
`tNext ~ Exp(lambda)` (inverse-transform sampling) drawn from a per-atom PCG32.
An atom's state (nuclide id + RNG state, 12 bytes) is designed to port 1:1 to a
compute shader (planned: the same code on the GPU). Stiffness does not matter:
an atom lives from event to event.

**Validation.** `bateman_test` runs synthetic chains (linear, branching, stiff
with lambdas from 1e-9 to 1e8) and the real chains at characteristic times;
comparisons are scored in MC sigmas (tolerance 4.5 sigma). Current status:
**pass 115, fail 0, skip 60** (skips are species with an expectation below 50
atoms — MC noise — and points with poor analytic conditioning).

## Next (plan steps 4–6)

1. The decay chain graph on screen: nodes, half-life labels, click to highlight.
2. Port the MC to a compute shader (same PCG32, per-atom state 1:1) —
   10^7–10^8 atoms in real time.
3. Branching / atom history: a particle's path highlighted along the chain.

## License

This project is licensed under the [MIT License](LICENSE). Third-party
libraries vendored in `external/` keep their own licenses: ImGui — MIT
(`external/imgui/LICENSE.txt`); glad (generated loader) — see the SPDX header
in `external/glad/src/gl.c`; Khronos platform header — its own permissive
notice (`external/glad/include/KHR/khrplatform.h`); GLM —
`external/glm/copying.txt`.