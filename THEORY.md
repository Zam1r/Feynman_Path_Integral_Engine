# What the path integral is actually testing

This is the short version of the theory behind the program, and the one result I want to take to my supervisor. The milestone numbers stay in `WRITEUP.md`. The Gap page of `PathIntegral.exe` is this note, drawn.

## The idea, in one paragraph

Ordinary quantum mechanics gives me a wavefunction. Feynman gives me a sum over paths. In imaginary time the sum is a weighted average, and the weight of a path is `exp(-ε E)`, with `E` the energy along the path and `ε = T / N` the time step. A closed path that starts and ends at `x` is one sample of the density at `x`. Average enough of them and I get `|φ(x)|²`. That is the whole of chapter 14, and it is what the milestone checks.

## What I can read off a path

The histogram of where the path sits is the ground-state density. Two other things are sitting in the same path, and I had not been using them.

The slope of `⟨x(τ) x(0)⟩` is the gap, `E1 - E0`. For the oscillator that gap is 1, and the correlator follows it.

The energy itself can be read off the links. Each link contributes a kinetic piece and `V` at the midpoint, with a subtraction that removes the divergence of a short time step:

`E = 1/(2ε) - ⟨ m (Δx)² / (2 ℏ ε²) ⟩ + ⟨ V / ℏ ⟩`

On a coarse lattice this comes out a little under the Schrödinger energy. It is the lattice energy, and it moves toward the true one as `ε` shrinks. It is a check, not a replacement for the eigenvalue.

How often the path changes sign is the third reading. An oscillator path crosses zero many times. A path sat in one side of a double well hardly ever does.

## The result

I solve the Schrödinger equation twice, for the ground state and the first excited state, on the same grid the rest of the program uses.

For the oscillator, `E0 = 0.500` and `E1 = 1.500`, so the gap is 1. The long chain crosses zero 10.6 times per path. The energy read off its links is `0.43 ± 0.01`, a bit under `0.50`, which is the coarse time step and not a bug in the eigenvalue. Its correlator follows the gap.

For the double well `V = 0.8 (x² - 2.25)²`, `E0 = 1.731` and `E1 = 1.798`, so the gap is `0.067`. That gap is the tunnel splitting: the two wells talk to each other only by leaking through the barrier. The chain on the same potential crosses zero 0.09 times per path. Its correlator drops a little at the first time step, which is the fast motion inside one well, and then it sits flat. A log-slope fit of that flat line does not give `0.067`.

So the splitting is in the spectrum, and a local Metropolis move does not see it. The histogram still shows two peaks, because half the chains are started on each side. That is not the same as a path that tunnels.

The Gap page plots the splitting against `ℏ`. It is not a slow drift.

| ℏ | gap |
| --- | --- |
| 0.45 | 0.00005 |
| 0.60 | 0.0013 |
| 0.80 | 0.015 |
| 1.00 | 0.067 |
| 1.30 | 0.26 |
| 1.70 | 0.71 |
| 2.20 | 1.52 |

From `ℏ = 2.2` down to `ℏ = 0.45` the gap falls by about five orders of magnitude. That is the instanton effect the textbooks describe, measured here as two eigenvalues of the same Hamiltonian the paths are sampling. I have not assumed a WKB formula. The curve is the solver.

## What I want to do next

The missing object is one update that proposes a kink: flip a whole stretch of the path from one well to the other, accept or reject it with the same action, and only then fit the gap. If that fit comes back near `0.067`, the path integral and the spectrum agree on the thing local updates cannot reach. If it does not, the action or the fit is wrong and I will know which.

That is a small change to `metropolis.hpp` and a number I can defend. It is also the step that makes the method worth having, because a matrix already gives me `0.067` in a second. The path integral earns its place when the same splitting can be read off a path.
