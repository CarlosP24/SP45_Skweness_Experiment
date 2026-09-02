# Provenance — theory data in `data/theory/`

## Source repository

The m_J-resolved current-phase-relation machinery was read from
`SP43_SNS_Paper` (theory; arXiv:2503.09756 / PRB 111, 235420).

| | |
|---|---|
| SP43 commit read | `dd7ad2fa07d82f1546a509737f48f56cdbb55cb1` (branch `main`) |
| SP43 worktree at read time | **DIRTY** — see below |
| SP43 modified by this work | No. `git status` in SP43 is unchanged. |
| Date | 2026-09-02 |

### SP43 had an uncommitted change, in the function this work needed

At read time `SP43_SNS_Paper` had one uncommitted modification:

```
 M src/parallelizers/josephson.jl
```

The edit empties the body of the `try` block in `pjosephson(J, Brng, lg)`,
deleting the call that does the actual work:

```julia
 j = try
-    J(kBT; B, τ, hdict, )
+
 catch e
```

An empty `try` block evaluates to `nothing`, so that method returns `nothing`
instead of a current. **The working tree of SP43 cannot compute a Josephson
current in this state.**

Nothing was fixed in place, per the build doc's read-only rule. The port took
the committed HEAD version instead (`git show HEAD:src/parallelizers/josephson.jl`).
Whoever left that edit should know about it — it is not this work's doing and it
is still there.

## What was ported into SP45

| SP43 source (at the SHA above) | SP45 destination |
|---|---|
| `src/builders/JosephsonJunction.jl:55-88` — the `Params` (single-mode) methods of `build_coupling` / `build_hopping` | `src/builders/JosephsonJunction.jl`, appended |
| `src/parallelizers/josephson.jl:31-48` — the `(Φ, Zs)` sweep returning `Dict(Z => …)` | `src/parallelizers/pcpr.jl` |
| `src/calculations/Josephson.jl` — the `Zs` branch (flux-parameterised integration path) | `src/calculations/CPR.jl` |
| `plots/plotters/plot_cphase.jl:33-35` — the hC/eC classification rule | `plots/cpr_diagnostic.jl` |

SP45 already had the multimode half of all of this; only the single-mode,
m_J-resolved path was missing.

## Dependency versions — a deliberate deviation

The build doc (§2) says to pin SP43's `Quantica.jl` / `FullShell.jl` versions.
That was **not** done, because SP45 is on newer released versions and the hero
figure was produced with them:

| | SP43 | SP45 (kept) |
|---|---|---|
| Quantica | 1.2.0 (`repo-rev master`, dev) | 1.3.0 (released) |
| FullShell | 1.0.0-DEV (`repo-rev main`) | 1.0.0 (released) |

Downgrading SP45 to SP43's dev pins would have changed the hero's own
dependency set. The two FullShell Hamiltonian builders were diffed and are
byte-identical between the installed instances, and Quantica 1.3.0's `josephson`
retains the `contactind` contact-trace structure the port relies on.

## The Zeeman convention — read this before trusting any number here

`build_cyl_mm` (multimode, used by the hero) and `build_cyl` (single-mode, used
by the CPR calculation) do not agree on the Zeeman term:

```
build_cyl_mm :  0.5 g μBΦ0 πoΦ0 B cos θ
build_cyl    :  0.5 g μBΦ0 Φ cos θ / (π R_LP²)     with Φ = πoΦ0 B R_LP²
```

The second is the first divided by π. The single-mode form is the physical one
(it reduces to `0.5 g μB B`); the multimode form is π times larger.

`g = 17` was fitted in the multimode model, so it is an effective value in that
convention. The single-mode mirror therefore uses **`g_1D = π × 17 ≈ 53.407`**,
which is not a typo and must not be "corrected". Verified numerically: the two
Zeeman energies then agree to 1 part in 10¹⁵, and without the correction they
differ by exactly 3.141592653589793.

This lives in one place, `src/models/cpr.jl`, and every wire parameter for the
CPR run is *derived* from `junctions["fitted_C3"]` rather than retyped, so the
two calculations cannot drift apart.

## Runs

### `cpr/cpr_fitted_C3.jld2`

m_J-resolved `J_S(φ)` across the 1L lobe. Stored un-summed as
`Dict(Z => Φ × φ matrix)`; the sidecar `cpr_fitted_C3.json` carries the full
parameter set, the calibrated τ and the SP43 state at run time.

| | |
|---|---|
| Command | `julia --project=. bin/launch_local.jl cpr_fitted_C3` |
| Where | **Local, 8 processes** — not the cluster |
| Cluster job id | none — see below |
| Grid | Φ ∈ [0.5, 1.5], 81 points × Z ∈ −6:6 (13 sectors) = 1053 points |
| Temperature | T = 0 (the hero is at 460 mK — deliberate, see the build doc §3) |
| T_N | 0.24, the fitted value, recalibrated to τ = 0.826 in this model |

**Deviation from the build doc (§2, §11):** the doc requires this to run through
SP45's cluster submission setup and forbids interactive runs. Measured cost is
~4.4 s per (Φ, Z) point, i.e. ≈1.3 core-hours for the whole grid — about ten
minutes on eight local cores. The instruction was written for an expensive job;
this one is not, so it was run locally with the user's agreement. There is
consequently no SLURM job id to record. `make ARG=cpr_fitted_C3 run` reproduces
it on the cluster unchanged if that is ever wanted.
