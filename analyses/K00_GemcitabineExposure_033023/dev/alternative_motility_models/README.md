# Alternative Motility Model Interim Notes

This directory contains a prototype comparison of simple generative motility
models for accepted yellow reconstructed tracks in the low-dose gemcitabine
range. The current report is:

- `alternative_motility_model_report.html`
- `alternative_motility_model_report.Rmd`
- `artifacts/alternative_motility_model_fits.rds`

Run from the repository root:

```bash
scripts/agentRrunner.sh analyses/K00_GemcitabineExposure_033023/dev/alternative_motility_models/scripts/fit_alternative_motility_models.R --force=true
scripts/agentRrunner.sh analyses/K00_GemcitabineExposure_033023/dev/alternative_motility_models/scripts/render_report.R
```

The current analysis uses all selected low-dose conditions from 0 through
25 nM gemcitabine: 0, 3.125, 6.25, 12.5, and 25 nM. Models are fit on training
sites when a site holdout is available and evaluated on the held-out track
skeletons.

## Interim Conclusions After Pixel Rounding

The most important update is that all simulated model coordinates are now
rounded to nearest-pixel observations before predictive checks are computed.
Latent movement remains continuous inside the simulators, but the comparison
layer treats every simulated track as an integer-pixel observed trajectory,
matching the observed coordinate grid.

This rounding step substantially improves the explanation of the lag-1 step
cosine autocorrelation distribution. In particular, the observed orthogonal
lag-1 pairs are strongly concentrated among very small pixel steps, mostly near
0 or 1 pixel. Once model observations are rounded to pixels, the sharp cosine
features near -1, 0, and 1 look much less like evidence for a special biological
turning process and much more like a consequence of pixel-grid observation,
short steps, centroid jitter, and integer displacement geometry.

The one-frame speed distribution remains reasonably well described after
rounding. The best reduced models still capture the broad one-frame step scale,
so the short-lag speed mismatch is no longer the dominant concern.

The main remaining failure is the MSD curve. None of the retained alternatives
describes MSD well beyond roughly lag 2 or 3 frames. This suggests that the
models can reproduce short-step marginal behavior and much of the lag-1
autocorrelation shape, but they still miss longer-lag spatial structure. The
missing component may involve track-level heterogeneity, time-varying motility,
local confinement over longer windows, unmodeled drift/crowding, or additional
observation and linking structure not captured by the current simple generators.

The comparable OU refit is especially concerning as a workflow reference. In
the current rounded low-dose comparison, `OU refit comparable` is last by the
aggregate composite rank score. It is also weak on the one-frame speed and
lagged-displacement distribution checks, while not resolving the longer-lag MSD
problem. This strongly argues for re-evaluating the OU model's placement in the
main workflow: it should not be treated as the default mechanistic explanation
for these yellow-track motility summaries without further revision or a clear
justification for the specific quantities it is intended to model.

For now, the rounded-observation results support a narrower interpretation:
pixel-level observation effects plausibly explain much of the lag-1
autocorrelation distribution, but the longer-lag MSD discrepancy remains open
and should drive the next model-development step.
