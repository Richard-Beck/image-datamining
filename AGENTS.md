# Repository Notes for Codex

## General development guidelines
The user really hates waiting an extended period of time for the LLM to run endless "smoke tests", 
inspections of git diffs, and other lengthy tasks. In general, unless explicitly asked otherwise,
the user wants a more-or-less working prototype in their hands as quickly as possible. 


## R Execution

Use `scripts/agentRrunner.sh` as the default wrapper for all R execution by the
LLM agent. This includes `Rscript`, `R CMD`, `rmarkdown::render()`, and any other
R commands run for inspection, analysis, plotting, or report rendering.

## Slurm QOS

Current QOS settings from
`sacctmgr show qos format=Name,Priority,MaxTRESPU,MaxJobsPU,MaxSubmitJobsPU,GrpTRES -P`:

```text
Name|Priority|MaxTRESPU|MaxJobsPU|MaxSubmitPU|GrpTRES
normal|0|cpu=64|1000001||
small|0|cpu=475|1000001||
large|0|cpu=1425|1000001||
xlarge|0|cpu=1800|1000001||
partsmall|0||1000001||cpu=250
medium|0|cpu=950|1000001||
xxlarge|0|cpu=3000|1000001||
```

For future Slurm submissions, choose the QOS that is likely to get the fastest
deployment for the requested job shape instead of defaulting to `normal`. For
large arrays of small single-CPU jobs, `small` or `large` may be more appropriate
than `normal` depending on the desired concurrency and current limits.

Do not add Slurm array throttles such as `--array=1-100%10` by default. A `%`
limit should only be used when there is a clear resource, filesystem, scheduler,
or user-requested reason to cap concurrency. When a higher-concurrency QOS is
chosen for many small jobs, leaving the array unthrottled is usually the intended
behavior.

## CellposeSAM / CPSAM

This repo has access to GPU execution, but Codex does not currently have permission
to submit or run GPU jobs. Only the user should run and test CellposeSAM/CPSAM jobs
for now.

The conda cpose environment is a useful default for running python scripts.

When adding or editing CellposeSAM scripts:

- Make GPU execution the default behavior.
- Provide an explicit CPU override only when useful for debugging or portability.
- Do not test CellposeSAM/CPSAM execution from Codex. Prepare scripts, commands,
  configs, and expected outputs, then leave execution and validation to the user.
- It is fine for Codex to inspect code, prepare inputs, and make non-CPSAM plots or
  summaries that do not invoke CellposeSAM.
