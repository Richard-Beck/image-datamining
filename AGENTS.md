# Repository Notes for Codex

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
