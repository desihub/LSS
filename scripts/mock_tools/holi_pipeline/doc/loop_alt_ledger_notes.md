# Notes on `loop_alt_ledger()` and the TileTracker file

Summary of a chat session (2026-10-29/30) about `LSS.SV3.altmtltools.loop_alt_ledger()`
and the `TileTracker.ecsv` file it consumes.

## Call tree of `loop_alt_ledger()`

**Callers:**
- [bin/SingleNode-runAltMTLParallel.py](../../../bin/SingleNode-runAltMTLParallel.py#L112)
- [bin/UpdateAltMTLParallel.py](../../../bin/UpdateAltMTLParallel.py#L36)
- [bin/runAltMTLParallel.py](../../../bin/runAltMTLParallel.py#L121)
- [bin/runAltMTLRealizations.py](../../../bin/runAltMTLRealizations.py#L306)
- [scripts/mock_tools/holi_pipeline/runAltMTL.py](runAltMTL.py#L212)

**Callees** (directly invoked from the function body in
[py/LSS/SV3/altmtltools.py](../../../py/LSS/SV3/altmtltools.py#L2375)):
- `do_fiberassignment(...)`
- `make_fibermaps(...)`
- `update_alt_ledger(...)`
- `reprocess_alt_ledger(...)`
- `update_lya_1b(...)` (from `desitarget.mtl`)
- `add_new_ledgers(...)`
- `process_vetoes_altmtl(...)`
- `write_amtl_tile_tracker(...)` (called inside the branches above)

## Loading the observation plan

The "observation plan" (ordered list of tile pointings with a date/time) is loaded here:

```python
altMTLTileTrackerFN = makeTileTrackerFN(altmtldir, survey=survey, obscon=obscon)
altMTLTileTracker = Table.read(altMTLTileTrackerFN)

actionList = altMTLTileTracker[np.invert(altMTLTileTracker["DONEFLAG"])]
actionList.sort(["ACTIONTIME"])
```
([py/LSS/SV3/altmtltools.py](../../../py/LSS/SV3/altmtltools.py#L2551-L2558))


The TileTracker file (`{survey}survey-{obscon}obscon-TileTracker.ecsv`) is a table with
columns `TILEID`, `ACTIONTYPE`, `ACTIONTIME`, `DONEFLAG`, `ARCHIVEDATE`. `loop_alt_ledger`
filters out already-done actions (`DONEFLAG` False) and processes the remainder in
chronological order (`for action in actionList:`).

## Resuming after an interruption

The resume mechanism relies entirely on the persistent `DONEFLAG` column of the
TileTracker file, combined with writing it back to disk immediately after each
processed action:

1. At every call, `loop_alt_ledger()` re-reads the TileTracker file from disk and keeps
   only the actions with `DONEFLAG == False`.
2. Each action is marked done (via `write_amtl_tile_tracker(...)`) as soon as it is
   processed, so a restarted run automatically skips everything already completed and
   resumes with the first remaining action, sorted by `ACTIONTIME`.

## TileTracker file location

Created by `makeTileTracker()` / `makeTileTrackerFN()`
([py/LSS/SV3/altmtltools.py](../../../py/LSS/SV3/altmtltools.py#L663-L665)) with the
naming template:

```
{dirName}/{survey}survey-{OBSCON}obscon-TileTracker.ecsv
```

For the *holi* pipeline (`mock_dir = /pscratch/sd/j/jcolley/test_fits3`), the real path is:

```
{mock_dir}/altmtl{SEED}/Univ000/{survey}survey-{OBSCON}obscon-TileTracker.ecsv
```
It's the same for all simulation.

## `ACTIONTYPE` values

No dedicated formal documentation exists; the following table was reconstructed from the
code in [altmtltools.py](../../../py/LSS/SV3/altmtltools.py):

| ACTIONTYPE | Triggered by | Description |
|---|---|---|
| `fa` | `makeTileTracker()` reads the `MTLTIME`/`RUNDATE` header of `fiberassign-{tileid}.fits` ([L829-L836](../../../py/LSS/SV3/altmtltools.py#L829-L836)) | Fiber assignment for a tile. Calls `do_fiberassignment()` + `make_fibermaps()` ([L2581-L2614](../../../py/LSS/SV3/altmtltools.py#L2581-L2614)). |
| `update` | First entry in `mtl-done-tiles.ecsv` for a given tile ([L892-L902](../../../py/LSS/SV3/altmtltools.py#L892-L902)) | Updates the MTL ledger with redshifts (`zcat`) from the tile's first observation. Calls `update_alt_ledger()` ([L2615-L2632](../../../py/LSS/SV3/altmtltools.py#L2615-L2632)). |
| `reproc` | Subsequent entries in `mtl-done-tiles.ecsv` for the same tile ([L899-L902](../../../py/LSS/SV3/altmtltools.py#L899-L902)) | Reprocesses an already-observed tile with new redshift data. Calls `reprocess_alt_ledger()` ([L2633-L2642](../../../py/LSS/SV3/altmtltools.py#L2633-L2642)). |
| `lya1b` | Special action auto-added if obscon=dark and an action exists after 2025-07-21 ([L973-L985](../../../py/LSS/SV3/altmtltools.py#L973-L985)) | Artificially increments `NUMOBS` for Lyman-alpha QSOs to reproduce the real survey decision (DESI PR #845). Calls `update_lya_1b()` ([L2644-L2658](../../../py/LSS/SV3/altmtltools.py#L2644-L2658)). |
| `addnew` | Dates listed in the `{OBSCON}-ledgers.yaml` file ([L987-L997](../../../py/LSS/SV3/altmtltools.py#L987-L997)) | Creates new MTL ledgers for HEALPixels added after the fact. Calls `add_new_ledgers()` ([L2659-L2677](../../../py/LSS/SV3/altmtltools.py#L2659-L2677)). |
| `veto` | Entries in `mtl-done-vetoes.ecsv` within the date window (notably BRIGHT1B obscon) ([L949-L963](../../../py/LSS/SV3/altmtltools.py#L949-L963)) | Vetoes specific targets in the ledgers after the fact. Calls `process_vetoes_altmtl()` ([L2679-L2688](../../../py/LSS/SV3/altmtltools.py#L2679-L2688)). |

The exhaustive list of valid values is enforced by the final `raise ValueError` in the
`if/elif` chain inside `loop_alt_ledger()`
([L2689-L2691](../../../py/LSS/SV3/altmtltools.py#L2689-L2691)).

This documentation was also added as an English "Notes" section in the docstring of
`loop_alt_ledger()` itself
([py/LSS/SV3/altmtltools.py](../../../py/LSS/SV3/altmtltools.py#L2468-L2493)).

## Ordering logic of `fa` / `update` / `reproc` in the TileTracker file

The file follows a double logic:

1. **Global chronological sort**: all rows, across all `TILEID`s, are ordered by
   increasing `ACTIONTIME` (this is explicitly done in `loop_alt_ledger()` via
   `actionList.sort(["ACTIONTIME"])`).
2. **Per-tile causal constraint**: for a given tile, the order is necessarily
   `fa` -> `update` -> `reproc` (`reproc` only if there are multiple
   observations/redshifts), because each action functionally depends on the previous one:
   - `fa` (fiber assignment) must exist before the tile can be observed.
   - `update` corresponds to the **first** redshift entry in `mtl-done-tiles.ecsv`
     for that tile.
   - `reproc` corresponds to **subsequent** entries (reprocessing) for the same tile,
     so it can only happen after an `update`.

Example verified on tile 3414:
```
3414 fa      2021-09-23T21:14:24+00:00
3414 update  2021-09-27T19:43:47+00:00
3414 reproc  2021-12-14T18:42:25+00:00
```

