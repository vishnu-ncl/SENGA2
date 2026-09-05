# Archer2 8-node 512³ execute/tiling sweep

Compares tiling (`OPS_TILING_MAXDEPTH`, `OPS_CACHE_SIZE`, large tiles) against
the per-species `ops_execute` (`SENGA_OPS_EXEC=+ispec` / `all`) on 8 nodes
(1024 ranks), 512³, 10 steps.

Grid size is a **compile-time** parameter. `cont.dat` alone cannot change it.

## 1. Build on Archer2

```bash
source /work/e647/e647/ireguly5/source_archer2_gcc
cd /work/e647/e647/ireguly5/SENGA2_OPS/src
make senga2_mpi_tiled SENGA_NX=512 -j8
```

Default `SENGA_NX` is 256. The job aborts after the first run if the binary
still prints `nxglbl=256`.

## 2. Submit

```bash
cd /work/e647/e647/ireguly5/SENGA2_OPS
mkdir -p results
sbatch scripts/submit_archer2_8node_512.sh
```

The script copies `input/cont_512_10.dat` over `input/cont.dat` for the job
and restores the original on exit. Last-step HDF5 dumps stay off.

## 3. Configs (one `srun` each, same allocation)

| name | tiling | `SENGA_OPS_EXEC` |
|---|---|---|
| `base` | `-OPS_DIAGS=2` only | default (`ispec` off) |
| `largetile` | `MAXDEPTH=6 TILESIZE=1000` | default |
| `cache1` | `MAXDEPTH=6 CACHE_SIZE=1` | default |
| `cache2` | `MAXDEPTH=6 CACHE_SIZE=2` | default |
| `cache4` | `MAXDEPTH=6 CACHE_SIZE=4` | default |
| `cache2_ispec` | `MAXDEPTH=6 CACHE_SIZE=2` | `+ispec` |
| `cache2_all` | `MAXDEPTH=6 CACHE_SIZE=2` | `all` |

Default exec (unset `SENGA_OPS_EXEC`) is the comms-optimised set: keepers on,
**`ispec` off**. `+ispec` restores the per-species flush (more MPI, better
cache tiles). `all` turns every manual execute site on.

## 4. Bring results back

```bash
scp <login>:/work/e647/e647/ireguly5/SENGA2_OPS/results/archer2_8node_512_results.tgz .
tar xzf archer2_8node_512_results.tgz
python3 scripts/collate_senga_logs.py results
```

The tarball already contains `summary.csv` and `summary.md` if the job's
collate step ran. Re-running the collator locally is safe.

Columns that matter: `kernel`, `tiled_halo`, `kernel_plus_halo`, `msgs`, `mb`,
`total`. Do not compare `rhscal_ttime` across `ispec` on/off — the fused
species plan is flushed after that timer stops.
