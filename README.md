# Differential Abundance with sccomp

Demo for the Compute and Conquer blog: testing cell-type composition changes in *Tal1*⁻/⁻ mouse chimeric embryos with [sccomp](https://bioconductor.org/packages/sccomp/).

## Data

`MouseGastrulationData::Tal1ChimeraData()` (Pijuan-Sala *et al.*, 2019). E8.5 chimeric embryos in which tdTomato⁺ *Tal1*⁻/⁻ ESCs were injected into wild-type blastocysts. Four samples, all from one pool:

| Sample | Cells |
|---|---|
| 1, 2 | tdTomato⁺ (*Tal1*⁻/⁻) |
| 3, 4 | tdTomato⁻ (wild-type host) |

Cells labelled `Doublet` or `Stripped` are removed before analysis.

## Analysis

`code/01-Differential-Abundance.R`:

1. Loads the data, log-normalises, converts to Seurat.
2. Fits `sccomp_estimate(~ tomato)` on `celltype.mapped` counts per sample, using:
   - `inference_method = "hmc"` (full MCMC instead of the default Pathfinder approximation)
   - `prior_mean$coefficients = c(0, 3)` (wider than the default `c(0, 1)`, which shrinks large effects too strongly with only 2 replicates per group)
3. Tests effects with `sccomp_test()` (default threshold ±0.1 logit).
4. Plots the `tomatoTRUE` effect per cell type with 95% credible intervals. Red = `c_FDR < 0.05`.

## Outputs

| File | Contents |
|---|---|
| `outputs/figures/sccomp_output.png` | Composition effect plot |

## Requirements

R with Seurat, scuttle, dplyr, ggplot2, sccomp (≥ 1.10.0) and MouseGastrulationData. The script expects sccomp in `~/Rlib-sccomp` and redirects its Stan model cache to `~/.sccomp_models/`. Edit `wd` at the top of the script to match your path.

## Notes

- Positive effects mean a larger share among tdTomato⁺ (*Tal1*⁻/⁻) cells.
- Effects are **relative**. Blood is about 33% of wild-type cells and absent from knockout cells, so most other cell types show a small positive effect (≈ +0.4 logit) from that redistribution alone.
- With n = 2 per group, results are sensitive to the prior. Report the settings used.

## References

- Pijuan-Sala B, *et al.* (2019). A single-cell molecular map of mouse gastrulation and early organogenesis. *Nature* 566, 490–495. [doi:10.1038/s41586-019-0933-9](https://doi.org/10.1038/s41586-019-0933-9)
- Mangiola S, *et al.* (2023). sccomp: Robust differential composition and variability analysis for single-cell data. *PNAS* 120(33), e2203828120. [doi:10.1073/pnas.2203828120](https://doi.org/10.1073/pnas.2203828120)
