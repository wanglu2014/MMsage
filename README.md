# compar_plot_0513

Minimal reproducible folder for 8-metrics supplementary PDFs.

## Structure

- `code/compare_methods_by_x_0513.py`
- `code/requirements.txt`
- `output/figures/` (target PDFs)
- `output/input_cache/` (auto-copied soil source file for traceability)

## Data Source Mapping

Paths are relative to this folder and resolved by the script from
`DATA_ROOT = ROOT_DIR.parent / "data"`, i.e. the package-level `data/`:

- Soil source (used directly):  
  `../data/soil_0512/soilori-ROC_0412result_replaced.csv`
- Simulation release source (used directly):  
  `../data/sim_0512/Sim release ROC_0206resmodi.csv`
- Simulation consumption source (used directly):  
  `../data/sim_0512/Sim cons_ROC_0303_resneg.csv`

`supple_compar_plot_0513/data` is not used by the script; it is a legacy copy kept
for fidelity with the 2026-05-13 package.

## Run

From the package root:

```bash
python supple_compar_plot_0513/code/compare_methods_by_x_0513.py
```

## Expected Target Outputs

- `output/figures/fig_soil_release_8metrics_nejm_0504.pdf`
- `output/figures/fig_simulation_release_8metrics_nejm_0504.pdf`
- `output/figures/fig_simulation_consumption_8metrics_nejm_0504.pdf`
