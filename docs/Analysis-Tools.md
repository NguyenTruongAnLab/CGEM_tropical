# Analysis Tools

The repository ships Python helpers in `tools/` for three common post-run tasks:

1. **Check the forcing data** before or after a run.
2. **Plot model outputs** without any observation file.
3. **Overlay validation data** when observation tables are available.

These tools are output-format agnostic: they read the default CSV outputs and
the faster CGEM binary (`.bin`) outputs through `tools/cgem_read.py`.

## Install dependencies

```bash
pip install -r tools/requirements.txt
```

Run all commands from the repository root.

## Before the model run

```bash
python tools/check_inputs.py
```

This writes:

- `OUT/figures/fig_input_diagnostics.png`
- `OUT/tables/check_inputs_report.csv`

## After the model run

### Model-only figures

```bash
python tools/plot_profiles.py --output-dir OUT
python tools/plot_seasonal.py --output-dir OUT
python tools/plot_hydrodynamics.py --output-dir OUT
python tools/plot_reactions.py --output-dir OUT
```

### With bundled validation observations

```bash
python tools/plot_profiles.py --output-dir OUT --with-obs
python tools/plot_seasonal.py --output-dir OUT --with-obs
python tools/plot_hydrodynamics.py --output-dir OUT --with-obs
```

The bundled observation tables live in `INPUT/Validation/` and are already in
model units.

## Converting new field data

If your new observation tables are still in field units, create converted copies
with:

```bash
python tools/convert_obs_to_model_units.py
```

Converted files are written to `INPUT/Validation/converted/` so the source files
remain untouched.

## Where to read more

- `tools/README.md` — command-line details for every analysis script
- `INPUT/Validation/README.md` — expected observation schema and units
- `Documentation/Getting-Started.md` — broader run workflow for the model itself