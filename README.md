# scPortrait-to-IMC helper

Batch-load Imaging Mass Cytometry (IMC) datasets into [scPortrait](https://github.com/MannLabs/scPortrait) using a single, whole-cell segmentation mask per ROI.

## Layout

```
C:\GitHub\scPortrait_to_IMC
├── imc_to_single_cells.py  # helper script
├── config_imc.yml          # sample extraction config
├── run_imc.slurm           # example SLURM batch script
└── README.md
```

Keep your data elsewhere, e.g.:

```
IMC_RUN/
├── images/
│   ├── ROI_001/
│   │   ├── 47_001_Pt198_PanCK.tiff
│   │   ├── 48_001_Ir191_DNA1.tiff
│   │   └── 12_001_Yb174_CD45.tiff
│   └── ROI_002/
│       ├── 47_002_Pt198_PanCK.tiff
│       ├── 48_002_Ir191_DNA1.tiff
│       └── 12_002_Yb174_CD45.tiff
└── masks/
    ├── ROI_001.tiff
    └── ROI_002.tiff
```

Each ROI has its own folder containing all channel TIFFs for that ROI. Filenames should follow `{chan#}_{roi#}_{str}_{channel name}.tiff` so the helper can infer the clean channel name from the final segment (e.g., `..._PanCK`). The mask directory contains one whole-cell mask per ROI (same filename stem as the ROI directory).

## Configuration

Edit `config_imc.yml` to suit your hardware (thread count, output image size, normalization). The config is passed directly to scPortrait's `Project`, so any extraction parameters supported by `HDF5CellExtraction` can be added here.

## Usage

1. Install scPortrait in your environment (`pip install scportrait tifffile`).
2. Use the helper either via CLI or by importing `process_imc_rois`.

### Command line

```powershell
python C:\GitHub\scPortrait_to_IMC\imc_to_single_cells.py `
  --channels-dir D:\IMC_RUN\images `
  --mask-dir D:\IMC_RUN\masks `
  --projects-root D:\IMC_projects `
  --config C:\GitHub\scPortrait_to_IMC\config_imc.yml `
  --mask-expand-px 2
```

Key behaviour:

- Every ROI with both a folder in `--channels-dir` and a mask in `--mask-dir` is processed automatically.
- Use `--roi ROI_001 --roi ROI_010` to limit the run to specific ROIs.
- Use `--overwrite` to force regeneration of existing scPortrait project folders.
- Use `--image-ext` to extend/override the default search extensions (`.tif .tiff .ome.tif .ome.tiff`).
- Use `--mask-expand-px N` to dilate each labelled cell mask by `N` pixels prior to extraction (helps capture thin membranes).

### From Python

```python
from pathlib import Path
from imc_to_single_cells import process_imc_rois

successes, failures = process_imc_rois(
    channels_dir=Path(r"D:\IMC_RUN\images"),
    mask_dir=Path(r"D:\IMC_RUN\masks"),
    projects_root=Path(r"D:\IMC_projects"),
    config_path=Path(r"C:\GitHub\scPortrait_to_IMC\config_imc.yml"),
    mask_expand_px=2,
)

print("Created:", successes)
print("Failed:", failures)
```

`successes` maps each ROI name to the generated `single_cells.h5sc` path, while `failures` (if any) maps ROI names to the raised exception for further inspection.

### SLURM batch example

Submit `run_imc.slurm` to execute the helper on an HPC cluster:

```bash
sbatch C:\GitHub\scPortrait_to_IMC\run_imc.slurm
```

Adjust the SBATCH resources, module/conda commands, and file paths inside the script to match your cluster environment. Logs land in `run_imc.slurm`'s `logs/` folder (created automatically by SLURM).

For each ROI `<ID>`, the script creates `--projects-root/<ID>/.../single_cells.h5sc`. Inspect the result with:

```python
from scportrait.io.h5sc import read_h5sc
adata = read_h5sc(r"D:\IMC_projects\ROI_001\extraction\data\single_cells.h5sc")
print(adata.obsm["single_cell_images"].shape)
```

## Verification checklist

- Confirm that `D:\IMC_projects\<ROI>\extraction\data\single_cells.h5sc` exists for every ROI.
- Load the file via `read_h5sc` to ensure the AnnData object contains the expected number of cells/channels.
- Re-run with `--debug` to view detailed logs if a ROI fails.
