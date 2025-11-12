# scPortrait-to-IMC helper

Batch-load Imaging Mass Cytometry (IMC) datasets into [scPortrait](https://github.com/MannLabs/scPortrait) using a single, whole-cell segmentation mask per ROI.

## Layout

```
C:\GitHub\scPortrait_to_IMC
├── imc_to_single_cells.py  # helper script
├── config_imc.yml          # sample extraction config
└── README.md
```

Keep your data elsewhere, e.g.:

```
IMC_RUN/
├── images/
│   ├── ROI_001/
│   │   ├── CD45.tiff
│   │   ├── PanCK.tiff
│   │   └── DNA1.tiff
│   └── ROI_002/
│       ├── CD45.tiff
│       ├── PanCK.tiff
│       └── DNA1.tiff
└── masks/
    ├── ROI_001.tiff
    └── ROI_002.tiff
```

Each ROI has its own folder containing all channel TIFFs for that ROI. The mask directory contains one whole-cell mask per ROI (same filename stem as the ROI directory).

## Configuration

Edit `config_imc.yml` to suit your hardware (thread count, output image size, normalization). The config is passed directly to scPortrait's `Project`, so any extraction parameters supported by `HDF5CellExtraction` can be added here.

## Usage

1. Install scPortrait in your environment (`pip install scportrait tifffile`).
2. Run the helper, pointing it at your IMC folders and config:

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
