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
│   ├── CD45/
│   │   ├── ROI_001.tiff
│   │   └── ROI_002.tiff
│   └── PanCK/
│       ├── ROI_001.tiff
│       └── ROI_002.tiff
└── masks/
    ├── ROI_001.tiff
    └── ROI_002.tiff
```

Each channel lives in its own subfolder and uses consistent ROI filenames. The mask directory contains a matching TIFF per ROI (single whole-cell segmentation).

## Configuration

Edit `config_imc.yml` to suit your hardware (thread count, output image size, normalization). The config is passed directly to scPortrait’s `Project`, so additional extraction parameters can be added under the `HDF5CellExtraction` block.

## Usage

1. Install scPortrait in your environment (`pip install scportrait tifffile`).
2. Run the helper, pointing it at your IMC folders and config:

```powershell
python C:\GitHub\scPortrait_to_IMC\imc_to_single_cells.py `
  --channels-dir D:\IMC_RUN\images `
  --mask-dir D:\IMC_RUN\masks `
  --projects-root D:\IMC_projects `
  --config C:\GitHub\scPortrait_to_IMC\config_imc.yml
```

Key behaviour:

- All ROIs present in both `--mask-dir` and the channel folders are processed automatically.
- Use `--roi ROI_001 --roi ROI_010` to limit the run to specific ROIs.
- Use `--overwrite` to force regeneration of existing scPortrait project folders.
- Use `--image-ext` to extend/override the default search extensions (`.tif .tiff .ome.tif .ome.tiff`).

For each ROI `<ID>`, the script creates `--projects-root/<ID>/.../single_cells.h5sc`. Inspect the result with:

```python
from scportrait.io.h5sc import read_h5sc
adata = read_h5sc(r"D:\IMC_projects\ROI_001\extraction\data\single_cells.h5sc")
print(adata.obsm["single_cell_images"].shape)
```

## Verification checklist

- Confirm that `D:\IMC_projects\<ROI>\extraction\data\single_cells.h5sc` exists for every ROI.
- Run `read_h5sc` as shown above to ensure the AnnData object loads and contains the expected number of cells/channels.
- Optionally re-run with `--debug` to view scPortrait’s log output if a ROI fails.
