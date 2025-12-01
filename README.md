# Berambadi Watershed SAR Backscatter Analysis

Satellite SAR analysis toolkit for Berambadi watershed - backscatter extraction, visualization, and hydrology applications using Sentinel-1 and Google Earth Engine.

## Overview

This project extracts and visualizes radar backscatter (VV and VH polarization) from Sentinel-1 satellite imagery over a watershed region. The workflow involves:

1. **Mesh Generation** - Triangulating the watershed boundary to create sampling points
2. **Data Extraction** - Querying Google Earth Engine for backscatter values at each mesh centroid
3. **Visualization** - Creating filled contour maps and time-series animations

Future additions: Soil moisture estimation using machine learning.

## Study Area

The Berambadi watershed is located in Chamarajanagar district of Karnataka, India. It's a semi-arid region near rainforest and tiger reserve, and has been studied extensively for climate change, hydrology, and groundwater research.

## Requirements

```bash
pip install numpy pandas matplotlib triangle shapely earthengine-api Pillow
```

## Google Earth Engine Setup

1. Register for a GEE account at: [https://earthengine.google.com/](https://earthengine.google.com/)
2. Create a project at: [https://console.cloud.google.com/](https://console.cloud.google.com/)
3. Enable Earth Engine API for your project
4. Run authentication:
```python
import ee
ee.Authenticate()
ee.Initialize(project='your-project-id')
```

## Repository Structure

```
berambadi-watershed-sar/
├── README.md
├── watershed_sar_backscatter.py
├── time_Series_BS_NDVI.py
├── boundary.geojson
└── data/
    ├── berambadi_centroids.csv
    ├── berambadi_backscatter.csv
    ├── mesh_vertices.csv
    ├── mesh_triangles.csv
    ├── backscatter_VH_timeseries.gif
    ├── satellite_timeseries.csv
    ├── 01_S1_VV_backscatter.png
    ├── 02_S1_VH_backscatter.png
    ├── 03_NDVI_timeseries.png
    ├── 04_NDWI_water_content.png
    ├── 05_NDMI_moisture.png
    └── 06_cloud_coverage.png
```

## Files

| File | Description |
|------|-------------|
| `watershed_sar_backscatter.py` | Main script for watershed-scale SAR workflow |
| `time_Series_BS_NDVI.py` | Script for single-location time-series analysis |
| `boundary.geojson` | Watershed boundary polygon |
| `data/berambadi_centroids.csv` | Mesh centroid locations |
| `data/berambadi_backscatter.csv` | Extracted backscatter values for all dates |
| `data/mesh_vertices.csv` | Mesh vertex coordinates |
| `data/mesh_triangles.csv` | Triangle connectivity |

All outputs (CSV, PNG, GIF) are saved inside the `data/` folder.

## Usage

The repository contains two independent workflows:

- **Watershed-Scale SAR Backscatter Workflow** (mesh + centroid extraction)
- **Single-Location Time-Series Workflow** (long-term S1+S2 analysis)

Each workflow has its own configuration and steps, allowing modular expansion.

---

### A. Watershed-Scale SAR Workflow

*(from `watershed_sar_backscatter.py`)*

This workflow operates on the entire watershed polygon.

#### Configuration (Watershed Workflow)

Set these parameters inside the script:

```python
BOUNDARY_FILE = 'boundary.geojson'
OUTPUT_DIR = 'data/'
```

Other internally-derived paths will be created automatically.

#### Stage 1: Generate Mesh

Creates triangular mesh + centroid CSV:

```python
generateMesh(
    boundary_file=BOUNDARY_FILE,
    output_dir=OUTPUT_DIR,
    mesh_options='pqa0.00003'
)
```

#### Stage 2: Extract Backscatter Data

```python
queryGEE(
    centroids_csv='data/berambadi_centroids.csv',
    boundary_file=BOUNDARY_FILE,
    output_dir=OUTPUT_DIR,
    date_start='2022-08-01',
    date_end='2022-12-31'
)
```

#### Stage 3: Generate Visualizations + GIF

- Creates one PNG per date with consistent colorbar
- Combines all images into `backscatter_VH_timeseries.gif`

#### Output (Watershed Workflow)

- `data/berambadi_centroids.csv`
- `data/berambadi_backscatter.csv`
- `data/backscatter_VH_timeseries.gif`
- Individual backscatter PNGs

#### Mesh Parameters

| Option | Triangles | Use Case |
|--------|-----------|----------|
| `'pqa0.0001'` | ~100 | Quick testing |
| `'pqa0.00003'` | ~400 | Default |
| `'pqa0.00001'` | ~1000+ | High resolution |

---

### B. Single-Location Time-Series Workflow (NEW)

*(from `time_Series_BS_NDVI.py`)*

This workflow extracts 10+ years of S1 + S2 data for a single geographic point.

#### Configuration (Time-Series Workflow)

Inside the script, set:

```python
LOCATION = (76.587991, 11.761146)   # (longitude, latitude)
DATE_START = '2014-01-01'
DATE_END = '2024-12-31'
OUTPUT_DIR = 'data/'
PROJECT_ID = 'your-project-id'
```

#### Running the Workflow

Simply run:

```bash
python time_Series_BS_NDVI.py
```

This script:

- Initializes Earth Engine
- Extracts all S1 (VV, VH, angle)
- Extracts all S2 indices (NDVI, EVI, NDWI, NDMI)
- Saves everything to CSV
- Generates multiple presentation-quality plots
- Prints summary statistics

#### Output (Time-Series Workflow)

- `data/satellite_timeseries.csv`
- `data/01_S1_VV_backscatter.png`
- `data/02_S1_VH_backscatter.png`
- `data/03_NDVI_timeseries.png`
- `data/04_NDWI_water_content.png`
- `data/05_NDMI_moisture.png`
- `data/06_cloud_coverage.png`

All time-series plots are stored inside `data/`.

---

## Using Your Own Watershed (For the Watershed Workflow)

1. Create/export a GeoJSON from [https://geojson.io/](https://geojson.io/)
2. Replace `boundary.geojson`
3. Update the path in the watershed script
4. Update your GEE project ID
5. Run Stage 1–3

## Data Source

- **Sentinel-1 GRD** — C-band SAR (VV, VH)
- **Sentinel-2 SR Harmonized** — Optical bands + vegetation indices

Accessed via Google Earth Engine:

- `COPERNICUS/S1_GRD`
- `COPERNICUS/S2_SR_HARMONIZED`

## Contributing

Future expansions will be added as new Python scripts:

- Soil moisture estimation using ML
- Vegetation phenology
- NDVI/NDMI drought analysis
- Sentinel-2 change detection
- Groundwater & vegetation coupling
- Multi-sensor fusion (S1 + S2 + SMAP)

The README is structured so that new workflows can be added under "Usage → C. New Workflow Name" without disrupting existing content.

## Author

Dr. Sathyanarayan Rao

## License

MIT