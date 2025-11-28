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

1. Register for a GEE account at: https://earthengine.google.com/
2. Create a project at: https://console.cloud.google.com/
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
├── watershed_analysis.py
├── boundary.geojson
└── data/
    ├── berambadi_centroids.csv
    ├── berambadi_backscatter.csv
    ├── mesh_vertices.csv
    ├── mesh_triangles.csv
    └── backscatter_VH_timeseries.gif
```

## Files

| File | Description |
|------|-------------|
| `watershed_analysis.py` | Main script with all processing functions |
| `boundary.geojson` | Watershed boundary polygon |
| `data/berambadi_centroids.csv` | Mesh centroid locations |
| `data/berambadi_backscatter.csv` | Extracted backscatter values for all dates |
| `data/mesh_vertices.csv` | Mesh vertex coordinates |
| `data/mesh_triangles.csv` | Triangle connectivity |

## Configuration

Before running, update these paths in `watershed_analysis.py`:

```python
# ============================================
# EDIT THESE PATHS FOR YOUR SYSTEM
# ============================================
BOUNDARY_FILE = 'boundary.geojson'
OUTPUT_DIR = 'data/'
```

And update your GEE project ID:

```python
ee.Initialize(project='your-project-id')  # Replace with your project ID
```

## Usage

The script is modular - run each stage by uncommenting the relevant section:

### Stage 1: Generate Mesh

Creates triangular mesh from boundary and saves centroid locations.

```python
generateMesh(
    boundary_file=BOUNDARY_FILE,
    output_dir=OUTPUT_DIR,
    mesh_options='pqa0.00003'
)
```

### Stage 2: Authenticate Google Earth Engine

```python
import ee
ee.Authenticate()  # Run once, opens browser for authentication
ee.Initialize(project='your-project-id')
```

### Stage 3: Extract Backscatter Data

Queries Sentinel-1 backscatter for all available dates in the date range.

```python
queryGEE(
    centroids_csv=CENTROIDS_CSV,
    boundary_file=BOUNDARY_FILE,
    output_dir=OUTPUT_DIR,
    date_start='2022-08-01',
    date_end='2022-12-31'
)
```

### Stage 4 & 5: Generate Visualizations and GIF

Generates PNG for each date with fixed colorbar and combines them into an animated GIF.

## Output

- Individual backscatter maps for each Sentinel-1 acquisition date
- Animated GIF showing temporal variation

![VH Backscatter Time Series](data/backscatter_VH_timeseries.gif)

## Mesh Parameters

Adjust `mesh_options` to control triangle density:

| Option | Triangles | Use Case |
|--------|-----------|----------|
| `'pqa0.0001'` | ~100 | Quick testing |
| `'pqa0.00003'` | ~400 | Default, balanced |
| `'pqa0.00001'` | ~1000+ | High resolution |

## Using Your Own Watershed

1. Create a GeoJSON boundary file for your region using https://geojson.io/
2. Replace `boundary.geojson` with your file
3. Update `BOUNDARY_FILE` path in the script
4. Update GEE project ID
5. Run Stage 1 to generate mesh, then Stage 2-5

## Data Source

- **Sentinel-1 GRD** - C-band SAR imagery from ESA/Copernicus
- Accessed via Google Earth Engine: `COPERNICUS/S1_GRD`
- Temporal resolution: ~12 days
- Spatial resolution: 10m

## Author

Dr. Sathyanarayan Rao

## License

MIT

## Contributing

Contributions welcome! Future plans include:
- Soil moisture estimation using ML
- Additional satellite data sources (Sentinel-2 optical)
- NDVI analysis
- Groundwater correlation studies
