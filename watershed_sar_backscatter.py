"""
Berambadi Watershed Processing
==============================
Functions: generateMesh(), queryGEE(), plotBackscatter()
Run stage by stage - comment/uncomment as needed.
"""

import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import triangle
from shapely.geometry import Polygon


def generateMesh(boundary_file, output_dir, mesh_options='pqa0.00003'):
    """Generate mesh from boundary, save centroids and mesh to CSV."""
    
    print("=" * 60)
    print("BERAMBADI WATERSHED MESH GENERATOR")
    print("=" * 60)
    
    # 1. Read boundary
    with open(boundary_file, 'r') as f:
        watershed = json.load(f)
    coords = watershed['features'][0]['geometry']['coordinates'][0]
    polygon = Polygon(coords)
    
    print(f"\n[1] Boundary loaded from: {boundary_file}")
    print(f"    Number of boundary points: {len(coords)}")
    print(f"    Polygon area: {polygon.area:.6f} sq degrees")
    print(f"    Bounds: {polygon.bounds}")
    
    # 2. Prepare for triangulation
    vertices = np.array(coords[:-1])
    n_vertices = len(vertices)
    segments = np.array([(i, (i + 1) % n_vertices) for i in range(n_vertices)])
    mesh_input = {'vertices': vertices, 'segments': segments}
    
    print(f"\n[2] Prepared mesh input:")
    print(f"    Vertices: {len(vertices)}")
    print(f"    Segments: {len(segments)}")
    
    # 3. Triangulate
    print(f"\n[3] Triangulating with options: '{mesh_options}'")
    mesh = triangle.triangulate(mesh_input, mesh_options)
    triangles_idx = mesh['triangles']
    mesh_vertices = mesh['vertices']
    print(f"    Triangles created: {len(triangles_idx)}")
    print(f"    Mesh vertices: {len(mesh_vertices)}")
    
    # 4. Compute centroids
    centroids = []
    triangle_areas = []
    for tri in triangles_idx:
        pts = mesh_vertices[tri]
        centroids.append(pts.mean(axis=0))
        x, y = pts[:, 0], pts[:, 1]
        triangle_areas.append(0.5 * abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1))))
    centroids = np.array(centroids)
    triangle_areas = np.array(triangle_areas)
    
    print(f"\n[4] Computed centroids:")
    print(f"    Number of centroids: {len(centroids)}")
    print(f"    Average triangle area: {triangle_areas.mean():.8f} sq degrees")
    
    # 5. Visualize
    print(f"\n[5] Creating visualization...")
    fig, axes = plt.subplots(1, 2, figsize=(16, 7), dpi=100)
    
    ax1 = axes[0]
    for tri in triangles_idx:
        pts = mesh_vertices[tri]
        pts_closed = np.vstack([pts, pts[0]])
        ax1.plot(pts_closed[:, 0], pts_closed[:, 1], 'b-', linewidth=0.3, alpha=0.7)
    ax1.scatter(centroids[:, 0], centroids[:, 1], c='red', s=8, zorder=5, label=f'Centroids ({len(centroids)})')
    boundary = np.array(coords)
    ax1.plot(boundary[:, 0], boundary[:, 1], 'k-', linewidth=2, label='Boundary')
    ax1.scatter(vertices[:, 0], vertices[:, 1], c='green', s=20, zorder=6, label=f'Boundary pts ({len(vertices)})')
    ax1.set_xlabel('Longitude')
    ax1.set_ylabel('Latitude')
    ax1.set_title(f'Berambadi Watershed - Triangular Mesh\n{len(triangles_idx)} triangles')
    ax1.legend(loc='upper right')
    ax1.set_aspect('equal')
    ax1.grid(True, alpha=0.3)
    
    ax2 = axes[1]
    sc = ax2.scatter(centroids[:, 0], centroids[:, 1], c=np.arange(len(centroids)), cmap='viridis', s=15, alpha=0.8)
    ax2.plot(boundary[:, 0], boundary[:, 1], 'k-', linewidth=2)
    ax2.set_xlabel('Longitude')
    ax2.set_ylabel('Latitude')
    ax2.set_title(f'Centroid Locations\n(Ready for GEE extraction)')
    ax2.set_aspect('equal')
    ax2.grid(True, alpha=0.3)
    plt.colorbar(sc, ax=ax2, label='Centroid Index')
    
    plt.tight_layout()
    plt.savefig(output_dir + 'berambadi_mesh.png', dpi=150, bbox_inches='tight')
    print(f"    Saved: {output_dir}berambadi_mesh.png")
    plt.show()
    
    # 6. Save centroids
    centroid_df = pd.DataFrame({
        'id': range(len(centroids)),
        'lon': centroids[:, 0],
        'lat': centroids[:, 1],
        'triangle_area': triangle_areas
    })
    centroid_df.to_csv(output_dir + 'berambadi_centroids.csv', index=False)
    print(f"\n[6] Saved centroid coordinates: {output_dir}berambadi_centroids.csv")
    
    # 7. Save mesh vertices and triangles for plotting later
    np.savetxt(output_dir + 'mesh_vertices.csv', mesh_vertices, delimiter=',', header='lon,lat', comments='')
    np.savetxt(output_dir + 'mesh_triangles.csv', triangles_idx, delimiter=',', fmt='%d', header='v0,v1,v2', comments='')
    print(f"[7] Saved mesh data:")
    print(f"    - {output_dir}mesh_vertices.csv")
    print(f"    - {output_dir}mesh_triangles.csv")
    
    print("\n" + "=" * 60)
    print("MESH GENERATION COMPLETE!")
    print("=" * 60)


def queryGEE(centroids_csv, boundary_file, output_dir, date_start='2022-08-01', date_end='2022-12-31'):
    """Query Sentinel-1 backscatter at each centroid location for ALL available dates."""
    
    print("\n" + "=" * 60)
    print("GOOGLE EARTH ENGINE EXTRACTION")
    print("=" * 60)
    
    # 1. Load centroids
    centroid_df = pd.read_csv(centroids_csv)
    print(f"\n[1] Loading centroids from: {centroids_csv}")
    print(f"    Loaded {len(centroid_df)} centroids")
    
    # 2. Load boundary for AOI
    with open(boundary_file, 'r') as f:
        watershed = json.load(f)
    coords = watershed['features'][0]['geometry']['coordinates'][0]
    aoi = ee.Geometry.Polygon(coords)
    
    # 3. Get all available Sentinel-1 images in date range
    print(f"\n[2] Loading Sentinel-1 imagery...")
    print(f"    Date range: {date_start} to {date_end}")
    
    collection = ee.ImageCollection('COPERNICUS/S1_GRD') \
        .filterBounds(aoi) \
        .filterDate(date_start, date_end)
    
    # Get list of dates
    dates = collection.aggregate_array('system:time_start').getInfo()
    dates = sorted(set([ee.Date(d).format('YYYY-MM-dd').getInfo() for d in dates]))
    
    print(f"    Found {len(dates)} available dates:")
    for d in dates:
        print(f"      - {d}")
    
    # 4. Query each date and each centroid
    print(f"\n[3] Querying backscatter...")
    print(f"    Total queries: {len(dates)} dates × {len(centroid_df)} centroids = {len(dates) * len(centroid_df)}")
    print("    This may take a while...\n")
    
    results = []
    
    for date_idx, date_str in enumerate(dates):
        print(f"[Date {date_idx+1}/{len(dates)}] Querying: {date_str}")
        
        # Get image for this date
        image = collection.filterDate(date_str, ee.Date(date_str).advance(1, 'day')).first().clip(aoi)
        
        for i, row in centroid_df.iterrows():
            lon, lat = row['lon'], row['lat']
            
            if (i + 1) % 50 == 0:
                print(f"    Processing centroid: {i+1}/{len(centroid_df)}")
            
            point = ee.Geometry.Point([lon, lat])
            value = image.reduceRegion(reducer=ee.Reducer.first(), geometry=point, scale=10).getInfo()
            
            results.append({
                'date': date_str,
                'id': row['id'],
                'lon': lon,
                'lat': lat,
                'VV': value.get('VV'),
                'VH': value.get('VH'),
                'angle': value.get('angle')
            })
        
        print(f"    ✓ Completed {date_str}")
    
    print(f"\n    Done! Queried {len(results)} total points.")
    
    # 5. Save results
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_dir + 'berambadi_backscatter.csv', index=False)
    
    print(f"\n[4] Saved: {output_dir}berambadi_backscatter.csv")
    print("\nResults preview:")
    print(results_df.head(15))
    print(f"\nDates in file: {results_df['date'].nunique()}")
    print(f"Centroids per date: {len(centroid_df)}")
    print(f"Total rows: {len(results_df)}")
    
    print("\n" + "=" * 60)
    print("GEE EXTRACTION COMPLETE!")
    print("=" * 60)
    
    return results_df


def plotBackscatter(backscatter_csv, mesh_vertices_csv, mesh_triangles_csv, output_dir,
                    day=1, quantity='VV', cmap='viridis'):
    """
    Plot backscatter as filled triangles (FEM style).
    """
    
    # Load data
    df = pd.read_csv(backscatter_csv)
    vertices = np.loadtxt(mesh_vertices_csv, delimiter=',', skiprows=1)
    triangles = np.loadtxt(mesh_triangles_csv, delimiter=',', skiprows=1, dtype=int)
    
    # Get available dates
    dates = sorted(df['date'].unique())
    print(f"Available dates ({len(dates)} total):")
    for i, d in enumerate(dates):
        print(f"  Day {i+1}: {d}")
    
    # Select date
    if day < 1 or day > len(dates):
        print(f"\nError: day must be between 1 and {len(dates)}")
        return
    
    selected_date = dates[day - 1]
    print(f"\n→ Plotting Day {day}: {selected_date}, Quantity: {quantity}")
    
    # Get backscatter values for selected date
    data = df[df['date'] == selected_date].sort_values('id')
    values = data[quantity].values
    
    # Compute color limits using percentiles
    vmin = np.percentile(values, 5)
    vmax = np.percentile(values, 95)
    
    # Create filled triangle polygons
    polys = [vertices[tri] for tri in triangles]
    
    # Plot
    fig, ax = plt.subplots(figsize=(10, 6), dpi=100)
    
    collection = PolyCollection(polys, array=values, cmap=cmap, edgecolors='none')
    collection.set_clim(vmin, vmax)
    ax.add_collection(collection)
    ax.autoscale()
    
    ax.set_title(f'{quantity} Backscatter (dB) - {selected_date}')
    ax.set_aspect('equal')
    ax.axis('off')
    
    plt.colorbar(collection, ax=ax, label=f'{quantity} (dB)', shrink=0.7)
    
    plt.tight_layout()
    save_path = output_dir + f'backscatter_{quantity}_{selected_date}.png'
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    print(f"Saved: {save_path}")
    plt.show()

# ============================================
# MAIN - RUN STAGE BY STAGE
# ============================================
if __name__ == "__main__":
    
    # Paths
    BOUNDARY_FILE = 'boundary.geojson'
    OUTPUT_DIR = 'data/'
    GEE_PROJECT = 'your-gee-project-id'  # Replace with your GEE project ID
    DATE_START = '2022-08-01'
    DATE_END = '2022-12-31'
    
    # Derived paths
    CENTROIDS_CSV = OUTPUT_DIR + 'centroids.csv'
    BACKSCATTER_CSV = OUTPUT_DIR + 'backscatter.csv'
    MESH_VERTICES_CSV = OUTPUT_DIR + 'mesh_vertices.csv'
    MESH_TRIANGLES_CSV = OUTPUT_DIR + 'mesh_triangles.csv'
    
    # ------------------------------------------
    # STAGE 1: Generate Mesh (re-run to save mesh files)
    # ------------------------------------------
    # generateMesh(
    #     boundary_file=BOUNDARY_FILE,
    #     output_dir=OUTPUT_DIR,
    #     mesh_options='pqa0.00003'
    # )
    
    # ------------------------------------------
    # STAGE 2: Authenticate GEE (skip - already done)
    # ------------------------------------------
    # import ee
    # ee.Authenticate()
    # ee.Initialize(project='ee-berambadi')
    # print("\n✓ Earth Engine initialized!")
    
    # ------------------------------------------
    # STAGE 3: Query GEE (skip - already done)
    # ------------------------------------------
    # queryGEE(
    #     centroids_csv=CENTROIDS_CSV,
    #     boundary_file=BOUNDARY_FILE,
    #     output_dir=OUTPUT_DIR,
    #     date_start='2022-08-01',
    #     date_end='2022-12-31'
    # )
    
    # ------------------------------------------
    # STAGE 4: Generate all plots with FIXED colorbar
    # ------------------------------------------
    import glob
    from PIL import Image
    
    df = pd.read_csv(BACKSCATTER_CSV)
    vertices = np.loadtxt(MESH_VERTICES_CSV, delimiter=',', skiprows=1)
    triangles = np.loadtxt(MESH_TRIANGLES_CSV, delimiter=',', skiprows=1, dtype=int)
    
    quantity = 'VH'
    dates = sorted(df['date'].unique())
    num_dates = len(dates)
    print(f"Total dates: {num_dates}")
    
    # Fixed colorbar limits (global across all dates)
    all_values = df[quantity].values
    vmin = np.percentile(all_values, 5)
    vmax = np.percentile(all_values, 95)
    print(f"Fixed color range: {vmin:.2f} to {vmax:.2f} dB")
    
    # Create polygons once
    polys = [vertices[tri] for tri in triangles]
    
    # Generate all PNGs
    for day_idx, selected_date in enumerate(dates):
        print(f"Plotting [{day_idx+1}/{num_dates}]: {selected_date}")
        
        data = df[df['date'] == selected_date].sort_values('id')
        values = data[quantity].values
        
        fig, ax = plt.subplots(figsize=(10, 6), dpi=100)
        collection = PolyCollection(polys, array=values, cmap='viridis', edgecolors='none')
        collection.set_clim(vmin, vmax)
        ax.add_collection(collection)
        ax.autoscale()
        ax.set_title(f'{quantity} Backscatter (dB) - {selected_date}')
        ax.set_aspect('equal')
        ax.axis('off')
        plt.colorbar(collection, ax=ax, label=f'{quantity} (dB)', shrink=0.7)
        
        plt.tight_layout()
        plt.savefig(OUTPUT_DIR + f'backscatter_{quantity}_{selected_date}.png', dpi=150, bbox_inches='tight')
        plt.close('all')
    
    print("All PNGs generated!")
    
    # ------------------------------------------
    # STAGE 5: Create GIF
    # ------------------------------------------
    png_files = sorted(glob.glob(OUTPUT_DIR + f'backscatter_{quantity}_*.png'))
    print(f"Found {len(png_files)} images")
    
    images = [Image.open(f) for f in png_files]
    gif_path = OUTPUT_DIR + f'backscatter_{quantity}_timeseries.gif'
    images[0].save(
        gif_path,
        save_all=True,
        append_images=images[1:],
        duration=500,
        loop=0
    )
    print(f"Saved: {gif_path}")
