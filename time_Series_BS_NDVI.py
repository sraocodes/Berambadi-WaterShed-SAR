"""
Single Location Time Series Extraction
==================================================
Extract ALL Sentinel-1 and Sentinel-2 data for a single location.
No filtering - keeps all observations including cloudy ones.
"""

import ee
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from scipy import stats
from datetime import datetime
import os


def initialize_gee(project_id='ee-berambadi'):
    """Initialize Google Earth Engine."""
    try:
        ee.Initialize(project=project_id)
        print("✓ Earth Engine initialized!")
    except:
        print("Authenticating Earth Engine...")
        ee.Authenticate()
        ee.Initialize(project=project_id)
        print("✓ Earth Engine authenticated and initialized!")


def extract_all_data(location, date_start, date_end, buffer_size=10):
    """
    Extract ALL Sentinel-1 and Sentinel-2 data for a location.
    Keeps cloudy observations - no filtering!
    
    Parameters:
    -----------
    location : tuple
        (longitude, latitude) of the point
    date_start : str
        Start date in 'YYYY-MM-DD' format
    date_end : str
        End date in 'YYYY-MM-DD' format
    buffer_size : int
        Buffer radius around point in meters
    
    Returns:
    --------
    DataFrame with all observations
    """
    
    print(f"\nExtracting data from {date_start} to {date_end}")
    print(f"Location: {location[1]:.6f}°N, {location[0]:.6f}°E")
    print("Keeping ALL observations (including cloudy ones)\n")
    
    # Create point geometry
    point = ee.Geometry.Point(location).buffer(buffer_size)
    
    all_records = []
    
    # =========================================
    # PART 1: SENTINEL-1 (RADAR - works in all weather)
    # =========================================
    print("[1/2] Processing Sentinel-1 (Radar)...")
    
    s1_collection = ee.ImageCollection('COPERNICUS/S1_GRD') \
        .filterBounds(point) \
        .filterDate(date_start, date_end) \
        .select(['VV', 'VH', 'angle'])
    
    s1_count = s1_collection.size().getInfo()
    print(f"      Found {s1_count} Sentinel-1 images")
    
    if s1_count > 0:
        # Process each S1 image
        def process_s1_image(image):
            values = image.reduceRegion(
                reducer=ee.Reducer.mean(),
                geometry=point,
                scale=10,
                maxPixels=1e9
            )
            
            return ee.Feature(None, {
                'date': ee.Date(image.get('system:time_start')).format('YYYY-MM-dd'),
                'timestamp': image.get('system:time_start'),
                'satellite': 'Sentinel-1',
                'VV_dB': values.get('VV'),
                'VH_dB': values.get('VH'),
                'angle': values.get('angle'),
                'VH_VV_ratio': ee.Number(values.get('VH')).divide(values.get('VV'))
            })
        
        s1_features = s1_collection.map(process_s1_image)
        s1_data = ee.FeatureCollection(s1_features).getInfo()
        
        # Convert to records
        for feature in s1_data['features']:
            props = feature['properties']
            if props.get('VV_dB') is not None:
                all_records.append({
                    'date': props['date'],
                    'timestamp': pd.to_datetime(props['timestamp'], unit='ms'),
                    'satellite': 'S1',
                    'VV_dB': props['VV_dB'],
                    'VH_dB': props['VH_dB'],
                    'angle': props['angle'],
                    'VH_VV_ratio': props['VH_VV_ratio'],
                    # S2 fields will be empty for S1 records
                    'cloud_percentage': np.nan,
                    'NDVI': np.nan,
                    'EVI': np.nan,
                    'NDWI': np.nan,
                    'NDMI': np.nan,
                    'Blue': np.nan,
                    'Green': np.nan,
                    'Red': np.nan,
                    'NIR': np.nan,
                    'SWIR1': np.nan,
                    'SWIR2': np.nan
                })
        
        print(f"      Processed {len([r for r in all_records if r['satellite'] == 'S1'])} S1 observations")
    
    # =========================================
    # PART 2: SENTINEL-2 (OPTICAL - affected by clouds)
    # =========================================
    print("\n[2/2] Processing Sentinel-2 (Optical)...")
    
    # NO CLOUD FILTER - we keep everything!
    s2_collection = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED') \
        .filterBounds(point) \
        .filterDate(date_start, date_end)
    
    s2_count = s2_collection.size().getInfo()
    print(f"      Found {s2_count} Sentinel-2 images (including cloudy)")
    
    if s2_count > 0:
        # Process each S2 image
        def process_s2_image(image):
            # Scale the bands (S2 values need to be divided by 10000)
            scaled = image.select(['B2', 'B3', 'B4', 'B8', 'B11', 'B12']).divide(10000)
            
            # Calculate indices
            nir = scaled.select('B8')
            red = scaled.select('B4')
            green = scaled.select('B3')
            blue = scaled.select('B2')
            swir1 = scaled.select('B11')
            
            # NDVI = (NIR - Red) / (NIR + Red)
            ndvi = nir.subtract(red).divide(nir.add(red)).rename('NDVI')
            
            # EVI = 2.5 * (NIR - Red) / (NIR + 6*Red - 7.5*Blue + 1)
            evi = nir.subtract(red).multiply(2.5).divide(
                nir.add(red.multiply(6)).subtract(blue.multiply(7.5)).add(1)
            ).rename('EVI')
            
            # NDWI = (Green - NIR) / (Green + NIR)
            ndwi = green.subtract(nir).divide(green.add(nir)).rename('NDWI')
            
            # NDMI = (NIR - SWIR1) / (NIR + SWIR1)
            ndmi = nir.subtract(swir1).divide(nir.add(swir1)).rename('NDMI')
            
            # Combine all bands
            combined = scaled.addBands([ndvi, evi, ndwi, ndmi])
            
            # Get values
            values = combined.reduceRegion(
                reducer=ee.Reducer.mean(),
                geometry=point,
                scale=10,
                maxPixels=1e9
            )
            
            return ee.Feature(None, {
                'date': ee.Date(image.get('system:time_start')).format('YYYY-MM-dd'),
                'timestamp': image.get('system:time_start'),
                'satellite': 'Sentinel-2',
                'cloud_percentage': image.get('CLOUDY_PIXEL_PERCENTAGE'),
                'NDVI': values.get('NDVI'),
                'EVI': values.get('EVI'),
                'NDWI': values.get('NDWI'),
                'NDMI': values.get('NDMI'),
                'Blue': values.get('B2'),
                'Green': values.get('B3'),
                'Red': values.get('B4'),
                'NIR': values.get('B8'),
                'SWIR1': values.get('B11'),
                'SWIR2': values.get('B12')
            })
        
        s2_features = s2_collection.map(process_s2_image)
        s2_data = ee.FeatureCollection(s2_features).getInfo()
        
        # Convert to records
        for feature in s2_data['features']:
            props = feature['properties']
            if props.get('NDVI') is not None:
                all_records.append({
                    'date': props['date'],
                    'timestamp': pd.to_datetime(props['timestamp'], unit='ms'),
                    'satellite': 'S2',
                    # S1 fields will be empty for S2 records
                    'VV_dB': np.nan,
                    'VH_dB': np.nan,
                    'angle': np.nan,
                    'VH_VV_ratio': np.nan,
                    # S2 fields
                    'cloud_percentage': props['cloud_percentage'],
                    'NDVI': props['NDVI'],
                    'EVI': props['EVI'],
                    'NDWI': props['NDWI'],
                    'NDMI': props['NDMI'],
                    'Blue': props['Blue'],
                    'Green': props['Green'],
                    'Red': props['Red'],
                    'NIR': props['NIR'],
                    'SWIR1': props['SWIR1'],
                    'SWIR2': props['SWIR2']
                })
        
        print(f"      Processed {len([r for r in all_records if r['satellite'] == 'S2'])} S2 observations")
        
        # Show cloud statistics
        cloud_stats = [r['cloud_percentage'] for r in all_records if r['satellite'] == 'S2']
        if cloud_stats:
            print(f"      Cloud cover: min={min(cloud_stats):.1f}%, max={max(cloud_stats):.1f}%, mean={np.mean(cloud_stats):.1f}%")
    
    # Convert to DataFrame and sort by date
    df = pd.DataFrame(all_records)
    if not df.empty:
        df = df.sort_values('timestamp').reset_index(drop=True)
        df['latitude'] = location[1]
        df['longitude'] = location[0]
        
        # Reorder columns for better readability
        column_order = ['date', 'timestamp', 'satellite', 'latitude', 'longitude',
                       'cloud_percentage', 'VV_dB', 'VH_dB', 'angle', 'VH_VV_ratio',
                       'NDVI', 'EVI', 'NDWI', 'NDMI', 
                       'Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2']
        df = df[column_order]
    
    print(f"\n✓ Total observations: {len(df)}")
    return df


def remove_outliers(data, column, z_threshold=3):
    """
    Remove outliers using Z-score method.
    Returns a copy of the data with outliers replaced by NaN.
    """
    import scipy.stats as stats
    
    data_copy = data.copy()
    
    if column in data_copy.columns and not data_copy[column].isna().all():
        # Remove impossible values based on physical limits
        if 'NDVI' in column or 'EVI' in column or 'NDWI' in column or 'NDMI' in column:
            # These indices should be between -1 and 1
            data_copy.loc[(data_copy[column] > 1) | (data_copy[column] < -1), column] = np.nan
        
        # Remove statistical outliers
        valid_data = data_copy[column].dropna()
        if len(valid_data) > 0:
            z_scores = np.abs(stats.zscore(valid_data))
            outlier_mask = z_scores > z_threshold
            outlier_indices = valid_data[outlier_mask].index
            data_copy.loc[outlier_indices, column] = np.nan
    
    return data_copy


def create_presentation_plots(df, output_dir):
    """
    Create individual presentation-ready plots for slideshow.
    Each plot is saved as a separate PNG file with consistent size and formatting.
    """
    
    if df.empty:
        print("No data to plot")
        return
    
    # Set up matplotlib for presentation graphics
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    
    plt.rcParams['figure.dpi'] = 100
    plt.rcParams['savefig.dpi'] = 200
    plt.rcParams['font.size'] = 18
    plt.rcParams['axes.labelsize'] = 20
    plt.rcParams['axes.titlesize'] = 24
    plt.rcParams['xtick.labelsize'] = 16
    plt.rcParams['ytick.labelsize'] = 16
    plt.rcParams['legend.fontsize'] = 18
    
    # Standard figure size for all plots (16:9 aspect ratio for presentations)
    fig_size = (16, 9)
    
    # Separate S1 and S2 data
    s1_data = df[df['satellite'] == 'S1'].copy()
    s2_data = df[df['satellite'] == 'S2'].copy()
    
    print("\nCreating presentation plots...")
    
    # =========================================
    # PLOT 1: Sentinel-1 VV Backscatter
    # =========================================
    if not s1_data.empty:
        s1_clean = remove_outliers(s1_data, 'VV_dB')
        
        fig, ax = plt.subplots(figsize=fig_size)
        ax.scatter(s1_clean['timestamp'], s1_clean['VV_dB'], 
                  color='navy', s=40, alpha=0.7, edgecolors='black', linewidth=0.5)
        
        ax.set_xlabel('Date', fontsize=20, fontweight='bold')
        ax.set_ylabel('VV Backscatter (dB)', fontsize=20, fontweight='bold')
        ax.set_title('Sentinel-1 VV Polarization', fontsize=26, fontweight='bold', pad=20)
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.xaxis.set_major_locator(mdates.YearLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
        
        # Add statistics box
        vv_mean = s1_clean['VV_dB'].dropna().mean()
        vv_std = s1_clean['VV_dB'].dropna().std()
        textstr = f'Mean: {vv_mean:.1f} dB\nStd: {vv_std:.1f} dB\nn = {s1_clean["VV_dB"].notna().sum()}'
        ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=16,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, '01_S1_VV_backscatter.png'), 
                   bbox_inches='tight', facecolor='white')
        plt.close()
        print("  ✓ Saved: 01_S1_VV_backscatter.png")
    
    # =========================================
    # PLOT 2: Sentinel-1 VH Backscatter
    # =========================================
    if not s1_data.empty:
        s1_clean = remove_outliers(s1_data, 'VH_dB')
        
        fig, ax = plt.subplots(figsize=fig_size)
        ax.scatter(s1_clean['timestamp'], s1_clean['VH_dB'], 
                  color='darkred', s=40, alpha=0.7, edgecolors='black', linewidth=0.5)
        
        ax.set_xlabel('Date', fontsize=20, fontweight='bold')
        ax.set_ylabel('VH Backscatter (dB)', fontsize=20, fontweight='bold')
        ax.set_title('Sentinel-1 VH Polarization', fontsize=26, fontweight='bold', pad=20)
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.xaxis.set_major_locator(mdates.YearLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
        
        # Add statistics box
        vh_mean = s1_clean['VH_dB'].dropna().mean()
        vh_std = s1_clean['VH_dB'].dropna().std()
        textstr = f'Mean: {vh_mean:.1f} dB\nStd: {vh_std:.1f} dB\nn = {s1_clean["VH_dB"].notna().sum()}'
        ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=16,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.8))
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, '02_S1_VH_backscatter.png'), 
                   bbox_inches='tight', facecolor='white')
        plt.close()
        print("  ✓ Saved: 02_S1_VH_backscatter.png")
    
    # =========================================
    # PLOT 3: NDVI with Cloud Information
    # =========================================
    if not s2_data.empty:
        s2_clean = remove_outliers(s2_data, 'NDVI')
        
        fig, ax = plt.subplots(figsize=fig_size)
        
        # Color by cloud percentage
        scatter = ax.scatter(s2_clean['timestamp'], s2_clean['NDVI'], 
                           c=s2_clean['cloud_percentage'], cmap='RdYlGn_r',
                           s=50, alpha=0.8, vmin=0, vmax=100,
                           edgecolors='black', linewidth=0.5)
        
        cbar = plt.colorbar(scatter, ax=ax, pad=0.02)
        cbar.set_label('Cloud Coverage (%)', fontsize=18)
        
        ax.set_xlabel('Date', fontsize=20, fontweight='bold')
        ax.set_ylabel('NDVI', fontsize=20, fontweight='bold')
        ax.set_title('Vegetation Index (NDVI)', fontsize=26, fontweight='bold', pad=20)
        ax.set_ylim(-0.1, 1.0)
        
        # Add reference lines
        ax.axhline(y=0.2, color='orange', linestyle='--', alpha=0.5, linewidth=2)
        ax.axhline(y=0.5, color='green', linestyle='--', alpha=0.5, linewidth=2)
        ax.axhline(y=0.7, color='darkgreen', linestyle='--', alpha=0.5, linewidth=2)
        
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.xaxis.set_major_locator(mdates.YearLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
        
        # Statistics for clear days
        clear_data = s2_clean[s2_clean['cloud_percentage'] < 20]
        if not clear_data.empty:
            ndvi_mean = clear_data['NDVI'].dropna().mean()
            textstr = f'Clear Sky Mean: {ndvi_mean:.3f}\nn_clear = {len(clear_data)}\nn_total = {len(s2_clean)}'
            ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=16,
                    verticalalignment='top', bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, '03_NDVI_timeseries.png'), 
                   bbox_inches='tight', facecolor='white')
        plt.close()
        print("  ✓ Saved: 03_NDVI_timeseries.png")
    
    # =========================================
    # PLOT 4: Water Content Index (NDWI)
    # =========================================
    if not s2_data.empty:
        s2_clean = remove_outliers(s2_data, 'NDWI')
        # Only use clear observations for water indices
        s2_clear = s2_clean[s2_clean['cloud_percentage'] < 30].copy()
        
        fig, ax = plt.subplots(figsize=fig_size)
        ax.scatter(s2_clear['timestamp'], s2_clear['NDWI'], 
                  color='dodgerblue', s=40, alpha=0.7, edgecolors='black', linewidth=0.5)
        
        ax.set_xlabel('Date', fontsize=20, fontweight='bold')
        ax.set_ylabel('NDWI', fontsize=20, fontweight='bold')
        ax.set_title('Water Content Index (NDWI)', fontsize=26, fontweight='bold', pad=20)
        ax.set_ylim(-0.8, 0.8)
        ax.axhline(y=0, color='gray', linestyle='-', alpha=0.3, linewidth=1)
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.xaxis.set_major_locator(mdates.YearLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
        
        # Add statistics
        ndwi_mean = s2_clear['NDWI'].dropna().mean()
        ndwi_std = s2_clear['NDWI'].dropna().std()
        textstr = f'Mean: {ndwi_mean:.3f}\nStd: {ndwi_std:.3f}\n(cloud < 30%)'
        ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=16,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, '04_NDWI_water_content.png'), 
                   bbox_inches='tight', facecolor='white')
        plt.close()
        print("  ✓ Saved: 04_NDWI_water_content.png")
    
    # =========================================
    # PLOT 5: Moisture Index (NDMI)
    # =========================================
    if not s2_data.empty:
        s2_clean = remove_outliers(s2_data, 'NDMI')
        s2_clear = s2_clean[s2_clean['cloud_percentage'] < 30].copy()
        
        fig, ax = plt.subplots(figsize=fig_size)
        ax.scatter(s2_clear['timestamp'], s2_clear['NDMI'], 
                  color='teal', s=40, alpha=0.7, edgecolors='black', linewidth=0.5)
        
        ax.set_xlabel('Date', fontsize=20, fontweight='bold')
        ax.set_ylabel('NDMI', fontsize=20, fontweight='bold')
        ax.set_title('Vegetation Moisture Index (NDMI)', fontsize=26, fontweight='bold', pad=20)
        ax.set_ylim(-0.8, 0.8)
        ax.axhline(y=0, color='gray', linestyle='-', alpha=0.3, linewidth=1)
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.xaxis.set_major_locator(mdates.YearLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
        
        # Add statistics
        ndmi_mean = s2_clear['NDMI'].dropna().mean()
        ndmi_std = s2_clear['NDMI'].dropna().std()
        textstr = f'Mean: {ndmi_mean:.3f}\nStd: {ndmi_std:.3f}\n(cloud < 30%)'
        ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=16,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='lightcyan', alpha=0.8))
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, '05_NDMI_moisture.png'), 
                   bbox_inches='tight', facecolor='white')
        plt.close()
        print("  ✓ Saved: 05_NDMI_moisture.png")
    
    # =========================================
    # PLOT 6: Cloud Coverage Distribution
    # =========================================
    if not s2_data.empty:
        fig, ax = plt.subplots(figsize=fig_size)
        
        # Create histogram
        ax.hist(s2_data['cloud_percentage'], bins=20, range=(0, 100), 
               color='skyblue', edgecolor='black', alpha=0.7)
        
        ax.set_xlabel('Cloud Coverage (%)', fontsize=20, fontweight='bold')
        ax.set_ylabel('Number of Observations', fontsize=20, fontweight='bold')
        ax.set_title('Cloud Coverage Distribution', fontsize=26, fontweight='bold', pad=20)
        ax.grid(True, alpha=0.3, linestyle='--', axis='y')
        
        # Add vertical lines for categories
        ax.axvline(x=20, color='green', linestyle='--', linewidth=2, label='Clear (<20%)')
        ax.axvline(x=50, color='orange', linestyle='--', linewidth=2, label='Partial (20-50%)')
        ax.axvline(x=80, color='red', linestyle='--', linewidth=2, label='Cloudy (>50%)')
        
        # Add statistics
        clear_count = len(s2_data[s2_data['cloud_percentage'] < 20])
        partial_count = len(s2_data[(s2_data['cloud_percentage'] >= 20) & (s2_data['cloud_percentage'] < 50)])
        cloudy_count = len(s2_data[s2_data['cloud_percentage'] >= 50])
        
        textstr = f'Clear: {clear_count} ({100*clear_count/len(s2_data):.1f}%)\nPartial: {partial_count} ({100*partial_count/len(s2_data):.1f}%)\nCloudy: {cloudy_count} ({100*cloudy_count/len(s2_data):.1f}%)'
        ax.text(0.98, 0.98, textstr, transform=ax.transAxes, fontsize=16,
                verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        ax.legend(loc='upper left', fontsize=16)
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, '06_cloud_coverage.png'), 
                   bbox_inches='tight', facecolor='white')
        plt.close()
        print("  ✓ Saved: 06_cloud_coverage.png")
    
    print(f"\n✓ Created {6} presentation-ready plots in {output_dir}")


def print_summary_statistics(df):
    """
    Print summary statistics.
    """
    
    print("\n" + "="*60)
    print("SUMMARY STATISTICS")
    print("="*60)
    
    s1_data = df[df['satellite'] == 'S1']
    s2_data = df[df['satellite'] == 'S2']
    
    if not s1_data.empty:
        print("\n[Sentinel-1 Radar]")
        print(f"  Observations: {len(s1_data)}")
        print(f"  Date range: {s1_data['date'].min()} to {s1_data['date'].max()}")
        print(f"  VV mean: {s1_data['VV_dB'].mean():.2f} dB")
        print(f"  VH mean: {s1_data['VH_dB'].mean():.2f} dB")
    
    if not s2_data.empty:
        print("\n[Sentinel-2 Optical]")
        print(f"  Observations: {len(s2_data)}")
        print(f"  Date range: {s2_data['date'].min()} to {s2_data['date'].max()}")
        
        # Cloud statistics
        print(f"\n  Cloud coverage:")
        print(f"    < 10%: {len(s2_data[s2_data['cloud_percentage'] < 10])} observations")
        print(f"    10-30%: {len(s2_data[(s2_data['cloud_percentage'] >= 10) & (s2_data['cloud_percentage'] < 30)])} observations")
        print(f"    30-50%: {len(s2_data[(s2_data['cloud_percentage'] >= 30) & (s2_data['cloud_percentage'] < 50)])} observations")
        print(f"    > 50%: {len(s2_data[s2_data['cloud_percentage'] >= 50])} observations")
        
        # NDVI statistics (for cloud-free observations)
        clear_data = s2_data[s2_data['cloud_percentage'] < 20]
        if not clear_data.empty:
            print(f"\n  NDVI (clear days, <20% cloud):")
            print(f"    Mean: {clear_data['NDVI'].mean():.3f}")
            print(f"    Max: {clear_data['NDVI'].max():.3f}")
            print(f"    Min: {clear_data['NDVI'].min():.3f}")


# ============================================
# MAIN EXECUTION
# ============================================
if __name__ == "__main__":
    
    # =========================================
    # USER CONFIGURATION!
    # =========================================
    
    LOCATION = (76.587991, 11.761146)  # Exact Berambadi site from data paper
    
    # Time period (for last decade, use 2014-2024)
    DATE_START = '2014-01-01'
    DATE_END = '2024-12-31'
    
    # Output directory
    OUTPUT_DIR = './data/'
    
    # GEE Project ID
    PROJECT_ID = 'ee'  # Change to your project
    
    # =========================================
    # EXECUTION
    # =========================================
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    print("="*60)
    print("SIMPLIFIED SATELLITE TIME SERIES EXTRACTION")
    print("="*60)
    
    # Step 1: Initialize GEE
    initialize_gee(project_id=PROJECT_ID)
    
    # Step 2: Extract ALL data (no filtering)
    df = extract_all_data(
        location=LOCATION,
        date_start=DATE_START,
        date_end=DATE_END,
        buffer_size=10  # 10m buffer
    )
    
    # Step 3: Save to single CSV
    if not df.empty:
        output_file = os.path.join(OUTPUT_DIR, 'satellite_timeseries.csv')
        df.to_csv(output_file, index=False)
        print(f"\n✓ Data saved: {output_file}")
        
        # Step 4: Create presentation plots
        create_presentation_plots(df, OUTPUT_DIR)
        
        # Step 5: Print summary
        print_summary_statistics(df)
    else:
        print("\n❌ No data found for the specified location and time period")
    
    print("\n" + "="*60)
    print("EXTRACTION COMPLETE!")
    print("="*60)