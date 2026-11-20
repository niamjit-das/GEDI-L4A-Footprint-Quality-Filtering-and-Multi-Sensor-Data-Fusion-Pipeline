# GEDI-L4A-Footprint-Quality-Filtering-and-Multi-Sensor-Data-Fusion-Pipeline
 GEDI L4A Footprint Quality Filtering and Multi-Sensor Data Fusion Pipeline 
 // Study Area Matrix: 18,439.58 ha Landscape Boundary
// Author: Niamjit Das
// UPDATED: per-image sampling, optional buffer, configurable filters

// 0. CONFIGURATION (change these to trade sample size vs. data quality)

var BUFFER_M        = 3000;          // buffer around MNP in metres (0 = strict boundary)
var SE_RATIO_MAX    = 0.5;           // relative SE threshold: try 0.75 or 1.0 for more samples
var REQUIRE_NON_DEGRADED = false;    // true = also require degrade_flag == 0 (stricter)
var GEDI_START      = '2019-04-01';
var GEDI_END        = '2025-12-31';  // later months are used only if present in the collection
var S2_START        = '2021-01-01';  // dry-season composite, closer to GEDI period
var S2_END          = '2021-03-31';  // (set to 2026 dates if you prefer your original composite)
var FOREST_NDVI     = 0.68;          // forest / non-forest threshold
var EDGE_MAX_M      = 2000;          // max distance-to-edge radius


// 1. Load and Delineate the Landscape (with optional buffer)

var MNP = ee.FeatureCollection('projects/silent-window-479006-n4/assets/MNP_ShapeFile');
var roi = BUFFER_M > 0 ? MNP.geometry().buffer(BUFFER_M) : MNP.geometry();

Map.centerObject(MNP, 11);
Map.addLayer(roi, {color: 'orange'}, 'Sampling ROI (buffered)', false);
Map.addLayer(MNP, {color: 'darkgreen'}, 'Madhupur Sal Forest Matrix');


// 2. Load and Cloud-Mask Sentinel-2 Level-2A Surface Reflectance

function maskS2clouds(image) {
  var qa = image.select('QA60');
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
    .and(qa.bitwiseAnd(cirrusBitMask).eq(0));
  return image.updateMask(mask).divide(10000)
    .copyProperties(image, ['system:time_start']);
}

// Clipped to the buffered ROI (not MNP) so buffer footprints are not masked out
var sentinel = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(roi)
  .filterDate(S2_START, S2_END)
  .map(maskS2clouds)
  .median()
  .clip(roi);


// 3. Vegetation Predictors

var ndvi = sentinel.normalizedDifference(['B8', 'B4']).rename('NDVI');

var evi = sentinel.expression(
  '2.5 * ((B8 - B4) / (B8 + 6 * B4 - 7.5 * B2 + 1))', {
    'B8': sentinel.select('B8'),
    'B4': sentinel.select('B4'),
    'B2': sentinel.select('B2')
  }).rename('EVI');

var ndre  = sentinel.normalizedDifference(['B8', 'B5']).rename('NDRE');
var swir1 = sentinel.select('B11').rename('B11');
var swir2 = sentinel.select('B12').rename('B12');


// 4. Continuous Euclidean Distance-to-Edge Predictor

var forestMask = ndvi.gt(FOREST_NDVI);
var nonForest  = forestMask.not();
var distanceToEdge = nonForest
  .distance(ee.Kernel.euclidean(EDGE_MAX_M, 'meters'))
  .unmask(EDGE_MAX_M)               // interior forest beyond max radius = max distance, not dropped
  .rename('distance_to_edge');

var predictorStack = ndvi.addBands([evi, ndre, swir1, swir2, distanceToEdge]);


// 5. GEDI L4A Monthly Rasters (25 m footprints) + Quality Mask

var gediCol = ee.ImageCollection('LARSE/GEDI/GEDI04_A_002_MONTHLY')
  .filterBounds(roi)
  .filterDate(GEDI_START, GEDI_END);

print('GEDI monthly images in range:', gediCol.size());

var qualityMask = function(img) {
  var agbd = img.select('agbd');
  var se   = img.select('agbd_se');
  var mask = img.select('l4_quality_flag').eq(1)
    .and(img.select('algorithm_run_flag').eq(1))
    .and(se.lte(agbd.multiply(SE_RATIO_MAX)));   // relative SE (avoids divide-by-zero)
  if (REQUIRE_NON_DEGRADED) {
    mask = mask.and(img.select('degrade_flag').eq(0));
  }
  return img.select(['agbd', 'agbd_se']).updateMask(mask);
};


// 6. Multi-Sensor Fusion: sample EACH monthly image separately
//    (mosaic() would collapse repeat hits on the same pixel into one row)

var perMonth = gediCol.map(function(img) {
  var g = qualityMask(img);
  return predictorStack.addBands(g)
    .sample({
      region: roi,
      scale: 25,                    // match GEDI footprint size
      projection: img.projection(),
      geometries: true,
      tileScale: 4
    })
    .map(function(f) {
      var c = f.geometry().coordinates();
      return f.set({
        longitude: c.get(0),
        latitude: c.get(1),
        date: img.date().format('YYYY-MM')
      });
    });
}).flatten();

var fusedDataset = perMonth.filter(ee.Filter.notNull(['agbd', 'NDVI']));
print('Fused sample count', fusedDataset.size());
print('Sample preview', fusedDataset.limit(5));


// 7. Export CSV to Google Drive

Export.table.toDrive({
  collection: fusedDataset,
  description: 'Gedi_Sentinel2_Landscape_Matrix_Fused',
  fileFormat: 'CSV',
  selectors: ['longitude', 'latitude', 'date', 'NDVI', 'EVI', 'NDRE',
              'B11', 'B12', 'distance_to_edge', 'agbd', 'agbd_se']
});


// 8. Map Diagnostics

var cleanGediDisplay = gediCol.map(qualityMask).mosaic().clip(roi);

Map.addLayer(ndvi, {min: 0.2, max: 0.9, palette: ['yellow', 'green', 'darkgreen']}, 'NDVI Composite');
Map.addLayer(distanceToEdge, {min: 0, max: 500, palette: ['red', 'orange', 'white']}, 'Distance to Edge');
Map.addLayer(cleanGediDisplay.select('agbd'),
  {min: 0, max: 300, palette: ['yellow', 'green', 'darkgreen']}, 'Filtered GEDI AGBD (25 m)');
