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



