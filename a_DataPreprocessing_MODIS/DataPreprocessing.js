// AUTOMATICALLY GENERATED: location from saved link.
Map.setCenter(264.8, 34.8, 4)

/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var WaterBoundary = ee.FeatureCollection("users/blackwhitezou/shapefile"),
    MODIS = ee.ImageCollection("MODIS/006/MOD09A1"),
    StudyArea = ee.FeatureCollection("users/blackwhitezou/StudyAreaV4");
/***** End of imports. If edited, may not auto-convert in the playground. *****/
//************** Minor thesis GEE script **************//
// Author: Yuting Zou
// Date: June 2020
//*****************************************************//

// **************** Step1: Preprocessing for MODIS data ********************//
// Define dates
var iniDate = ee.Date.fromYMD(2000,1,1);
var endDate = ee.Date.fromYMD(2017,12,31);
// bands
var modisBands = ['sur_refl_b03','sur_refl_b04','sur_refl_b01','sur_refl_b02','sur_refl_b05', 'sur_refl_b06','sur_refl_b07'];
var lsBands = ['blue','green','red','nir','swir1','swir2','swir3'];
 
// helper function to extract the QA bits
function getQABits(image, start, end, newName) {
// Compute the bits we need to extract.
var pattern = 0;
for (var i = start; i <= end; i++) {
pattern += Math.pow(2, i);
}
// Return a single band image of the extracted QA bits, giving the band
// a new name.
return image.select([0], [newName])
.bitwiseAnd(pattern)
.rightShift(start);
}
 
// A function to mask out cloudy pixels.
function maskQuality(image) {
// Select the QA band.
var QA = image.select('StateQA');
// Get the internal_cloud_algorithm_flag bit.
var internalQuality = getQABits(QA,8, 13, 'internal_quality_flag');
// Return an image masking out cloudy areas.
  return image.updateMask(internalQuality.eq(0)).clipToCollection(StudyArea);
}
 
// A function to mask out cloudy pixels.
function maskQuality1(image) {
// Select the QA band.
var QA = image.select('StateQA');
// Get the internal_cloud_algorithm_flag bit.
var internalQuality1 = getQABits(QA,0, 2, 'internal_quality_flag1');
// Return an image masking out cloudy areas.
  return image.updateMask(internalQuality1.eq(0)).clipToCollection(StudyArea);
}

// create cloud free composite
var noCloud = MODIS.filterDate(iniDate,endDate)
.map(maskQuality)
.map(maskQuality1)
.select(modisBands,lsBands);

// **************** Step2: Water multispectral indices ********************//
var Indices = function(image) {
  // Compute the normalized spectral indices: MNDWI, LSWI, NDPI using an expression.
  var mndwi = image.normalizedDifference(['green', 'swir2']).rename('MNDWI');
  var lswi = image.normalizedDifference(['nir', 'swir2']).rename('LSWI');
  var ndpi = image.normalizedDifference(['swir2', 'green']).rename('NDPI');
  var ndmi = image.normalizedDifference(['nir', 'swir1']).rename('NDMI');
  
  // Compute the AWEI using an expression.
    var awei = image.expression(
    '4*(green - swir1) - (nir * 0.25 + swir2 * 2.75)',{
      nir: image.select('nir'),
      green: image.select('green'),
      swir1: image.select('swir1'),
      swir2: image.select('swir2')
    }).rename('AWEI');
    
  // Compute the TCBI using an expression.
  var tcbi = image.expression(
    'red*0.3956 + nir * 0.4718 + blue * 0.3354 + green * 0.3834 + swir1 * 0.3946 + swir2 * 0.3434+ swir3 * 0.2964',{
      red: image.select('red'),
      nir: image.select('nir'),
      blue: image.select('blue'),
      green: image.select('green'),
      swir1: image.select('swir1'),
      swir2: image.select('swir2'),
      swir3: image.select('swir3')
    }).rename('TCBI');
    
  // Compute the CWI using an expression.
  var cwi = image.expression(
    '((nir-red)/(nir+red) + swir3 + 0.4) * 100',{
      red: image.select('red'),
      nir: image.select('nir'),
      swir3: image.select('swir3')
    }).rename('CWI');
  return image.addBands(mndwi).addBands(lswi).addBands(ndpi).addBands(ndmi).addBands(awei).addBands(tcbi).addBands(cwi);
};

var MODIS_Indices = noCloud.map(Indices);

// ********************** Step2: Export results ***********************//

// Export MNDWI time series
var MNDWI_EXPORT = MODIS_Indices.select('MNDWI').map(function(image) {
  return image.select('MNDWI').reduceRegions({
    collection: ee.FeatureCollection(StudyArea), 
    reducer: ee.Reducer.mean().setOutputs(['MNDWI']), 
    scale: 500,
  })// reduceRegion doesn't return any output if the image doesn't intersect
    // with the point or if the image is masked out due to cloud
    // If there was no ndvi value found, we set the ndvi to a NoData value -9999
    .map(function(feature) {
    var MNDWI = ee.List([feature.get('MNDWI'), -9999])
      .reduce(ee.Reducer.firstNonNull())
    return feature.set({'MNDWI': MNDWI, 'imageID': image.id()})
    })
  }).flatten();
var format = function(table, rowId, colId) {
  var rows = table.distinct(rowId); 
  var joined = ee.Join.saveAll('matches').apply({
    primary: rows, 
    secondary: table, 
    condition: ee.Filter.equals({
      leftField: rowId, 
      rightField: rowId
    })
  });
         
  return joined.map(function(row) {
      var values = ee.List(row.get('matches'))
        .map(function(feature) {
          feature = ee.Feature(feature);
          return [feature.get(colId), feature.get('MNDWI')];
        });
      return row.select([rowId]).set(ee.Dictionary(values.flatten()));
    });
};
var sentinelResults = format(MNDWI_EXPORT, 'imageID', 'LakeName');

Export.table.toDrive({
  collection: sentinelResults,
  folder: 'MinorThesis',
  description: 'MNDWI',
  fileFormat: 'CSV'
});

// Export LSWI time series
var LSWI_EXPORT = MODIS_Indices.select('LSWI').map(function(image) {
  return image.select('LSWI').reduceRegions({
    collection: ee.FeatureCollection(StudyArea), 
    reducer: ee.Reducer.mean().setOutputs(['LSWI']), 
    scale: 500,
  })// reduceRegion doesn't return any output if the image doesn't intersect
    // with the point or if the image is masked out due to cloud
    // If there was no ndvi value found, we set the ndvi to a NoData value -9999
    .map(function(feature) {
    var LSWI = ee.List([feature.get('LSWI'), -9999])
      .reduce(ee.Reducer.firstNonNull());
    return feature.set({'LSWI': LSWI, 'imageID': image.id()});
    });
  }).flatten();
var format = function(table, rowId, colId) {
  var rows = table.distinct(rowId); 
  var joined = ee.Join.saveAll('matches').apply({
    primary: rows, 
    secondary: table, 
    condition: ee.Filter.equals({
      leftField: rowId, 
      rightField: rowId
    })
  });
         
  return joined.map(function(row) {
      var values = ee.List(row.get('matches'))
        .map(function(feature) {
          feature = ee.Feature(feature);
          return [feature.get(colId), feature.get('LSWI')];
        });
      return row.select([rowId]).set(ee.Dictionary(values.flatten()));
    });
};
var sentinelResults = format(LSWI_EXPORT, 'imageID', 'LakeName');

Export.table.toDrive({
  collection: sentinelResults,
  folder: 'MinorThesis',
  description: 'LSWI',
  fileFormat: 'CSV'
});

// Export AWEI time series
var AWEI_EXPORT = MODIS_Indices.select('AWEI').map(function(image) {
  return image.select('AWEI').reduceRegions({
    collection: ee.FeatureCollection(StudyArea), 
    reducer: ee.Reducer.mean().setOutputs(['AWEI']), 
    scale: 500,
  })// reduceRegion doesn't return any output if the image doesn't intersect
    // with the point or if the image is masked out due to cloud
    // If there was no ndvi value found, we set the ndvi to a NoData value -9999
    .map(function(feature) {
    var AWEI = ee.List([feature.get('AWEI'), -9999])
      .reduce(ee.Reducer.firstNonNull());
    return feature.set({'AWEI': AWEI, 'imageID': image.id()});
    });
  }).flatten();
var format = function(table, rowId, colId) {
  var rows = table.distinct(rowId); 
  var joined = ee.Join.saveAll('matches').apply({
    primary: rows, 
    secondary: table, 
    condition: ee.Filter.equals({
      leftField: rowId, 
      rightField: rowId
    })
  });
         
  return joined.map(function(row) {
      var values = ee.List(row.get('matches'))
        .map(function(feature) {
          feature = ee.Feature(feature);
          return [feature.get(colId), feature.get('AWEI')];
        });
      return row.select([rowId]).set(ee.Dictionary(values.flatten()));
    });
};
var sentinelResults = format(AWEI_EXPORT, 'imageID', 'LakeName');

Export.table.toDrive({
  collection: sentinelResults,
  folder: 'MinorThesis',
  description: 'AWEI',
  fileFormat: 'CSV'
});

// Export NDMI time series
var NDMI_EXPORT = MODIS_Indices.select('NDMI').map(function(image) {
  return image.select('NDMI').reduceRegions({
    collection: ee.FeatureCollection(StudyArea), 
    reducer: ee.Reducer.mean().setOutputs(['NDMI']), 
    scale: 500,
  })// reduceRegion doesn't return any output if the image doesn't intersect
    // with the point or if the image is masked out due to cloud
    // If there was no ndvi value found, we set the ndvi to a NoData value -9999
    .map(function(feature) {
    var NDMI = ee.List([feature.get('NDMI'), -9999])
      .reduce(ee.Reducer.firstNonNull());
    return feature.set({'NDMI': NDMI, 'imageID': image.id()});
    });
  }).flatten();
var format = function(table, rowId, colId) {
  var rows = table.distinct(rowId); 
  var joined = ee.Join.saveAll('matches').apply({
    primary: rows, 
    secondary: table, 
    condition: ee.Filter.equals({
      leftField: rowId, 
      rightField: rowId
    })
  });
         
  return joined.map(function(row) {
      var values = ee.List(row.get('matches'))
        .map(function(feature) {
          feature = ee.Feature(feature);
          return [feature.get(colId), feature.get('NDMI')];
        });
      return row.select([rowId]).set(ee.Dictionary(values.flatten()));
    });
};
var sentinelResults = format(NDMI_EXPORT, 'imageID', 'LakeName');

Export.table.toDrive({
  collection: sentinelResults,
  folder: 'MinorThesis',
  description: 'NDMI',
  fileFormat: 'CSV'
});


// Export CWI time series
var CWI_EXPORT = MODIS_Indices.select('CWI').map(function(image) {
  return image.select('CWI').reduceRegions({
    collection: ee.FeatureCollection(StudyArea), 
    reducer: ee.Reducer.mean().setOutputs(['CWI']), 
    scale: 500,
  })// reduceRegion doesn't return any output if the image doesn't intersect
    // with the point or if the image is masked out due to cloud
    // If there was no ndvi value found, we set the ndvi to a NoData value -9999
    .map(function(feature) {
    var CWI = ee.List([feature.get('CWI'), -9999])
      .reduce(ee.Reducer.firstNonNull());
    return feature.set({'CWI': CWI, 'imageID': image.id()});
    });
  }).flatten();
var format = function(table, rowId, colId) {
  var rows = table.distinct(rowId); 
  var joined = ee.Join.saveAll('matches').apply({
    primary: rows, 
    secondary: table, 
    condition: ee.Filter.equals({
      leftField: rowId, 
      rightField: rowId
    })
  });
         
  return joined.map(function(row) {
      var values = ee.List(row.get('matches'))
        .map(function(feature) {
          feature = ee.Feature(feature);
          return [feature.get(colId), feature.get('CWI')];
        });
      return row.select([rowId]).set(ee.Dictionary(values.flatten()));
    });
};
var sentinelResults = format(CWI_EXPORT, 'imageID', 'LakeName');

Export.table.toDrive({
  collection: sentinelResults,
  folder: 'MinorThesis',
  description: 'CWI',
  fileFormat: 'CSV'
});

// Export TCBI time series
var TCBI_EXPORT = MODIS_Indices.select('TCBI').map(function(image) {
  return image.select('TCBI').reduceRegions({
    collection: ee.FeatureCollection(StudyArea), 
    reducer: ee.Reducer.mean().setOutputs(['TCBI']), 
    scale: 500,
  })// reduceRegion doesn't return any output if the image doesn't intersect
    // with the point or if the image is masked out due to cloud
    // If there was no ndvi value found, we set the ndvi to a NoData value -9999
    .map(function(feature) {
    var TCBI = ee.List([feature.get('TCBI'), -9999])
      .reduce(ee.Reducer.firstNonNull());
    return feature.set({'TCBI': TCBI, 'imageID': image.id()});
    });
  }).flatten();
var format = function(table, rowId, colId) {
  var rows = table.distinct(rowId); 
  var joined = ee.Join.saveAll('matches').apply({
    primary: rows, 
    secondary: table, 
    condition: ee.Filter.equals({
      leftField: rowId, 
      rightField: rowId
    })
  });
         
  return joined.map(function(row) {
      var values = ee.List(row.get('matches'))
        .map(function(feature) {
          feature = ee.Feature(feature);
          return [feature.get(colId), feature.get('TCBI')];
        });
      return row.select([rowId]).set(ee.Dictionary(values.flatten()));
    });
};
var sentinelResults = format(TCBI_EXPORT, 'imageID', 'LakeName');

Export.table.toDrive({
  collection: sentinelResults,
  folder: 'MinorThesis',
  description: 'TCBI',
  fileFormat: 'CSV'
});