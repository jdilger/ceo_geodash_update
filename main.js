//////////////////////////////////////////////////////////////////////////
// UTILS
//////////////////////////////////////////////////////////////////////////
 
var visualization = {
  landsat7 : {min: 0.007, max: 0.11, gamma: [1.3, 1.3, 1.3] },
  landsat5 : {min: 0, max: 0.325, gamma: [1.3, 1.3, 1.3]} ,
  sentinel : {min: 0, max: 0.3 } ,
  sentinelFalse_20: { bands :['B8','B11','B4'], min: 0, max: 3000},
  sentinelFalse: { bands :['B8','B11','B4'], min: 0, max: 3000},
  landsat8 : {min: 0.007, max: 0.11, gamma: [1.3, 1.3, 1.3]},
  //landsat8 : {min: 0, max: 0.3, gamma: [1.2, 1.2, 1.2]},
  asterFalse : { bands: ["nir_median","swir1_median",  "red_median" ], min:0, max:0.5},
  vizParams2000_gfcFirst : {"opacity":1,"bands":["first_b40","first_b50","first_b30"],"min":7.04,"max":104.96,"gamma":1},
  vizParams2018_gfcLast : {"opacity":1,"bands":["last_b40","last_b50","last_b30"],"min":10.72,"max":126.28,"gamma":1},

  
};

function panSharpenLandsat(landsatCollectionName, bandsToSharpen, bandsCloudScore, startPeriod, endPeriod ){

	// Load a Landsat collection. (GREENEST PIXEL, NO NEED TO CLEAN CLOUDS!! )
	var image= ee.ImageCollection( landsatCollectionName )
	// Select the bands of interest to avoid taking up memory.
	//.select(bandsCloudScore)
	// Filter to get only six months of data.
	.filterDate( startPeriod , endPeriod )
	.mosaic();
	image = ee.Image(image );


	// Convert the RGB bands to the HSV color space.
	var hsv = image.select( bandsToSharpen ).rgbToHsv();

	// Swap in the panchromatic band and convert back to RGB.
	// The pan band is B8 for both Landsat 7 and 8
	var sharpened = ee.Image.cat([
		hsv.select('hue'), hsv.select('saturation'), image.select('B8')
	]).hsvToRgb();

	return sharpened;
}


function cloudMask(im) {
  // Opaque and cirrus cloud masks cause bits 10 and 11 in QA60 to be set,
  // so values less than 1024 are cloud-free
  var mask = ee.Image(0).where(im.select('QA60').gte(1024), 1).not();
  return im.updateMask(mask);
}

function sharpenSentinel(image){
	// sharpen see e.g. http://www.cse.psu.edu/~rtc12/CSE486/lecture11_6pp.pdf
	var log = image
    .convolve(ee.Kernel.gaussian(10, 7, 'meters')) // G
    .convolve(ee.Kernel.laplacian8(0.5)); // L of G

	var sharpened = image.subtract(log);
	return sharpened;
}

// From the GEE forum. Thread : https://groups.google.com/forum/#!searchin/google-earth-engine-developers/sentinel$20cloud|sort:date/google-earth-engine-developers/i63DS-Dg8Sg/LLk9ZM0mAgAJ
function sentinelCloudScore(toa) {
  // authors: Matt Hancher, Chris Hewig and Ian Housman
  
  function rescale(img, thresholds) { 
    return img.subtract(thresholds[0]).divide(thresholds[1] - thresholds[0]);
  }
  
  // Compute several indicators of cloudyness and take the minimum of them.
  var score = ee.Image(1);
  
  //Clouds are reasonably bright
  score = score.min(rescale(toa.select(['blue']), [0.1, 0.5]));
  score = score.min(rescale(toa.select(['aerosol']), [0.1, 0.3]));
  score = score.min(rescale(toa.select(['aerosol']).add(toa.select(['cirrus'])), [0.15, 0.2]));
  score = score.min(rescale(toa.select(['red']).add(toa.select(['green'])).add(toa.select('blue')), [0.2, 0.8]));

  //Clouds are moist
  var ndmi = toa.normalizedDifference(['red4','swir1']);
  score=score.min(rescale(ndmi, [-0.1, 0.1]));
  
  // However, clouds are not snow.
  var ndsi = toa.normalizedDifference(['green', 'swir1']);
  score=score.min(rescale(ndsi, [0.8, 0.6]));
  
  // a (somewhat arbitrary) threshold 
  var cloudScoreThreshold = 0.5;
  var cloud = score.gt(cloudScoreThreshold);
  
  return cloud;
}


function ESAcloud(toa) {
  // author: Nick Clinton
  
  var qa = toa.select('QA60');
  
  // Bits 10 and 11 are clouds and cirrus, respectively.
  var cloudBitMask = Math.pow(2, 10);
  var cirrusBitMask = Math.pow(2, 11);
  
  // clear if both flags set to zero.
  var clear = qa.bitwiseAnd(cloudBitMask).eq(0).and(
             qa.bitwiseAnd(cirrusBitMask).eq(0));
  
  var cloud = clear.eq(0);
  
  return cloud;
}

// do we want this? I think we could use smething else
function shadowMask(toa,cloud){
  // Author: Gennadii Donchyts
  // License: Apache 2.0
  
  // solar geometry (radians)
  var azimuth =ee.Number(toa.get('solar_azimuth')).multiply(Math.PI).divide(180.0).add(ee.Number(0.5).multiply(Math.PI));
  var zenith  =ee.Number(0.5).multiply(Math.PI ).subtract(ee.Number(toa.get('solar_zenith')).multiply(Math.PI).divide(180.0));

  // find where cloud shadows should be based on solar geometry
  var nominalScale = cloud.projection().nominalScale();
  var cloudHeights = ee.List.sequence(200,10000,500);
  var shadows = cloudHeights.map(function(cloudHeight){
    cloudHeight = ee.Number(cloudHeight);
    var shadowVector = zenith.tan().multiply(cloudHeight);
    var x = azimuth.cos().multiply(shadowVector).divide(nominalScale).round();
    var y = azimuth.sin().multiply(shadowVector).divide(nominalScale).round();
    return cloud.changeProj(cloud.projection(), cloud.projection().translate(x, y));
  });
  var potentialShadow = ee.ImageCollection.fromImages(shadows).max();
  
  // shadows are not clouds
  var potentialShadow1 = potentialShadow.and(cloud.not());
  
  // ['B1','B2','B3','B4','B6','B8','B8A','B9','B10', 'B11','B12']
  // ['aerosol', 'blue', 'green', 'red', 'red2','nir','red4','h2o', 'cirrus','swir1', 'swir2']
                      
  // (modified by Sam Murphy) dark pixel detection 
  var darkPixels = toa.normalizedDifference(['green', 'swir2']).gt(0.25).rename(['dark_pixels']);
  
  // shadows are dark
  var shadow = potentialShadow1.and(darkPixels).rename('shadows');
  
  return shadow;
}

function imgRename( img){
  return img.select(
                      ['B1','B2','B3','B4','B6','B8','B8A','B9','B10', 'B11','B12','QA60']
                      ,['aerosol', 'blue', 'green', 'red', 'red2','nir','red4','h2o', 'cirrus','swir1', 'swir2','QA60']
                    )
                    //.divide(10000)
                    ;
}

function sentinel2toa(img) {
  img = imgRename( img );
  img = img
          .set('solar_azimuth',img.get('MEAN_SOLAR_AZIMUTH_ANGLE'))
          .set('solar_zenith',img.get('MEAN_SOLAR_ZENITH_ANGLE') )
          .set('system:time_start',img.get('system:time_start'));
  return img;
}

function cloud_and_shadow_mask(img) {
  
  var toa = sentinel2toa(img);
  
  var cloud = ESAcloud(toa);
  // Thi method "sentinelCloudScore" gives a better result but is slower
  // var cloud = sentinelCloudScore(toa);
  
  var shadow = shadowMask(toa,cloud);
  
  var mask = cloud.or(shadow).eq(0);
  
  return toa.updateMask(mask);
  
}
//////////////////////////////////////////////////////////////////////////
// SENTINEL-2
//////////////////////////////////////////////////////////////////////////
/***
 * Sentinel-2 produces multiple images, resultsing sometimes 4x more images than the actual size.
 * This is bad for any statistical analysis.
 * 
 * by Genna
 *
 * This function mosaics images by time.
 */

function mosaicByTime (images) {
  var TIME_FIELD = 'system:time_start';

  var distinct = images.distinct([TIME_FIELD]);

  var filter = ee.Filter.equals({ leftField: TIME_FIELD, rightField: TIME_FIELD });
  var results = ee.Join.saveAll('matches').apply(distinct, images, filter);

  // mosaic
  results = results.map(
    function(i) {
      var mosaic = ee.ImageCollection.fromImages(i.get('matches'))
          .sort('system:index')
          .mosaic();
  
      return mosaic.copyProperties(i).set(TIME_FIELD, i.get(TIME_FIELD));
    }
  );

  return ee.ImageCollection(results);
};

function getSentinel2CloudFreeImage ( geometry, startDate, endDate ){
  var images =  getSentinel2Images(startDate, endDate, geometry, true );
  var out = images.median().divide(10000);
  return out;
}

function getSentinel2Images (start, end, geometry, mask ){
      var sentinelImages = ee.ImageCollection('COPERNICUS/S2')
      .filterBounds(geometry)
      .filterMetadata( "CLOUDY_PIXEL_PERCENTAGE", 'less_than', 70 ) // Filter out tiles with too much tree cover
      .filterDate(start, end);
      
      if( mask ){
        sentinelImages = sentinelImages.map(cloud_and_shadow_mask);
      }
      
      return sentinelImages;
};

function getImageTimeStart( geometry, timeStart ){
    var equalDate = ee.Filter.equals('system:time_start', timeStart);
    return ee.Image(
      ee.ImageCollection('COPERNICUS/S2')
      .filterBounds(geometry).filter(equalDate).first());
}
//////////////////////////////////////////////////////////////////////////
// OTHER INDICES
//////////////////////////////////////////////////////////////////////////
// function getSaviIndices( imageCollection, region, nirBand, redBand ){
  
//   print( imageCollection.bandNames() );  
  
//   // Calculate the SAVI index for each of the pixels within the clipped images
//   var imageCollectionWithSAVI = 
//     imageCollection.map(
//       function(image){
        
//       if( 
//         image.bandNames().contains( redBand) 
//         ){
//         var saviImage = ee.Image(0).expression(
//             '(1 + L) * float(nir - red)/ (nir + red + L)',
//             {
//               'nir': image.select(nirBand),
//               'red': image.select(redBand),
//               'L': 0.5
//             }).select([0],['SAVI']);
//           return saviImage;
       
//       }else{
//         return image;
//       }
//     }
      
//   );

//   return  imageCollectionWithSAVI;
// }

// function getOsaviIndices( imageCollection, region, nirBand, redBand ){
    
//   // Calculate the OSAVI index for each of the pixels within the clipped images
//   var imageCollectionWithOSAVI = 
//     imageCollection.map(
//       function(image){
      
//         var osaviImage = ee.Image(0).expression(
//           '(1 + L) * float(nir - red)/ (nir + red + L)',
//           {
//             'nir': image.select(nirBand),
//             'red': image.select(redBand),
//             'L': 0.16
//           }).select([0],['OSAVI']);
          
       
//         return ee.Image.cat(image,osaviImage);
//       }
      
//   );

//   return  imageCollectionWithOSAVI;
// }


// function getMsaviIndices( imageCollection, region, nirBand, redBand ){
    
//   // Calculate the MSAVI index for each of the pixels within the clipped images
//   var imageCollectionWithMSAVI = 
//     imageCollection.map(
//       function(image){
      
//         var msaviImage = ee.Image(0).expression(
//           '( 2*nir + 1 - sqrt( pow( (2*nir+1), 2 ) - 8*(nir - red) ) )/2',
//           {
//             'nir': image.select(nirBand), 
//             'red': image.select(redBand)
//           }).select([0],['MSAVI']);
          
       
//         return ee.Image.cat(image,msaviImage);
//       }
      
//   );

//   return  imageCollectionWithMSAVI;
// }

//////////////////////////////////////////////////////////////////////////
// ASTER
//////////////////////////////////////////////////////////////////////////
/***
* Aster TOA reflectance and temperature from DN. This example is for the ASTER VNIR subsystem (i.e. bands 1,2,3N)
* 
* Algorithm: 
*  http://www.pancroma.com/downloads/ASTER%20Temperature%20and%20Reflectance.pdf
*  https://lpdaac.usgs.gov/sites/default/files/public/product_documentation/aster_l1t_users_guide.pdf
* 
* Authors: 
*  Gennadii Donchyts (GD), gennadiy.donchyts@gmail.com
*  Sam Murphy (SM), samsammurphy@gmail.com - DN > TOA
* 
* Changelog:
*  GD 2016-06: DN > temperature
*  SM 2016-06: DN > TOA
*  GD 2017-06: added swir1
* 
* License: MIT, https://opensource.org/licenses/MIT
*/

/***
* Rescales to given ranges
*/
function rescale (img, exp, thresholds) {
  return img.expression(exp, {img: img}).subtract(thresholds[0]).divide(thresholds[1] - thresholds[0]);
}

/*** 
* Aster radiometric correction algorithms
*/
var Aster = {
  radiance: {
    fromDN: function(image) {
      // Gain coefficients are dynamic (i.e. can be high, normal, low_1 or low_2)
      var multiplier = ee.Image([
        ee.Number(image.get('GAIN_COEFFICIENT_B01')),
        ee.Number(image.get('GAIN_COEFFICIENT_B02')),
        ee.Number(image.get('GAIN_COEFFICIENT_B3N')),
        ee.Number(image.get('GAIN_COEFFICIENT_B04'))
        ]).float();
      
      // Apply correction
      var radiance = image.select(['B01', 'B02', 'B3N', 'B04'], ['green','red','nir','swir1']).subtract(1).multiply(multiplier);
      
      // Define properties required for reflectance calculation
      var solar_z = ee.Number(90).subtract(image.get('SOLAR_ELEVATION'));
      
      return radiance.set({
        'system:time_start': image.get('system:time_start'),
        'solar_zenith':solar_z
      });
    }
  },
  
  reflectance: {
    fromRad: function(rad) {
      // calculate day of year from time stamp
      var date = ee.Date(rad.get('system:time_start'));
      var jan01 = ee.Date.fromYMD(date.get('year'),1,1);
      var doy = date.difference(jan01,'day').add(1);

      // Earth-Sun distance squared (d2) 
      var d = ee.Number(doy).subtract(4).multiply(0.017202).cos().multiply(-0.01672).add(1); // http://physics.stackexchange.com/questions/177949/earth-sun-distance-on-a-given-day-of-the-year
      var d2 = d.multiply(d);
      
      // mean exoatmospheric solar irradiance (ESUN)
      var ESUN = [1847, 1553, 1118, 232.5]; // from Thome et al (A) see http://www.pancroma.com/downloads/ASTER%20Temperature%20and%20Reflectance.pdf
      
      // cosine of solar zenith angle (cosz)
      var solar_z = ee.Number(rad.get('solar_zenith'));
      var cosz = solar_z.multiply(Math.PI).divide(180).cos();

      // calculate reflectance
      var scalarFactors = ee.Number(Math.PI).multiply(d2).divide(cosz);
      var scalarApplied = rad.multiply(scalarFactors);
      var reflectance = scalarApplied.divide(ESUN);
      
      return reflectance;
    }
  },
  
  temperature: {
    fromDN: function(image) {
      var bands = ['B10', 'B11', 'B12', 'B13', 'B14'];
      var multiplier = ee.Image([0.006822, 0.006780, 0.006590, 0.005693, 0.005225]);
      var k1 = ee.Image([3040.136402, 2482.375199, 1935.060183, 866.468575, 641.326517]);
      var k2 = ee.Image([1735.337945, 1666.398761, 1585.420044, 1350.069147, 1271.221673]);
  
      var radiance = image.select(bands).subtract(1).multiply(multiplier);
      var t = k2.divide(k1.divide(radiance).add(1).log()).rename(bands);
      
      return t;
    }
  },
  
  cloudScore: function(image) {
    // Compute several indicators of cloudyness and take the minimum of them.
    var score = ee.Image(1.0);

    // Snow is reasonably bright in all visible bands.
    score = score.min(rescale(image, 'img.red + img.green', [0.2, 0.8]));

    // Excluded this for snow reasonably bright in all infrared bands.
    score = score.min(rescale(image, 'img.nir + img.swir1', [0.2, 0.4]));

    // Clouds are reasonably cool in temperature.
    score = score.min(rescale(image.resample('bicubic'), '(img.B10 + img.B12 + img.B14) / 3.0', [293, 280]));

    // However, clouds are not snow.
    //let ndsi = img.normalizedDifference(['red', 'swir']);
    //score = score.min(rescale(ndsi, 'img', [0.8, 0.6])).aside(show, 'score ndsi')

    return score;
  },
  
  TOA: function(image) {
    var radiance = Aster.radiance.fromDN(image);
    var reflectance = Aster.reflectance.fromRad(radiance);
    var temperature = Aster.temperature.fromDN(image);

    var result = reflectance.addBands(temperature);
    result = result.set('system:time_start', image.get('system:time_start'));
    result = result.copyProperties(image);
    result = ee.Image(result);
    
    return result;
  }
};


function getAsterForYear(year){
  // ASTER image collection with VNIR filter in the ROI
  var aster = ee.ImageCollection('ASTER/AST_L1T_003').filterDate(year+'-01-01', year+'-12-31').filter(ee.Filter.and(
      ee.Filter.listContains('ORIGINAL_BANDS_PRESENT', 'B01'),
      ee.Filter.listContains('ORIGINAL_BANDS_PRESENT', 'B02'),
      ee.Filter.listContains('ORIGINAL_BANDS_PRESENT', 'B3N'),
      ee.Filter.listContains('ORIGINAL_BANDS_PRESENT', 'B04')
      ));
  var refcol = aster.map(Aster.TOA);
  var median = refcol.reduce(ee.Reducer.median());
  return median;
}

//////////////////////////////////////////////////////////////////////////
// UI
//////////////////////////////////////////////////////////////////////////

var PROPERTY_IMAGE_NAME = "image_name"; 

function addCheckBox( map, visualization, images, name){
  
  var onChangeFunction = function( checked ){
      var image;
      if( checked ){
       image = images[1];
      }else{
        image = images[0];
      }
      
      image = image.visualize(visualization );
      map.layers().set( 0, image);
  };
    
  var checkBox = ui.Checkbox({
    label : name,
    onChange: onChangeFunction,
    style: {stretch: 'horizontal'}
  });
  

  map.add(checkBox);
  
  return map;
}

// function addCheckBoxMult( map, visualizationArray, images, name){
  
//   var onChangeFunction = function( checked ){
//       var image;
//       if( checked ){
//       image = images[1];
//       image = image.visualize(visualizationArray[1] );
//       }else{
//         image = images[0];
//         image = image.visualize(visualizationArray[0] );
//       }
      
      
//       map.layers().set( 0, image);
//   };
    
//   var checkBox = ui.Checkbox({
//     label : name,
//     onChange: onChangeFunction,
//     style: {stretch: 'horizontal'}
//   });
//   onChangeFunction( false );
  
//   map.add(checkBox);
  
//   return map;
// }

// function getLabelImage(image){
//   return ee.Algorithms.If(
//       image.get("system:time_start"),
//       ee.Date(image.get("system:time_start") ).format( "dd, MMMM, yyyy"),
//       image.get(PROPERTY_IMAGE_NAME)
//     );
// }

// function getImageryDates( images ){  

//   var getImageryDate = function(current, aggregate){
//     var label = getLabelImage(current);
//     return ee.Dictionary( { date: label });  
//   };

//   var dateValues = images.iterate(getImageryDate );
//   return dateValues;
// }

function addPlotToMap(map, plot, subplots){
    map.addLayer(plot, {}, "Plot");
    map.centerObject( ee.Geometry(plot) , 15 );
    if( subplots ){
      map.addLayer(subplots, {}, "Subplots");
    }
    // Center Object not working   
 /*   
    var center = ee.Geometry(plot).centroid();
    center.coordinates().evaluate(
      function( coords ){
        map.setCenter( coords[0], coords[1], 15 );
      }
    );
*/

    map.setControlVisibility(false);
    return map;
}

function createL7SliderMap( image, visualization, year, name, plot, subplots){
    var map = ui.Map();
    addL7Slider(map, year);
    map.addLayer(image, visualization, name);
    // Add the plot polygon to the map
    return addPlotToMap(map, plot, subplots);
}

function createL8SliderMap( image, visualization, year, name, plot, subplots){
    var map = ui.Map();
    addL8Slider(map, year);
    map.addLayer(image, visualization, name);
    // Add the plot polygon to the map
    return addPlotToMap(map, plot, subplots);
}

function createL5SliderMap( image, visualization, year, name, plot, subplots){
    var map = ui.Map();
    addL5Slider(map, year);
    map.addLayer(image, visualization, name);
    // Add the plot polygon to the map
    return addPlotToMap(map, plot, subplots);
}

// TODO: see if these need to be global or not
var mapL7,mapL8,mapL5,mapAster;

function showLandsat7ForYearWrapper(year){
  showLandsat7ForYear( mapL7, year);
}

function showLandsat8ForYearWrapper(year){
  showLandsat8ForYear( mapL8, year);
}

function showLandsat5ForYearWrapper( year){
  showLandsat5ForYear( mapL5, year);
}

function showAsterForYearWrapper( year){
  showAsterForYear( mapAster, year);
}


function showLandsat5ForYear(map, year){
  var landsat = (ee.Image) ( ee.ImageCollection( "LANDSAT/LT5_L1T_ANNUAL_GREENEST_TOA" ).select(['B4', 'B5', 'B3']).filterDate( year + "-01-01", year + "-12-31" ).first() );
  landsat = landsat.visualize(visualization.landsat5 );
  map.layers().set( 0, landsat);
}

function showLandsat7ForYear (map, year){
  var landsat;
  if( year == 2000 ){
    landsat = gfc.select(["first_b40","first_b50","first_b30"]).divide(1000);
  }else{
    landsat = panSharpenLandsat('LANDSAT/LE7_L1T_ANNUAL_GREENEST_TOA' , ['B4', 'B5', 'B3'], ['B1', 'B2', 'B3', 'B4', 'B5', 'B6_VCID_1' , 'B6_VCID_2' , 'B7', 'B8'], year + "-01-01", year + "-12-31"  );
  }
  landsat = landsat.visualize(visualization.landsat7 );
  map.layers().set( 0, landsat);
}

function showLandsat8ForYear (map, year){
  var landsat;
  if( year == 2018 ){
    landsat = gfc.select(["last_b40","last_b50","last_b30"]).divide(1000);
  }else{
    landsat =  panSharpenLandsat('LANDSAT/LC8_L1T_ANNUAL_GREENEST_TOA' , ['B5', 'B6', 'B4'], ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B10', 'B11'], year + "-01-01", year + "-12-31"  );
  }
  landsat = landsat.visualize( visualization.landsat8 );
  map.layers().set( 0, landsat);
}

function showAsterForYear(map, year){  
  var asterImage = getAsterForYear(year );
  asterImage = asterImage.visualize( visualization.asterFalse );
  map.layers().set( 0, asterImage);
}


function addL7Slider(mapWithSlider , selectedYear ){
  // Create a label and slider.
  var label = ui.Label('Landsat 7 False Color Yearly Mosaic');
  var onChangeFunction = showLandsat7ForYearWrapper;
  
  mapL7 = mapWithSlider;
  
  var slider = ui.Slider({
    min: 2000,
    max: 2017,
    step: 1,
    onChange: onChangeFunction,
    style: {stretch: 'horizontal'}
  });
  
  // Create a panel that contains both the slider and the label.
  var panelSlider = ui.Panel({
    widgets: [label, slider],
    layout: ui.Panel.Layout.flow('vertical'),
    style: {
      position: 'top-center',
      padding: '7px'
    }
  });
 
  mapWithSlider.add(panelSlider);
  slider.setValue(selectedYear,false);
}

function addL8Slider(mapWithSlider , selectedYear ){
  // Create a label and slider.
  var label = ui.Label('Landsat 8 False Color Yearly Mosaic');
  var onChangeFunction = showLandsat8ForYearWrapper;
  
  mapL8 = mapWithSlider;
  
  var slider = ui.Slider({
    min: 2013,
    max: 2018,
    step: 1,
    onChange: onChangeFunction,
    style: {stretch: 'horizontal'}
  });
  
  // Create a panel that contains both the slider and the label.
  var panelSlider = ui.Panel({
    widgets: [label, slider],
    layout: ui.Panel.Layout.flow('vertical'),
    style: {
      position: 'top-center',
      padding: '7px'
    }
  });
 
  mapWithSlider.add(panelSlider);
  slider.setValue(selectedYear,false);
}

function addL5Slider(mapWithSlider , selectedYear ){
  // Create a label and slider.
  var label = ui.Label('Landsat 5 False Color');
 
  mapL5 = mapWithSlider;
 
  var onChangeFunction = showLandsat5ForYearWrapper;
    
  var slider = ui.Slider({
    min: 1984,
    max: 2000,
    step: 1,
    onChange: onChangeFunction,
    style: {stretch: 'horizontal'}
  });
  
  // Create a panel that contains both the slider and the label.
  var panelSlider = ui.Panel({
    widgets: [label, slider],
    layout: ui.Panel.Layout.flow('vertical'),
    style: {
      position: 'top-center',
      padding: '7px'
    }
  });
 
  mapWithSlider.add(panelSlider);
  slider.setValue(selectedYear,false);
}


//////////////////////////////////////////////////////////////////////////
// CHARTS

// Auxiliary Sentinel Map variable, do not remove
var sentinelAux = 0;

function maskModisQuality(x){
  var Q = x.select(['DetailedQA']);
  var shadow = Q.bitwiseAnd(ee.Image.constant(32768));
  var cloud = Q.bitwiseAnd(ee.Image.constant(1024));
  var aerosol = Q.bitwiseAnd(ee.Image.constant(192));
  var filter = shadow.add(cloud).add(aerosol);
  return filter.lte(64)
}

function getModis16NdviCollection(){
  
  var modisNDVI = ee.ImageCollection('MODIS/006/MOD13Q1');
  
  modisNDVI = modisNDVI.map(function(img){
    var mask = maskModisQuality(img)
    var NDVI = img.select(['NDVI']).multiply(0.0001)
    var masked = NDVI.updateMask(mask)
    return masked.copyProperties(img, ['system:time_start', 'system:time_end'])
  });
  
  return modisNDVI;
}

function getName(name){
  name = name || "Plot";
  return name + " - ";
}

function getModis16DayNDVI(plot, start, end, name ){
  // Load the MODIS  Vegetation Index composite. Select the NDVI band. Resolution of the pixels is 250 meters.
  var modisNoaaNdvi = ee.ImageCollection( getModis16NdviCollection() ).filterDate(start, end).select('NDVI');
    
  var modisNoaaTimeSeries = ui.Chart.image.series(modisNoaaNdvi, plot, ee.Reducer.mean(), 30);
  modisNoaaTimeSeries = modisNoaaTimeSeries
  .setOptions({  
    title: getName(name) + 'MOD13Q1 Modis Vegetation Indices 16-Day Global 250m',
    hAxis: {
      title: 'Date',gridlines: {count: 10,}
    },
    vAxis: {
      title: 'NDVI',viewWindowMode: 'explicit', viewWindow: {max: 1,min: -0.25,},gridlines: {count: 5}
    },
    series: [
      {color: 'green', visibleInLegend: true}
    ]
    ,
    trendlines: {
      1: {
        visibleInLegend: false,
        color: 'blue',
        labelInLegend : "Trend",
        opacity : 0.7
      }
    }
  });
  
  // Show the MODIS NDVI chart on the console
  return modisNoaaTimeSeries;
};

function ndviCalculcation(image){
    var ndvi = image.normalizedDifference(["nir","red"]).rename("NDVI");
    ndvi = ndvi.set("system:time_start", image.get("system:time_start"));
    return ndvi
}

function getMergedLandsat(plot, start, end){
  var LANDSAT_7_THRESHOLD_DATE = "2015-01-01";
  var LANDSAT_THRESHOLD_CLOUDCOVER = 40;

  var landsat7Images = ee.ImageCollection("LANDSAT/LE07/C01/T1_TOA")
          .filterBounds( plot )
          .filterDate(start, LANDSAT_7_THRESHOLD_DATE)
          .filter(ee.Filter.lt('CLOUD_COVER', LANDSAT_THRESHOLD_CLOUDCOVER ))
          .select( [ 'B4', 'B5', 'B3'], ['nir' ,'swir', 'red'] );
  var landsat8Images = ee.ImageCollection("LANDSAT/LC08/C01/T1_TOA")
          .filterBounds( plot )
          .filterDate(LANDSAT_7_THRESHOLD_DATE, end)
          .filter(ee.Filter.lt('CLOUD_COVER', LANDSAT_THRESHOLD_CLOUDCOVER ))
          .select( ['B5','B6','B4'], ['nir','swir','red'] );
  return landsat7Images.merge( landsat8Images );
}

function getImageLandsat( plot, start, end, timeStart ){
    var equalDate = ee.Filter.equals('system:time_start', timeStart);
    return ee.Image(
      getMergedLandsat(plot, start, end)
      .filter(equalDate).first());
}

function getLandsatNDVI(plot, start, end, name, landsat7Map, landsat8Map ){


  var landsatNDVI = getMergedLandsat(plot, start, end).map( ndviCalculcation );
  
  
  var landsatNDVITimeSeries = ui.Chart.image.series(landsatNDVI, plot, ee.Reducer.mean(), 30);
  
  landsatNDVITimeSeries = landsatNDVITimeSeries
  .setOptions({  
    title: getName(name) + 'Landsat 7/8 NDVI - Click on graph to load specific image in visor',
    hAxis: {
      title: 'Date',gridlines: {count: 10,}
    },
    vAxis: {
      title: 'NDVI',viewWindowMode: 'explicit', viewWindow: {max: 1,min: -1,},gridlines: {count: 5}
    },
    series: [
      {color: 'green', visibleInLegend: true}
    ]
    ,
    trendlines: {
      1: {
        visibleInLegend: false,
        color: 'blue',
        labelInLegend : "Trend",
        opacity : 0.7
      }
    }
  });
  
  // When the chart is clicked, update the map and label.
  landsatNDVITimeSeries.onClick(function(xValue, yValue, seriesName) {
    if (!xValue) return;  // Selection was cleared.
  
    // Show the image for the clicked date.
    var imageForDate = getImageLandsat( plot, start, end, xValue );
    
    var landsatLayer = ui.Map.Layer(imageForDate, { bands:['nir' ,'swir', 'red'], min:0.05, max:0.65 , gama:1.3} );
    var landsatSensor = imageForDate.get('LANDSAT_PRODUCT_ID').getInfo() + "";
   
    var mapToChange;
    if( landsatSensor.slice(0,4) === "LC08"  ){
      mapToChange = landsat8Map;
    }else{
      mapToChange = landsat7Map;
    }

    var plotlayer = mapToChange.layers().get(1);
    mapToChange.clear();
    var label = ui.Label((new Date(xValue)).toUTCString());
    mapToChange.add(label);
    mapToChange.layers().reset([landsatLayer, plotlayer]);
  });

  // Show the MODIS NDVI chart on the console
  return landsatNDVITimeSeries;
  
}

function getSentinel2NDVI(plot, start, end, name, sentinelMap, sentinelMap2 ){
  // Load the MODIS  Vegetation Index composite. Select the NDVI band. Resolution of the pixels is 250 meters.
  var s2 =  getSentinel2Images(start, end, plot, true);
  
  s2 =  mosaicByTime(s2);
    // ['B1','B2','B3','B4','B6','B8','B8A','B9','B10', 'B11','B12']
  // ['aerosol', 'blue', 'green', 'red', 'red2','nir','red4','h2o', 'cirrus','swir1', 'swir2']
  s2 = s2.map( ndviCalculcation ); 
   
  var sentinel2TimeSeries = ui.Chart.image.series(s2, plot, ee.Reducer.mean(), 20);
  sentinel2TimeSeries = sentinel2TimeSeries
  .setOptions({  
    title: getName(name) + 'Sentinel-2 NDVI - Click on graph to load specific image in visor',
    hAxis: {
      title: 'Date',gridlines: {count: 10,}
    },
    vAxis: {
      title: 'NDVI',viewWindowMode: 'explicit', viewWindow: {max: 1,min: -1,},gridlines: {count: 5}
    },
    series: [
      {color: 'green', visibleInLegend: true}
    ]
    ,
    trendlines: {
      1: {
        visibleInLegend: false,
        color: 'blue',
        labelInLegend : "Trend",
        opacity : 0.7
      }
    }
  });
  
  // When the chart is clicked, update the map and label.
  sentinel2TimeSeries.onClick(function(xValue, yValue, seriesName) {
    if (!xValue) return;  // Selection was cleared.
  
    // Show the image for the clicked date.
    var imageForDate =  getImageTimeStart( plot, xValue );
    
    var sentinelLayer = ui.Map.Layer(imageForDate,  visualization.sentinelFalse );
    
    var mapToChange = sentinelMap;
    if( sentinelMap2 ){
      sentinelAux++;
      if( sentinelAux % 2 === 0 ){ 
        mapToChange = sentinelMap2;
      }
    }
    
    //sentinelMap.add(sentinelLayer);
    
    var plotlayer = mapToChange.layers().get(1);
    mapToChange.clear();
    var label = ui.Label((new Date(xValue)).toUTCString());
    mapToChange.add(label);
    mapToChange.layers().reset([sentinelLayer, plotlayer]);
    // Show a label with the date on the map.
    //label.setValue((new Date(xValue)).toUTCString());
  });

  // Show the MODIS NDVI chart on the console
  return sentinel2TimeSeries;
  

}

//////////////////////////////////////////////////////////////////////////
// EXPORTS
//////////////////////////////////////////////////////////////////////////
exports = { 

          // getLandsat7MonthlyNDVI : getLandsat7MonthlyNDVI,
          // getLandsat8MonthlyNDVI : getLandsat8MonthlyNDVI,
          getSentinel2NDVI : getSentinel2NDVI,
          getModis16DayNDVI : getModis16DayNDVI,
          getLandsatNDVI : getLandsatNDVI,
          createL7SliderMap : createL7SliderMap,
          createL8SliderMap : createL8SliderMap,
          createL5SliderMap : createL5SliderMap,
          addPlotToMap : addPlotToMap,
          PROPERTY_IMAGE_NAME : PROPERTY_IMAGE_NAME,
          getSentinel2CloudFreeImage : getSentinel2CloudFreeImage,
          sentinel2toa : sentinel2toa,
          mosaicByTime : mosaicByTime,
          getSentinel2Images : getSentinel2Images,
          sharpenSentinel : sharpenSentinel,
          cloudMask : cloudMask,
          visualization : visualization
          };
