var main = require("users/TEST/ceo_dev:main");
// var charts = require("users/collectearth/ce_scripts:common/charts_v2");
// var utils = require("users/collectearth/ce_scripts:common/utils");
// var ui_methods = require("users/collectearth/ce_scripts:common/ui");
// var s2 = require("users/TEST/ceo_dev:sentinel2");
var gfc = ee.Image("UMD/hansen/global_forest_change_2018_v1_6");

 

  
// Script automatically produced by Collect Earth for the plot that has been clicked on Google Earht. See bottom of the script for more info on customization.
// This script will show graphs of NDVI, EVI and NDWI mean-values for the pixels contained within the plot.
var addCharts = function( plot, start, end, name, sentinelMap, landsat7Map, landsat8Map, panel ){
    panel.widgets().add( main.getModis16DayNDVI( plot, start, end, name) );
    panel.widgets().add( main.getLandsatNDVI( plot, start, end, name, landsat7Map, landsat8Map  ) );
    panel.widgets().add( main.getSentinel2NDVI( plot, start, end, name,sentinelMap ) );
};

var processPlotInfo = function( plot, start, end, lastYearsDate, subplots, plotId ){
    var startYear = "-01-01"; 
    var endYear = "-12-31";
    var lastYearsDateObject = new Date( lastYearsDate );
    var lastYearStart = lastYearsDateObject.getFullYear() + startYear;
    var lastYearEnd = lastYearsDateObject.getFullYear() +  endYear;
    // var geo = plot.buffer( 50000 );
    //Array sorting is useful for obtaining custom quality mosaics which involve reducing a subset of image bands according to the values in a different band. The following example sorts by a cloud index, then gets the mean of the least cloudy subset of images in the collection:
    // var sharpenedLandsat8_false =  main.panSharpenLandsat('LANDSAT/LC8_L1T_ANNUAL_GREENEST_TOA' , ['B5', 'B6', 'B4'], ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B10', 'B11'], '2015'+startYear, '2015' + endYear, geo);
    // var sharpenedLandsat7_2001 = utils.panSharpenLandsat('LANDSAT/LE7_L1T_ANNUAL_GREENEST_TOA' , ['B4', 'B5', 'B3'], ['B1', 'B2', 'B3', 'B4', 'B5', 'B6_VCID_1' , 'B6_VCID_2' , 'B7', 'B8'], '2001'+startYear, '2001' + endYear);
   
    var gfcFirst_2000 = gfc.select(["first_b40","first_b50","first_b30"]).divide(1000);
    var gfcLast_2018 = gfc.select(["last_b40","last_b50","last_b30"]).divide(1000);
   
    // Create a map for each visualization option.
    var maps = [];
    
    var landsat8Map =  main.createL8SliderMap(
      gfcLast_2018, 
       main.visualization.landsat8,
      2018,
      "Landsat 8 2018 GFC mosaic", 
      plot, 
      subplots);
    
    var landsat7Map =  main.createL7SliderMap(
      gfcFirst_2000, 
       main.visualization.landsat7, 
      2000, 
      "Landsat 7 2000 GFC mosaic", 
      plot, 
      subplots);

    // get s2 collection, no fancy processing, also not used?
    // var sentinelImages =  main.getSentinel2Images(start,end,plot);
    // gets median of sd and ed imagery, 
    var sentinelImage12MonthsComposite =  main.getSentinel2CloudFreeImage( plot, lastYearsDate, end );
    //sharpens image, shouldnt be to computationally expensive
    sentinelImage12MonthsComposite =  main.sharpenSentinel( sentinelImage12MonthsComposite)
      .select(['nir','red', 'swir1'])
      .rename(['B8','B4','B11']).multiply( ee.Image(10000));

    var labelS2 = "Composite "+  lastYearsDate +" - "+ end;

    sentinelImage12MonthsComposite = ee.Image(sentinelImage12MonthsComposite).set(
       main.PROPERTY_IMAGE_NAME,
      labelS2
    );
    
    var label = ui.Label('Sentinel 2 : Composite of last 12 months. To select single image click on Sentinel NDVI chart');
    
    var sentinelMap = ui.Map();
    sentinelMap.add(label);
    sentinelMap.addLayer(sentinelImage12MonthsComposite.visualize( main.visualization.sentinelFalse));
    sentinelMap =  main.addPlotToMap(sentinelMap, plot, subplots);
    
    
    maps.push(sentinelMap);
    maps.push(landsat8Map);
    maps.push(landsat7Map);
    
    var linker = ui.Map.Linker(maps);
    // Create a grid of maps.
    
    var panelLeft =  ui.Panel([sentinelMap], null, {stretch: 'both'});
    var panelRight =  ui.Panel([landsat8Map, landsat7Map], null, {stretch: 'both'});
    
    var mapGrid = ui.Panel([ panelLeft,panelRight],
      ui.Panel.Layout.Flow('horizontal'), {stretch: 'both'}
    );

    // Enable zooming on the top-left map.
    sentinelMap.setControlVisibility({zoomControl: true, fullscreenControl:true});
    
    // Show the scale (e.g. '500m') on the bottom-right map.
    landsat7Map.setControlVisibility({scaleControl: true, fullscreenControl:true});
    
    landsat8Map.setControlVisibility({ fullscreenControl:true});
    
    
    // Create a panel to hold title, intro text, chart and legend components.
    var inspectorPanel = ui.Panel({style: {width: '30%'}});
    addCharts( plot, start, end , plotId, sentinelMap,  landsat7Map , landsat8Map, inspectorPanel  );
    
    // Add the maps and title to the ui.root.
    ui.root.clear();

    ui.root.add(ui.SplitPanel(mapGrid, inspectorPanel));
    ui.root.setLayout(ui.Panel.Layout.Flow('vertical'));
    
    Map.centerObject( plot );
};

exports = { processPlotInfo:processPlotInfo};



// dev

var startTime = ui.url.get('startTime', "2011-01-01");
var now = new Date( Date.now() );
var today = now.getFullYear()+'-'+ ( now.getMonth() +1) +'-'+now.getDate();
var lastYearToday = ( now.getFullYear() -1 )+'-'+ ( now.getMonth() +1) +'-'+now.getDate();
var endTime = ui.url.get('endTime', today);
var plotId = ui.url.get('plotId', 'Plot');
var geometry = geo//ee.Geometry( obj );

processPlotInfo( geometry, startTime, endTime, lastYearToday, null, plotId  );