// Exected URL parameters:
// ee_script : one of "sentinelAndLandsat", ""sentinelAndLandsatWithFires"
// startTime : Defautls to 2000-01-01
// endTime : Defaults to todays date
// geoJson : This is the only parameter that MUST be filled
// plotId : The name of the plot being assessed, defaults to Plot

var ee_script = ui.url.get('ee_script', "sentinelAndLandsat");
var script;


switch(ee_script) {
  case 'sentinelAndLandsat':
    script = require("users/TEST/ceo_dev:main");
    break;
  default:
    script = require("users/TEST/ceo_dev:main");
}

// The date that is used as the start of the chart ( if the dataset is available )
// You can change the start date manually and hit the button "Run""again to reload the charts using the different time series
var startTime = ui.url.get('startTime', "2000-01-01");

var now = new Date( Date.now() );
var today = now.getFullYear()+'-'+ ( now.getMonth() +1) +'-'+now.getDate();
var lastYearToday = ( now.getFullYear() -1 )+'-'+ ( now.getMonth() +1) +'-'+now.getDate();

// The last date for which the chart is generated.
var endTime = ui.url.get('endTime', today);

var plotId = ui.url.get('plotId', 'Plot');
// var geoJson = ui.url.get('geoJson');
// var obj = JSON.parse(geoJson);
var geometry = geo//ee.Geometry( obj );

script.processPlotInfo( geometry, startTime, endTime, lastYearToday, null, plotId  );