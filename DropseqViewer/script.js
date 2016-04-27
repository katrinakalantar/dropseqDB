var app = angular.module('myapp', []);

app.service('sharedProperties', function ($q) {

    var self = this, defer = $q.defer();
    var property = 'First';
    
    return{
      getProperty: function(){
        return property;
      },
      observeProperty: function(){
        return defer.promise;
      },
      setProperty: function(value){
        property = value;
      }
    }
    
  
      /*
        var property = 'First';

        return {
            getProperty: function () {
                return property;
            },
            setProperty: function(value) {
                property = value;
            }
        };
        */
    });

//function MyCtrl($scope, sharedProperties) {
app.controller('MyCtrl', function($scope, sharedProperties){
    $scope.foo = null;
    $scope.doSomething = function () {
        sharedProperties.setProperty($scope.foo)
        console.log(sharedProperties.getProperty())
    }
});
//}


app.directive('summary', function(){
  return {
    restrict: 'A',
    replace: false,
    scope: {
      summary: '='
    },
    templateUrl: 'summary.html',
    link: function(scope, elem, attrs){
      console.log('directive: call to load pie1');
      scope.summary();
    }
  }
});



app.controller('SummaryController', function($scope, sharedProperties){
  console.log('controller glued');
  
    
  $scope.reloadTSNE = function () {
        drawTSNE('#TSNE',sharedProperties.getProperty()); 
        console.log("reload TSNE for "  + sharedProperties.getProperty())
    }
  
var drawTSNE = function(where, geneID) {
  
  if (where === '#TSNE') {
      d3.select(where).html('');
    }

var tooltip = d3.select("body")
  .append("div")
  .style("position", "absolute")
  .style("z-index", "10")
  .attr("class", "hoverinfo")
  .style("visibility", "hidden");

var color_vector=["#00b300","#008ae6","#000066","#ff0000","#993d00", "#330066","#800000"]
var color = d3.scale.linear()
    .domain([0, 5])
    .range(["blue","red"]);
var opacity = d3.scale.linear()
    .domain([0,5])
    .range([0,1])
var point_size = d3.scale.linear()
    .domain([0, 15])
    .range([5,8]);
var color_classif = d3.scale.linear()
    .domain([0,1])
    .range(["black","yellow"]);
var margin = {top: 20, right: 20, bottom: 30, left: 40},
    width = 760 - margin.left - margin.right,
    height = 400 - margin.top - margin.bottom;
var x = d3.scale.linear()
    .range([0, width]);
var y = d3.scale.linear()
    .range([height, 0]);
var svg = d3.select("#TSNE").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .call(d3.behavior.zoom().on("zoom", function () {
    svg.attr("transform", "translate(" + d3.event.translate + ")" + " scale(" + d3.event.scale + ")")
  }))
  .append("g")
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");
d3.tsv("lupus_sm.tsv", function(error, data) {
  if (error) throw error;
  // Coerce the data to numbers.
  data.forEach(function(d) {
    d.x = +d.tSNEx;
    d.y = +d.tSNEy;
    d.variance = +d.variance;
    d.classif = +d.classif;
    d.PDPN = +d.PDPN;
    a = geneID;
    d.gid = d[a]
    //d.CellID = +d.CellID;
  });
  // Compute the scales’ domains.
  x.domain(d3.extent(data, function(d) { return d.x; })).nice();
  y.domain(d3.extent(data, function(d) { return d.y; })).nice();
  // Add the x-axis.
  svg.append("g")
      .attr("class", "x axis")
      .attr("transform", "translate(0," + height + ")")
      .call(d3.svg.axis().scale(x).orient("bottom"));
  // Add the y-axis.
  svg.append("g")
      .attr("class", "y axis")
      .call(d3.svg.axis().scale(y).orient("left"));
  // Add the points!
  svg.selectAll(".point")
      .data(data)
    .enter().append("circle")
      .attr("class", "point")
      .style("stroke", function(d){return color_vector[d.classif]; })
      .style("stroke-width",1)
      .style("fill-opacity", function(d) {return opacity(d.gid);})
      .style("stroke-opacity", 1)
      .style("fill", function(d){return color_vector[d.classif];})
      .attr("r", function(d) {return point_size(d.variance); })
      .attr("cx", function(d) { return x(d.x); })
      .attr("cy", function(d) { return y(d.y); })
      //.on("mouseover", function(){return tooltip.style("visibility", "visible").text(data, function(d){return d.CellID;});})
      .on("mouseover", function(d,i){
        d3.select(this).attr("fill", "#F77E1C")
        tooltip.html(d.CellID+" : "+d.classif+" : " + d.variance);
        return tooltip.style("visibility", "visible");
    })
      .on("mousemove", function(){return tooltip.style("top",
        (d3.event.pageY-10)+"px").style("left",(d3.event.pageX+10)+"px");})
      .on("mouseout", function(){return tooltip.style("visibility", "hidden");})
      .on('click', function(d){
        var cellID=d.CellID;
        var classif=d.classif;
        console.log("hello - cell: " + d.cellID);
        console.log("variance: " + d.variance);
        console.log("classif: " + d.classif)
        $scope.visualizeCell(cellID, classif);
      });
      
      
      
      var legend = [0, 1, 2, 3, 4, 5, 6];

    svg.append("g")
       .attr("transform", "translate(600,10)")
       .selectAll("rect")
       .data(legend)
       .enter()
       .append("rect")
       .attr("x", function(d, i) {
           return i * 12;  //Bar width of 12 plus 1 for padding
       })
       .attr("y", function(d) {
           return 20;
       })
      .attr("width", 10)
      .attr("fill",function(d){
          return color_vector[d];
      })
      .attr("height", 10);
      
      svg.append("g")
         .attr("transform", "translate(600,10)")
         .selectAll("text")
            .data(legend)
            .enter()
            .append("text")
            .text(function(d) {
                    return d;
               })
           .attr("x", function(d, i) {
               return i * 12 + 5;

           })
           .attr("y", function(d, i) {
                return 40;
           })
           .attr("font-family", "sans-serif")
           .attr("font-size", "8px")
           .attr("fill", "black")
           .attr("text-anchor", "middle");      
           
       svg.append("text").attr("transform", "translate(600,15)")
       .text("Classification:")
       .attr("font-family", "sans-serif")
       .attr("font-size", "12px");             
      
      
      
  });

};



var drawCellHM = function(where, cellID, classif){
  
  if (where === '#cellHM') {
      d3.select(where).html('');
    }

console.log("inside drawCellHM method - " + cellID)
d3.text("lupus_sm.tsv", function(datasetText){

var parsedCSV = d3.tsv.parseRows(datasetText);
console.log("inside drawCellHM method - " + cellID)
var cellOI = cellID;
console.log("cellID: " + cellOI)
var classOI = classif;

// add tooltip
var tooltip = d3.select("body")
  .append("div")
  .style("position", "absolute")
  .style("z-index", "10")
  .attr("class", "hoverinfo")
  .style("visibility", "hidden");

// parse out x and y axis labels
var xlab = parsedCSV[0].slice(1,parsedCSV[0].length-4);
var ylab = ["cell_OI: " + cellOI, "avg expression within class", "overall avg expression"];

// parse the datamatrix
var data_matrix = [];
var single_cell = [];
var sum_total = [];
var sum_within_class = [];
var count_total = 0;
var count_within = 0;
for(var i = 1; i < parsedCSV.length;i++){//-4; i++){ // parsedCSV.length; i++){
    var temp = parsedCSV[i].slice(1,parsedCSV[i].length-4);
    if(parsedCSV[i][0] == cellOI){
      single_cell = temp;
      console.log("cellid is equal to expected")
      console.log("expected class: " + classOI)
      console.log("class : " + parsedCSV[i][parsedCSV[0].length-4])
    }
    console.log(parsedCSV[i][parsedCSV[0].length-4] )
    if (parsedCSV[i][parsedCSV[0].length-4] == classOI){
        console.log("classIO = class")
        sum_within_class = sum_within_class.concat([temp]);
        count_within = count_within + 1;
    } 
    sum_total = sum_total.concat([temp]);
    count_total = count_total + 1;
}
console.log("count_within: " + count_within)
console.log("count_total: " + count_total)

swc = [];
sumtot=[];
for(var l = 0; l < sum_within_class[0].length; l++){
  sum = 0;
  console.log("sum_within_class.length: " + sum_within_class.length)
  for(var j = 0; j< sum_within_class.length; j++){
    sum = sum + parseFloat(sum_within_class[j][l]);
  }
  swc.push((sum/count_within));
}
for(var l = 0; l < sum_total[0].length; l++){
  sum_tot= 0;
  for(var j = 0; j< sum_total.length; j++){
    sum_tot= sum_tot + parseFloat(sum_total[j][l]);
  }
  sumtot.push((sum_tot/count_total));
}
data_matrix = [single_cell, swc, sumtot];



// find marginal sums
var sum_rows = []; // sum of the rows
for(var i = 0; i < data_matrix.length; i++){
    var temp_sum = 0;
    for(var j = 0; j < data_matrix[i].length; j++){
        temp_sum = parseFloat(temp_sum) + parseFloat(data_matrix[i][j]);
    }
    sum_rows = sum_rows.concat([temp_sum]);
}

var sum_col = []; // sum of the columns
for(var i = 0; i < data_matrix[0].length; i++){
    var temp_sum = 0;
    for(var j = 0; j < data_matrix.length; j++){
        temp_sum = parseFloat(temp_sum) + parseFloat(data_matrix[j][i]);
    }
    sum_col = sum_col.concat([temp_sum]);
}                    

function colorPicker(d) {
	if( d <= 0 ){
		return "rgb(255, 247, 251)"; // 0
	} else if (d > 0 && d<= 2){
		return "rgb(236, 226, 240)"; // 1
	} else if (d > 2 && d<= 4) {
		return "rgb(208, 209, 230)"; // 2
	} else if (d > 4 && d<= 6) {
		return "rgb(166, 189, 219)"; // 3
	} else if (d > 6 && d<= 8){
		return "rgb(103, 169, 207)"; // 4
	} else if (d > 8 && d<= 10){
		return "rgb(54, 144, 192)"; // 5
	} else if (d > 10 && d<= 12){
		return "rgb(2, 129, 138)"; // 6
	} else if (d > 12 && d<= 14) {
		return "rgb(1, 108, 89)"; // 7
	} else if (d > 14){
		return "rgb(1, 70, 54)"; // 7+
	} else {
		return "pink"; // number not in range
	}
}


var margin = {top: 160, right: 160, bottom: 80, left: 160},
    width = data_matrix[1].length*12;
    height = data_matrix.length*12;                    
    
var svg = d3.select("#cellHM").append("svg");
    svg.attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .style("border-collapse", "collapse")
    .append("g")
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")")
    .selectAll("g")
    .data(data_matrix)
    .enter()
    .append("g")
    .attr("class", "row")
    .attr("transform", function(d, i){ return "translate(0," + (i*12) + ")"; })
    
    .selectAll("rect")
    .data(function(d){return d;})
    .enter().append("rect")
    .attr("width", "10px")
    .attr("height", "10px")
    .attr("fill", function(d) {
        return colorPicker(d);
    })
    .attr("x", function(d,i) {
        return (i * 12);
    })
    .on("mouseover", function(d,i){
        d3.select(this).attr("fill", "#F77E1C")
        tooltip.html(d);
        return tooltip.style("visibility", "visible");
    })
    .on("mouseout", function(){d3.select(this).attr("fill", function(d) {
        tooltip.style("visibility", "hidden");
        return colorPicker(d);
    })})
    .on("mousemove", function(){
        return tooltip.style("top",(d3.event.pageY+10)+"px").style("left",(d3.event.pageX+10)+"px");
        })
    .on("mousedown", function(){d3.select(this).attr("fill", "red")})
    .on("mouseup", function(){d3.select(this).attr("fill", "#F77E1C")})
    .on('click', function(d,i){
        $scope.visualizeGeneDistribution(xlab[i]);//d.geneID);
        //if (where === '#gene1') {
        //  $scope.visualizeCell(cellID);  
        //}
      })
    //.attr("onclick", "javascript:showBlast('MMGA_0001')");                        
    .attr("alt", function(d) {
        return d;
    })
    // axis
    var xScale = d3.scale.ordinal()
                          .domain(xlab)
                          .rangePoints([0+5, width-5]);
    var xAxis = d3.svg.axis()
                        .scale(xScale)
                        .orient("left")
                        
    var yScale = d3.scale.ordinal()
                         .domain(ylab)
                         .rangePoints([0+4, height-4]);
    var yAxis = d3.svg.axis()
                       .scale(yScale)
                       .orient("left")                    
    svg.append("g")
        .attr("class", "axis")
        .attr("transform", "translate("+ margin.left + "," + (margin.top + height) + ") rotate(270)")
        .call(xAxis);
        
    svg.append("g")
        .attr("class", "axis")
        .attr("transform", "translate("+ (margin.left -2) + "," + margin.top + ") rotate(0)")
        .call(yAxis);

    // 
    
    var col_hist_size = 0.5* margin.top;

    var col_hist_scale = d3.scale.linear()
                         .domain([0, d3.max(sum_col)])
                         .range([0, col_hist_size]);
        
    svg.append("g")
       //.attr("transform", "translate(" + margin.left + "," + margin.top + ")")
       .attr("transform", "translate(" + margin.left + ",0)")
       .selectAll("rect")
       .data(sum_col)
       .enter()
       .append("rect")
       .attr("x", function(d, i) {
           return i * 12;  //Bar width of 12 plus 1 for padding
       })
       .attr("y", function(d) {
              return margin.top - col_hist_scale(d) - 5; // spacing of 5
        })
       .attr("width", 10)
       .attr("fill","teal")
       .attr("height", function(d) {
           return col_hist_scale(d);
       });

      /*
       svg.append("g")
          .attr("transform", "translate(" + margin.left + ",0)")
          .selectAll("text")
          .data(sum_col)
          .enter()
          .append("text")
          .text(function(d) {
                  return d;
             })
         .attr("x", function(d, i) {
              return i * 12 + 2 + 3;
         })
         .attr("y", function(d) {
              return margin.top - col_hist_scale(d) - 5 - 3;
         })
         .attr("font-family", "sans-serif")
         .attr("font-size", "8px")
         .attr("fill", "black")
         .attr("text-anchor", "middle");*/
         
     svg.append("g").append("line")
        .attr("x1", margin.left).attr("x2", width + margin.right - 2).attr("y1", margin.top - 3).attr("y2", margin.top - 3).attr("stroke","black");
        
        
        // row histogram
        
        var row_hist_size = 0.8* margin.right;

        var row_hist_scale = d3.scale.linear()
                             .domain([0, d3.max(sum_rows)])
                             .range([0, row_hist_size]);
        
        
        svg.append("g")
        .attr("transform", "translate(" + (width + margin.right) + "," + margin.top + ")")
        .selectAll("rect")
        .data(sum_rows)
        .enter()
        .append("rect")
        .attr("x", function(d, i) {
            return 0; // spacing of 5
        })
        .attr("y", function(d, i) {
            return i * 12;  //Bar width of 12 plus 1 for padding
              
        })
        .attr("width",  function(d) {
               return row_hist_scale(d);
        })
        .attr("fill","rgb(54, 144, 192)")
        .attr("height", 10);
        
        svg.append("g")
           .attr("transform", "translate(" + (width + margin.right) + "," + margin.top + ")")
           .selectAll("text")
           .data(sum_rows)
           .enter()
           .append("text")
           .text(function(d) {
                   return d;//+"xxx";  //MODIFIED HERE TO SEE IF THIS IS THE GENE NAMES
              })
          .attr("x", function(d, i) {
              return row_hist_scale(d) + 5;

          })
          .attr("y", function(d, i) {
               return i * 12 + 8;
          })
          .attr("font-family", "sans-serif")
          .attr("font-size", "8px")
          .attr("fill", "black");
    
    var legend = [0, 2, 4, 6, 8, 10, 12, 14, 16];

    svg.append("g")
       .attr("transform", "translate(20,60)")
       .selectAll("rect")
       .data(legend)
       .enter()
       .append("rect")
       .attr("x", function(d, i) {
           return i * 12;  //Bar width of 12 plus 1 for padding
       })
       .attr("y", function(d) {
           return 20;
       })
      .attr("width", 10)
      .attr("fill",function(d){
          return colorPicker(d);
      })
      .attr("height", 10);
      
      svg.append("g")
         .attr("transform", "translate(20,60)")
         .selectAll("text")
            .data(legend)
            .enter()
            .append("text")
            .text(function(d) {
                    return d;
               })
           .attr("x", function(d, i) {
               return i * 12 + 5;

           })
           .attr("y", function(d, i) {
                return 40;
           })
           .attr("font-family", "sans-serif")
           .attr("font-size", "8px")
           .attr("fill", "black")
           .attr("text-anchor", "middle");      
           
       svg.append("text").attr("transform", "translate(20,65)")
       .text("Legend:")
       .attr("font-family", "sans-serif")
       .attr("font-size", "12px");                                       
   console.log("end of cell dispplay method")
                 
});
  
  
};







var drawGeneDist = function(where, geneID){
  
  if (where === '#geneDist') {
      d3.select(where).html('');
    }

console.log("inside drawgeneDistr method - " + geneID)

function parser(d) {
    d.classif = +d.classif
    d.gDFFB = +d.DFFB;
    d.gHES2 = +d.HES2;
    d.gH6PD = +d.H6PD;
    var a = geneID
    console.log("a: " + a)
    d.gene = d[a]
    //console.log("d.HES1: " + d.HES2)
    //var a = "HES2"
    //console.log("d.a: " + d.a)
    console.log("d[a]: " + d[a])
    //d.gid = +d.geneID;
    return d;
}

function genehist(tsvdata) {
    console.log(tsvdata)
    var binsize = 2;
    var minbin = 0;
    var maxbin = 100;
    var numbins = (maxbin - minbin) / binsize;

    // whitespace on either side of the bars in units of MPG
    var binmargin = .05; 
    var margin = {top: 10, right: 30, bottom: 50, left: 60};
    var width = 450 - margin.left - margin.right;
    var height = 250 - margin.top - margin.bottom;

    // Set the limits of the x axis
    var xmin = minbin - 1
    var xmax = maxbin + 1

    histdata = new Array(numbins);
    for (var i = 0; i < numbins; i++) {
		histdata[i] = { numfill: 0, meta: "" };
	}

	// Fill histdata with y-axis values and meta data
    tsvdata.forEach(function(d) {
		var bin = Math.floor((d.gene - minbin) / binsize);
		if ((bin.toString() != "NaN") && (bin < histdata.length)) {
			histdata[bin].numfill += 1;
		}
    });
    
    console.log(histdata);

    // This scale is for determining the widths of the histogram bars
    // Must start at 0 or else x(binsize a.k.a dx) will be negative
    var x = d3.scale.linear()
	  .domain([0, (xmax - xmin)])
	  .range([0, width]);

    // Scale for the placement of the bars
    var x2 = d3.scale.linear()
	  .domain([xmin, xmax])
	  .range([0, width]);
	
    var y = d3.scale.linear()
	  .domain([0, d3.max(histdata, function(d) { 
						return d.numfill; 
						})])
	  .range([height, 0]);

    var xAxis = d3.svg.axis()
	  .scale(x2)
	  .orient("bottom");
    var yAxis = d3.svg.axis()
	  .scale(y)
	  .ticks(8)
	  .orient("left");

/*
    var tip = d3.tip()
	  .attr('class', 'd3-tip')
	  .direction('e')
	  .offset([0, 20])
	  .html(function(d) {
	    return '<table id="tiptable">' + d.meta + "</table>";
	});
*/

var tooltip = d3.select("body")
  .append("div")
  .style("position", "absolute")
  .style("z-index", "10")
  .attr("class", "hoverinfo")
  .style("visibility", "hidden");

    // put the graph in the "mpg" div
    var svg = d3.select("#geneDist").append("svg")
	  .attr("width", width + margin.left + margin.right)
	  .attr("height", height + margin.top + margin.bottom)
	  .append("g")
	  .attr("transform", "translate(" + margin.left + "," + 
						margin.top + ")");
						
						
		svg.append("text")
        .attr("x", (width / 2))             
        .attr("y", 0 - (margin.top / 10))
        .attr("text-anchor", "middle")  
        .style("font-size", "10px")  
        .text("Distribution of Gene: " + geneID);

    //svg.call(tooltip);

    // set up the bars
    var bar = svg.selectAll(".bar")
	  .data(histdata)
	  .enter().append("g")
	  .attr("class", "bar")
	  .attr("transform", function(d, i) { return "translate(" + 
	       x2(i * binsize + minbin) + "," + y(d.numfill) + ")"; })
	  .on("mouseover", function(d,i){
        d3.select(this).attr("fill", "#F77E1C")
        tooltip.html(d.CellID);
        return tooltip.style("visibility", "visible");
      })
    .on("mousemove", function(){return tooltip.style("top",
        (d3.event.pageY-10)+"px").style("left",(d3.event.pageX+10)+"px");})
    .on("mouseout", function(){d3.select(this).attr("fill", function(d) {
        tooltip.style("visibility", "hidden");
        return colorPicker(d);
    })});
	  //.on('mouseover', tooltip.show)
	  //.on('mouseout', tooltip.hide);

    // add rectangles of correct size at correct location
    bar.append("rect")
	  .attr("x", x(binmargin))
	  .attr("width", x(binsize - 2 * binmargin))
	  .attr("height", function(d) { return height - y(d.numfill); });

    // add the x axis and x-label
    svg.append("g")
	  .attr("class", "x axis")
	  .attr("transform", "translate(0," + height + ")")
	  .call(xAxis);
    svg.append("text")
	  .attr("class", "xlabel")
	  .attr("text-anchor", "middle")
	  .attr("x", width / 2)
	  .attr("y", height + margin.bottom)
	  .text("expression level")
	  .attr("font-family", "sans-serif");

    // add the y axis and y-label
    svg.append("g")
	  .attr("class", "y axis")
	  .attr("transform", "translate(0,0)")
	  .call(yAxis);
    svg.append("text")
	  .attr("class", "ylabel")
	  .attr("y", 0 - margin.left) // x and y switched due to rotation
	  .attr("x", 0 - (height / 2))
	  .attr("dy", "1em")
	  .attr("transform", "rotate(-90)")
	  .style("text-anchor", "middle")
	  .text("cells")
	  .attr("font-family", "sans-serif");
}

// Read in .csv data and make graph
d3.tsv("lupus_sm.tsv",parser,//"lupus_Sm.tsv", parser,
       function(error, tsvdata) {
	   genehist(tsvdata);
}); 


}; //end of visualizeGeneDistr method

  $scope.visualizeTSNE = function(){//sharedProperties){
    console.log('conroller: will load genes1');
    drawTSNE('#TSNE',sharedProperties.getProperty());  
  };

  $scope.visualizeCell = function(cellID, classif) {
    console.log('loading cell for ' + cellID);
    drawCellHM('#cellHM', cellID, classif);
  };
  
  $scope.visualizeGeneDistribution = function(geneID) {
    console.log('loading gene distribution for ' + geneID);
    drawGeneDist('#geneDist', geneID);
  };


});