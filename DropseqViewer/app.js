var app = angular.module('myapp', []);

app.service('sharedProperties', function ($q) {

    var self = this, defer = $q.defer();
    var property = 'First';
    var analysisDir = 'C:\\'
    var analysisName = null
    
    return{
      getProperty: function(){
        return property;
      },
      getAnalysisDir: function(){
        return analysisDir;
      },
      getAnalysisName: function(){
        return analysisName;
      },
      observeProperty: function(){
        return defer.promise;
      },
      setProperty: function(value){
        property = value;
      },
      setAnalysisDir: function(value){
      	analysisDir = value;
      },
      setAnalysisName: function(value){
      	analysisName = value;
      },
    }
});

app.controller('MyCtrl', function($scope, sharedProperties){
    $scope.foo = null;
    $scope.doSomething = function () {
        sharedProperties.setProperty($scope.foo)
        //console.log(sharedProperties.getProperty())
    }
    $scope.foo2 = null;
    $scope.sendAnalysisDir = function () {
        sharedProperties.setAnalysisDir($scope.foo2)
        //console.log(sharedProperties.getAnalysisDir())
    }
    $scope.foo3 = null;
    $scope.sendAnalysisName = function () {
        sharedProperties.setAnalysisName($scope.foo3)
        //console.log(sharedProperties.getAnalysisName())
    }
});

app.directive('summary', function(){
  return {
    restrict: 'A',
    replace: false,
    scope: {
      summary: '='
    },
    templateUrl: 'summary.html',
    link: function(scope, elem, attrs){
      scope.summary();
    }
  }
});



app.controller('SummaryController', function($scope,$http, sharedProperties){
  console.log('controller glued');

//var exec = require("child_process").exec
   
  $scope.reloadTSNE = function () {
        drawTSNE('#TSNE',sharedProperties.getProperty(),sharedProperties.getAnalysisDir()+"\\d3data.tsv"); 
        //PIP('#TSNE');
        console.log("reload TSNE for "  + sharedProperties.getProperty())
    }

  var selected = []

  $scope.writeCells = function(){
  	//console.log("WRITING CELLS")
  	//console.log(selected)
  	var filename = sharedProperties.getAnalysisDir() + "\\\\" + sharedProperties.getAnalysisName()+".txt"
  	//console.log(filename)
  	var blob = new Blob([selected], {type: "text/plain;charset=utf-8"});
  	saveAs(blob, filename);
  }

  $scope.upregulate = function(){
  	exec("C:\\Users\\KATRINA\\Anaconda3\\python.exe \\scipts\\test.py", $output);
      /*$http({
         url: "/scripts/test.py",
         success: function(response){
           console.log("inside success function")
         }
      });*/
  }

    //var TSNE_file = sharedProperties.getAnalysisDir()+"\\lupus_sm.tsv";
  
	var drawTSNE = function(where, geneIDlist, TSNE_file) {
  
  
		if (where === '#TSNE') {
			d3.select(where).html('');
		}


		var coords = [];//,

	    var line = d3.svg.line();//,

	      // Set the behavior for each part
	      // of the drag.
	    var drag = d3.behavior.drag()
	                  .on("dragstart", function() {
	                    // Empty the coords array.
	                    coords = [];
	                    svg = d3.select(this);

	                    // If a selection line already exists,
	                    // remove it.
	                    svg.select(".selection").remove();

	                    // Add a new selection line.
	                    svg.append("path").attr({"class": "selection"});
	                  })
	                  .on("drag", function() {
	                    // Store the mouse's current position
	                    coords.push(d3.mouse(this));

	                    svg = d3.select(this);

	                    // Change the path of the selection line
	                    // to represent the area where the mouse
	                    // has been dragged.
	                    svg.select(".selection").attr({
	                      d: line(coords)
	                    });

	                    // Figure out which dots are inside the
	                    // drawn path and highlight them.
	                    selected = [];
	                    svg.selectAll("circle.dot").each(function(d, i) {
	                      point = [d3.select(this).attr("cx"), d3.select(this).attr("cy")];
	                      if (pointInPolygon(point, coords)) {
	                      	//console.log(d.CellID)
	                        selected.push(d.CellID);//id);
	                      }
	                    });
	                    highlight(selected);
	                    //console.log(selected)
	                  })
	                  .on("dragend", function() {
	                    svg = d3.select(this);
	                    // If the user clicks without having
	                    // drawn a path, remove any paths
	                    // that were drawn previously.
	                    if (coords.length === 0) {
	                      d3.selectAll("svg path").remove();
	                      unhighlight();
	                      return;
	                    }

	                    // Draw a path between the first point
	                    // and the last point, to close the path.
	                    svg.append("path").attr({
	                      "class": "terminator",
	                      d: line([coords[0], coords[coords.length-1]])
	                    });
	                  });



		var tooltip = d3.select("body")
		  .append("div")
		  .style("position", "absolute")
		  .style("z-index", "10")
		  .attr("class", "hoverinfo")
		  .style("visibility", "hidden");

		//var color_vector=["#00b300","#008ae6","#000066","#ff0000","#993d00", "#330066","#800000"]
		//var color_vector=["#ff5050","#52d4ff","#52ff52","#ffa852","#d452ff","#5252ff","#ffff52"]
		var color_vector=["#ffff52","#005ce6","#ff9900","#cc0066","#33ccff","#00e600","#a64dff","#ff0040","#ff0040","#ff0000","#ff0000","#0099ff","#0099ff","#0099ff"]
		var color = d3.scale.linear()
			.domain([0, 5])
			.range(["blue","red"]);
		var number_of_genes = 1;
		var opacity = d3.scale.linear()
			.domain([0,8000*parseInt(number_of_genes)])   ///THIS IS NOT WORKING (I DON"T THINK)...I don't think the domain is updating
			.range([0,1])
		var point_size = d3.scale.linear()//log()
			.domain([1, 10000000])
			.range([2,8])
		var color_classif = d3.scale.linear()
			.domain([0,1])
			.range(["black","yellow"]);
		var margin = {top: 20, right: 20, bottom: 30, left: 40},
			width = 960 - margin.left - margin.right,
			height = 500 - margin.top - margin.bottom;
		var x = d3.scale.linear()
			.range([0, width]);
		var y = d3.scale.linear()
			.range([height, 0]);
		var svg = d3.select("#TSNE").append("svg")
			.attr("width", width + margin.left + margin.right)
			.attr("height", height + margin.top + margin.bottom)
			.call(drag)
			/*.call(d3.behavior.zoom().on("zoom", function () {
			svg.attr("transform", "translate(" + d3.event.translate + ")" + " scale(" + d3.event.scale + ")")
		  }))*/
		  .append("g")
			.attr("transform", "translate(" + margin.left + "," + margin.top + ")");
		d3.tsv(TSNE_file, function(error,data){//"lupus_sm.tsv", function(error, data) {
		  if (error) throw error;
		  // Coerce the data to numbers.
		  data.forEach(function(d) {
			d.x = +d.tSNEx;
			d.y = +d.tSNEy;
			d.variance = +d.variance;
			d.classif = +d.classif;
			d.PDPN = +d.PDPN;
			//a = geneID;
			//d.gid = d[a]
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
			  .attr("class", "dot")
			  .style("stroke", function(d){return color_vector[d.classif/100]; })
			  .style("stroke-width",1)
			  .style("fill-opacity", function(d) {
			  	var GL = geneIDlist.split(',');
			  	var totalExpn = 0;
			  	var genes_count = 0
			  	for (var i = 0, j = GL.length - 1; i < GL.length; j = i++) {
			  		a = GL[i]
			  		if (d[a] != null){
			  			//console.log(d[a])
			  			totalExpn = totalExpn + parseInt(d[a])
			  			genes_count = genes_count + 1
			  		}
			  	}
			  	number_of_genes = genes_count
			  	return opacity(totalExpn);//d.gid);
			  })    ////MAKE THE GENE ID LIST UPDATE HERE.
			  .style("stroke-opacity", 1)
			  .style("fill", function(d){return color_vector[d.classif/100];})
			  .attr("r", function(d) {return point_size(d.variance); }) //{return 5;})// 
			  .attr("cx", function(d) { return x(d.x); })
			  .attr("cy", function(d) { return y(d.y); })
			  //.on("mouseover", function(){return tooltip.style("visibility", "visible").text(data, function(d){return d.CellID;});})
			  .on("mouseover", function(d,i){
				d3.select(this).attr("fill", "#F77E1C")
				tooltip.html(d.CellID);
				return tooltip.style("visibility", "visible");
			})
			  .on("mousemove", function(){return tooltip.style("top",
				(d3.event.pageY-10)+"px").style("left",(d3.event.pageX+10)+"px");})
			  //.on("mouseout", function(){return tooltip.style("visibility", "hidden");})
			  .on('click', function(d){
				var cellID=d.CellID;
				var classif=d.classif;
				//$scope.visualizeCell(cellID, classif);
			  })
			.on("mouseover", function(d, i) {
                    // Highlight circles on mouseover, but
                    // only if a path hasn't been drawn.
                    if (d3.selectAll("svg path").empty()) {
                      highlight([d.CellID]);
                    }

				d3.select(this).attr("fill", "#F77E1C")
				tooltip.html(d.CellID);//+" : "+d.classif+" : " + d.variance+" : " + d.x+" : " + d.y);
				return tooltip.style("visibility", "visible");

                  })
             .on("mouseout", function(d, i) {
                    // If a path hasn't been drawn,
                    // unhighlight the highlighted circles.
                    if (d3.selectAll("svg path").empty()) {
                      unhighlight();
                    }
                    return tooltip.style("visibility", "hidden");
             })
			  //);
			  
			  
			  /*
			  var legend = [0, 1, 2, 3, 4, 5, 6];

			svg.append("g")
			   .attr("transform", "translate(800,10)")
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
				 .attr("transform", "translate(800,10)")
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
				   .attr("text-anchor", "middle")
				   .attr({"class": "g-dot"});      
				   
				   
			   svg.append("text").attr("transform", "translate(800,15)")
			   .text("Classification:")
			   .attr("font-family", "sans-serif")
			   .attr("font-size", "12px");             
			  */
			  
			  
		  });



// from https://github.com/substack/point-in-polygon
  function pointInPolygon (point, vs) {
    var xi, xj, i, intersect,
        x = point[0],
        y = point[1],
        inside = false;
    for (var i = 0, j = vs.length - 1; i < vs.length; j = i++) {
      xi = vs[i][0]-margin.left,  //the SVG is off by a margin - must account for that in path-checking or else the points will be off from selection
      yi = vs[i][1]-margin.top,
      xj = vs[j][0]-margin.left,
      yj = vs[j][1]-margin.top,
      intersect = ((yi > y) != (yj > y))
          && (x < (xj - xi) * (y - yi) / (yj - yi) + xi);
      if (intersect) inside = !inside;
    }
    return inside;
  }

  function unhighlight() {
    d3.selectAll("circle.dot").classed("highlighted", false);
  }

  function highlight(ids) {
    // First unhighlight all the circles.
    unhighlight();

    // Find the circles that have an id
    // in the array of ids given, and 
    // highlight those.
    d3.selectAll("circle.dot").filter(function(d, i) {
      return ids.indexOf(d.id) > -1;
    })
    .classed("highlighted", true);
  }



	};
	
	//draw the heatmap for the cell - currently not instantiated 
	var drawCellHM = function(where, cellID, classif){
	  
		if (where === '#cellHM') {
			d3.select(where).html('');
		}

		//console.log("inside drawCellHM method - " + cellID)
		d3.text(TSNE_file, function(datasetText){

		var parsedCSV = d3.tsv.parseRows(datasetText);
		//console.log("inside drawCellHM method - " + cellID)
		var cellOI = cellID;
		//console.log("cellID: " + cellOI)
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
		var ylab = ["cell_OI: " + cellOI, "avg expression within class"+classOI, "overall avg expression"];

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
			  //console.log("cellid is equal to expected")
			  //console.log("expected class: " + classOI)
			  //console.log("class : " + parsedCSV[i][parsedCSV[0].length-4])
			}
			//console.log(parsedCSV[i][parsedCSV[0].length-4] )
			if (parsedCSV[i][parsedCSV[0].length-4] == classOI){
				//console.log("classIO = class")
				sum_within_class = sum_within_class.concat([temp]);
				count_within = count_within + 1;
			} 
			sum_total = sum_total.concat([temp]);
			count_total = count_total + 1;
		}
		//console.log("count_within: " + count_within)
		//console.log("count_total: " + count_total)

		swc = [];
		sumtot=[];
		for(var l = 0; l < sum_within_class[0].length; l++){
		  sum = 0;
		  //console.log("sum_within_class.length: " + sum_within_class.length)
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
		//0,500,1000,1500,2000,2500,3000,3500
		
				function colorPicker(d) {
			if( d <= 0 ){
				return "rgb(255, 247, 251)"; // 0
			} else if (d > 0 && d<= 500){
				return "rgb(236, 226, 240)"; // 1
			} else if (d > 500 && d<= 1000) {
				return "rgb(208, 209, 230)"; // 2
			} else if (d > 1000 && d<= 1500) {
				return "rgb(166, 189, 219)"; // 3
			} else if (d > 1500 && d<= 2000){
				return "rgb(103, 169, 207)"; // 4
			} else if (d > 2000 && d<= 2500){
				return "rgb(54, 144, 192)"; // 5
			} else if (d > 2500 && d<= 3000){
				return "rgb(2, 129, 138)"; // 6
			} else if (d > 3000 && d<= 3500) {
				return "rgb(1, 108, 89)"; // 7
			} else if (d > 3500){
				return "rgb(1, 70, 54)"; // 7+
			} else {
				return "pink"; // number not in range
			}
		}
		/*
		function colorPicker(d) {
			if( d <= 0 ){
				return "rgb(255, 247, 251)"; // 0
			} else if (d > 0 && d<= 2000){
				return "rgb(236, 226, 240)"; // 1
			} else if (d > 2000 && d<= 4000) {
				return "rgb(208, 209, 230)"; // 2
			} else if (d > 4000 && d<= 6000) {
				return "rgb(166, 189, 219)"; // 3
			} else if (d > 6000 && d<= 8000){
				return "rgb(103, 169, 207)"; // 4
			} else if (d > 8000 && d<= 10000){
				return "rgb(54, 144, 192)"; // 5
			} else if (d > 10000 && d<= 12000){
				return "rgb(2, 129, 138)"; // 6
			} else if (d > 12000 && d<= 14000) {
				return "rgb(1, 108, 89)"; // 7
			} else if (d > 14000){
				return "rgb(1, 70, 54)"; // 7+
			} else {
				return "pink"; // number not in range
			}
		}*/


		var margin ={top: 160, right: 160, bottom: 80, left: 160}, //{top: 160, right: 160, bottom: 80, left: 160},
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
						   return d;
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
			
			//var legend = [0,2000,4000,6000,8000,10000,12000,14000]; //[0, 2, 4, 6, 8, 10, 12, 14, 16,18,20,22];
			//var legend = [0,1000,2000,3000,4000,5000,6000,7000]; //[0, 2, 4, 6, 8, 10, 12, 14, 16,18,20,22];
			var legend = [0,500,1000,1500,2000,2500,3000,3500]; //[0, 2, 4, 6, 8, 10, 12, 14, 16,18,20,22];

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
		   //console.log("end of cell dispplay method")
						 
		});
		  
		  
	};

	//Draw the distribution of expression values for the selected gene - currently not instantiated
	var drawGeneDist = function(where, geneID){
	  
		/*
		if (where === '#geneDist') {
			d3.select(where).html('');
		}*/

		function parser(d) {
			d.classif = +d.classif
			d.gDFFB = +d.DFFB;
			d.gHES2 = +d.HES2;
			d.gH6PD = +d.H6PD;
			var a = geneID
			//console.log("a: " + a)
			d.gene = d[a]
			//console.log("d[a]: " + d[a])
			return d;
		}

		function genehist(tsvdata) {
			//console.log(tsvdata)
			var binsize = 2;
			var minbin = 0;
			var maxbin = 100;
			var numbins = (maxbin - minbin) / binsize;

			// whitespace on either side of the bars 
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
			
			//console.log(histdata);

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
			  .attr("border","1px solid black;")
			  .attr("style", "outline: thin solid black;")   //This will do the job
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
		d3.tsv(TSNE_file, parser,
			   function(error, tsvdata) {
			   genehist(tsvdata);
		}); 


	}; //end of visualizeGeneDistr method

  $scope.visualizeTSNE = function(){//sharedProperties){
    console.log('conroller: will load genes1');
    drawTSNE('#TSNE',sharedProperties.getProperty(),sharedProperties.getAnalysisDir()+"\\d3data.tsv");
	//PIP('#TSNE');
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