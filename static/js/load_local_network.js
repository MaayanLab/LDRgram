
function load_chik_cutoff(inst_cutoff){

		// switch view: remove class and zscore buttons 
		d3.select('#class_buttons').style('display','none');
		// d3.select('#zscore_buttons').style('display','none');

		// set up wait message before request is made 
		$.blockUI({ css: { 
		        border: 'none', 
		        padding: '15px', 
		        backgroundColor: '#000', 
		        '-webkit-border-radius': '10px', 
		        '-moz-border-radius': '10px', 
		        opacity: .8, 
		        color: '#fff' 
		    } });

	// callback function for click tile 
	click_tile = function click_tile(tile_info){
		console.log(tile_info);	
	}

	//!! temporary change to test similarity matrix 
	// d3.json('/LDRgram/static/networks/chik_cutoff_'+inst_cutoff+'.json', function(network_data){
	d3.json('static/networks/LDR_as_cl.json', function(network_data){

		console.log(network_data)

		// new way of making clustergram 
		////////////////////////////////////////
		// define the outer margins of the visualization 
		var outer_margins = {
		    'top':10,
		    'bottom':40,
		    'left':200,
		    'right':45
		  };

		// define arguments object 
		var arguments_obj = {
			'network_data': network_data,
			'svg_div_id': 'svg_div',
			'row_label':'Assays',
			'col_label':'Cell Lines',
		  'outer_margins': outer_margins,
		  // 'input_domain':10
		  'click_tile':click_tile,
		  'order':'rank'
		};

		//!! define mock res_color_dict
		res_color_dict = {};

		// make clustergram: pass network_data and the div name where the svg should be made 
		d3_clustergram.make_clust( arguments_obj );

	  // turn off the wait sign 
	  $.unblockUI();
	});

}


