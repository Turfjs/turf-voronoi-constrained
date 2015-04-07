var 
	fc = require('turf-featurecollection'),
	pol = require('turf-polygon'),
	tin = require('turf-tin-constrained'),
	linestring = require('turf-linestring'),
	point= require('turf-point'),
	midpoint = require('turf-midpoint');

/**
*	Function that process a {@link Polygon} and gives a constrained Voronoi diagram
* @param {Feature<(Polygon)>} poly - single Polygon Feature
* @param {[Feature<(Point)>]} steinerpoints - Array of Steiner Point features: {@linkhttp://www.iue.tuwien.ac.at/phd/fleischmann/node54.html|Steiner Points and Steiner Triangulation} 
* @return {FeatureCollection<(Linestring)>} - walls of the Voronoi cells
* @author	Abel Vázquez
* @version 0.0.6
*/
module.exports = function(poly, steiner){
	
	if (poly.geometry === void 0 || poly.geometry.type !== 'Polygon' ) throw('"turf-voronoi-constrained" only accepts polygon type input');
	
	var
		delaunay = tin(poly, steiner),
		d = delaunay.features,
		vertex =poly.geometry.coordinates.slice(0),
		sl,
		walls = [],
		rad, node, nodex,
		used =[], results =[],
		
		// same wall
		samew = function (w1, w2){
			var b = samep(w1[2], w2[2]) && samep(w1[3], w2[3]);
			return b || samep(w1[2], w2[3]) && samep(w1[3], w2[2]);
		},
		
		// same point
		samep = function (point1, point2){
			if (point1 == void 0 || point2 == void 0) debugger;
			return (point1[0]==point2[0] && point1[1]==point2[1]);
		},
		
		// voro radius, each one has [cc, mid] (circumcenter -> mid point of opposite segment)
		radius = function(tri){	
			var
				points = tri.geometry.coordinates[0].map(function(a){ return point(a)}),
				nor1 = midnormal(points[0], points[1]),
				nor2 = midnormal(points[1], points[2]),
				nor3 = midnormal(points[2], points[0]),
				cc = segmentX(nor1[0], nor1[1], nor2[0], nor2[1]);
			// distance factor from mid to cc
			nor1[4] = Math.abs(cc.on1);
			nor2[4] = Math.abs(cc.on2);
			nor3[4] = Math.abs(segmentX(nor1[0], nor1[1], nor3[0], nor3[1]).on2);
			if (cc.x == null) debugger;
			nor1[0] = nor2[0] = nor3[0] = [cc.x, cc.y];
			return [nor1,nor2,nor3] ;
		},
		
		// get normal line in middle point of segment
		midnormal = function(point1, point2){
			var
				midt = midpoint(point1, point2),
				mid = midt.geometry.coordinates,
				p1 = point1.geometry.coordinates,
				p2 = point2.geometry.coordinates,
				vector = [p2[1]-p1[1], p1[0]-p2[0]],		// TODO: normalize this vector to make distance ratios real
				cc = [mid[0]+vector[0], mid[1]+vector[1]];
			return [cc, mid,point1.geometry.coordinates, point2.geometry.coordinates];
		},
		
		// calculates the intersection between two segments and the factors of distance from each starting point
		segmentX =function (p11, p12, p21,p22) {
			var 
				a = p11[1] - p21[1],
				b = p11[0] - p21[0],
				c = ((p22[1] - p21[1]) * (p12[0]- p11[0])) - ((p22[0] - p21[0]) * (p12[1] - p11[1])), 
				d =  ((p22[0] - p21[0]) * a) - ((p22[1] -p21[1]) * b),
				e =  ((p12[0] - p11[0]) * a) - ((p12[1] -p11[1]) * b),
				f = {
					x: null,
					y: null,
					on1: false,
					on2: false
				};	
			if (c === 0) return f;
			a = d / c;
			b = e / c;
			f.x = p11[0] + (a * (p12[0]- p11[0]));
			f.y = p11[1] + (a * (p12[1] - p11[1]));
			f.on1 = a;
			f.on2 = b;
			return f;
		};
	
	
	var k = 1;
	steiner = steiner || [];	// init steiner
	sl = steiner.slice(0).length; 
	
	// Crescent accuracy loop
	while (k>0){
		k=0;
		d = delaunay.features;
		var La =[], Lb, Lc, Ld, coor, j, edge, ring, tmp = [];
		for (var i=0; i<d.length; i++){
			coor = d[i].geometry.coordinates[0];
			rad = radius(d[i]); 
			cc = rad[0][0];
			
			// cosine law, check for angles > 90º
			for (var j=0; j<3;j++){
				La[j] = Math.pow(coor[j][0]-coor[j+1][0],2) + Math.pow(coor[j][1]-coor[j+1][1],2);
			}
			for (j=0; j<3;j++){
				Lb = (j ==0) ? La[2] : La[j-1]; //-1
				Lc = La[j];	//+1
				Ld = (j == 2) ? La[0] : La[j+1]; // opposite
				if (Lb+Lc-Ld<0) break;	// angle > 90º
			}
			
			if (j>2) continue;
			
			if (d[i].properties.constrained_edge[j] == false){	// if the segment opposed to the obtuse angle is not a constrained edge => add its circumcenter to Steiner points array
				steiner.push(cc);
				k++;
			}else{	// if the segment opposed to the obtuse angle is a constrained edge =>  add its mid point to the polygon ring definition
				edge = [cc,midpoint(point(coor[(j==0)?2:j-1]),point(coor[(j==2)?0:j+1])).geometry.coordinates,coor[(j==0)?2:j-1],coor[(j==2)?0:j+1]];
				for (var m=0; m<vertex.length; m++) {	//rings
					ring = vertex[m];
					tmp[m] = [ring[0]];
					for (var n=0; n< ring.length-1;n++){
						if (samew(edge, [0,0,ring[n], ring[n+1]])==true)  {
								tmp[m].push(edge[1]);
								k++;
						}
						tmp[m].push(ring[n+1]);
					}
					if (samew(edge, [0,0,ring[ring.length-1], ring[0]])) { 
						tmp[m].push(edge[1]);
						k++;
					}
				}
				vertex = tmp.slice(0);
			}
				
		}
		
		// if there are new points, find Delaunay again with the new data
		if (k>0) {
			var tmp = tin(pol(vertex), steiner);
			if (tmp!=-1){
				delaunay = tmp;
				d = delaunay.features
			} else {
				k=0
			}
			
		} 

	}
	
	if (d.length-delaunay.features.length>0) console.warn((d.length-delaunay.features.length) + " triangles out of "+d.length+" don't meet Delaunay condition");
	
	d = delaunay.features;
	for (var i=0; i<d.length; i++){
		walls = walls.concat(radius(d[i]));
	}
	
	for (var i=0; i<walls.length; i++){
		if (used.indexOf(i)>-1) continue; // already used wall
		for (var j=0; j<walls.length; j++){ // for each wall, any other wall
			if (j == i || used.indexOf(j)>-1) continue;	// wall itself or already used
			if (samew(walls[i], walls[j])){	 //share triangle segment
				results.push([walls[i][0],walls[j][0]])
				used.push(i);
				used.push(j);
			}
			
		}
		
		// wall of a border triangle
		if  (used.indexOf(i)<0){
			results.push([walls[i][0],walls[i][1]])
			used.push(i);
		}	
		
	}
	return fc(results.map(function(w){return linestring(w)}));
	
}