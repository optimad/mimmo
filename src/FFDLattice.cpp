/*---------------------------------------------------------------------------*\
 * 
 *  MiMMO
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of MiMMO.
 *
 *  MiMMO is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  MiMMO is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with MiMMO. If not, see <http://www.gnu.org/licenses/>.
 *
 \ *---------------------------------------------------------------------------*/

#include "FFDLattice.hpp"

using namespace std;

// IMPLEMENTATION OF FFDLatticeBox***********************************************//
/*
 *	\date			09/feb/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief Free Form Deformation of a 3D surface and point clouds, with cartesian box structured lattice.
 *
 *	Free Form deformation tool for 3D geometries (surface and point clouds). Basically, it builds a 3D box around the geometry and 
 *  and set a structured cartesian mesh of control points on it (lattice). Displacements of each control point is linked to the geometry inside 
 *  the box by means of a NURBS volumetric parameterization. Deformation will be applied only to those portion of geometry
 *  encased into the lattice box.
 *
 */

/*! Basic Constructor. Doing nothing.*/
FFDLatticeBox::FFDLatticeBox(){
	
	m_knots.resize(3);
	m_multK.resize(3);
	m_deg.resize(3,1);
	
};
/*! Custom Constructor. Build the lattice mesh, but still not the Nurbs knots structure. 
 *  Link class to BaseManipulation parents and the target geometry. 
 *
 * \param[in] origin origin of the lattice box
 * \param[in] span   span of the lattice box
 * \param[in] dimension number of nodal points in each direction of the box  
 * \param[in] geometry  pointer to MimmoObject class containg target geometry
 * \param[in] parent	pointer to a BaseManipulation object parent of the current class
 * 
 */
FFDLatticeBox::FFDLatticeBox(darray3E origin, darray3E span, ivector1D dimension, 
							 MimmoObject* geometry, BaseManipulation * child):BaseManipulation(geometry, child){
								
	m_knots.resize(3);
	m_multK.resize(3);
	m_deg.resize(3,1);
	setMesh(origin, span[0],span[1],span[2],dimension[0],dimension[1],dimension[2]);
};

/*! Custom Constructor. Build the lattice mesh, but still not the Nurbs knots structure. 
 *  BaseManipulation parents and the target geometry are still unlinked. 
 *
 * \param[in] origin origin of the lattice box
 * \param[in] span   span of the lattice box
 * \param[in] dimension number of nodal points in each direction of the box  
 */
FFDLatticeBox::FFDLatticeBox(darray3E origin, darray3E span, ivector1D dimension){
	m_knots.resize(3);
	m_multK.resize(3);
	m_deg.resize(3,1);
	setMesh(origin, span[0],span[1],span[2],dimension[0],dimension[1],dimension[2]);
};

/*! Custom Constructor. Build the lattice mesh, assign degrees to Nurbs curve and automatically compile
 *  Nurbs knots structure. Link class to BaseManipulation parents and the target geometry. 
 *
 * \param[in] origin origin of the lattice box
 * \param[in] span   span of the lattice box
 * \param[in] dimension number of nodal points in each direction of the box  
 * \param[in] degrees   Nurbs curve degrees in each space direction 
 * \param[in] geometry  pointer to MimmoObject class containg target geometry
 * \param[in] parent	pointer to a BaseManipulation object parent of the current class
 * 
 */

FFDLatticeBox::FFDLatticeBox(darray3E origin, darray3E span, ivector1D dimension, ivector1D degrees, 
			 MimmoObject* geometry, BaseManipulation * child):BaseManipulation(geometry, child){
								 
	m_knots.resize(3);
	m_multK.resize(3);
	m_deg.resize(3,1);
	setMesh(origin, span[0],span[1],span[2],dimension[0],dimension[1],dimension[2]);
	setDimension(dimension, degrees);

};


/*! Custom Constructor. Build the lattice mesh, assign degrees to Nurbs curve and automatically compile
 *  Nurbs knots structure. BaseManipulation parents and the target geometry unlinked. 
 *
 * \param[in] origin origin of the lattice box
 * \param[in] span   span of the lattice box
 * \param[in] dimension number of nodal points in each direction of the box  
 * \param[in] degrees   Nurbs curve degrees in each space direction 
 */
FFDLatticeBox::FFDLatticeBox(darray3E origin, darray3E span, ivector1D dimension, ivector1D degrees){
	 m_knots.resize(3);
	 m_multK.resize(3);
	 m_deg.resize(3,1);
	setMesh(origin, span[0],span[1],span[2],dimension[0],dimension[1],dimension[2]);
	setDimension(dimension, degrees);
};

FFDLatticeBox::~FFDLatticeBox(){
	cleanAll();
	freeContainer(m_knots);
	freeContainer(m_multK);
	freeContainer(m_deg);
	freeContainer(m_weights);
};

/*! Copy Constructor
 *\param[in] other FFDLatticeBox where copy from
 */ 
FFDLatticeBox::FFDLatticeBox(const FFDLatticeBox & other){
	*this = other;
};

/*! Copy Operator
 * \param[in] other FFDLatticeBox where copy from
 */ 
FFDLatticeBox & FFDLatticeBox::operator=(const FFDLatticeBox & other){
	
	*(static_cast<BaseManipulation *>(this))  = *(static_cast<const BaseManipulation *>(&other));
	*(static_cast<UCubicMesh *>(this))  = *(static_cast<const UCubicMesh *>(&other));

	m_deg = other.m_deg;
	m_knots = other.m_knots;
	m_multK = other.m_multK;
	m_weights = other.m_weights;
	
	return(*this);
};

/*!Clean all stuff in your lattice */
void FFDLatticeBox::cleanAll(){
	
	
	freeContainer(m_knots);
	freeContainer(m_multK);
	freeContainer(m_deg);
	clearMesh();
	clear();

	m_knots.resize(3);
	m_multK.resize(3);
	m_deg.resize(3,1);
};

/*! Return a vector of six elements reporting the real number of knots effectively stored in the current class (first 3 elements)
 * and the theoretical number of knots (last 3 elements) for Nurbs representation (see Nurbs Books of Peigl)
 * \param[out] result six element vector
 */
ivector1D 	FFDLatticeBox::getKnotsDimension(){
	ivector1D res(6,0);
	
		//effectively stored
		res[0] = m_nx - m_deg[0] + 2;
		res[1] = m_ny - m_deg[1] + 2;
		res[2] = m_nz - m_deg[2] + 2;
		
		//theoretical value
		res[3] = m_nx + m_deg[0] + 2;
		res[4] = m_ny + m_deg[1] + 2;
		res[5] = m_nz + m_deg[2] + 2;
		return(res);
};

/*! Return weight actually set for each control node 
 * \param[out] result list of weights
 */
dvector1D	FFDLatticeBox::getWeights(){
	return(m_weights);
};

/*! Return knots structure and knot multiplicity vector for the current Lattice
 * \param[out] knots knots list 
 * \param[out] multK multiplicity knot list 
 */ 
void 		FFDLatticeBox::returnKnotsStructure(dvector2D & knots, ivector2D &multK){
	knots.resize(3);
	multK.resize(3);
	
	returnKnotsStructure("x", knots[0], multK[0]);
	returnKnotsStructure("y", knots[1], multK[1]);
	returnKnotsStructure("z", knots[2], multK[2]);
};

/*! Return knots structure and knot multiplicity vector for a specified Nurbs curve "dir"
 * \param[in] dir x,y,z identifies nurbs curve in x,y,and z direction respectively
 * \param[out] knots knots list
 * \param[out] multK multiplicity knot list
 */ 
void 		FFDLatticeBox::returnKnotsStructure( std::string dir, dvector1D & knots, ivector1D & multK){
	int k = (dir=="x") + 2*(dir=="y") + 3*(dir=="z") -1;    
	if(k<0){return;}
	
	knots.resize(m_knots[k].size());
	multK.resize(m_multK[k].size());
	
	knots = m_knots[k];   
	multK = m_multK[k];
};

/*! Get number of control nodes in each space direction
 * \param[out] result number of control nodes in x,y and z direction
 */
ivector1D   FFDLatticeBox::getDimension(){
			ivector1D result(3,0);
			result[0] = m_nx+1;
			result[1] = m_ny+1;
			result[2] = m_nz+1;
			return(result);
};

/*! Set number of control nodes in each space direction.Nurbs curves are treated as
 * Bezier curves, their degree is automatically set. Weights are reset to unitary value
 * \param[in] dimension vector of control nodes numbers in each direction
 */
void		FFDLatticeBox::setDimension(ivector1D dimension){
		
		if(dimension.size() < 3) return;
		//Pure Bezier
		m_nx = std::max(2,dimension[0]) -1;
		m_ny = std::max(2,dimension[1]) -1;
		m_nz = std::max(2,dimension[2]) -1;
	
		resizeDisplacements(m_nx+1, m_ny+1, m_nz+1);
		
		rebaseMesh();
		
		m_deg[0] = m_nx;
		m_deg[1] = m_ny;
		m_deg[2] = m_nz;
	
		//setting knots and eventually weights to non-rational B-Spline
		setKnotsStructure();
	
		freeContainer(m_weights);
		ivector1D dd = getDimension();
		int size= dd[0]*dd[1]*dd[2];
		m_weights.resize(size, 1.0);
};

/*! Set number of control nodes in each space direction and degrees of Nurbs curves. 
 *  Weights are reset to unitary value
 * \param[in] dimension vector of control nodes numbers in each direction
 * \param[in] curveDegrees vector of degree of nurbs curve in each direction
 */
void		FFDLatticeBox::setDimension(ivector1D dimension, ivector1D curveDegrees){
	if(dimension.size() < 3) return;
	//Pure Bezier
	m_nx = std::max(2,dimension[0]) -1;
	m_ny = std::max(2,dimension[1]) -1;
	m_nz = std::max(2,dimension[2]) -1;
	
	resizeDisplacements(m_nx+1, m_ny+1,m_nz+1);
	
	rebaseMesh();
	
	m_deg[0] = std::min(m_nx,std::max(1,curveDegrees[0]));
	m_deg[1] = std::min(m_ny,std::max(1,curveDegrees[1]));
	m_deg[2] = std::min(m_nz,std::max(1,curveDegrees[2]));
	
	//setting knots and eventually weights to non-rational B-Spline
	setKnotsStructure();
	
	freeContainer(m_weights);
	ivector1D dd = getDimension();
	int size= dd[0]*dd[1]*dd[2];
	m_weights.resize(size, 1.0);
};

/*! set Origin of your mesh
 * \param[in] point new origin
 */
void		FFDLatticeBox::setOrigin(darray3E origin){
	m_origin = origin;
};

/*! Set span of your lattice box
 * \param[in] span array of width, height and depth of your box;
 */
void		FFDLatticeBox::setSpan(darray3E span){
			m_span = span;
			rebaseMesh();
};


/*! Set lattice mesh origin, span size, dimension and nodal data structure. 
 *  Knots structure is built providing curve degrees as in case of a Pure Bezier Volumetric 
 *  Parameterization, that is degX = nx-1, degY = ny-1, degZ=nz-1.
 *   
 * \param[in] origin point origin in global reference system
 * \param[in] spanX width span -> must be > 0;
 * \param[in] spanY height span -> must be > 0;
 * \param[in] spanZ depth span -> must be > 0;
 * \param[in] nx   number of control points in x-direction,2 is default;
 * \param[in] ny   number of control points in y-direction,2 is default;
 * \param[in] nz   number of control points in z-direction,2 is default;
 */
void FFDLatticeBox::setMesh(darray3E origin, double spanX, double spanY, double spanZ, int nx, int ny, int nz){
	
	clearMesh();
	nx = std::max(2, nx); 
	ny = std::max(2, ny); 
	nz = std::max(2, nz); 
	//set Origin and Span
	m_origin = origin;
	m_span[0] = std::fmax(0, spanX); 
	m_span[1] = std::fmax(0, spanY); 
	m_span[2] = std::fmax(0, spanZ); 
	
	// Number of mesh cells
	m_nx = nx-1; m_ny = ny-1; m_nz = nz-1;
	// Resize mesh data structure
	resizeMesh();
	// Create it 
	
	// get mesh spacing
	m_dx = m_span[0]/((double) m_nx);
	m_dy = m_span[1]/((double) m_ny);	
	m_dz = m_span[2]/((double) m_nz);
	// get point distro;
	for (int i = 0; i < m_nx+1; i++) {m_xedge[i] = ((double) i) * m_dx;} 
	for (int i = 0; i < m_ny+1; i++) {m_yedge[i] = ((double) i) * m_dy;}
	for (int i = 0; i < m_nz+1; i++) {m_zedge[i] = ((double) i) * m_dz;}
	// get cell distro
	for (int i = 0; i < m_nx; i++) {m_xnode[i] = m_xedge[i] + 0.5 * m_dx;}
	for (int i = 0; i < m_ny; i++) {m_ynode[i] = m_yedge[i] + 0.5 * m_dy;}
	for (int i = 0; i < m_nz; i++) {m_znode[i] = m_zedge[i] + 0.5 * m_dz;}
	
	
	m_deg[0] = m_nx;
	m_deg[1] = m_ny;
	m_deg[2] = m_nz;
	
	//setting knots and eventually weights to non-rational B-Spline
	setKnotsStructure();
	
	freeContainer(m_weights);
	ivector1D dd = getDimension();
	int size= dd[0]*dd[1]*dd[2];
	m_weights.resize(size, 1.0);
	
	
	//reallocate your displacement node
	resizeDisplacements(nx,ny,nz);
};

/*! Modify a weight of a control node. Access to a node in global indexing
 * \param[in] val weight value
 * \param[in] index index of the control node -> gloab indexing
 */
void 		FFDLatticeBox::setNodalWeight(double val, int index){
			m_weights[index] =  val;
};

/*! Modify a weight of a control node. Access to a node in cartesian indexing
 * \param[in] val weight value
 * \param[in] i index of x coordinate
 * \param[in] j index of y coordinates 
 * \param[in] k index of z coordinate
 */
void 		FFDLatticeBox::setNodalWeight(double val, int i, int j, int k){
		int index = accessPointData(i,j,k);
		setNodalWeight(val, index);
};


/*! Plot your current lattice as a structured grid to *vtu file. Wrapped method of plotGrid of father class UCubicMesh.
 * \param[in] directory output directory
 * \param[in] filename  output filename w/out tag
 * \param[in] counter   integer identifier of the file
 * \param[in] ascii     boolean flag for 0-"ascii" or 1-"appended" writing 
 * \param[in] deformed  boolean flag for plotting 0-"original lattice", 1-"deformed lattice"
 */
void		FFDLatticeBox::plotGrid(std::string directory, std::string filename,int counter, bool ascii, bool deformed){
		
		if(deformed){
				ivector1D n =getDimension();
				const dvecarr3E * disp = getDisplacements(); 
				int size = n[0]*n[1]*n[2];
				dvecarr3E data(size);
				for(int i=0; i<size; ++i){
					data[i] = getGlobalPoint(i) + (*disp)[i];
				}
				BASE_UStructMesh::plotGrid(directory, filename, counter, ascii, &data);
		}else{
			dvecarr3E* pnull = NULL;
			BASE_UStructMesh::plotGrid(directory, filename, counter, ascii,  pnull);
			
		}
			
	
};

/*! Plot your current lattice as a point cloud to *vtu file.Wrapped method of plotCloud of father class UCubicMesh.
 * \param[in] directory output directory
 * \param[in] filename  output filename w/out tag
 * \param[in] counter   integer identifier of the file
 * \param[in] ascii     boolean flag for 0-"ascii" or 1-"appended" writing 
 * \param[in] deformed  boolean flag for plotting 0-"original lattice", 1-"deformed lattice"
 */
void		FFDLatticeBox::plotCloud(std::string directory, std::string filename, int counter, bool ascii, bool deformed){
	
	if(deformed){
		ivector1D n =getDimension();
		const dvecarr3E * disp = getDisplacements(); 
		int size = n[0]*n[1]*n[2];
		dvecarr3E data(size);
		for(int i=0; i<size; ++i){
			data[i] = getGlobalPoint(i) + (*disp)[i];
		}
		BASE_UStructMesh::plotCloud(directory, filename, counter, ascii, &data);
	}else{
		dvecarr3E* pnull = NULL;
		BASE_UStructMesh::plotCloud(directory, filename, counter, ascii, pnull);
		
	}
	
};

/*! Given pointer to a reference geometry and, execute deformation w/ the current setup */
void 		FFDLatticeBox::execute(){
			
	//TODO see todo note on "dvecarr3E 	FFDLatticeBox::apply(ivector1D & list)" method of the class. 
			MimmoObject * container = getGeometry();
			if(container == NULL ) return;
			
//			recoverDisplacementsOut();
			ivector1D map;
			dvecarr3E localdef = apply(map);
			
			//reset displacement in a unique vector
			bitpit::Patch * tri = container->getGeometry();
			int size = tri->getVertexCount();
			
			dvecarr3E result(size, darray3E{0,0,0});
			
			for(int i=0; i<map.size(); ++i){
				result[map[i]] = localdef[i];
			}
			for (int i=0; i<getNChild(); i++){
				setDisplacementsOut(i, result);
			}

};

/*! Apply current deformation setup to a single 3D point. If point is not included in lattice return zero
 * \param[in] point coordinate of the points 
 * \param[out] result point displacement 
 */
darray3E 	FFDLatticeBox::apply(darray3E & point){
	darray3E result;
	result.fill(0.0);
	if(!isPointIncluded(point)) return result;
	return(nurbsEvaluator(point));
};
/*! Apply current deformation setup to geometry linked as a MimmoObject container, member of the class 
 * (see method getGeometry).If MimmoObject member m_geometry is NULL,return void results. 
 * \param[out] result list of non-zero displacement of m_geometry vertices 
 * \param[out] map list of ids of non-zero displaced vertex belonging to geometry
 */
dvecarr3E 	FFDLatticeBox::apply(ivector1D & list){
	//TODO now geometry displacement points are recovered from a list of effectively displaced points +
	// an internal map called list. This list contained an absolute id (get_id, set_id methods) of the geometry vertex, according to the 
	// logic of bitpit bitpit::Patch object for geometry data structure. You can rethink this outputs when your knowledge on
	//how the object bitpit::Patch work is more mature.
	
	dvecarr3E result;
	MimmoObject * container = getGeometry();
	if(container == NULL ) return result;
	
	bitpit::Patch * tri = container->getGeometry();
	
	freeContainer(list);
	list = includeCloudPoints(tri);

	const long nVert = tri->getVertexCount();
	result.resize(list.size(), darray3E{0,0,0});
	
	for(int i=0; i<list.size(); ++i){
		long id  = list[i];
		darray3E target = tri->getVertex(id).getCoords();
		result[i] = nurbsEvaluator(target);
	}
	
	return(result);
};

/*! Apply current deformation setup to a custom list of 3D points. Only points contained into the lattice will be deformed, 
 *  diplacement of the others will be set to zero. If point list is NULL,return void results. 
 * \param[in] point pointer to a list of 3D points. 
 * \param[out] result displacements of points 
 */
dvecarr3E 	FFDLatticeBox::apply(dvecarr3E * point){
	
	dvecarr3E result;
	if(point ==NULL ) return result;
	
	result.resize(point->size(), darray3E{0,0,0});
	ivector1D list = includeCloudPoints(*point);
	
	for(int i=0; i<list.size(); ++i){
		darray3E target = (*point)[list[i]];
		result[list[i]] = nurbsEvaluator(target);
	}
	
	return(result);
};

/*! Apply knots structure after modifications to Nurbs curve degrees member */
void 		FFDLatticeBox::setKnotsStructure(){
	setKnotsStructure("x");
	setKnotsStructure("y");
	setKnotsStructure("z");
}; 
/*! Apply knots structure after modifications to Nurbs curve degrees member, to a specific curve knots structure.
 * Curve is identified by a string "x", "y" or "z" for Nurbs curves in x,y,z, directions respectively.
 * \param[in] dir string identifier
 */
void 		FFDLatticeBox::setKnotsStructure( std::string dir){
	
	int k = (dir=="x") + 2*(dir=="y") + 3*(dir=="z") -1;    
	if(k<0){return;}
	
	freeContainer(m_knots[k]);
	freeContainer(m_multK[k]);
	
	ivector1D n = getDimension();
	
	m_deg[k] = min(m_deg[k], n[k]-1);
	int kEff = n[k] - m_deg[k]+ 1;
	m_knots[k].resize(kEff);
	m_multK[k].resize(kEff,1);
	
	//set knots grid in the dir space
	darray3E spanLatt  = getSpan();
	darray3E orig = getOrigin();
	double dKn = spanLatt[k]/(kEff-1); 
	for(int i=0; i<kEff; ++i){
		m_knots[k][i] = orig[k] + i*dKn;
	}
	
	m_multK[k][0]=m_multK[k][kEff-1]= m_deg[k]+1;
	
}; 

/*! Return displacement of a given point, under the deformation effect of the whole Lattice. 
 * \param[in] coord 3D point
 * \param[out] result displacement  
 */
darray3E 	FFDLatticeBox::nurbsEvaluator(darray3E & point){
	
	const dvecarr3E *disp = getDisplacements();
	ivector1D n = getDimension();
	// get reference Interval int the knot matrix
	ivector1D knotInterval(3,0);
	dvector2D BSbasis(3);
	for(int i=0; i<3; i++){
		knotInterval[i] = getKnotInterval(point[i],i);
		BSbasis[i] = basisITS0(knotInterval[i], i, point[i]);
	}
	
	//get loads in homogeneous coordinate.
	// get local knot.
	
	dvector1D valH(4,0.0);
	dvector2D temp2((n[0])*(n[1]), dvector1D(4,0.0));
	dvector2D work;
	int counter, locT;
	// start of Dimensional reduction algorithm........................................................
	work.resize(n[2], dvector1D(4, 0.0));
	counter=0;
	//reducing dimension z
	for (int j=0; j<(n[0])*(n[1]); ++j){ //  counting deformation due to every load applied in lattice control points
		
		for(int k=counter; k<counter+n[2]; ++k){	
			for(int intv=0; intv<3; intv++){ work[k-counter][intv] = m_weights[k] * (*disp)[k][intv];}
			work[k-counter][3] = m_weights[k];
		}
		counter +=(n[2]);
		temp2[j] = getNurbsPoint(knotInterval[2], BSbasis[2], work);
	}// next j3
	
	freeContainer(work);
	work.resize(n[1],dvector1D(4,0));
	counter=0;
	for (int j=0; j<(n[0]); ++j){ //  counting deformation due to every load applied in lattice control points
		
		for(int k=counter; k<counter+n[1]; ++k){	
			work[k-counter] = temp2[k];
		}
		counter +=(n[1]);
		temp2[j] = getNurbsPoint(knotInterval[1],BSbasis[1], work);
	}// next j
	
	temp2.resize(n[0]);
	valH = getNurbsPoint(knotInterval[0],BSbasis[0], temp2);
	
	darray3E outres;
	for(int i=0; i<3; ++i){outres[i] = valH[i]/valH[3];}
	
	return(outres);
	
}; 

/*! Return a specified component of a displacement of a given point, under the deformation effect of the whole Lattice. 
 * \param[in] coord 3D point
 * \param[in] intV component of displacement vector (0,1,2)
 * \param[out] result return displacement disp[intV]
 */
double 		FFDLatticeBox::nurbsEvaluatorScalar(darray3E & coord, int intV){

	const dvecarr3E *disp = getDisplacements();
	ivector1D n = getDimension();
	// get reference Interval int the knot matrix
	ivector1D knotInterval(3,0);
	dvector2D BSbasis(3);
	for(int i=0; i<3; i++){
		knotInterval[i] = getKnotInterval(coord[i],i);
		BSbasis[i] = basisITS0(knotInterval[i], i, coord[i]);
	}
	
	//get loads in homogeneous coordinate.
	// get local knot.
	
	dvector1D valH(2,0.0);
	dvector2D temp2((n[0])*(n[1]), dvector1D(2,0.0));
	dvector2D work;
	int counter, locT;
	// start of Dimensional reduction algorithm........................................................
	work.resize(n[2], dvector1D(2, 0.0));
	counter=0;
	//reducing dimension z
	for (int j=0; j<(n[0])*(n[1]); ++j){ //  counting deformation due to every load applied in lattice control points
		
		for(int k=counter; k<counter+n[2]; ++k){	
			work[k-counter][0] = m_weights[k] * (*disp)[k][intV];
			work[k-counter][1] = m_weights[k];
		}
		counter +=(n[2]);
		temp2[j] = getNurbsPoint(knotInterval[2], BSbasis[2], work);
	}// next j3
	
	freeContainer(work);
	work.resize(n[1],dvector1D(2,0));
	counter=0;
	for (int j=0; j<(n[0]); ++j){ //  counting deformation due to every load applied in lattice control points
		
		for(int k=counter; k<counter+n[1]; ++k){	
			work[k-counter] = temp2[k];
		}
		counter +=(n[1]);
		temp2[j] = getNurbsPoint(knotInterval[1],BSbasis[1], work);
	}// next j
	
	temp2.resize(n[0]);
	valH = getNurbsPoint(knotInterval[0],BSbasis[0], temp2);
	return(valH[0]/valH[1]);
};

/*!Evaluate NURBS displacement of a point coord via ITS0 algorithm (called after basisITS0 method).
 * Note nodal displacements must be given in homogeneous 4-D coordinate (w*x, w*y, w*z, w), where w is the weight
 * and x,y,z are the components of the original displacement.
 *  \param[in] k local knot interval in which a point coord resides -> theoretical knot indexing
 *  \param[in] basis pre-calculated B-Spline basis of point coord
 *  \param[in] loads nodal displacement of the nurbs curve in exam
 *  \param[out] result point displacement in coord, in homogeneous coordinates
*/
dvector1D 	FFDLatticeBox::getNurbsPoint(int k, dvector1D & basis, dvector2D & loads){
	
	// provide theoretical interval index, basis functions for that point, and loads row-> complete.
	int p = basis.size() -1;
	dvector1D res(loads[0].size(),0.0);
	
	for(int i=0; i <= p; ++i){
		res = res + basis[i] * loads[k - p + i]; 
	}//next[i];
	
	return(res);
};

/*!Return the local basis function of a Nurbs Curve.Please refer to NURBS book of PEIGL 
 * for this Inverted Triangular Scheme Algorithm (pag 74);* 
 *\param[in] k  local knot interval in which coord resides -> theoretical knot indexing, 
 *\param[in] pos identifies which nurbs curve of lattice (3 curve for 3 box direction) you are pointing
 *\param[in] coord the evaluation point on the curve
 *\param[out] result local basis of ITS algorithm-> coefficient of interpolation of control nodes for position u of curve pos
 */
dvector1D 	FFDLatticeBox::basisITS0(int k, int pos, double coord){
	
	//return local basis function given the local interval in theoretical knot index,
	//local degree of the curve -> Please refer to NURBS book of PEIGL for this Inverted Triangular Scheme Algorithm (pag 74);
	int dd = m_deg[pos]; 
	dvector1D basis(dd+1,1);
	dvector1D left(dd+1,0), right(dd+1,0);
	
	for(int j = 1; j <= dd; ++j){
		double saved = 0.0;
		left[j] = coord - getKnotValue(k+1-j, pos);
		right[j]= getKnotValue(k+j, pos) - coord;
		
		for(int r = 0; r < j; ++r){
			double tmp = basis[r]/(right[r+1] + left[j-r]);
			basis[r] = saved + right[r+1] * tmp;
			saved = left[j-r] * tmp;
		}//next r	
		
		basis[j] = saved;  
	}//next j
	
	return(basis);
};	

/*! Rebuild your mesh after a modification is applied in origin, span and dimension parameters*/
void FFDLatticeBox::rebaseMesh(){
	
	m_nx = std::max(1, m_nx);
	m_ny = std::max(1, m_ny);
	m_nz = std::max(1, m_nz);

	m_span[0] = std::fmax(0, m_span[0]); 
	m_span[1] = std::fmax(0, m_span[1]); 
	m_span[2] = std::fmax(0, m_span[2]); 

	reshapeNodalStructure();
	m_dx = m_span[0]/((double) m_nx);
	m_dy = m_span[1]/((double) m_ny);	
	m_dz = m_span[2]/((double) m_nz);
	// get point distro;
	for (int i = 0; i < m_nx+1; i++) {m_xedge[i] = ((double) i) * m_dx;} 
	for (int i = 0; i < m_ny+1; i++) {m_yedge[i] = ((double) i) * m_dy;}
	for (int i = 0; i < m_nz+1; i++) {m_zedge[i] = ((double) i) * m_dz;}
	// get cell distro
	for (int i = 0; i < m_nx; i++) {m_xnode[i] = m_xedge[i] + 0.5 * m_dx;}
	for (int i = 0; i < m_ny; i++) {m_ynode[i] = m_yedge[i] + 0.5 * m_dy;}
	for (int i = 0; i < m_nz; i++) {m_znode[i] = m_zedge[i] + 0.5 * m_dz;}
}

/*! Given a knots distribution for one curve in direction "dir", return the index of 
 * interval which a coord belongs to
 * \param[in] coord target position
 * \param[in] dir   0,1,2 identifies the three direction x,y,z in space, and the relative knots distribution
 * \param[out] result return the interval index.
 */ 
int  		FFDLatticeBox::getKnotInterval(double coord, int dir){
	
	int size = m_knots[dir].size(); 
	if(coord< m_knots[dir][0] ){ return(getTheoreticalKnotIndex(0, dir));}
	if(coord >= m_knots[dir][size-1]){ return(getTheoreticalKnotIndex(size-2, dir));}
	
	int low = 0; 
	int high = size-1; 
	int mid = (low + high)/2;
	while( coord < m_knots[dir][mid] || coord>= m_knots[dir][mid+1]){
		if(coord < m_knots[dir][mid])	{high=mid;}
		else				{low=mid;}
		mid = (low+high)/2;
	}
	return(getTheoreticalKnotIndex(mid, dir));
};

/*! Return value of a knot for a given its theoretical index and a direction in space 
 * \param[in] index theoretical index of the knot
 * \param[in] dir 0,1,2 identifies the three direction x,y,z in space, and the relative knots distribution
 */
double 		FFDLatticeBox::getKnotValue(int index, int dir){
	int target = getKnotIndex(index, dir);
	if(target ==-1){return(-1.0);}
	return(m_knots[dir][target]);
};
/*! Return a knot real index given its theoretical index and a direction in space 
 * \param[in] index theoretical index of the knot
 * \param[in] dir 0,1,2 identifies the three direction x,y,z in space, and the relative knots distribution
 */
int 		FFDLatticeBox::getKnotIndex(int index ,int dir){
	int res = index+1;
	int counter=-1;
	while (res > 0){
		res = res - m_multK[dir][counter+1];
		counter++;
	}
	return(counter);
};

/*! Return a knot theoretical index given its real index and a direction in space 
 * \param[in] index theoretical index of the knot
 * \param[in] dir 0,1,2 identifies the three direction x,y,z in space, and the relative knots distribution
 */
int 		FFDLatticeBox::getTheoreticalKnotIndex(int locIndex,int dir){
	if(locIndex <0){return(-1);}
	return(m_multK[dir][locIndex] + getTheoreticalKnotIndex(locIndex-1, dir));
};

/*! Resize BaseManipulation class member m_displ to fit a total number od degree of freedom nx*ny*nz.
 * Old structure is deleted and reset to zero.
 *  \param[in] nx number of control nodes in x direction
 *  \param[in] ny number of control nodes in y direction 
 *  \param[in] nz number of control nodes in z direction
 */
void FFDLatticeBox::resizeDisplacements(int nx, int ny,int nz){
	//reallocate your displacement node
	dvecarr3E * displ = getDisplacements();
	freeContainer(*displ);
	int size = nx*ny*nz;
	displ->resize(size, darray3E{0,0,0});
	m_ndeg = size;
}

