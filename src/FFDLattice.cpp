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

// IMPLEMENTATION OF FFDLATTICE ***********************************************//
/*
 *	\date			09/feb/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief Free Form Deformation of a 3D surface and point clouds, with structured lattice.
 *
 *	Free Form deformation tool for 3D geometries (surface and point clouds). Basically, it builds an elemental 3D shape 
 *  (box, sphere, cylinder or part of them) around the geometry and set a structured cartesian mesh of control 
 *  points on it (lattice). Displacements of each control point is linked to the geometry inside 
 *  the shape by means of a NURBS volumetric parameterization. Deformation will be applied only to 
 *  those portion of geometry encased into the 3D shape.
 *
 */

/*! Basic Constructor. Doing nothing.*/
FFDLattice::FFDLattice(){
	m_knots.resize(3);
	m_mapEff.resize(3);
	m_deg.resize(3,1);
	
};
/*! Custom constructor.Set lattice mesh, dimensions and curve degree for Nurbs trivariate parameterization.
 *   
 * \param[in] origin point origin in global reference system
 * \param[in] span span for each shape coordinate in space (local r.s.)
 * \param[in] type BasicShape::ShapeType enum identifies the shape
 * \param[in] dimensions number of control nodes for each direction
 * \param[in] degrees   curve degrees for each direction;
 */
FFDLattice::FFDLattice(darray3E &origin, darray3E & span, BasicShape::ShapeType type, ivector1D &dimensions, 
						ivector1D & degrees):FFDLattice(){
	setMesh(origin, span, type, dimensions, degrees);
};

/*! Custom Constructor.Set lattice mesh, dimensions and curve degree for Rational Bezier trivariate parameterization.
 *  Knots structure is built with curve degrees as in case of a Pure Bezier Volumetric 
 *  Parameterization, that is degX = nx-1, degY = ny-1, degZ=nz-1.
 *   
 * \param[in] origin point origin in global reference system
 * \param[in] span span for each shape coordinate in space (local r.s.)
 * \param[in] type BasicShape::ShapeType enum identifies the shape
 * \param[in] dimensions number of control nodes for each direction
 */
FFDLattice::FFDLattice(darray3E &origin, darray3E & span, BasicShape::ShapeType type, ivector1D &dimensions 
					   ):FFDLattice(){
	   setMesh(origin, span, type, dimensions);
};

/*! Custom Constructor.Set lattice mesh, dimensions and curve degree for Nurbs trivariate parameterization.
 *   
 * \param[in] shape pointer to an external BasicShape object
 * \param[in] dimensions number of control nodes for each direction
 * \param[in] degrees   curve degrees for each direction;
 */
FFDLattice::FFDLattice(BasicShape * shape, ivector1D &dimensions, ivector1D & degrees):FFDLattice(){
	setMesh(shape, dimensions, degrees);
};

/*! Set lattice mesh, dimensions and curve degree for Rational Bezier trivariate parameterization.
 *  Knots structure is built with curve degrees as in case of a Pure Bezier Volumetric 
 *  Parameterization, that is degX = nx-1, degY = ny-1, degZ=nz-1.
 *   
 * \param[in] shape pointer to an external BasicShape object
 * \param[in] dimensions number of control nodes for each direction
 *
 */
FFDLattice::FFDLattice(BasicShape * shape, ivector1D &dimensions):FFDLattice(){
	setMesh(shape, dimensions);
};

/*! Destructor */
FFDLattice::~FFDLattice(){};

/*! Copy Constructor
 *\param[in] other FFDLattice where copy from
 */ 
FFDLattice::FFDLattice(const FFDLattice & other){
	*this = other;
};

/*! Copy Operator
 * \param[in] other FFDLattice where copy from
 */ 
FFDLattice & FFDLattice::operator=(const FFDLattice & other){
	
	*(static_cast<UStructMesh *>(this))  = *(static_cast<const UStructMesh *>(&other));

	m_deg = other.m_deg;
	m_knots = other.m_knots;
	m_mapEff = other.m_mapEff;
	m_weights = other.m_weights;
	return(*this);
};

/*!Clean all stuffs in your lattice */
void FFDLattice::clearLattice(){
		clear(); //base manipulation stuff clear
		clearMesh(); // structured mesh cleaned
		clearKnots(); //clear all knots stuff;
};

/*! Return a vector of six elements reporting the real number of knots effectively stored in the current class (first 3 elements)
 * and the theoretical number of knots (last 3 elements) for Nurbs representation (see Nurbs Books of Peigl)
 * \param[out] result six element vector
 */
ivector1D 	FFDLattice::getKnotsDimension(){
	ivector1D res(6,0);
	
		//effectively stored
		res[0] = m_knots[0].size();
		res[1] = m_knots[1].size();
		res[2] = m_knots[2].size();
		
		//theoretical stored
		res[3] = m_mapEff[0].size();
		res[4] = m_mapEff[1].size();
		res[5] = m_mapEff[2].size();
		return(res);
};

/*! Return weight actually set for each control node 
 * \param[out] result list of weights
 */
dvector1D	FFDLattice::getWeights(){
	return(m_weights);
};

/*! Return knots structure and theoretical map of knots distributions for the current Lattice
 * \param[out] knots knots list 
 * \param[out] mapTheo  map of knots theoretical distribution
 */ 
void 		FFDLattice::returnKnotsStructure(dvector2D & knots, ivector2D &mapTheo){
	knots.resize(3);
	mapTheo.resize(3);
	
	returnKnotsStructure(0, knots[0], mapTheo[0]);
	returnKnotsStructure(1, knots[1], mapTheo[1]);
	returnKnotsStructure(2, knots[2], mapTheo[2]);
};

/*! Return knots structure and knot multiplicity vector for a specified Nurbs curve "dir"
 * \param[in] dir integer (0,1,2) identifies nurbs curve in x,y,and z direction respectively
 * \param[out] knots knots list
 * \param[out] mapT theoretical knot map distribution
 */ 
void 		FFDLattice::returnKnotsStructure( int dir, dvector1D & knots, ivector1D & mapT){
	if(dir<0 || dir>2){return;}
	
	knots.resize(m_knots[dir].size());
	mapT.resize(m_mapEff[dir].size());
	
	knots = m_knots[dir];   
	mapT = m_mapEff[dir];
};

/*! Set number of control nodes in each space direction.Nurbs curves are treated as
 * Bezier curves, their degree is automatically set. Weights are reset to unitary value
 * \param[in] dimension vector of control nodes numbers in each direction
 */
void		FFDLattice::setDimension(ivector1D &dimensions){
		
		if(dimensions.size() < 3 || getShape() ==NULL) return;
		ivector1D dimLimit(3,2);
		switch(getShape()->getShapeType()){
			case BasicShape::ShapeType::CYLINDER :
				dimLimit[1] = 5;
				break;
			case BasicShape::ShapeType::SPHERE :
				dimLimit[1] = 5; dimLimit[2] = 3;
				break;
			default://CUBE
				break;
		}
		
		//check on dimensions and eventual closed loops on coordinates.
		m_nx = std::max(dimensions[0], dimLimit[0])-1;
		m_ny = std::max(dimensions[1], dimLimit[1])-1;
		m_nz = std::max(dimensions[2], dimLimit[2])-1;
		
		rebaseMesh();
		
		m_deg[0] = m_nx;
		m_deg[1] = m_ny;
		m_deg[2] = m_nz;
	
		//setting knots and eventually weights to non-rational B-Spline
		setKnotsStructure();
	
		freeContainer(m_weights);
		int size= (m_nx+1)*(m_ny+1)*(m_nz+1);
		m_weights.resize(size, 1.0);
		
		resizeDisplacements(m_nx+1, m_ny+1, m_nz+1);
};

/*! Set number of control nodes in each space direction and degrees of Nurbs curves. 
 *  Weights are reset to unitary value
 * \param[in] dimension vector of control nodes numbers in each direction
 * \param[in] degrees vector of degree of nurbs curve in each direction
 */
void		FFDLattice::setDimension(ivector1D &dimensions, ivector1D &degrees){
	
	if(dimensions.size() < 3 || degrees.size() <3 || getShape() ==NULL) return;
	
	ivector1D dimLimit(3,2);
	switch(getShape()->getShapeType()){
		case BasicShape::ShapeType::CYLINDER :
			dimLimit[1] = 5;
			break;
		case BasicShape::ShapeType::SPHERE :
			dimLimit[1] = 5; dimLimit[2] = 3;
			break;
		default://CUBE
			break;
	}
	
	//check on dimensions and eventual closed loops on coordinates.
	m_nx = std::max(dimensions[0], dimLimit[0])-1;
	m_ny = std::max(dimensions[1], dimLimit[1])-1;
	m_nz = std::max(dimensions[2], dimLimit[2])-1;
	
	rebaseMesh();
	
	m_deg[0] = std::min(m_nx,std::max(1,degrees[0]));
	m_deg[1] = std::min(m_ny,std::max(1,degrees[1]));
	m_deg[2] = std::min(m_nz,std::max(1,degrees[2]));
	
	//setting knots and eventually weights to non-rational B-Spline
	setKnotsStructure();
	
	freeContainer(m_weights);
	int size= (m_nx+1)*(m_ny+1)*(m_nz+1);
	m_weights.resize(size, 1.0);
	
	resizeDisplacements(m_nx+1, m_ny+1, m_nz+1);
};

/*! Set lattice mesh, dimensions and curve degree for Nurbs trivariate parameterization.
 *   
 * \param[in] origin point origin in global reference system
 * \param[in] span span for each shape coordinate in space (local r.s.)
 * \param[in] type BasicShape::ShapeType enum identifies the shape
 * \param[in] dimensions number of control nodes for each direction
 * \param[in] degrees   curve degrees for each direction;
 */
void FFDLattice::setMesh(darray3E &origin,darray3E & span, BasicShape::ShapeType type, ivector1D & dimensions, ivector1D & degrees){
	
	clearMesh();
	UStructMesh::setMesh(origin,span,type,dimensions);
	
	m_deg[0] = std::min(m_nx,std::max(1,degrees[0]));
	m_deg[1] = std::min(m_ny,std::max(1,degrees[1]));
	m_deg[2] = std::min(m_nz,std::max(1,degrees[2]));
	
	//setting knots and eventually weights to non-rational B-Spline
	setKnotsStructure();
	
	//reset your weights
	freeContainer(m_weights);
	ivector1D dd = getDimension();
	int size= dd[0]*dd[1]*dd[2];
	m_weights.resize(size, 1.0);
	
	//reallocate your displacement node
	resizeDisplacements(dd[0],dd[1],dd[2]);
};

/*! Set lattice mesh, dimensions and curve degree for Rational Bezier trivariate parameterization.
 *  Knots structure is built with curve degrees as in case of a Pure Bezier Volumetric 
 *  Parameterization, that is degX = nx-1, degY = ny-1, degZ=nz-1.
 *   
 * \param[in] origin point origin in global reference system
 * \param[in] span span for each shape coordinate in space (local r.s.)
 * \param[in] type BasicShape::ShapeType enum identifies the shape
 * \param[in] dimensions number of control nodes for each direction
 */
void FFDLattice::setMesh(darray3E &origin,darray3E & span, BasicShape::ShapeType type, ivector1D & dimensions){
	
	clearMesh();
	UStructMesh::setMesh(origin,span,type,dimensions);
	
	m_deg[0] = m_nx;
	m_deg[1] = m_ny;
	m_deg[2] = m_nz;
	
	//setting knots and eventually weights to non-rational B-Spline
	setKnotsStructure();
	
	//reset your weights
	freeContainer(m_weights);
	ivector1D dd = getDimension();
	int size= dd[0]*dd[1]*dd[2];
	m_weights.resize(size, 1.0);
	
	//reallocate your displacement node
	resizeDisplacements(dd[0],dd[1],dd[2]);
};

/*! Set lattice mesh, dimensions and curve degree for Nurbs trivariate parameterization.
 *   
 * \param[in] shape pointer to an external BasicShape object
 * \param[in] dimensions number of control nodes for each direction
 * \param[in] degrees   curve degrees for each direction;
 */
void FFDLattice::setMesh(BasicShape * shape, ivector1D & dimensions, ivector1D & degrees){
	
	clearMesh();
	UStructMesh::setMesh(shape,dimensions);
	
	m_deg[0] = std::min(m_nx,std::max(1,degrees[0]));
	m_deg[1] = std::min(m_ny,std::max(1,degrees[1]));
	m_deg[2] = std::min(m_nz,std::max(1,degrees[2]));
	
	//setting knots and eventually weights to non-rational B-Spline
	setKnotsStructure();
	
	//reset your weights
	freeContainer(m_weights);
	ivector1D dd = getDimension();
	int size= dd[0]*dd[1]*dd[2];
	m_weights.resize(size, 1.0);
	
	//reallocate your displacement node
	resizeDisplacements(dd[0],dd[1],dd[2]);
};

/*! Set lattice mesh, dimensions and curve degree for Rational Bezier trivariate parameterization.
 *  Knots structure is built with curve degrees as in case of a Pure Bezier Volumetric 
 *  Parameterization, that is degX = nx-1, degY = ny-1, degZ=nz-1.
 *   
 * \param[in] shape pointer to an external BasicShape object
 * \param[in] dimensions number of control nodes for each direction
 *
 */
void FFDLattice::setMesh(BasicShape * shape, ivector1D & dimensions){
	
	clearMesh();
	UStructMesh::setMesh(shape,dimensions);
	
	m_deg[0] = m_nx;
	m_deg[1] = m_ny;
	m_deg[2] = m_nz;
	
	//setting knots and eventually weights to non-rational B-Spline
	setKnotsStructure();
	
	//reset your weights
	freeContainer(m_weights);
	ivector1D dd = getDimension();
	int size= dd[0]*dd[1]*dd[2];
	m_weights.resize(size, 1.0);
	
	//reallocate your displacement node
	resizeDisplacements(dd[0],dd[1],dd[2]);
};

/*! Modify a weight of a control node. Access to a node in global indexing
 * \param[in] val weight value
 * \param[in] index index of the control node -> gloab indexing
 */
void 		FFDLattice::setNodalWeight(double val, int index){
			m_weights[index] =  val;
};

/*! Modify a weight of a control node. Access to a node in cartesian indexing
 * \param[in] val weight value
 * \param[in] i index of x coordinate
 * \param[in] j index of y coordinates 
 * \param[in] k index of z coordinate
 */
void 		FFDLattice::setNodalWeight(double val, int i, int j, int k){
		int index = accessPointIndex(i,j,k);
		setNodalWeight(val, index);
};


/*! Plot your current lattice as a structured grid to *vtu file. Wrapped method of plotGrid of father class UCubicMesh.
 * \param[in] directory output directory
 * \param[in] filename  output filename w/out tag
 * \param[in] counter   integer identifier of the file
 * \param[in] binary     boolean flag for 0-"ascii" or 1-"appended" writing 
 * \param[in] deformed  boolean flag for plotting 0-"original lattice", 1-"deformed lattice"
 */
void		FFDLattice::plotGrid(std::string directory, std::string filename,int counter, bool binary, bool deformed){
		
		if(deformed){
				ivector1D n =getDimension();
				const dvecarr3E * disp = getDisplacements(); 
				int size = n[0]*n[1]*n[2];
				dvecarr3E data(size);
				for(int i=0; i<size; ++i){
					data[i] = getGlobalPoint(i) + (*disp)[i];
				}
			UStructMesh::plotGrid(directory, filename, counter, binary, &data);
		}else{
			dvecarr3E* pnull = NULL;
			UStructMesh::plotGrid(directory, filename, counter, binary,  pnull);
			
		}
			
	
};

/*! Plot your current lattice as a point cloud to *vtu file.Wrapped method of plotCloud of father class UCubicMesh.
 * \param[in] directory output directory
 * \param[in] filename  output filename w/out tag
 * \param[in] counter   integer identifier of the file
 * \param[in] binary     boolean flag for 0-"ascii" or 1-"appended" writing 
 * \param[in] deformed  boolean flag for plotting 0-"original lattice", 1-"deformed lattice"
 */
void		FFDLattice::plotCloud(std::string directory, std::string filename, int counter, bool binary, bool deformed){
	
	if(deformed){
		ivector1D n = getDimension();
		const dvecarr3E * disp = getDisplacements(); 
		int size = n[0]*n[1]*n[2];
		dvecarr3E data(size);
		for(int i=0; i<size; ++i){
			data[i] = getGlobalPoint(i) + (*disp)[i];
		}
		UStructMesh::plotCloud(directory, filename, counter, binary, &data);
	}else{
		dvecarr3E* pnull = NULL;
		UStructMesh::plotCloud(directory, filename, counter, binary, pnull);
		
	}
	
};

/*! TODO per Edo: ma che fa sto metodo me lo potevi pure lasciar detto no? :-D */
void
FFDLattice::setInfo(){
	m_info->m_naxes = 3;

	ivector1D n =getDimension();
	int size = n[0]*n[1]*n[2];
	m_info->m_coords.resize(size, dvector1D(3));
	for(int i=0; i<size; ++i){
		darray3E data = getGlobalPoint(i);
		for (int j=0; j<3; j++){
			m_info->m_coords[i][j] = data[j];
		}
	}

}

/*! Given pointer to a reference geometry and, execute deformation w/ the current setup */
void 		FFDLattice::execute(){
			
	//TODO see todo note on "dvecarr3E 	FFDLatticeBox::apply(ivector1D & list)" method of the class. 
			MimmoObject * container = getGeometry();
			if(container == NULL ) return;
			
			ivector1D map;
			dvecarr3E localdef = apply(map);
			
			//reset displacement in a unique vector
			bitpit::Patch * tri = container->getGeometry();
			int size = tri->getVertexCount();
			
			dvecarr3E result(size, darray3E{0,0,0});
			
			for(int i=0; i<map.size(); ++i){
				result[map[i]] = localdef[i];
			}
	//debug
// 			for(int i=0; i<size; ++i){
// 				darray3E pV = container->getVertex(i);
// 				pV = pV + result[i];
// 				container->modifyVertex(pV, i);
// 			}	
			
			
			for (int i=0; i<getNChild(); i++){
				setDisplacementsOut(i, result);
			}
};

/*! Apply current deformation setup to a single 3D point. If point is not included in lattice return zero
 * \param[in] point coordinate of the points 
 * \param[out] result point displacement 
 */
darray3E 	FFDLattice::apply(darray3E & point){
	darray3E result;
	result.fill(0.0);
	if(!getShape()->isPointIncluded(point)) return result;
	return(nurbsEvaluator(point));
};

/*! Apply current deformation setup to geometry linked as a MimmoObject container, member of the class 
 * (see method getGeometry).If MimmoObject member m_geometry is NULL,return void results. 
 * \param[out] result list of non-zero displacement of m_geometry vertices 
 * \param[out] map list of ids of non-zero displaced vertex belonging to geometry
 */
dvecarr3E 	FFDLattice::apply(ivector1D & list){
	//TODO now geometry displacement points are recovered from a list of effectively displaced points +
	// an internal map called list. This list contained an absolute id (get_id, set_id methods) of the geometry vertex, according to the 
	// logic of bitpit bitpit::Patch object for geometry data structure. You can rethink this outputs when your knowledge on
	//how the object bitpit::Patch work is more mature.
	
	dvecarr3E result;
	MimmoObject * container = getGeometry();
	if(container == NULL ) return result;
	
	bitpit::Patch * tri = container->getGeometry();
	
	freeContainer(list);
	list = getShape()->includeCloudPoints(tri);

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
dvecarr3E 	FFDLattice::apply(dvecarr3E * point){
	
	dvecarr3E result;
	if(point ==NULL ) return result;
	
	result.resize(point->size(), darray3E{0,0,0});
	ivector1D list = getShape()->includeCloudPoints(*point);
	
	for(int i=0; i<list.size(); ++i){
		darray3E target = (*point)[list[i]];
		result[list[i]] = nurbsEvaluator(target);
	}
	
	return(result);
};

/*! Return displacement of a given point, under the deformation effect of the whole Lattice. 
 * \param[in] coord 3D point
 * \param[out] result displacement  
 */
darray3E 	FFDLattice::nurbsEvaluator(darray3E & pointOr){
	
	darray3E point = getShape()->toLocalCoord(pointOr);
	
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
	dvector2D temp2(n[0]*n[1], dvector1D(4,0.0));
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
		counter +=n[2];
		//dvector2D work2 = getWorkLoad(2, work);
		dvector2D work2;
		getWorkLoad(2, work, work2);
		temp2[j] = getNurbsPoint(knotInterval[2], BSbasis[2], work2);
	}// next j3
	
	freeContainer(work);
	work.resize(n[1],dvector1D(4,0));
	counter=0;
	for (int j=0; j<(n[0]); ++j){ //  counting deformation due to every load applied in lattice control points
		
		for(int k=counter; k<counter+n[1]; ++k){	
			work[k-counter] = temp2[k];
		}
		counter +=(n[1]);
		//dvector2D work2 = getWorkLoad(1, work);
		dvector2D work2;
		getWorkLoad(2, work, work2);
		temp2[j] = getNurbsPoint(knotInterval[1],BSbasis[1], work2);
	}// next j
	
	temp2.resize(n[0]);
	//dvector2D work2 = getWorkLoad(0, temp2);
	dvector2D work2;
	getWorkLoad(2, temp2, work2);
	valH = getNurbsPoint(knotInterval[0],BSbasis[0], work2);
	
	darray3E outres;
	for(int i=0; i<3; ++i){outres[i] = valH[i]/valH[3];}
	
	return(outres);
	
}; 

/*! Return a specified component of a displacement of a given point, under the deformation effect of the whole Lattice. 
 * \param[in] coord 3D point
 * \param[in] intV component of displacement vector (0,1,2)
 * \param[out] result return displacement disp[intV]
 */
double 		FFDLattice::nurbsEvaluatorScalar(darray3E & coordOr, int intV){
	
	darray3E coord = getShape()->toLocalCoord(coordOr);
	
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
		//dvector2D work2 = getWorkLoad(2, work);
		dvector2D work2;
		getWorkLoad(2, work, work2);
		temp2[j] = getNurbsPoint(knotInterval[2], BSbasis[2], work2);
	}// next j3
	
	freeContainer(work);
	work.resize(n[1],dvector1D(2,0));
	counter=0;
	for (int j=0; j<(n[0]); ++j){ //  counting deformation due to every load applied in lattice control points
		
		for(int k=counter; k<counter+n[1]; ++k){	
			work[k-counter] = temp2[k];
		}
		counter +=(n[1]);
		//dvector2D work2 = getWorkLoad(1, work);
		dvector2D work2;
		getWorkLoad(2, work, work2);
		temp2[j] = getNurbsPoint(knotInterval[1],BSbasis[1], work2);
	}// next j
	
	temp2.resize(n[0]);
	//dvector2D work2 = getWorkLoad(1, temp2);
	dvector2D work2;
	getWorkLoad(2, temp2, work2);
	valH = getNurbsPoint(knotInterval[0],BSbasis[0], work2);
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
dvector1D 	FFDLattice::getNurbsPoint(int k, dvector1D & basis, dvector2D & loads){
	
	// provide theoretical interval index, basis functions for that point, and loads row-> complete.
	int p = basis.size() -1;
	int rs = loads[0].size();
	dvector1D res(rs,0.0);
	
	for(int i=0; i <= p; ++i){
		for (int j=0; j<rs; ++j){
			res[j] += basis[i] * loads[k - p + i][j];
		}
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
dvector1D 	FFDLattice::basisITS0(int k, int pos, double coord){
	
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

/*! Given a vector of nodal diplacement loads, associated to a curve dir(0,1,2)
 *  wrap the right list of loads, according to the type of structure (periodic or clamped)
 * associated to the curve
 * \param[in] dir int id 0,1,2 for x,y,z direction Nurbs curve
 * \param[in] loads list of homogeneous nodal displacements associated to dir Nurb curve
 */
dvector2D	FFDLattice::getWorkLoad(int dir, dvector2D & loads){
	
	dvector2D result;
	bool loop = getShape()->areClosedLoops(dir);
	
	if(loop){
		ivector1D dim = getDimension();
		int nn = dim[dir]+m_deg[dir];
		result.resize(nn, dvector1D(loads[0].size(),0));
		int preNNumb = (m_deg[dir]-1)/2 + (m_deg[dir]-1)%2;
		int postNNumb = (m_deg[dir]-1) - preNNumb;

		// loads on extrema get averaged
		result[preNNumb] = 0.5*(loads[0] + loads[dimdir -1]);
		result[preNNumb + dimdir -1] = result[preNNumb];

		// set the other internal loads
		for(int i=1; i<dimdir-1; ++i){
			result[i+preNNumb] = loads[i];
		}
		//postpend the first preNNumb loads
		int pInd = loads.size() - preNNumb -1;
		for(int i=0; i<preNNumb; ++i){
			result[i] = loads[pInd + i];
		}
		//prepend the last postNNumb loads.
		pInd = 1;
		for(int i=0; i<=postNNumb; ++i){
			result[i+preNNumb+dimdir] = loads[pInd+i];
		}
	}else{
		return(loads);
	}

	return(result);
};

/*! Given a vector of nodal diplacement loads, associated to a curve dir(0,1,2)
 *  wrap the right list of loads, according to the type of structure (periodic or clamped)
 * associated to the curve
 * \param[in] dir int id 0,1,2 for x,y,z direction Nurbs curve
 * \param[in] loads list of homogeneous nodal displacements associated to dir Nurb curve
 */
void	FFDLattice::getWorkLoad(int dir, dvector2D & loads, dvector2D & result){

//	dvector2D result;
	bvector1D loop = getShape()->areClosedLoops();

	if(loop[dir]){
		//ivector1D dim = getDimension();
		int dimdir = getDimension()[dir];
		int nn = dimdir+m_deg[dir];
		int ls = loads[0].size();
		result.resize(nn, dvector1D(ls,0));

		int preNNumb = (m_deg[dir]-1)/2 + (m_deg[dir]-1)%2;
		int postNNumb = (m_deg[dir]-1) - preNNumb;

		// loads on extrema get averaged
		for (int j=0; j<ls; ++j){
			result[preNNumb][j] = 0.5*(loads[0][j] + loads[dimdir -1][j]);
			result[preNNumb + dimdir -1][j] = result[preNNumb][j];
		}
		// set the other internal loads
		for(int i=1; i<dimdir-1; ++i){
			for (int j=0; j<ls; ++j){
				result[i+preNNumb][j] = loads[i][j];
			}
		}
		//postpend the first preNNumb loads
		int pInd = loads.size() - preNNumb -1;
		for(int i=0; i<preNNumb; ++i){
			for (int j=0; j<ls; ++j){
				result[i][j] = loads[pInd + i][j];
			}
		}
		//prepend the last postNNumb loads.
		pInd = 1;
		for(int i=0; i<=postNNumb; ++i){
			for (int j=0; j<ls; ++j){
				result[i+preNNumb+dimdir][j] = loads[pInd+i][j];
			}
		}
	}else{
		return(loads);
	}

};

/*!Return list of equispaced knots for the Nurbs curve in a specific lattice direction
 * \param[in] dir 0,1,2 int identifier of Lattice's Nurbs Curve.
 */
dvector1D	FFDLattice::getNodeSpacing(int dir){
	
	dvector1D result;
	ivector1D dim = getDimension();
	darray3E span = getShape()->getLocalSpan();
	bool loop = getShape()->areClosedLoops(dir);
	
	if(loop){
		int nn = dim[dir]+m_deg[dir]-1;
		result.resize(nn);
		double dKn = span[dir]/(dim[dir]-1);
		
		int retroOrigin = (m_deg[dir]-1)/2 + (m_deg[dir]-1)%2;
		double origin = -1.0 * retroOrigin * dKn;
		
		for(int i=0; i<nn; ++i){
			result[i] = origin + i*dKn;
		}
		
	}else{
		
		int nn = dim[dir];
		result.resize(nn);
		double dKn = span[dir]/(dim[dir]-1);
		for(int i=0; i<nn; ++i){
			result[i] = i*dKn;
		}
	}
	
	return(result);
};

/*!Clean all knots stuff in your lattice */
void FFDLattice::clearKnots(){
	
	freeContainer(m_knots);
	freeContainer(m_mapEff);
	freeContainer(m_deg);
	freeContainer(m_weights);
	m_knots.resize(3);
	m_mapEff.resize(3);
	m_deg.resize(3,1);
	
};

/*! Apply knots structure after modifications to Nurbs curve degrees member */
void 		FFDLattice::setKnotsStructure(){
		for(int i=0; i<3; i++){
			setKnotsStructure(i,getShape()->areClosedLoops(i));
		}
}

/*! Apply knots structure after modifications to Nurbs curve degrees member, to a specific curve knots structure.
 * Curve is identified by a int 0,1,2 for Nurbs curves in x,y,z, directions respectively.
 * \param[in] dir int identifier
 * \param[in] flag identifies a closed periodic curve (true) or a clamped one (false)
 */
void 		FFDLattice::setKnotsStructure(int dir, bool flag){
	
	//recover number of node for direction dir;
	ivector1D dim = getDimension();
	int nn = dim[dir];
	
	// free necessary knot structures
	freeContainer(m_knots[dir]);
	freeContainer(m_mapEff[dir]);
	
	dvector1D equinode = getNodeSpacing(dir);
	
	if(!flag){ //clamped curve structure 
		
		m_deg[dir] = min(m_deg[dir], nn-1);
		int kEff = nn - m_deg[dir] + 1;
		int kTheo = nn + m_deg[dir] + 1;
		m_knots[dir].resize(kEff);
		m_mapEff[dir].resize(kTheo,0);
		
		//set knots grid in the dir space considered 1 for this example
		m_knots[dir][0] = equinode[0];
		m_knots[dir][kEff-1] = equinode[equinode.size()-1];
		
		for(int i=1; i<kEff-1; ++i){
			m_knots[dir][i] = 0.0;
			for(int j=i; j<=i+m_deg[dir]-1; ++j){
				m_knots[dir][i] += equinode[j]/((double) m_deg[dir]);
			}
		}
		
		for(int i=m_deg[dir]; i<(m_deg[dir]+kEff); i++){
			m_mapEff[dir][i]=i-m_deg[dir];
		}
		
		for(int i=(m_deg[dir]+kEff); i<kTheo; i++){
			m_mapEff[dir][i]=kEff-1;
		}
		
	}else{ //periodic curve structure
		m_deg[dir] = min(m_deg[dir], nn-1);
		int nEff = nn + m_deg[dir] - 1;
		int kEff = nEff -m_deg[dir] + 1;
		int kTheo = nEff +m_deg[dir] + 1;
		m_knots[dir].resize(kTheo);
		m_mapEff[dir].resize(kTheo,0);
		
		//set knots grid in the dir space considered 1 for this example
		
		m_knots[dir][m_deg[dir]] = 0.0;
		for(int j=0; j<=m_deg[dir]-1; ++j){
			m_knots[dir][m_deg[dir]] += equinode[j]/((double) m_deg[dir]);
		}
		if(abs(m_knots[dir][m_deg[dir]])<1.e-12) m_knots[dir][m_deg[dir]] = 0.0;
		
		for(int i=1; i<kEff; ++i){
			m_knots[dir][i+m_deg[dir]] = 0.0;
			for(int j=i; j<=i+m_deg[dir]-1; ++j){
				m_knots[dir][i+m_deg[dir]] += equinode[j]/((double) m_deg[dir]);
			}
		}
		
		int kend = m_deg[dir] + kEff-1;
		//unclamp the other nodes
		for(int i=0; i<m_deg[dir]; ++i){
			m_knots[dir][m_deg[dir]-i-1] = m_knots[dir][m_deg[dir]-i] - (m_knots[dir][kend-i] - m_knots[dir][kend-i-1]);
			m_knots[dir][kend +1 +i] = m_knots[dir][kend + i] + (m_knots[dir][m_deg[dir]+i+1] - m_knots[dir][m_deg[dir]+i]);
		}
		
		for(int i=0; i<kTheo; i++){
			m_mapEff[dir][i]=i;
		}
	}	
	
};

/*! Given a knots distribution for one curve in direction "dir", return the index of 
 * interval which a coord belongs to
 * \param[in] coord target position
 * \param[in] dir   0,1,2 identifies the three direction x,y,z in space, and the relative knots distribution
 * \param[out] result return the interval index.
 */ 
int  		FFDLattice::getKnotInterval(double coord, int dir){
	
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
double 		FFDLattice::getKnotValue(int index, int dir){
	int target = getKnotIndex(index, dir);
	if(target ==-1){return(-1.0);}
	return(m_knots[dir][target]);
};
/*! Return a knot real index given its theoretical index and a direction in space 
 * \param[in] index theoretical index of the knot
 * \param[in] dir 0,1,2 identifies the three direction x,y,z in space, and the relative knots distribution
 */
int 		FFDLattice::getKnotIndex(int index ,int dir){
	if(index < 0 || index>=m_mapEff[dir].size()) return -1;
	return(m_mapEff[dir][index]);
};

/*! Return a knot theoretical index given its real index and a direction in space 
 * \param[in] index theoretical index of the knot
 * \param[in] dir 0,1,2 identifies the three direction x,y,z in space, and the relative knots distribution
 */
int 		FFDLattice::getTheoreticalKnotIndex(int locIndex,int dir){
	if(locIndex <0 || locIndex >= m_knots[dir].size()){return(-1);}
	
	// search from the end your m_mapEff vector
	ivector1D::reverse_iterator it = find(m_mapEff[dir].rbegin(), m_mapEff[dir].rend(), locIndex);
	int result = std::distance(m_mapEff[dir].begin(), (it.base()-1));
	return(result);
};

/*! Resize BaseManipulation class member m_displ to fit a total number od degree of freedom nx*ny*nz.
 * Old structure is deleted and reset to zero.
 *  \param[in] nx number of control nodes in x direction
 *  \param[in] ny number of control nodes in y direction 
 *  \param[in] nz number of control nodes in z direction
 */
void 		FFDLattice::resizeDisplacements(int nx, int ny,int nz){
	//reallocate your displacement node
	dvecarr3E * displ = getDisplacements();
	freeContainer(*displ);
	int size = nx*ny*nz;
	displ->resize(size, darray3E{0,0,0});
	m_ndeg = size;
}

