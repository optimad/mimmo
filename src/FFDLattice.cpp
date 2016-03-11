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
	m_mapNodes.resize(3);
	m_globalDispl = false;
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
	m_globalDispl = other.m_globalDispl;
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
		
		if(getShape() ==NULL) return;
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

		//reorder dimensions
		orderDimension();

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

	//reorder dimensions
	orderDimension();

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

	//reorder dimensions
	orderDimension();
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

	//reorder dimensions
	orderDimension();

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

	//reorder dimensions
	orderDimension();

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

	//reorder dimensions
	orderDimension();

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

/*!Set the displacement of the degree of freedom currently stored in the object.
 * \param[in] displacements Displacements of lattice nodal degrees of freedom.
 * \param[in] flag set displacements as Global(true) or Local (false)
 */
void
FFDLattice::setDisplacements(dvecarr3E & displacements){
	
	if(m_ndeg != displacements.size() || getShape() == NULL) return;
	BaseManipulation::setDisplacements(displacements);
	checkPeriodicDirections();
};

/*! Check if displacements are meant as global-true or local-false*/
bool
FFDLattice::isDisplGlobal(){return(m_globalDispl);}

/*! Set if displacements are meant as global-true or local-false*/
void
FFDLattice::setDisplGlobal(bool flag){m_globalDispl = flag;}

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
				dvecarr3E dispXYZ;
				if(isDisplGlobal()){
					dispXYZ = *(getDisplacements());
				}else{
					dispXYZ = convertDisplToXYZ(); 
				}
				int size = n[0]*n[1]*n[2];
				dvecarr3E data(size);
				for(int i=0; i<size; ++i){
					data[i] = getGlobalPoint(i) + dispXYZ[i];
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
		
		dvecarr3E dispXYZ;
		if(isDisplGlobal()){
			dispXYZ = *(getDisplacements());
		}else{
			dispXYZ = convertDisplToXYZ(); 
		}
		
		int size = n[0]*n[1]*n[2];
		dvecarr3E data(size);
		for(int i=0; i<size; ++i){
			data[i] = getGlobalPoint(i) + dispXYZ[i];
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

	darray3E origin = getShape()->getOrigin();
	m_info->m_origin.resize(3);
	for (int i=0; i<3; i++){
		m_info->m_origin[i] = origin[i];
	}

	dmatrix33E sdr = getShape()->getRefSystem();
	m_info->m_axes.resize(3, dvector1D(3));
	for (int i=0; i<3; i++){
		for (int j=0; j<3; j++){
		m_info->m_axes[i][j] = sdr[i][j];
		}
	}

	ivector1D n = getDimension();
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
			
			MimmoObject * container = getGeometry();
			if(container == NULL ) return;
			
			ivector1D map;
			dvecarr3E localdef = apply(map);
			
			//reset displacement in a unique vector
			int size = container->getNVertex();
			
			dvecarr3E result(size, darray3E{0,0,0});
			
			for(int i=0; i<map.size(); ++i){
				result[container->getMapDataInv(map[i])] = localdef[i];
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
	
	MimmoObject * container = getGeometry();
	if(container == NULL ) return dvecarr3E(0);
	
	bitpit::PatchKernel * tri = container->getGeometry();
	
	freeContainer(list);
	list = getShape()->includeCloudPoints(tri);
	
	return(nurbsEvaluator(list));

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

/*! Convert a target displacement (expressed in local shape ref frame) in XYZ frame
 *	\param[in] target  target displacement
 * 	\param[out] result displacement in xyz ref frame
 */
darray3E FFDLattice::convertDisplToXYZ(darray3E & target, int i){
	
	darray3E scaling = getShape()->getScaling();
	darray3E work;
	for(int i=0; i<3; ++i){
		work[i] = target[i]/scaling[i];
	}
	work += getLocalPoint(i);
	darray3E result = transfToGlobal(work) -  getGlobalPoint(i);
	return(result);
};

/*! Convert and return all target displacements (expressed in local shape ref frame) in XYZ frame
 *	\param[in] target  target displacement
 * 	\param[out] result displacement in xyz ref frame
 */
dvecarr3E FFDLattice::convertDisplToXYZ(){
	
	dvecarr3E * displ = getDisplacements();
	int sizeD = displ->size();
	
	dvecarr3E result(sizeD);
	for(int i=0; i<sizeD; ++i){
		result[i] = convertDisplToXYZ((*displ)[i],i);
	}
	return(result);
};

/*! Return displacement of a given point, under the deformation effect of the whole Lattice. 
 * \param[in] coord 3D point
 * \param[out] result displacement  
 */
darray3E 	FFDLattice::nurbsEvaluator(darray3E & pointOr){
	
	darray3E point = transfToLocal(pointOr);
	darray3E scaling = getShape()->getScaling();
	
	ivector1D knotInterval(3,0);
	dvector2D BSbasis(3);
	dvector1D valH(4,0), temp1(4,0),temp2(4,0), zeros(4,0);
	int uind, vind, wind, index;
	
	dvecarr3E *displ = getDisplacements();

	int i0 = m_mapdeg[0];
	int i1 = m_mapdeg[1];
	int i2 = m_mapdeg[2];

	iarray3E mappedIndex;
	
	for(int i=0; i<3; i++){
		knotInterval[i] = getKnotInterval(point[i],i);
		BSbasis[i] = basisITS0(knotInterval[i], i, point[i]);
	}
	
	uind = knotInterval[i0] - m_deg[i0];
	vind = knotInterval[i1] - m_deg[i1];
	wind = knotInterval[i2] - m_deg[i2];

	for(int i=0; i<=m_deg[i0]; ++i){
		
		mappedIndex[0] = uind+i;
		temp1 = zeros;
		
		for(int j=0; j<=m_deg[i1]; ++j){
			
			mappedIndex[1] = vind+j;
			temp2 = zeros;
			
			for(int k=0; k<=m_deg[i2]; ++k){
				
				mappedIndex[2] = wind+k;
				index = accessMapNodes(mappedIndex[0],mappedIndex[1],mappedIndex[2]);

				for(int intv=0; intv<3; ++intv){
					temp2[intv] += BSbasis[i2][k]*m_weights[index]*(*displ)[index][intv];
				}	
				temp2[3] += BSbasis[i2][k]*m_weights[index];
			}
			for(int intv=0; intv<4; ++intv){
				temp1[intv] += BSbasis[i1][j]*temp2[intv];
			}	
			
		}
		for(int intv=0; intv<4; ++intv){
			valH[intv] += BSbasis[i0][i]*temp1[intv];
		}	
	}

	darray3E outres;
	if(isDisplGlobal()){
		for(int i=0; i<3; ++i){
			outres[i] =  valH[i]/valH[3];
		}
		
	}else{
		//summing scaled displ in local ref frame; 
		for(int i=0; i<3; ++i){
			point[i] +=  valH[i]/(valH[3]*scaling[i]);
		}
	
		//get final displ in global ref frame:
		outres = transfToGlobal(point) - pointOr;
	}
	return(outres);
	
}; 

/*! Return displacement of a given point, under the deformation effect of the whole Lattice.
 * \param[in] coord 3D point
 * \param[out] result displacement
 */
dvecarr3E 	FFDLattice::nurbsEvaluator(ivector1D & list){

	bitpit::PatchKernel * tri = getGeometry()->getGeometry();
	long id;
	int lsize = list.size();
	ivector1D::iterator it, itend = list.end();
	darray3E target, point;
	ivector1D knotInterval(3,0);
	dvector1D BSbasisi0, BSbasisi1, BSbasisi2;
	dvector1D valH(4,0), temp1(4,0),temp2(4,0);

	dvecarr3E *displ = getDisplacements();

	int uind, vind, wind, index;

	int intv, i, j, k;

	int i0 = m_mapdeg[0];
	int i1 = m_mapdeg[1];
	int i2 = m_mapdeg[2];
	iarray3E mappedIndex;

	int md0 = m_deg[i0];
	int md1 = m_deg[i1];
	int md2 = m_deg[i2];

	dvecarr3E outres(lsize);
	dvecarr3E::iterator itout = outres.begin();
	darray3E scaling = getShape()->getScaling();
	
	double bbasisw2;
	double bbasis1, bbasis0;

	for(it = list.begin(); it != itend; ++it){

		for(intv=0; intv<4; ++intv){
			valH[intv] = 0.0;
		}
		id  = *it;
		target = tri->getVertex(id).getCoords();
		point = transfToLocal(target);

		// get reference Interval int the knot matrix
		for(i=0; i<3; i++){
			knotInterval[i] = getKnotInterval(point[i],i);
		}
		BSbasisi0 = basisITS0(knotInterval[i0], i0, point[i0]);
		BSbasisi1 = basisITS0(knotInterval[i1], i1, point[i1]);
		BSbasisi2 = basisITS0(knotInterval[i2], i2, point[i2]);

		uind = knotInterval[i0] - md0;
		vind = knotInterval[i1] - md1;
		wind = knotInterval[i2] - md2;


		for(i=0; i<=md0; ++i){

			mappedIndex[i0] = uind + i;

			for(int intv=0; intv<4; ++intv){
				temp1[intv] = 0.0;
			}

			for(j=0; j<=md1; ++j){

				mappedIndex[i1] = vind + j;

				for(intv=0; intv<4; ++intv){
					temp2[intv] = 0.0;
				}

				for(k=0; k<=md2; ++k){

					mappedIndex[i2] = wind + k;

					index = accessMapNodes(mappedIndex[0], mappedIndex[1], mappedIndex[2]);

					bbasisw2 = BSbasisi2[k]*m_weights[index];

					for(intv=0; intv<3; ++intv){
						temp2[intv] += bbasisw2*(*displ)[index][intv];
					}
					temp2[3] += bbasisw2;

				}
				bbasis1 = BSbasisi1[j];
				for(intv=0; intv<4; ++intv){
					temp1[intv] += bbasis1*temp2[intv];
				}

			}
			bbasis0 = BSbasisi0[i];
			for(intv=0; intv<4; ++intv){
				valH[intv] += bbasis0*temp1[intv];
			}
		}


		if(isDisplGlobal()){

			//adding to local point displ rescaled
			for(i=0; i<3; ++i){
				(*itout)[i] = valH[i]/valH[3];
			}

		}else{

			//adding to local point displ rescaled
			for(i=0; i<3; ++i){
				 point[i]+= valH[i]/(valH[3]*scaling[i]);
			}

			//get absolute displ as difference of
			for(i=0; i<3; ++i){
				(*itout)[i] = transfToGlobal(point)[i] - target[i];
			}

		}

		itout++;

	}//next list id

	displ = NULL;
	itout = outres.end();


	return(outres);

};

/*! Return a specified component of a displacement of a given point, under the deformation effect of the whole Lattice. 
 * \param[in] coord 3D point
 * \param[in] intV component of displacement vector (0,1,2)
 * \param[out] result return displacement disp[intV]
 */
double 		FFDLattice::nurbsEvaluatorScalar(darray3E & coordOr, int targ){
	
	darray3E point = transfToLocal(coordOr);
	double scaling = getShape()->getScaling()[targ];
	
	ivector1D knotInterval(3,0);
	dvector2D BSbasis(3);
	dvector1D valH(2,0), temp1(2,0),temp2(2,0), zeros(2,0);
	int uind, vind, wind, index;
	
	dvecarr3E *displ = getDisplacements();
	
	int i0 = m_mapdeg[0];
	int i1 = m_mapdeg[1];
	int i2 = m_mapdeg[2];
	
	int md0 = m_deg[i0];
	int md1 = m_deg[i1];
	int md2 = m_deg[i2];

	iarray3E mappedIndex;
	
	for(int i=0; i<3; i++){
		knotInterval[i] = getKnotInterval(point[i],i);
		BSbasis[i] = basisITS0(knotInterval[i], i, point[i]);
	}
	
	uind = knotInterval[i0] - md0;
	vind = knotInterval[i1] - md1;
	wind = knotInterval[i2] - md2;
	
	for(int i=0; i<=md0; ++i){
		
		mappedIndex[0] = uind+i;
		temp1 = zeros;
		
		for(int j=0; j<=md1; ++j){
			
			mappedIndex[1] = vind+j;
			temp2 = zeros;
			
			for(int k=0; k<=md2; ++k){
				
				mappedIndex[2] = wind+k;
				index = accessMapNodes(mappedIndex[0],mappedIndex[1],mappedIndex[2]);
				
				temp2[0] += BSbasis[i2][k]*m_weights[index]*(*displ)[index][targ];
				temp2[1] += BSbasis[i2][k]*m_weights[index];
			}
			
			for(int intv=0; intv<2; ++intv){
				temp1[intv] += BSbasis[i1][j]*temp2[intv];
			}	
			
		}
		
		for(int intv=0; intv<2; ++intv){
			valH[intv] += BSbasis[i0][i]*temp1[intv];
		}	
	}

	darray3E res;
	if(isDisplGlobal()){
		res[targ] = valH[0]/valH[1];
	}else{
		point[targ] += valH[0]/(valH[1]*scaling);
		res = transfToGlobal(point)- coordOr;
	}
	return(res[targ]);
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
	int dd1 = m_deg[pos]+1;
	dvector1D basis(dd1,1);
	dvector1D left(dd1,0), right(dd1,0);
	double saved, tmp;
	
	for(int j = 1; j < dd1; ++j){
		saved = 0.0;
		left[j] = coord - getKnotValue(k+1-j, pos);
		right[j]= getKnotValue(k+j, pos) - coord;
		
		for(int r = 0; r < j; ++r){
			tmp = basis[r]/(right[r+1] + left[j-r]);
			basis[r] = saved + right[r+1] * tmp;
			saved = left[j-r] * tmp;
		}//next r	
		
		basis[j] = saved;  
	}//next j
	
	return(basis);
};	

/*!Return list of equispaced knots for the Nurbs curve in a specific lattice direction
 * \param[in] dir 0,1,2 int identifier of Lattice's Nurbs Curve.
 */
dvector1D	FFDLattice::getNodeSpacing(int dir){
	
	dvector1D result;
	int dim = getDimension()[dir];
	double span = getShape()->getLocalSpan()[dir];
	double locOr= getShape()->getLocalOrigin()[dir];
	bool loop = getShape()->areClosedLoops(dir);
	
	if(loop){
		int nn = dim+m_deg[dir]-1;
		result.resize(nn);
		double dKn = span/(dim-1);
		
		int retroOrigin = (m_deg[dir]-1)/2 + (m_deg[dir]-1)%2;
		double origin = locOr-1.0 * retroOrigin * dKn;
		
		for(int i=0; i<nn; ++i){
			result[i] = origin + i*dKn;
		}
		
	}else{
		
		int nn = dim;
		result.resize(nn);
		double dKn = span/(dim-1);
		for(int i=0; i<nn; ++i){
			result[i] =locOr+ i*dKn;
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
	freeContainer(m_mapNodes);
	m_knots.resize(3);
	m_mapEff.resize(3);
	m_deg.resize(3,1);
	m_mapNodes.resize(3);
	
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
	
	setMapNodes(dir);
	
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
	m_displ.clear();
	int size = nx*ny*nz;
	m_displ.resize(size, darray3E{0,0,0});
	m_ndeg = size;
}

/*! Check periodic displacements distribution for each of lattice coordinates.
 */
void FFDLattice::checkPeriodicDirections(){

	checkPeriodicDirections(0);
	checkPeriodicDirections(1);
	checkPeriodicDirections(2);
};

/*! Given a certain distribution of nodes in a specified direction, force displacements of nodes on extremes
 * to be an average of their pre-existent values 
 */
void FFDLattice::checkPeriodicDirections(int dir){
	
	if(!getShape()->areClosedLoops(dir)) return;
	
	iarray3E map, u;
	int index0, index1;
	darray3E val;
	map[0] = dir;
	map[1] = (dir+1)%3;
	map[2] = (dir+2)%3;
	
	ivector1D dim = getDimension();
	
	for(int k=0; k<dim[map[2]]; ++k){
		u[map[2]] = k;
		for(int j=0; j<dim[map[1]]; ++j){
			u[map[1]] = j;
		
			u[map[0]] = 0;
			index0 = accessPointIndex(u[0],u[1],u[2]);
			
			u[map[0]] = dim[map[0]]-1;
			index1 = accessPointIndex(u[0],u[1],u[2]);
			
			val = 0.5*(m_displ[index0] + m_displ[index1]);
			m_displ[index0] = val;
			m_displ[index1] = val;
		}
	}
};

/*! Fill m_mapnodes, to access correct displacement w knots structure 
 * theoretical knot indexing*/
void FFDLattice::setMapNodes( int ind){

		bool loop = getShape()->areClosedLoops(ind);
		int dimdir = getDimension()[ind];
		
		if(loop){
			int nn = dimdir+m_deg[ind];
			m_mapNodes[ind].resize(nn);
			
			int preNNumb = (m_deg[ind]-1)/2 + (m_deg[ind]-1)%2;
			int postNNumb = (m_deg[ind]-1) - preNNumb;
			
			// set the other internal loads
			for(int i=0; i<dimdir; ++i){
				m_mapNodes[ind][i+preNNumb] = i;
			}
			
			//postpend the first preNNumb loads
			int pInd = dimdir - preNNumb -1;
			for(int i=0; i<preNNumb; ++i){
				m_mapNodes[ind][i] = pInd + i;
			}
			//prepend the last postNNumb loads.
			pInd = 1;
			for(int i=0; i<=postNNumb; ++i){
				m_mapNodes[ind][i+preNNumb+dimdir] = pInd+i;
			}
			
		}else{
			m_mapNodes[ind].resize(dimdir);
			for (int i=0; i<dimdir; ++i){
				m_mapNodes[ind][i] = i;
			}
		}
};

///*! Return real global index of a nodal displacement, given its position i,j,k in knots indexing logic*/
//inline int FFDLattice::accessMapNodes(int i, int j, int k){
//	return(accessPointIndex(m_mapNodes[0][i], m_mapNodes[1][j], m_mapNodes[2][k]));
//};

/*! Fill m_mapdeg with the ordered indices of dimensions.
*/
void FFDLattice::orderDimension(){

	map<pair<int,int>, int > mapsort;
	mapsort[make_pair(m_nx,0)] = 0;
	mapsort[make_pair(m_ny,1)] = 1;
	mapsort[make_pair(m_nz,2)] = 2;

	int i=0;
	for (map<pair<int,int>, int >::iterator it = mapsort.begin(); it != mapsort.end(); ++it){
		m_mapdeg[i] = it->second;
		i++;
	}

};
