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

#include "Lattice.hpp"

using namespace std;

// IMPLEMENTATION OF FFDLATTICE ***********************************************//
/*
 *	\date			24/mar/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief Structured lattice.
 *
 *	Basically, it builds an elemental 3D shape
 *  (box, sphere, cylinder or part of them) around the geometry and set a structured cartesian mesh of control
 *  points on it (lattice). NO displacements for control points and NO NURBS parameters for FFD are present
 *  in this structure, only geometrical information are stored in the object.
 *
 *
 */

/*! Basic Constructor. Doing nothing.*/
Lattice::Lattice(){
	m_name = "MiMMO.Lattice";
};

///*! Custom constructor.Set lattice mesh, dimensions and curve degree for Nurbs trivariate parameterization.
// *
// * \param[in] origin point origin in global reference system
// * \param[in] span span for each shape coordinate in space (local r.s.)
// * \param[in] type BasicShape::ShapeType enum identifies the shape
// * \param[in] dimensions number of control nodes for each direction
// * \param[in] degrees   curve degrees for each direction;
// */
//Lattice::Lattice(darray3E &origin, darray3E & span, BasicShape::ShapeType type, ivector1D &dimensions,
//						ivector1D & degrees):Lattice(){
//	setMesh(origin, span, type, dimensions, degrees);
//};
//
///*! Custom Constructor.Set lattice mesh, dimensions and curve degree for Rational Bezier trivariate parameterization.
// *  Knots structure is built with curve degrees as in case of a Pure Bezier Volumetric
// *  Parameterization, that is degX = nx-1, degY = ny-1, degZ=nz-1.
// *
// * \param[in] origin point origin in global reference system
// * \param[in] span span for each shape coordinate in space (local r.s.)
// * \param[in] type BasicShape::ShapeType enum identifies the shape
// * \param[in] dimensions number of control nodes for each direction
// */
//Lattice::Lattice(darray3E &origin, darray3E & span, BasicShape::ShapeType type, ivector1D &dimensions
//					   ):Lattice(){
//	   setMesh(origin, span, type, dimensions);
//};
//
///*! Custom Constructor.Set lattice mesh, dimensions and curve degree for Nurbs trivariate parameterization.
// *
// * \param[in] shape pointer to an external BasicShape object
// * \param[in] dimensions number of control nodes for each direction
// * \param[in] degrees   curve degrees for each direction;
// */
//Lattice::Lattice(BasicShape * shape, ivector1D &dimensions, ivector1D & degrees):Lattice(){
//	setMesh(shape, dimensions, degrees);
//};
//
///*! Set lattice mesh, dimensions and curve degree for Rational Bezier trivariate parameterization.
// *  Knots structure is built with curve degrees as in case of a Pure Bezier Volumetric
// *  Parameterization, that is degX = nx-1, degY = ny-1, degZ=nz-1.
// *
// * \param[in] shape pointer to an external BasicShape object
// * \param[in] dimensions number of control nodes for each direction
// *
// */
//Lattice::Lattice(BasicShape * shape, ivector1D &dimensions):Lattice(){
//	setMesh(shape, dimensions);
//};

/*! Destructor */
Lattice::~Lattice(){};

/*! Copy Constructor
 *\param[in] other Lattice where copy from
 */
Lattice::Lattice(const Lattice & other){
	*this = other;
};

/*! Copy Operator
 * \param[in] other Lattice where copy from
 */
Lattice & Lattice::operator=(const Lattice & other){

	*(static_cast<UStructMesh *>(this))  = *(static_cast<const UStructMesh *>(&other));
	m_intMapDOF = other.m_intMapDOF;
	m_np = other.m_np;
	return(*this);
};

/*!Clean all stuffs in your lattice */
void Lattice::clearLattice(){
	clear(); //base manipulation stuff clear
	clearMesh(); // structured mesh cleaned
};

/*! Get number of control nodes in each space direction.
 * \return Array of control nodes numbers in each direction
 */
iarray3E		Lattice::getDimension(){
	iarray3E dim;
	dim[0] = m_nx + 1;
	dim[1] = m_ny + 1;
	dim[2] = m_nz + 1;
	return(dim);
}

/*! Get the total number of control nodes.
 * \return Number of control nodes
 */
double		Lattice::getNNodes(){
	return(m_np);
}



dvecarr3E
Lattice::getGlobalCoords(){
	int np = (getNNodes());
	dvecarr3E coords(np);
	int index, i0, i1, i2;
	for (int i=0; i<np; i++){
		index = accessGridFromDOF(i);
		accessPointIndex(index,i0,i1,i2);
		coords[i] = getGlobalPoint(i0,i1,i2);
	}
	return(coords);
};

dvecarr3E
Lattice::getLocalCoords(){
	int np = (getNNodes());
	dvecarr3E coords(np);
	int index, i0, i1, i2;
	for (int i=0; i<np; i++){
		index = accessGridFromDOF(i);
		accessPointIndex(index,i0,i1,i2);
		coords[i] = getLocalPoint(i0,i1,i2);
	}
	return(coords);
};



/*! Set number of control nodes in each space direction.Nurbs curves are treated as
 * Bezier curves, their degree is automatically set. Weights are reset to unitary value
 * \param[in] dimension vector of control nodes numbers in each direction
 */
void		Lattice::setDimension(ivector1D dimensions){

		if(getShape() ==NULL) return;
		if(dimensions.size() < 3 || getShape() ==NULL) return;
		ivector1D dimLimit(3,2);
		switch(getShapeType()){
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
//		m_nx = std::max(dimensions[0], dimLimit[0])-1;
//		m_ny = std::max(dimensions[1], dimLimit[1])-1;
//		m_nz = std::max(dimensions[2], dimLimit[2])-1;
		//TODO RIGHT IN THIS WAY?? WITH -1 NODES AND DISPL ARE NOT COHERENT
		m_nx = std::max(dimensions[0], dimLimit[0]);
		m_ny = std::max(dimensions[1], dimLimit[1]);
		m_nz = std::max(dimensions[2], dimLimit[2]);

		rebaseMesh();

		//setting knots and eventually weights to non-rational B-Spline
		resizeDisplacements(m_nx+1, m_ny+1, m_nz+1);

};

/*! Set number of control nodes in each space direction and degrees of Nurbs curves.
 *  Weights are reset to unitary value
 * \param[in] dimension vector of control nodes numbers in each direction
 * \param[in] degrees vector of degree of nurbs curve in each direction
 */
void		Lattice::setDimension(ivector1D &dimensions, ivector1D &degrees){

	if(dimensions.size() < 3 || degrees.size() <3 || getShape() ==NULL) return;

	ivector1D dimLimit(3,2);
	switch(getShapeType()){
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
//	m_nx = std::max(dimensions[0], dimLimit[0])-1;
//	m_ny = std::max(dimensions[1], dimLimit[1])-1;
//	m_nz = std::max(dimensions[2], dimLimit[2])-1;
	//TODO RIGHT IN THIS WAY?? WITH -1 NODES AND DISPL ARE NOT COHERENT
	m_nx = std::max(dimensions[0], dimLimit[0]);
	m_ny = std::max(dimensions[1], dimLimit[1]);
	m_nz = std::max(dimensions[2], dimLimit[2]);

	rebaseMesh();

	resizeDisplacements(m_nx+1, m_ny+1, m_nz+1);

};

/*! Set number of control nodes in each space direction.Nurbs curves are treated as
 * Bezier curves, their degree is automatically set. Weights are reset to unitary value
 * \param[in] dimension vector of control nodes numbers in each direction
 */
void		Lattice::setDimension(iarray3E dimensions){

		if(getShape() ==NULL) return;
		ivector1D dimLimit(3,2);
		switch(getShapeType()){
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
//		m_nx = std::max(dimensions[0], dimLimit[0])-1;
//		m_ny = std::max(dimensions[1], dimLimit[1])-1;
//		m_nz = std::max(dimensions[2], dimLimit[2])-1;
		//TODO RIGHT IN THIS WAY?? WITH -1 NODES AND DISPL ARE NOT COHERENT
		m_nx = std::max(dimensions[0], dimLimit[0]);
		m_ny = std::max(dimensions[1], dimLimit[1]);
		m_nz = std::max(dimensions[2], dimLimit[2]);

		rebaseMesh();

		resizeDisplacements(m_nx+1, m_ny+1, m_nz+1);

};

/*! Set span of your shape, according to its local reference system
 * \param[in] s0 first coordinate span
 * \param[in] s1 second coordinate span
 * \param[in] s2 third coordinate span
 * \param[in] flag if true, lattice is rebuilt according to the new input.TRUE is default
 */
void Lattice::setSpan(double s0, double s1, double s2, bool flag){
	UStructMesh::setSpan(s0,s1,s2, flag);
}

/*! Set span of your shape, according to its local reference system
 * Lattice is rebuilt according to the new input.
 * \param[in] Coordinates span
 */
void Lattice::setSpan(darray3E s){
	getShape()->setSpan( s[0], s[1], s[2]);
	UStructMesh::setSpan(s);
}

/*! Set coordinate's origin of your shape, according to its local reference system
 * \param[in] orig first coordinate origin
 * \param[in] dir 0,1,2 int flag identifying coordinate
 * \param[in] flag if true, lattice is rebuilt according to the new input.TRUE is default
 */
void Lattice::setInfLimits(double orig, int dir, bool flag){
	UStructMesh::setInfLimits( orig, dir, flag);
}

/*! Set coordinates' origin of your shape, according to its local reference system.
 *  Lattice is rebuilt according to the new input.
 * \param[in] orig coordinates origin
 */
void Lattice::setInfLimits(darray3E orig){
	UStructMesh::setInfLimits(orig);
}

/*! Set coordinate type of Lattice core shape. See BasicShape::CoordType enum
 * \param[in] type coordinate type
 * \param[in] dir  0,1,2 flag for coordinate
 * \param[in] flag if true, force lattice nodal structure to be updated.TRUE is default
 */
void Lattice::setCoordType(BasicShape::CoordType type, int dir, bool flag){
	getShape()->setCoordinateType(type,dir);
	if(flag){
		iarray3E dim = getDimension();
		resizeDisplacements(dim[0],dim[1],dim[2]);
	}
}

/*! Set x-coordinate type of Lattice core shape. See BasicShape::CoordType enum
 * Force lattice nodal structure to be updated.
 * \param[in] type coordinate type
 */
void Lattice::setCoordTypex(BasicShape::CoordType type){
	getShape()->setCoordinateType(type,0);
	iarray3E dim = getDimension();
	resizeDisplacements(dim[0],dim[1],dim[2]);
}

/*! Set y-coordinate type of Lattice core shape. See BasicShape::CoordType enum
 * Force lattice nodal structure to be updated.
 * \param[in] type coordinate type
 */
void Lattice::setCoordTypey(BasicShape::CoordType type){
	getShape()->setCoordinateType(type,1);
	iarray3E dim = getDimension();
	resizeDisplacements(dim[0],dim[1],dim[2]);
}

/*! Set z-coordinate type of Lattice core shape. See BasicShape::CoordType enum
 * Force lattice nodal structure to be updated.
 * \param[in] type coordinate type
 */
void Lattice::setCoordTypez(BasicShape::CoordType type){
	getShape()->setCoordinateType(type,2);
	iarray3E dim = getDimension();
	resizeDisplacements(dim[0],dim[1],dim[2]);
}

/*! Set coordinates type of Lattice core shape. See BasicShape::CoordType enum
 * Force lattice nodal structure to be updated.
 * \param[in] type coordinates type
 */
void Lattice::setCoordType(array<BasicShape::CoordType,3> type){
	for (int i=0; i<3; i++){
		getShape()->setCoordinateType(type[i],i);
	}
	iarray3E dim = getDimension();
	resizeDisplacements(dim[0],dim[1],dim[2]);
}

/*! Set lattice mesh, dimensions and curve degree for Nurbs trivariate parameterization.
 *
 * \param[in] origin point origin in global reference system
 * \param[in] span span for each shape coordinate in space (local r.s.)
 * \param[in] type BasicShape::ShapeType enum identifies the shape
 * \param[in] dimensions number of control nodes for each direction
 * \param[in] degrees   curve degrees for each direction;
 */
void Lattice::setMesh(darray3E &origin,darray3E & span, BasicShape::ShapeType type, ivector1D & dimensions, ivector1D & degrees){

	clearMesh();
	UStructMesh::setMesh(origin,span,type,dimensions);

	//reallocate your displacement node
	iarray3E dd = getDimension();
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
void Lattice::setMesh(darray3E &origin,darray3E & span, BasicShape::ShapeType type, ivector1D & dimensions){

	clearMesh();
	UStructMesh::setMesh(origin,span,type,dimensions);

	//reallocate your displacement node
	iarray3E dd = getDimension();
	resizeDisplacements(dd[0],dd[1],dd[2]);

};

/*! Set lattice mesh, dimensions and curve degree for Nurbs trivariate parameterization.
 *
 * \param[in] shape pointer to an external BasicShape object
 * \param[in] dimensions number of control nodes for each direction
 * \param[in] degrees   curve degrees for each direction;
 */
void Lattice::setMesh(BasicShape * shape, ivector1D & dimensions, ivector1D & degrees){

	clearMesh();
	UStructMesh::setMesh(shape,dimensions);

	//reallocate your displacement node
	iarray3E dd = getDimension();
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
void Lattice::setMesh(BasicShape * shape, ivector1D & dimensions){

	clearMesh();
	UStructMesh::setMesh(shape,dimensions);

	//reallocate your displacement node
	iarray3E dd = getDimension();
	resizeDisplacements(dd[0],dd[1],dd[2]);

};


/*! Find a corrispondent degree of freedom index of a lattice grid node
 * \param[in] index lattice grid global index
 * \param[out] result corrispondent DOF global index
 */
int Lattice::accessDOFFromGrid(int index){
	return(m_intMapDOF[index]);
}

/*! Find a corrispondent lattice grid index of a degree of freedom node
 * \param[in] index DOF global index
 * \param[out] result corrispondent lattice grid global index
 */
int Lattice::accessGridFromDOF(int index){
	return(posVectorFind(m_intMapDOF, index));
}

/*! Plot your current lattice as a structured grid to *vtu file. Wrapped method of plotGrid of father class UCubicMesh.
 * \param[in] directory output directory
 * \param[in] filename  output filename w/out tag
 * \param[in] counter   integer identifier of the file
 * \param[in] binary     boolean flag for 0-"ascii" or 1-"appended" writing
 * \param[in] deformed  boolean flag for plotting 0-"original lattice", 1-"deformed lattice"
 */
void		Lattice::plotGrid(std::string directory, std::string filename,int counter, bool binary){
	dvecarr3E* pnull = NULL;
	UStructMesh::plotGrid(directory, filename, counter, binary,  pnull);
};

/*! Plot your current lattice as a point cloud to *vtu file.Wrapped method of plotCloud of father class UCubicMesh.
 * \param[in] directory output directory
 * \param[in] filename  output filename w/out tag
 * \param[in] counter   integer identifier of the file
 * \param[in] binary     boolean flag for 0-"ascii" or 1-"appended" writing
 * \param[in] deformed  boolean flag for plotting 0-"original lattice", 1-"deformed lattice"
 */
void		Lattice::plotCloud(std::string directory, std::string filename, int counter, bool binary){
	dvecarr3E* pnull = NULL;
	UStructMesh::plotCloud(directory, filename, counter, binary, pnull);
};

/*! Given pointer to a reference geometry and, execute deformation w/ the current setup.
 * Result is stored in BaseManipulation IOData member m_result.
 */
void 		Lattice::execute(){
	rebaseMesh();
};

/*! Resize BaseManipulation class member m_displ to fit a total number od degree of freedom nx*ny*nz.
 * Old structure is deleted and reset to zero.
 *  \param[in] nx number of control nodes in x direction
 *  \param[in] ny number of control nodes in y direction
 *  \param[in] nz number of control nodes in z direction
 */
void 		Lattice::resizeDisplacements(int nx, int ny,int nz){
	//reallocate your displacement node
	m_intMapDOF.clear();
	m_intMapDOF.resize(nx*ny*nz, -1);
	ivector1D::iterator itMapBegin = m_intMapDOF.begin();
	ivector1D::iterator itMap = itMapBegin;
	ivector1D::iterator itMapEnd = m_intMapDOF.end();
	bvector1D info;
	m_np = reduceDimToDOF(nx,ny,nz, info);

	//set m_intMapDOF

	int target;
	int index;
	ivector1D dummy;

	int i0,i1,i2;
	switch(getShapeType()){

		case BasicShape::ShapeType::CYLINDER :
			target=0;
			while(itMap != itMapEnd){

				*itMap = target;
				index = std::distance(itMapBegin, itMap);
				accessPointIndex(index,i0,i1,i2);

				if(info[0] && i0 == 0){
					for(int k=0; k<ny;++k){
						m_intMapDOF[accessPointIndex(i0,k,i2)] = target;
					}
				}
				if(info[1] && i1 == 0){
					m_intMapDOF[accessPointIndex(i0,ny-1,i2)] = target;
				}

				itMap = find(m_intMapDOF.begin(), itMapEnd,-1);
				target++;
			}
			break;

		case BasicShape::ShapeType::SPHERE :

			target = 0;
			while(itMap != itMapEnd){

				*itMap = target;
				index = std::distance(itMapBegin, itMap);
				accessPointIndex(index,i0,i1,i2);

				if(info[0] && i0 == 0){
					for(int k1=0; k1<ny;++k1){
						for(int k2=0; k2<nz; ++k2){
							m_intMapDOF[accessPointIndex(i0,k1,k2)] = target;
						}
					}
				}

				if(info[1] && i1 == 0){
					m_intMapDOF[accessPointIndex(i0,ny-1,i2)] = target;
				}

				if(info[2] && i2 == 0){
					for(int k1=0; k1<ny; ++k1){
							m_intMapDOF[accessPointIndex(i0,k1,i2)] = target;
						}
				}

				if(info[3] && i2 == (nz-1)){
					for(int k1=0; k1<ny;++k1){
						m_intMapDOF[accessPointIndex(i0,k1,i2)] = target;
					}
				}

				itMap = find(m_intMapDOF.begin(), itMapEnd,-1);
				target++;
			}
			break;


		case BasicShape::ShapeType::CUBE :
			target = 0;
			while(itMap != itMapEnd){

				*itMap = target;
				itMap = find(m_intMapDOF.begin(), itMapEnd,-1);
				target++;
			}
			break;

		default: //doing nothing
			break;
	}//end switch

}

/*!Get the effective dof size of the lattice according to its shape. Return info
 * to build successfully m_intMapDOF
 */
int
Lattice::reduceDimToDOF(int nx, int ny, int nz, bvector1D & info){

	int delta = 0;
	int dum = 0;
	double dval;
	switch(getShapeType()){

		case BasicShape::ShapeType::CYLINDER :
			delta += nz;
			nx--;
			if(getCoordType(1) == BasicShape::CoordType::PERIODIC)	ny--;

			info.push_back(true);
			info.push_back(getCoordType(1) == BasicShape::CoordType::PERIODIC);
			break;

		case BasicShape::ShapeType::SPHERE :
			delta ++;
			nx--;
			if(getCoordType(1) == BasicShape::CoordType::PERIODIC)	ny--;
			dval = getInfLimits()[2];
			if(dval == 0.0)	{
				nz--;
				delta += nx;
			}
			if((dval + getLocalSpan()[2]) == M_PI){
				nz--;
				delta += nx;
			}

			info.push_back(true);
			info.push_back(getCoordType(1) == BasicShape::CoordType::PERIODIC);
			info.push_back(dval==0.0);
			info.push_back((dval + getLocalSpan()[2]) == M_PI);
			break;

		default:
			//doing nothing
			break;
	}

	int result = nx*ny*nz + delta;
	return(result);
};
