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
#include "Operators.hpp"
#include "customOperators.hpp"


using namespace std;
using namespace mimmo;

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
	m_np = 0;
	m_intMapDOF.clear();
	m_name = "MiMMO.Lattice";
};

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
	m_intMapDOF.clear();
};

/*! Get the total number of control nodes.
 * \return Number of control nodes
 */
int		Lattice::getNNodes(){
	return(m_np);
}

/*! Return position of effective mesh nodes in the lattice, in absolute ref system.
 *  Reimplemented from UstructMesh::getGlobalCoords 
 * \return effective mesh nodes position 
 */
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

/*! Return position of effective mesh nodes in the lattice, in local shape ref system.
 *  Reimplemented from UstructMesh::getLocalCoords. 
 * \return effective mesh nodes position 
 */
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

/*! Find a corrispondent degree of freedom index of a lattice grid node
 * \param[in] index lattice grid global index
 * \return corrispondent DOF global index
 */
int Lattice::accessDOFFromGrid(int index){
	return(m_intMapDOF[index]);
}

/*! Find a corrispondent lattice grid index of a degree of freedom node
 * \param[in] index DOF global index
 * \return corrispondent lattice grid global index
 */
int Lattice::accessGridFromDOF(int index){
	return(posVectorFind(m_intMapDOF, index));
}

/*! Plot your current lattice as a structured grid to *vtu file. Wrapped method of plotGrid of father class UCubicMesh.
 * \param[in] directory output directory
 * \param[in] filename  output filename w/out tag
 * \param[in] counter   integer identifier of the file
 * \param[in] binary     boolean flag for 0-"ascii" or 1-"appended" writing
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
 */
void		Lattice::plotCloud(std::string directory, std::string filename, int counter, bool binary){
	dvecarr3E* pnull = NULL;
	UStructMesh::plotCloud(directory, filename, counter, binary, pnull);
};

/*! Build your mesh and create wrapped map of effective degree of freedom of the current lattice mesh.
 *  Map is stored in internal member m_intMapDOF. Reimplemented from UstructMesh::build();
 */
void 		Lattice::build(){
	UStructMesh::build();
	resizeMapDof();
};

/*!Execute your object, that is recall method Lattice::build(). Implementation of pure virtual BaseManipulation::execute
 */
void 		Lattice::execute(){
	build();
};

/*! Resize map of effective nodes of the lattice grid to fit a total number od degree of freedom nx*ny*nz.
 * Old structure is deleted and reset to zero.
 */
void 		Lattice::resizeMapDof(){
	//reallocate your displacement node
	m_intMapDOF.clear();
	m_intMapDOF.resize((m_nx+1)*(m_ny+1)*(m_nz+1), -1);
	ivector1D::iterator itMapBegin = m_intMapDOF.begin();
	ivector1D::iterator itMap = itMapBegin;
	ivector1D::iterator itMapEnd = m_intMapDOF.end();
	bvector1D info;
	m_np = reduceDimToDOF(m_nx+1,m_ny+1,m_nz+1, info);

	//set m_intMapDOF

	int target;
	int index;
	ivector1D dummy;

	int i0,i1,i2;
	switch(getShapeType()){

		case ShapeType::CYLINDER :
			target=0;
			while(itMap != itMapEnd){

				*itMap = target;
				index = std::distance(itMapBegin, itMap);
				accessPointIndex(index,i0,i1,i2);

				if(info[0] && i0 == 0){
					for(int k=0; k<=m_ny;++k){
						m_intMapDOF[accessPointIndex(i0,k,i2)] = target;
					}
				}
				if(info[1] && i1 == 0){
					m_intMapDOF[accessPointIndex(i0,m_ny,i2)] = target;
				}

				itMap = find(m_intMapDOF.begin(), itMapEnd,-1);
				target++;
			}
			break;

		case ShapeType::SPHERE :

			target = 0;
			while(itMap != itMapEnd){

				*itMap = target;
				index = std::distance(itMapBegin, itMap);
				accessPointIndex(index,i0,i1,i2);

				if(info[0] && i0 == 0){
					for(int k1=0; k1<=m_ny;++k1){
						for(int k2=0; k2<=m_nz; ++k2){
							m_intMapDOF[accessPointIndex(i0,k1,k2)] = target;
						}
					}
				}

				if(info[1] && i1 == 0){
					m_intMapDOF[accessPointIndex(i0,m_ny-1,i2)] = target;
				}

				if(info[2] && i2 == 0){
					for(int k1=0; k1<=m_ny; ++k1){
							m_intMapDOF[accessPointIndex(i0,k1,i2)] = target;
						}
				}

				if(info[3] && i2 == (m_nz)){
					for(int k1=0; k1<=m_ny;++k1){
						m_intMapDOF[accessPointIndex(i0,k1,i2)] = target;
					}
				}

				itMap = find(m_intMapDOF.begin(), itMapEnd,-1);
				target++;
			}
			break;


		case ShapeType::CUBE :
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

		case ShapeType::CYLINDER :
			delta += nz;
			nx--;
			if(getCoordType(1) == CoordType::PERIODIC)	ny--;

			info.push_back(true);
			info.push_back(getCoordType(1) == CoordType::PERIODIC);
			break;

		case ShapeType::SPHERE :
			delta ++;
			nx--;
			if(getCoordType(1) == CoordType::PERIODIC)	ny--;
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
			info.push_back(getCoordType(1) == CoordType::PERIODIC);
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
