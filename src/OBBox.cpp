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

#include "OBBox.hpp"
#include "LinearAlgebra.hpp"
// #include <Eigen/Eigenvalues>
#include "lapacke.h"

using namespace std;
using namespace mimmo;

/*
 *	\date			24/mar/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief Oriented Bounding Box calculator.
 *
 *	Builds the oriented bounding box of a 3D object (Point Clouds or superficial tessellations), passed as a MimmoObject;
 *
 *	=========================================================
 * ~~~
 *	|-------------------------------------------------------------------------------------|
 *	|                    Port Input                                                       |
 *	|-------|-------------|---------------------------------------|-----------------------|
 *	|PortID | PortType    | variable/function                     | compatibilities       |
 *	|-------|-------------|---------------------------------------|-----------------------|
 *	| 99    | M_GEOM      | m_geometry                            | 			          |
 *	|-------|-------------|---------------------------------------|-----------------------|
 * 
 *
 *	|-----------------------------------------|
 *	|               Port Output               |
 *	|-------|-------------|-------------------|
 *	|PortID | PortType    | variable/function |
 *	|-------|-------------|-------------------|
 *	| 1     | M_POINT     | getOrigin         |
 *	| 22    | M_AXES      | getRefSystem      |
 *	| 23    | M_SPAN      | getSpan           |
 *	|-------|-------------|-------------------|
 * ~~~
 *	=========================================================
 *
 */

/*! Basic Constructor. Doing nothing.*/
OBBox::OBBox(){
	m_name = "MiMMO.OBBox";
	m_origin.fill(0.0);
	m_span.fill(1.0);
	int counter = 0;
	for(auto &val : m_axes)	{
		val.fill(0.0);
		val[counter] = 1.0;
		++counter;
	}	
};

/*! Destructor */
OBBox::~OBBox(){};

/*! Copy Constructor
 *\param[in] other OBBox where copy from
 */
OBBox::OBBox(const OBBox & other){
	*this = other;
};

/*! Copy Operator
 * \param[in] other OBBox where copy from
 */
OBBox & OBBox::operator=(const OBBox & other){

	*(static_cast<BaseManipulation *>(this))  = *(static_cast<const BaseManipulation *>(&other));
	m_origin = other.m_origin;
	m_span   = other.m_span;
	m_axes = other.m_axes;	
	return(*this);
};

/*! It builds the input/output ports of the object
 */
void OBBox::buildPorts(){

	bool built = true;
//creating input ports	
	built = (built && createPortIn<MimmoObject*, OBBox>(&m_geometry, M_GEOM));
// creating output ports
	built = (built && createPortOut<darray3E, OBBox>(this, &mimmo::OBBox::getOrigin, M_POINT));
	built = (built && createPortOut<dmatrix33E, OBBox>(this, &mimmo::OBBox::getAxes, M_AXES));
	built = (built && createPortOut<darray3E, OBBox>(this, &mimmo::OBBox::getSpan, M_SPAN));

	m_arePortsBuilt = built;
};

/*!Clean all stuffs in your class */
void OBBox::clearOBBox(){
	clear(); //base manipulation stuff clear
	m_origin.fill(0.0);
	m_span.fill(1.0);
	int counter = 0;
	for(auto &val : m_axes)	{
		val.fill(0.0);
		val[counter] = 1.0;
		++counter;
	}	
	
};

/*! 
 * Return the origin of the OBB.
 * \return Number of control nodes
 */
darray3E	OBBox::getOrigin(){
	return(m_origin);
}

/*! 
 * Return the span of the OBB.
 * \return Number of control nodes
 */
darray3E	OBBox::getSpan(){
	return(m_span);
}


/*! 
 * Return the oriented axes of the OBB.
 * \return Number of control nodes
 */
dmatrix33E	OBBox::getAxes(){
	return(m_axes);
}

/*!
 * Set your target geometry. Not supported volumetric tessellations(type =2).
 * Reimplemented from BaseManipulation::setGeometry().
 */
void	OBBox::setGeometry(MimmoObject * geo){
	if (geo->getType() == 2 )	{
		std::cout<<"WARNING: "<<m_name<<" does not support volumetric tessellation. Geometry not set"<<std::endl;
		return;
	}
	BaseManipulation::setGeometry(geo);
};

/*! Plot the OBB as a structured grid to *vtu file.
 * \param[in] directory output directory
 * \param[in] filename  output filename w/out tag
 * \param[in] counter   integer identifier of the file
 * \param[in] binary    boolean flag for 0-"ascii" or 1-"appended" writing
 */
void		OBBox::plot(std::string directory, std::string filename,int counter, bool binary){
	

	dvecarr3E activeP(8);
	
	activeP[0] =  - 0.5 * m_span;
	activeP[6] =    0.5 * m_span;
	
	activeP[1] = activeP[0]; activeP[1][0] += m_span[0];
	activeP[3] = activeP[0]; activeP[3][1] += m_span[1];
	activeP[2] = activeP[6]; activeP[2][2] += -1.0*m_span[2];
	
	activeP[7] = activeP[6]; activeP[7][0] += -1.0*m_span[0];
	activeP[5] = activeP[6]; activeP[5][1] += -1.0*m_span[1];
	activeP[4] = activeP[0]; activeP[4][2] += m_span[2];
	
	darray3E temp;
	dmatrix33E trasp = bitpit::linearalgebra::transpose(m_axes);
	
	for(auto &val : activeP){
		
		for(int i=0; i<3; ++i){
			temp[i] = dotProduct(val, trasp[i]);
		}
		val = temp + m_origin;
	}
	
	ivector2D activeConn(1);
	for(int i=0; i<8; ++i)	activeConn[0].push_back(i);

	bitpit::VTKFormat codex = bitpit::VTKFormat::ASCII;
	if(binary){codex=bitpit::VTKFormat::APPENDED;}
	bitpit::VTKElementType elDM = bitpit::VTKElementType::HEXAHEDRON;
	bitpit::VTKUnstructuredGrid vtk(directory, filename, elDM);
	vtk.setGeomData( bitpit::VTKUnstructuredField::POINTS, activeP) ;
	vtk.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, activeConn) ;
	vtk.setDimensions(1, 8);
	vtk.setCodex(codex);
	if(counter>=0){vtk.setCounter(counter);}
	
	vtk.write();
};


/*!Execute your object, calculate the OBBox of your geometry. Implementation of pure virtual BaseManipulation::execute
 */
void 		OBBox::execute(){
	if(m_geometry == NULL)	return;
	
	dmatrix33E covariance;
	darray3E etaPoint;
	
	evaluateCovarianceMatrix(covariance, etaPoint);
	
	std::cout<<"mass center"<<std::endl;
	std::cout<<etaPoint<<std::endl;
	std::cout<<"----------------"<<std::endl;
    m_axes = eigenVectors(covariance);	
	
	darray3E pmin, pmax;
	pmin.fill(1.e18);
	pmax.fill(-1.e18);
	double val;
	
	for(auto & vert: getGeometry()->getVertices()){
		darray3E coord = vert.getCoords(); 
		for(int i=0;i<3; ++i){
			val = dotProduct(coord, m_axes[i]);
			pmin[i] = std::fmin(pmin[i], val);
			pmax[i] = std::fmax(pmax[i], val);
		}
	}
	
	m_span = pmax - pmin;
	darray3E originLoc = 0.5*(pmin+pmax);
	dmatrix33E trasp = bitpit::linearalgebra::transpose(m_axes);
	for(int i=0; i<3; ++i){
		m_origin[i] = dotProduct(originLoc, trasp[i]);
	}
	std::cout<<"origin"<<std::endl;
	std::cout<<m_origin<<std::endl;
	std::cout<<"----------------"<<std::endl;
	std::cout<<"span"<<std::endl;
	std::cout<<m_span<<std::endl;
	std::cout<<"----------------"<<std::endl;
	
	
};

/*! 
 * Calculates and returns the eigenVectors of a 3x3 matrix.
 * \param[in] matrix	target matrix
 * \return	matrix of eigenvectors by column
 */
dmatrix33E 		OBBox::eigenVectors( dmatrix33E & matrix){
	
	dmatrix33E result;

// 	matrix[0][0] = 4.0;
// 	matrix[1][0]  = -3.0;
// 	matrix[2][0]  = -3.0;
// 	matrix[0][1]  = 6.0;
// 	matrix[1][1]  = -5.0;
// 	matrix[2][1]  = -6.0;
// 	matrix[0][2]  = 0.0;
// 	matrix[1][2]  = 0.0;
// 	matrix[2][2]  = -5.0;
	
	
		double * a = new double [9];
		double * u = new double [9];
		double * vt = new double [9];
		double * s = new double[3];
		double * superb = new double[3];
		int info;
		
		int k=0;
		for(int i=0; i<3; i++){
			for(int j=0; j<3; j++){
				//transpose
				a[k] = matrix[j][i];
				k++;
			}
		}

		
		info = LAPACKE_dgesvd( LAPACK_COL_MAJOR, 'N','S', 3, 3, a, 3, s, u, 3, vt, 3, superb);
		
		//solution norm
		for (int i=0; i<9; i++){		
			result[i/3][i%3] = vt[i];
		}	
		for(int i=0; i<3; ++i)	result[i] /= norm2(result[i]);

		//result[2] = crossProduct(result[0],result[1]);
		std::cout<<"lapack result"<<std::endl;
		
		std::cout<< "verify normal 0 1 "<< dotProduct(result[0],result[1]) <<std::endl;
		std::cout<< "verify normal 1 2 "<< dotProduct(result[1],result[2]) <<std::endl;
		std::cout<< "verify normal 2 0 "<< dotProduct(result[2],result[0]) <<std::endl;
	
		std::cout<<s[0] <<'\t' << result[0]<<std::endl;
		std::cout<<s[1] << '\t'<< result[1]<<std::endl;
		std::cout<<s[2] << '\t'<< result[2]<<std::endl;
		
		delete [] a; a = NULL;
		delete [] u; u = NULL;
		delete [] vt; vt = NULL;
		delete [] s; s = NULL;
		delete [] superb; superb = NULL;

// 	Eigen::Matrix<double,3,3> A;
// 	for(int i=0; i<3; i++){
// 		for(int j=0; j<3; j++){
// 			A(i,j) = matrix[i][j];
// 		}
// 	}
	
// 	Eigen::EigenSolver<Eigen::Matrix<double,3,3> > sol(A);
// 	
// 	auto v0 = sol.eigenvectors().col(0);
// 	auto v1 = sol.eigenvectors().col(1);
// 	auto v2 = sol.eigenvectors().col(2);
// 	
// 	std::cout<<"eigen result"<<std::endl;
// 	std::cout<<v0<<std::endl;
// 	std::cout<<v1<<std::endl;
// 	std::cout<<v2<<std::endl;
// 	

		return result;
}


/*! 
 * Calculates and returns covariance matrix of the target geometry and its average eta point.
 * The method intrinsecally distinguish between cloud point and tessellation, according to target MimmoObject geometry type.
 * \param[out] covariance matrix
 * \param[out] eta	average Point on the shape
 * \param[out] devP set of the target mesh point, translated by the average eta  of the covariance matrix(OPTIONAL, if null does not evaluate them)
 */
void 		OBBox::evaluateCovarianceMatrix( dmatrix33E & covariance, darray3E & eta){
	
	if(getGeometry()==NULL)	{return;}

	eta.fill(0.0);
	for(auto & val:covariance)	val.fill(0.0);
	int size, counter;
	dvector1D area;
	darray3E temp;
	double areaTot = 0.0;
	double maxVal= 1.e-18; 
	
	if(getGeometry()->getType() == 3){
		size = getGeometry()->getNVertex(); 
		//evaluate eta;
		for(auto & vert: getGeometry()->getVertices()){
			eta += vert.getCoords();
		}
		eta /= (double)size;
		
		//evaluate covariance
		for(auto & vert: getGeometry()->getVertices()){
			temp = vert.getCoords() - eta;
			
			for(int j=0; j<3; ++j){
				for(int k=j; k<3; ++k){
					covariance[j][k] += temp[j]*temp[k];
				}
			}

			for(int j=0; j<3; ++j){
				for(int k=j; k<3; ++k){
					covariance[j][k] /= (double) size;
					maxVal = std::fmax(maxVal, std::abs(covariance[j][k]));
				}
			}
			
		}
		
	}else{
		size = getGeometry()->getNCells(); 
		
		bitpit::SurfUnstructured * tri = static_cast<bitpit::SurfUnstructured * >(getGeometry()->getPatch());
		area.resize(size);
		counter = 0;
		for(auto &cell : tri->getCells()){
			area[counter] = tri->evalCellArea(cell.getId());
			areaTot +=area[counter];
			++counter;
		}
		for(auto & val: area)	val /= areaTot;
		
		//evaluate eta;
		counter = 0;
		darray3E pp;
		long vCount;
		for(auto & cell: tri->getCells()){
			vCount = cell.getVertexCount();
			pp.fill(0.0);
			for(int kk=0; kk<vCount; ++kk){
				pp += tri->getVertexCoords(cell.getVertex(kk)); 
			}
			
			eta += pp*area[counter];
			++counter;		
		}
		
		eta /= (3.0);
		
		
		counter = 0;
		dvecarr3E p2;
		long id;
		for(auto & cell: tri->getCells()){
			vCount = cell.getVertexCount();
			p2.resize(vCount);
			for(int kk=0; kk<vCount; ++kk){
				p2[kk] = tri->getVertexCoords(cell.getVertex(kk)) - eta; 
			}
			
			dmatrix33E tempM = evalCovTriangle(p2);
			for(int j=0; j<3; ++j){
				for(int k=j; k<3; ++k){
					covariance[j][k] += area[counter]*tempM[j][k];
				}
			}
			
			++counter;
		}
		
		
		for(int j=0; j<3; ++j){
			for(int k=j; k<3; ++k){
				covariance[j][k] /= (12.0*size); 
				maxVal = std::fmax(maxVal, std::abs(covariance[j][k]));
			}
		}
	}
	
	for(int j=0; j<3; ++j){
		for(int k=j; k<3; ++k){
			covariance[j][k] /= maxVal; 
		}	
	}
	
	covariance[1][0] = covariance[0][1];
	covariance[2][0] = covariance[0][2];
	covariance[2][1] = covariance[1][2];
	
};	

/*!
 * Evaluate covariance matrix of a triangle, given its 3 vertices.
 * Vertices are supposed to be already expressed w.r.t a target mass center.
 * \param[in]	vv list of 3 triangle vertices
 *\return 	local covariance matrix 
 */

dmatrix33E OBBox::evalCovTriangle(dvecarr3E & vv){
	
	dmatrix33E result;
	for(auto & val: result)	val.fill(0.0);
	
	for(int i=0; i<3; ++i){
		for(int j=0; j<3; ++j){
			dmatrix33E temp  = createMatrix(vv[i],vv[j]);
			
			for(int ki=0; ki<3; ++ki){
				for(int kj=0; kj<3; ++kj){
					result[ki][kj] += temp[ki][kj];
				}
			}
			
		}
	}

	return result;	
};

/*!
 * Perform multiplication of vector V1  and transposed vector V2
 * \return 3x3 matrix
 */
dmatrix33E OBBox::createMatrix(darray3E v1, darray3E v2){
	
	dmatrix33E result;

	for(int i=0; i<3; ++i){
		for(int j=0; j<3; ++j){
			result[i][j] = v1[i]*v2[j];
		}
	}
	
	return result;
};