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
#include "lapacke.h"

#include <chrono>

using namespace std::chrono;

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
 *	|PortID | PortType    | variable/function                     | DataType		      |
 *	|-------|-------------|---------------------------------------|-----------------------|
 *	| 99    | M_GEOM      | m_geometry                            |(SCALAR, MIMMO_)       |
 *	|-------|-------------|---------------------------------------|-----------------------|
 * 
 *
 *  |---------------------------------------------------------------|
 *	|               Port Output               						|
 *	|-------|-------------|-------------------|---------------------|
 *	|PortID | PortType    | variable/function | DataType		  	|
 *	|-------|-------------|-------------------|---------------------|
 *	| 20    | M_POINT     | getOrigin         |	(ARRAY3, FLOAT)		|
 *	| 22    | M_AXES      | getAxes 	      |	(ARR3ARR3, FLOAT)	|
 *	| 23    | M_SPAN      | getSpan           |	(ARRAY3, FLOAT)		|
 *  | 99    | M_GEOM      | getGeometry       |	(SCALAR, MIMMO_)	|
 *	|-------|-------------|-------------------|---------------------|
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
	built = (built && createPortIn<MimmoObject*, OBBox>(&m_geometry, M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
// creating output ports
	built = (built && createPortOut<darray3E, OBBox>(this, &mimmo::OBBox::getOrigin, M_POINT, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortOut<dmatrix33E, OBBox>(this, &mimmo::OBBox::getAxes, M_AXES, mimmo::pin::containerTAG::ARR3ARR3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortOut<darray3E, OBBox>(this, &mimmo::OBBox::getSpan, M_SPAN, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::FLOAT));
	
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
	dmatrix33E	trasp = bitpit::linearalgebra::transpose(m_axes);
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
	darray3E spectrum;
	darray3E etaPoint;
	
// 	steady_clock::time_point t1,t2,t3,t4,t5;
// 	duration<double> time_span;
	
// 	t1 = steady_clock::now();

	evaluateCovarianceMatrix(covariance, etaPoint);

// 	t2 = steady_clock::now();
	
    m_axes = eigenVectors(covariance, spectrum);	
	
// 	t3 = steady_clock::now();
	
	adjustBasis(m_axes, spectrum);
	
// 	t4 = steady_clock::now();
	
	darray3E pmin, pmax;
	pmin.fill(1.e18);
	pmax.fill(-1.e18);
	double val;
	
	dmatrix33E trasp = bitpit::linearalgebra::transpose(m_axes);
	
	for(auto & vert: getGeometry()->getVertices()){
		darray3E coord = vert.getCoords(); 
		for(int i=0;i<3; ++i){
			val = dotProduct(coord, m_axes[i]);
			pmin[i] = std::fmin(pmin[i], val);
			pmax[i] = std::fmax(pmax[i], val);
		}
	}
	
	m_span = pmax - pmin;
	//check if one of the span goes to 0;
	double avg_span = 0.0;
	for(auto & val: m_span)	avg_span+=val;
	avg_span /= 3.0;
	
	for(auto &val : m_span)	{
		val = std::fmax(val, 1.E-04*avg_span);
	}
	
	darray3E originLoc = 0.5*(pmin+pmax);
	for(int i=0; i<3; ++i){
		m_origin[i] = dotProduct(originLoc, trasp[i]);
	}

	double volOBB = m_span[0]*m_span[1]*m_span[2];	
	
	
	//check the axis aligned bounding box if volume is lesser then obb
	//take it instead of the oriented.
	
	pmin.fill(1.e18);
	pmax.fill(-1.e18);
	for(int i=0; i<3; ++i){
		trasp[i].fill(0.0);
		trasp[i][i] = 1.0;
	}
	
	for(auto & vert: getGeometry()->getVertices()){
		darray3E coord = vert.getCoords(); 
		for(int i=0;i<3; ++i){
			val = dotProduct(coord, trasp[i]);
			pmin[i] = std::fmin(pmin[i], val);
			pmax[i] = std::fmax(pmax[i], val);
		}
	}
	
	darray3E span2 = pmax - pmin;
	darray3E orig = 0.5*(pmin+pmax);
	//check if one of the span goes to 0;
	avg_span = 0.0;
	for(auto & val: span2)	avg_span+=val;
	avg_span /= 3.0;
	
	for(auto &val : span2)	{
		val = std::fmax(val, 1.E-04*avg_span);
	}
	
	double volAAA =span2[0]*span2[1]*span2[2];

	
	if(volAAA <= volOBB){
		for(int i=0; i<3; ++i)	m_axes[i] = trasp[i];
		m_span = span2;
		m_origin = orig;
	}

//	t5 = steady_clock::now();
/*
	time_span = duration_cast<duration<double>>(t2 - t1);
	std::cout << "eval covariance execution took me " << time_span.count() << " seconds."<<std::endl;
	time_span = duration_cast<duration<double>>(t3 - t2);
	std::cout << "eigenvectors execution took me " << time_span.count() << " seconds."<<std::endl;
	time_span = duration_cast<duration<double>>(t4 - t3);
	std::cout << "adjust basis execution took me " << time_span.count() << " seconds."<<std::endl;
	time_span = duration_cast<duration<double>>(t5 - t4);
	std::cout << "compute bbox execution took me " << time_span.count() << " seconds."<<std::endl;
	*/

	if(isPlotInExecution())		plotOptionalResults();

};

/*!
* Plot Optional results of the class, that is the oriented bounding box as *.vtu mesh
*/
void 	OBBox::plotOptionalResults(){
	
	std::string dir = m_outputPlot + "/";
	std::string nameGrid  = m_name;
	plot(dir, nameGrid, getClassCounter(), true );
}



/*! 
 * Calculates and returns the eigenVectors and eigenvalues of a 3x3 matrix.
 * \param[in] matrix	target matrix
 * \param[out]	eigenvalues eigenvalues of the matrix
 * \return	matrix of eigenvectors by column
 */
dmatrix33E 		OBBox::eigenVectors( dmatrix33E & matrix, darray3E & eigenvalues){
	
	dmatrix33E result;
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
		for(int i=0; i<2; ++i)	{
			result[i] /= norm2(result[i]);
			eigenvalues[i] = s[i];
		}	
		result[2] = crossProduct(result[0],result[1]);
		eigenvalues[2] = s[2];
		
		delete [] a; a = NULL;
		delete [] u; u = NULL;
		delete [] vt; vt = NULL;
		delete [] s; s = NULL;
		delete [] superb; superb = NULL;

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

/*!
 * Adjust basis of eigenvectors, given is eigenvalues spectrum. In order to best fit the 
 * 3D shape.  The rules are the following
 * 1) three real coincident eigenvalues, return fundamental axis as eigenvectors
 * 2) three real different eigenvalues, do nothing
 * 3) two coincident eigenvalues, one not, find best fit to shape for eigVec associated to coincident eigenvalues.
 * \param[in,out]	eigVec basis of eigenvectors by rows
 * \param[in]	eigVal eigenvalues associated to eigVec
 */
void 	OBBox::adjustBasis(dmatrix33E & eigVec, darray3E & eigVal){
	
	darray3E diff;
	double tol = 1.0E-3;
	
	diff[0] = std::abs(eigVal[1]-eigVal[0]);
	diff[1] = std::abs(eigVal[2]-eigVal[1]);
	diff[2] = std::abs(eigVal[0]-eigVal[2]);
	
	if(diff[0]>tol && diff[1]>tol && diff[2]>tol)	return;
	if(diff[0]<=tol && diff[1]<=tol && diff[2]<= tol){
		int counter=0;
		for(auto & val : eigVec){
			val.fill(0.0);
			val[counter] = 1.0;
			++counter;
			}
		return;
	}
	
	int guess = 0, third =1, stable=2;
	
	if(diff[1] <= tol)	{
		guess = 1;
		third = 2;
		stable= 0;
	}
	
	if(diff[2] <= tol)	{
		guess = 2;
		third = 0;
		stable= 1;
	}
		
	dmatrix33E axes, trasp;
	int counter=0;
	for(auto & val: eigVec){
		axes[counter] = val;
		++counter;
	}	
	
	darray3E	guessVec = axes[guess];
	darray3E	refVec   = axes[stable];
	
	int nstage = 15;
	int niterations = 8;
	
	double distance = M_PI/(2.0*(nstage));
	double volume, theta, val;
	std::map<double, double>	mapval;
	darray3E pmin, pmax, temp;
	
	for(int i=0; i<nstage; ++i){
		
		theta = i*distance;
		
		axes[guess] = guessVec*std::cos(theta) + std::sin(theta)*crossProduct(refVec, guessVec);
		axes[third] = crossProduct(axes[stable],axes[guess]);
		
		pmin.fill(1.e18);
		pmax.fill(-1.e18);
		
		//trasp = bitpit::linearalgebra::transpose(axes);
		
		for(auto & vert: getGeometry()->getVertices()){
			darray3E coord = vert.getCoords(); 
			for(int i=0;i<3; ++i){
				val = dotProduct(coord, axes[i]); //trasp[i]);
				pmin[i] = std::fmin(pmin[i], val);
				pmax[i] = std::fmax(pmax[i], val);
			}
		}
		
		temp = pmax - pmin;
		volume = temp[0]*temp[1]*temp[2];
		
		mapval[volume] = theta;
	}
	
	int it =0;
	while (it < niterations){
		
		
		distance /= 2.0;
		
		double thetadx = mapval.begin()->second + distance;
		if(thetadx >= M_PI/2.0)	thetadx += -1.0*M_PI/2.0;
		
		double thetasx = mapval.begin()->second - distance;
		if(thetasx <= 0.0)	thetasx += M_PI/2.0;
		
		//evaluate on the right
		axes[guess] = guessVec*std::cos(thetadx) + std::sin(thetadx)*crossProduct(refVec, guessVec);
		axes[third] = crossProduct(axes[stable],axes[guess]);
		
		
		//trasp = bitpit::linearalgebra::transpose(axes);
		pmin.fill(1.e18);
		pmax.fill(-1.e18);
		for(auto & vert: getGeometry()->getVertices()){
			darray3E coord = vert.getCoords(); 
			for(int i=0;i<3; ++i){
				val = dotProduct(coord,axes[i]); // trasp[i]);
				pmin[i] = std::fmin(pmin[i], val);
				pmax[i] = std::fmax(pmax[i], val);
			}
		}
		
		temp = pmax - pmin;
		volume = temp[0]*temp[1]*temp[2];
		
		mapval[volume] = thetadx;

		//evaluate on the left
		axes[guess] = guessVec*std::cos(thetasx) + std::sin(thetasx)*crossProduct(refVec, guessVec);
		axes[third] = crossProduct(axes[stable],axes[guess]);
		
		trasp = bitpit::linearalgebra::transpose(axes);
		pmin.fill(1.e18);
		pmax.fill(-1.e18);
		for(auto & vert: getGeometry()->getVertices()){
			darray3E coord = vert.getCoords(); 
			for(int i=0;i<3; ++i){
				val = dotProduct(coord,axes[i]); // trasp[i]);
				pmin[i] = std::fmin(pmin[i], val);
				pmax[i] = std::fmax(pmax[i], val);
			}
		}
		
		temp = pmax - pmin;
		volume = temp[0]*temp[1]*temp[2];
		
		mapval[volume] = thetasx;
		++it;
	}
	
	theta = mapval.begin()->second;
	eigVec[guess] = guessVec*std::cos(theta) + std::sin(theta)*crossProduct(refVec, guessVec);
	eigVec[third] = crossProduct(eigVec[stable],eigVec[guess]);

	
	
	return;
};