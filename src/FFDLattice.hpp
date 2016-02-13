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
\*---------------------------------------------------------------------------*/
#ifndef __FFDLATTICE_HPP__
#define __FFDLATTICE_HPP__

#include "Operators.hpp"
#include "customOperators.hpp"
#include "CoreSelectionPatches.hpp"
#include "BaseManipulation.hpp"
/*!
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
class FFDLatticeBox: public BaseManipulation, public HullCube {

protected:
	ivector1D m_deg;		/**< Nurbs curve degree for each of the possible 3 direction in space*/
	dvector2D m_knots;		/**< Nurbs curve knots for each of the possible 3 direction in space*/
	ivector2D m_multK;		/**< Nurbs curve knots multiplicity vector, for each direction */
	dvector1D m_weights;	/**< Weights of each control node*/

public:
	FFDLatticeBox();
	FFDLatticeBox(darray3E origin, darray3E span, ivector1D dimension, MimmoObject* geometry, BaseManipulation * parent=NULL);
	FFDLatticeBox(darray3E origin, darray3E span, ivector1D dimension);
	~FFDLatticeBox();

	//copy operators/constructors
	FFDLatticeBox(const FFDLatticeBox & other);
	FFDLatticeBox & operator=(const FFDLatticeBox & other);
	
	//clean structure;
	void cleanAll();
	void cleanKnots();
	
	//internal methods
	ivector1D 	getKnotsDimension();
	dvector1D   getWeights();
	void 		returnKnotsStructure(dvector2D &, ivector2D &);
	void 		returnKnotsStructure( std::string, dvector1D &, ivector1D &);
	ivector1D 	getDimension();
	
	void		setDimension(ivector1D dimension);
	void		setDimension(ivector1D dimension, ivector1D curveDegrees);
	void		setOrigin(darray3E origin);
	void		setSpan(darray3E span);
	void 		setMesh(darray3E origin, double spanX, double spanY, double spanZ, int nx, int ny, int nz);
	
	void 		setNodalWeight(double , int );
	void 		setNodalWeight(double , int, int, int);
	
	//plotting wrappers
	void		plotGrid(std::string directory, std::string filename, int counter, bool ascii, bool deformed);
	void		plotCloud(std::string directory, std::string filename, int counter, bool ascii, bool deformed);
	
	//execute deformation methods
	void 		execute();
	darray3E 	apply(darray3E & point);
	dvecarr3E 	apply(dvecarr3E * point);
	dvecarr3E 	apply(ivector1D & map);
	
	//relationship methods
protected:
	void 		recoverDisplacements();

private:

	//Nurbs Evaluators
	darray3E nurbsEvaluator(darray3E &); 
	double nurbsEvaluatorScalar(darray3E &, int);
	dvector1D getNurbsPoint(int k, dvector1D & basis, dvector2D & loads);
	dvector1D basisITS0(int k, int pos, double coord);	

	//Mesh utility
	void rebaseMesh();
	
	//knots mantenaince utilities
	void setKnotsStructure(); 
	void setKnotsStructure( std::string); 
	int  	getKnotInterval(double, int);
	double 	getKnotValue(int, int);
	int 	getKnotIndex(int,int);
	int 	getTheoreticalKnotIndex(int,int);

	//nodal displacement utility
	void resizeDisplacements(int, int, int);
	
};

#endif /* __FFDLATTICE_HPP__ */
