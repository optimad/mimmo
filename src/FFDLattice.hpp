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
#include "BasicMeshes.hpp"
#include "BaseManipulation.hpp"

/*!
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
class FFDLattice: public BaseManipulation, public UStructMesh {

protected:
	ivector1D m_deg;		/**< Nurbs curve degree for each of the possible 3 direction in space*/
	dvector2D m_knots;		/**< Nurbs curve knots for each of the possible 3 direction in space*/
	ivector2D m_mapEff;		/**< Nurbs map of theoretical node distribution */
	dvector1D m_weights;	/**< Weights of each control node*/
	ivector2D m_mapNodes;	/**< Internal map to access node index w/ knots structure theoretical indexing */
private:
	iarray3E  m_mapdeg;		/**< Map of curves degrees. Increasing order of curves degrees. */
	
public:
	FFDLattice();
	FFDLattice(darray3E &origin, darray3E & span, BasicShape::ShapeType type, ivector1D & dimension);
	FFDLattice(darray3E &origin, darray3E & span, BasicShape::ShapeType type, ivector1D & dimension, ivector1D & degrees);
	FFDLattice(BasicShape * shape, ivector1D & dimension);
	FFDLattice(BasicShape * shape, ivector1D & dimension, ivector1D & degrees);	
	virtual ~FFDLattice();

	//copy operators/constructors
	FFDLattice(const FFDLattice & other);
	FFDLattice & operator=(const FFDLattice & other);
	
	//clean structure;
	void clearLattice();
	
	//internal methods
	ivector1D 	getKnotsDimension();
	dvector1D   getWeights();
	void 		returnKnotsStructure(dvector2D &, ivector2D &);
	void 		returnKnotsStructure( int, dvector1D &, ivector1D &);
		
	void		setDimension(ivector1D &dimensions);
	void		setDimension(ivector1D &dimensions, ivector1D &curveDegrees);
	
	void 		setMesh(darray3E &origin, darray3E & span, BasicShape::ShapeType type, ivector1D & dimensions);
	void 		setMesh(darray3E &origin, darray3E & span, BasicShape::ShapeType type, ivector1D & dimensions, ivector1D & degrees);
	void 		setMesh(BasicShape * shape, ivector1D & dimension);
	void		setMesh(BasicShape * shape, ivector1D & dimension, ivector1D & degrees);
	
	void 		setNodalWeight(double , int );
	void 		setNodalWeight(double , int, int, int);
	
	//plotting wrappers
	void		plotGrid(std::string directory, std::string filename, int counter, bool binary, bool deformed);
	void		plotCloud(std::string directory, std::string filename, int counter, bool binary, bool deformed);
	
	//execute deformation methods
	void		setInfo();
	void 		execute();
	darray3E 	apply(darray3E & point);
	dvecarr3E 	apply(dvecarr3E * point);
	dvecarr3E 	apply(ivector1D & map);
	
protected:
	darray3E convertDisplToXYZ(darray3E &, int i);
	dvecarr3E convertDisplToXYZ();
private:

	//Nurbs Evaluators
	darray3E nurbsEvaluator(darray3E &); 
	dvecarr3E nurbsEvaluator(ivector1D &);
	double nurbsEvaluatorScalar(darray3E &, int);

	//Nurbs utilities
	dvector1D basisITS0(int k, int pos, double coord);	
	void getWorkLoad(int dir, dvector2D & loads, dvector2D & result);
	dvector1D getNodeSpacing(int dir);
	
	//knots mantenaince utilities
	void clearKnots();
	void setKnotsStructure(); 
	void setKnotsStructure(int dir, bool flag); 
	int  	getKnotInterval(double, int);
	double 	getKnotValue(int, int);
	int 	getKnotIndex(int,int);
	int 	getTheoreticalKnotIndex(int,int);

	//nodal displacement utility
	void resizeDisplacements(int, int, int);
	void setMapNodes(int ind);
	int  accessMapNodes(int,int,int);
	//dimension utilities
	void orderDimension();
	
};

/*! Return real global index of a nodal displacement, given its position i,j,k in knots indexing logic*/
inline int FFDLattice::accessMapNodes(int i, int j, int k){
	return(accessPointIndex(m_mapNodes[0][i], m_mapNodes[1][j], m_mapNodes[2][k]));
};

#endif /* __FFDLATTICE_HPP__ */
