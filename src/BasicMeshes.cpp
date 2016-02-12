#include "BasicMeshes.hpp"

using namespace std;

//*****************************************************************************************************************************
// BASE_USTRUCTMESH IMPLEMENTATION 
/*
 *	\date			31/12/2015
 *	\authors		Gallizio Federico
 * 	\authors 		Amos Tori
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief Virtual class for 3D uniform structured mesh --> based on ucartmesh paradigm of Bitpit 
 *
 *	Base class for uniform structured grids, suitable for derivation of meshes on pure cartesian, cylindrical or spherical
 *	coordinates system.
 */

/*! Basic Constructor */
BASE_UStructMesh::BASE_UStructMesh(){
	m_nx = 0; m_ny=0; m_nz=0;
	m_dx = 0; m_dy=0; m_dz=0;
};


/*! Basic destructor */
BASE_UStructMesh::~BASE_UStructMesh(){
	m_dx = m_dy = m_dz =0.0;
	m_nx = m_ny = m_nz =0;
	m_xnode.clear(); m_ynode.clear(); m_znode.clear();
	m_xedge.clear(); m_yedge.clear(); m_zedge.clear();
}

/*! Copy Constructor.
 * \param[in] other BASE_UStructMesh object where copy from
 */
BASE_UStructMesh::BASE_UStructMesh(const BASE_UStructMesh & other){
	// Number of cells
	m_nx = other.m_nx;
	m_ny = other.m_ny;
	m_nz = other.m_nz;
	
	// Cell Spacing
	m_dx = other.m_dx;
	m_dy = other.m_dy;
	m_dz = other.m_dz;
	
	// Mesh Origin & span	
	m_origin = other.m_origin;
	m_span = other.m_span;
	// Resize mesh data structure ----------------------------------------------- //
	resizeMesh();
	
	// Copy cell edges and cell centers ----------------------------------------- //
	m_xnode = other.m_xnode;
	m_ynode = other.m_ynode;
	m_znode = other.m_znode;
	m_xedge = other.m_xedge;
	m_yedge = other.m_yedge;
	m_zedge = other.m_zedge;
};

/*! Copy Operator.
 * \param[in] other BASE_UStructMesh object where copy from
 */
BASE_UStructMesh & BASE_UStructMesh::operator=(const BASE_UStructMesh & other){
	
	// Number of cells
	m_nx = other.m_nx;
	m_ny = other.m_ny;
	m_nz = other.m_nz;
	
	// Cell Spacing
	m_dx = other.m_dx;
	m_dy = other.m_dy;
	m_dz = other.m_dz;
	
	// Mesh Origin & span	
	m_origin = other.m_origin;
	m_span = other.m_span;
	// Resize mesh data structure ----------------------------------------------- //
	resizeMesh();
	
	// Copy cell edges and cell centers ----------------------------------------- //
	m_xnode = other.m_xnode;
	m_ynode = other.m_ynode;
	m_znode = other.m_znode;
	m_xedge = other.m_xedge;
	m_yedge = other.m_yedge;
	m_zedge = other.m_zedge;
	
	return(*this); 
};

/*!Clear the Mesh structure. Resize your nodal vector list to zero, but not destroy them.*/
void BASE_UStructMesh::clearMesh(){
	m_nx = 0; m_ny=0; m_nz =0;
	m_dx=0; m_dy=0; m_dz=0;
	
	m_origin.fill(0.0);
	m_span.fill(0.0);
	
	reshapeNodalStructure();
};  

/*! Resize the nodal vector lists to current mesh dimensions, but don't destroy them. Please use reshape and destroy methods to handle them */
void BASE_UStructMesh::resizeMesh(){
	// Cell centers
	m_xnode.resize(m_nx, 0.0);
	m_ynode.resize(m_ny, 0.0);
	m_znode.resize(m_nz, 0.0);
	
	// Points
	m_xedge.resize(m_nx+1, 0.0);
	m_yedge.resize(m_ny+1, 0.0);
	m_zedge.resize(m_nz+1, 0.0);
};

/*! Destroy the all nodal structures of the mesh. */
void BASE_UStructMesh::destroyNodalStructure(){
	freeContainer(m_xnode);
	freeContainer(m_ynode);
	freeContainer(m_znode);
	freeContainer(m_xedge);
	freeContainer(m_yedge);
	freeContainer(m_zedge);  
};

/*! Destroy the all nodal structures of the mesh, and reinitialize them to current mesh dimensions. */
void BASE_UStructMesh::reshapeNodalStructure(){
	destroyNodalStructure();
	resizeMesh();
};


/*! Translate your mesh in the 3D space 
 * \param[in] tx x coordinate translation
 * \param[in] ty y coordinate translation
 * \param[in] tz z coordinate translation
 */
void BASE_UStructMesh::translateMesh(double tx, double ty, double tz){
	m_origin[0]+= tx;
	m_origin[1]+= ty;
	m_origin[2]+= tz;
};


/*! Return current mesh origin */
darray3E BASE_UStructMesh::getOrigin(){return(m_origin);};

/*! Return current mesh span */
darray3E BASE_UStructMesh::getSpan(){return(m_span);};

/*! Return current mesh spacing */
darray3E BASE_UStructMesh::getSpacing(){
	darray3E res;
	res[0] = m_dx; res[1] =m_dy; res[2] = m_dz;
	return(res); 
};

/*! Return current dimension of the mesh (number of cells in each direction) */
ivector1D BASE_UStructMesh::getDimension(){
	ivector1D res(3,0);
	res[0] = m_nx; res[1] =m_ny; res[2] = m_nz;
	return(res); 
};

/*! Get n-th center cell coordinates in local mesh reference frame, 
 * given its global cell index on the mesh. 
 * \param[in] index cell index in the global nodal list.
 */
darray3E BASE_UStructMesh::getGridCCell(int index){
	
	int i,j,k;
	accessCellData(index, i,j,k);
	darray3E res = getGridCCell(i,j,k);
	
	return(res);
};
/*! Get n-th center cell coordinates in local mesh reference frame, 
 * given its cartesian cell indices on the mesh. 
 * \param[in] i_ x cartesian index.
 * \param[in] j_ y cartesian index.
 * \param[in] k_ z cartesian index.
 */
darray3E BASE_UStructMesh::getGridCCell(int i_, int j_, int k_){
	
	darray3E res{0,0,0};
	res[0] = m_xnode[i_];
	res[1] = m_ynode[j_];
	res[2] = m_znode[k_];
	
	return(res);
};

/*! Get n-th nodal vertex coordinates in local mesh reference frame, 
 * given its global point index on the mesh. 
 * \param[in] index point index in the global nodal list.
 */
darray3E BASE_UStructMesh::getGridPoint(int index){
	
	int i,j,k;
	accessPointData(index, i,j,k);
	darray3E res = getGridPoint(i,j,k);
	
	return(res);
};
/*! Get n-th nodal vertex coordinates in local reference frame, 
 * given its cartesian point indices on the mesh. 
 * \param[in] i_ x cartesian index.
 * \param[in] j_ y cartesian index.
 * \param[in] k_ z cartesian index.
 */
darray3E BASE_UStructMesh::getGridPoint(int i_, int j_, int k_){
	
	darray3E res{0,0,0};
	res[0] = m_xedge[i_];
	res[1] = m_yedge[j_];
	res[2] = m_zedge[k_];
	
	return(res);
};    

/*! Get n-th center cell coordinates in global absolute reference frame, 
 * given its global cell index on the mesh. 
 * \param[in] index cell index in the global nodal list.
 */
darray3E BASE_UStructMesh::getGlobalCCell(int index){
	
	darray3E res = getGridCCell(index);
	return(transfToGlobal(res));
};
/*! Get n-th center cell coordinates in global absolute reference frame, 
 * given its cartesian cell indices on the mesh. 
 * \param[in] i_ x cartesian index.
 * \param[in] j_ y cartesian index.
 * \param[in] k_ z cartesian index.
 */
darray3E BASE_UStructMesh::getGlobalCCell(int i_, int j_, int k_){
	
	darray3E res = getGridCCell(i_,j_,k_);
	return(transfToGlobal(res));
};

/*! Get n-th nodal vertex coordinates in global absolute reference frame, 
 * given its global point index on the mesh. 
 * \param[in] index point index in the global nodal list.
 */
darray3E BASE_UStructMesh::getGlobalPoint(int index){
	
	darray3E res = getGridPoint(index);
	return(transfToGlobal(res));
};
/*! Get n-th nodal vertex coordinates in global absolute reference frame, 
 * given its cartesian point indices on the mesh. 
 * \param[in] i_ x cartesian index.
 * \param[in] j_ y cartesian index.
 * \param[in] k_ z cartesian index.
 */
darray3E BASE_UStructMesh::getGlobalPoint(int i_, int j_, int k_){
	
	darray3E res = getGridPoint(i_,j_,k_);
	return(transfToGlobal(res));
};    

/*! Get neighbor vertices (by their global indices and in VTK-hexahedra ordered) of a given cell
 * \param[in] i_ x coordinate cell index
 * \param[in] j_ y coordinate cell index
 * \param[in] k_ z coordinate cell index
 */ 
ivector1D BASE_UStructMesh::getCellNeighs(int i_, int j_, int k_){
	
	ivector1D result(8);
	
	result[0] = accessPointData(i_, j_, k_);
	result[1] = accessPointData(i_+1, j_, k_);
	result[2] = accessPointData(i_+1, j_+1, k_);
	result[3] = accessPointData(i_, j_+1, k_);
	result[4] = accessPointData(i_, j_, k_+1);
	result[5] = accessPointData(i_+1, j_, k_+1);
	result[6] = accessPointData(i_+1, j_+1, k_+1);
	result[7] = accessPointData(i_, j_+1, k_+1);
	
	return(result);  
}

/*! Get neighbor vertices (by their global indices and in VTK-hexahedra ordered) of a given cell
 * \param[in] index cell global index
 */ 
ivector1D BASE_UStructMesh::getCellNeighs(int index){
	
	ivector1D pp(3,0);
	accessCellData(index, pp[0],pp[1], pp[2]);
	return(getCellNeighs(pp[0],pp[1],pp[2]));
}


/*! Return cartesian indices of the cell containing the target point in global reference frame
 * \param[in] point 3D coordinate of target point
 * \param[out] i x cell index
 * \param[out] j y cell index
 * \param[out] k z cell index 
 */ 
void BASE_UStructMesh::returnCellID(darray3E & point, int &i, int &j, int &k){
	
	darray3E P = transfToLocal(point);
	i = min(m_nx-1, max(0, (int) floor((P[0])/m_dx)));
	j = min(m_ny-1, max(0, (int) floor((P[1])/m_dy)));
	k = min(m_nz-1, max(0, (int) floor((P[2])/m_dz)));
};

/*! Return cartesian indices of the cell containing the target point in global reference framne
 * \param[in] point 3D coordinate of target point
 * \param[out] i x cell index
 * \param[out] j y cell index
 * \param[out] k z cell index 
 */ 
void BASE_UStructMesh::returnCellID(dvector1D & point, int &i, int &j, int &k){
	
	dvector1D P = transfToLocal(point);
	i = min(m_nx-1, max(0, (int) floor((P[0])/m_dx)));
	j = min(m_ny-1, max(0, (int) floor((P[1])/m_dy)));
	k = min(m_nz-1, max(0, (int) floor((P[2])/m_dz)));
};


/*! Return global index of the cell given its cartesian indices. Follows the ordering sequences z-y-x
 * \param[in] i x cartesian index
 *\param[in] j y cartesian index
 *\param[in] k z cartesian index
 *\param[out] result global index 
 */
int BASE_UStructMesh::accessCellData(int i, int j, int k){
	
	int index = m_ny * m_nz * i + m_nz * j + k;
	return(index);
};

/*! Return cartesian indices of the cell given its global index. Follows the ordering sequences z-y-x
 * \param[in] N_ global index 
 *\param[out] i x cartesian index
 *\param[out] j y cartesian index
 *\param[out] k z cartesian index
 */
void BASE_UStructMesh::accessCellData(int N_, int & i, int & j, int & k){
	k = N_ % m_nz;
	int index = N_ / m_nz;
	j = index % m_ny;
	i = index / m_ny; 
};

/*! Return global index of the point given its cartesian indices. Follows the ordering sequences z-y-x
 * \param[in] i x cartesian index
 *\param[in] j y cartesian index
 *\param[in] k z cartesian index
 *\param[out] result global index 
 */
int  BASE_UStructMesh::accessPointData(int i, int j, int k){
	
	int index = (m_ny+1) * (m_nz+1) * i + (m_nz+1) * j + k;
	return(index);
};

/*! Return cartesian indices of the point given its global index. Follows the ordering sequences z-y-x
 * \param[in] N_ global index 
 *\param[out] i x cartesian index
 *\param[out] j y cartesian index
 *\param[out] k z cartesian index
 */
void BASE_UStructMesh::accessPointData(int N_, int &i, int &j, int &k){
	
	k = N_ % (m_nz+1);
	int index = N_ / (m_nz+1);
	j = index % (m_ny+1);
	i = index / (m_ny+1); 
};


/*! Interpolate value of a given data field on a target point inside the mesh
 * \param[in] point target point
 * \param[in] celldata data field defined on centercells
 * \param[out] interpResult interpolated value
 */
double BASE_UStructMesh::interpolateCellData(darray3E & point, dvector1D & celldata){
	
	int i0, j0, k0, ip, jp, kp;
	double wx0,wx1,wy0,wy1,wz0,wz1;
	
	returnCellID(point, i0, j0, k0);
	darray3E P = transfToLocal(point);
	if (P[0] > m_xnode[i0]) 	{ ip = min(i0+1, m_nx-1);   }
	else                  	{ ip = max(0, i0-1);      }
	if (P[1] > m_ynode[j0]) 	{ jp = min(j0+1, m_ny-1);   }
	else                  	{ jp = max(0, j0-1);      }
	if (P[2] > m_znode[k0]) 	{ kp = min(k0+1, m_nz-1);   }
	else                  	{ kp = max(0, k0-1);      }
	
	// Interpolation weights
	wx1 = max(0.0, min(1.0, abs((P[0] - m_xnode[i0])/m_dx)));     wx0 = 1.0 - wx1;
	wy1 = max(0.0, min(1.0, abs((P[1] - m_ynode[j0])/m_dy)));     wy0 = 1.0 - wy1;
	wz1 = max(0.0, min(1.0, abs((P[2] - m_znode[k0])/m_dz)));     wz0 = 1.0 - wz1;
	
	// Interpolation
	double result  = wz0 * wx0 * wy0 * celldata[accessCellData(i0,j0,k0)]
	+ wz0 * wx0 * wy1 * celldata[accessCellData(i0,jp,k0)]
	+ wz0 * wx1 * wy0 * celldata[accessCellData(ip,j0,k0)]
	+ wz0 * wx1 * wy1 * celldata[accessCellData(ip,jp,k0)]
	+ wz1 * wx0 * wy0 * celldata[accessCellData(i0,j0,kp)]
	+ wz1 * wx0 * wy1 * celldata[accessCellData(i0,jp,kp)]
	+ wz1 * wx1 * wy0 * celldata[accessCellData(ip,j0,kp)]
	+ wz1 * wx1 * wy1 * celldata[accessCellData(ip,jp,kp)];
	
	return(result);      
};

/*! Interpolate value of a given data field on a target point inside the mesh
 * \param[in] point target point
 * \param[in] celldata data field defined on centercells
 * \param[out] interpResult interpolated value
 */
int BASE_UStructMesh::interpolateCellData(darray3E & point, ivector1D & celldata){
	
	int i0, j0, k0, ip, jp, kp;
	double wx0,wx1,wy0,wy1,wz0,wz1;
	
	returnCellID(point, i0, j0, k0);
	darray3E P = transfToLocal(point);
	if (P[0] > m_xnode[i0]) 	{ ip = min(i0+1, m_nx-1);   }
	else                  	{ ip = max(0, i0-1);      }
	if (P[1] > m_ynode[j0]) 	{ jp = min(j0+1, m_ny-1);   }
	else                  	{ jp = max(0, j0-1);      }
	if (P[2] > m_znode[k0]) 	{ kp = min(k0+1, m_nz-1);   }
	else                  	{ kp = max(0, k0-1);      }
	
	// Interpolation weights
	wx1 = max(0.0, min(1.0, abs((P[0] - m_xnode[i0])/m_dx)));     wx0 = 1.0 - wx1;
	wy1 = max(0.0, min(1.0, abs((P[1] - m_ynode[j0])/m_dy)));     wy0 = 1.0 - wy1;
	wz1 = max(0.0, min(1.0, abs((P[2] - m_znode[k0])/m_dz)));     wz0 = 1.0 - wz1;
	
	// Interpolation
	double result  = 	wz0 * wx0 * wy0 * (double)celldata[accessCellData(i0,j0,k0)]
	+ wz0 * wx0 * wy1 * (double)celldata[accessCellData(i0,jp,k0)]
	+ wz0 * wx1 * wy0 *(double) celldata[accessCellData(ip,j0,k0)]
	+ wz0 * wx1 * wy1 *(double) celldata[accessCellData(ip,jp,k0)]
	+ wz1 * wx0 * wy0 *(double) celldata[accessCellData(i0,j0,kp)]
	+ wz1 * wx0 * wy1 *(double) celldata[accessCellData(i0,jp,kp)]
	+ wz1 * wx1 * wy0 *(double) celldata[accessCellData(ip,j0,kp)]
	+ wz1 * wx1 * wy1 *(double) celldata[accessCellData(ip,jp,kp)]; 
	
	int result2 =std::floor(result+0.5); 		
	return(result2);
};

/*! Interpolate value of a given data field on a target point inside the mesh
 * \param[in] point target point
 * \param[in] celldata data field defined on centercells
 * \param[out] interpResult interpolated value
 */
darray3E BASE_UStructMesh::interpolateCellData(darray3E & point, dvecarr3E & celldata){
	
	int i0, j0, k0, ip, jp, kp;
	double wx0,wx1,wy0,wy1,wz0,wz1;
	
	returnCellID(point, i0, j0, k0);
	darray3E P = transfToLocal(point);
	if (P[0] > m_xnode[i0]) 	{ ip = min(i0+1, m_nx-1);   }
	else                  	{ ip = max(0, i0-1);      }
	if (P[1] > m_ynode[j0]) 	{ jp = min(j0+1, m_ny-1);   }
	else                  	{ jp = max(0, j0-1);      }
	if (P[2] > m_znode[k0]) 	{ kp = min(k0+1, m_nz-1);   }
	else                  	{ kp = max(0, k0-1);      }
	
	// Interpolation weights
	wx1 = max(0.0, min(1.0, abs((P[0] - m_xnode[i0])/m_dx)));     wx0 = 1.0 - wx1;
	wy1 = max(0.0, min(1.0, abs((P[1] - m_ynode[j0])/m_dy)));     wy0 = 1.0 - wy1;
	wz1 = max(0.0, min(1.0, abs((P[2] - m_znode[k0])/m_dz)));     wz0 = 1.0 - wz1;
	
	// Interpolation
	darray3E result = wz0 * wx0 * wy0 * celldata[accessCellData(i0,j0,k0)]
	+ wz0 * wx0 * wy1 * celldata[accessCellData(i0,jp,k0)]
	+ wz0 * wx1 * wy0 * celldata[accessCellData(ip,j0,k0)]
	+ wz0 * wx1 * wy1 * celldata[accessCellData(ip,jp,k0)]
	+ wz1 * wx0 * wy0 * celldata[accessCellData(i0,j0,kp)]
	+ wz1 * wx0 * wy1 * celldata[accessCellData(i0,jp,kp)]
	+ wz1 * wx1 * wy0 * celldata[accessCellData(ip,j0,kp)]
	+ wz1 * wx1 * wy1 * celldata[accessCellData(ip,jp,kp)];
	
	return(result);	      
};

/*! Interpolate value of a given data field on a target point inside the mesh
 * \param[in] point target point
 * \param[in] pointdata data field defined on grid points
 * \param[out] interpResult interpolated value
 */
double BASE_UStructMesh::interpolatePointData(darray3E & point, dvector1D & pointdata){
	
	int i0, j0, k0, ip, jp, kp;
	double wx0,wx1,wy0,wy1,wz0,wz1;
	
	darray3E P = transfToLocal(point);
	i0 = max(0, min(m_nx, (int) floor((P[0])/m_dx)));
	j0 = max(0, min(m_ny, (int) floor((P[1])/m_dy)));
	k0 = max(0, min(m_nz, (int) floor((P[2])/m_dz)));
	
	if (P[0] >= m_xedge[i0]) 	{ ip = min(i0+1, m_nx);   }
	else                  	{ ip = max(0, i0-1);      }
	if (P[1] >= m_yedge[j0]) 	{ jp = min(j0+1, m_ny);   }
	else                  	{ jp = max(0, j0-1);      }
	if (P[2] >= m_zedge[k0]) 	{ kp = min(k0+1, m_nz);   }
	else                  	{ kp = max(0, k0-1);      }
	
	// Interpolation weights
	wx1 = max(0.0, min(1.0, abs((P[0] - m_xedge[i0])/m_dx)));     wx0 = 1.0 - wx1;
	wy1 = max(0.0, min(1.0, abs((P[1] - m_yedge[j0])/m_dy)));     wy0 = 1.0 - wy1;
	wz1 = max(0.0, min(1.0, abs((P[2] - m_zedge[k0])/m_dz)));     wz0 = 1.0 - wz1;
	
	// Interpolation
	double result =  wz0 * wx0 * wy0 * pointdata[accessPointData(i0,j0,k0)]
	+ wz0 * wx0 * wy1 * pointdata[accessPointData(i0,jp,k0)]
	+ wz0 * wx1 * wy0 * pointdata[accessPointData(ip,j0,k0)]
	+ wz0 * wx1 * wy1 * pointdata[accessPointData(ip,jp,k0)]
	+ wz1 * wx0 * wy0 * pointdata[accessPointData(i0,j0,kp)]
	+ wz1 * wx0 * wy1 * pointdata[accessPointData(i0,jp,kp)]
	+ wz1 * wx1 * wy0 * pointdata[accessPointData(ip,j0,kp)]
	+ wz1 * wx1 * wy1 * pointdata[accessPointData(ip,jp,kp)];
	return(result);
};

/*! Interpolate value of a given data field on a target point inside the mesh
 * \param[in] point target point
 * \param[in] celldata data field defined on grid points
 * \param[out] interpResult interpolated value
 */
int BASE_UStructMesh::interpolatePointData(darray3E & point, ivector1D & pointdata){
	
	int i0, j0, k0, ip, jp, kp;
	double wx0,wx1,wy0,wy1,wz0,wz1;
	
	darray3E P = transfToLocal(point);
	i0 = max(0, min(m_nx, (int) floor((P[0] - m_origin[0])/m_dx)));
	j0 = max(0, min(m_ny, (int) floor((P[1] - m_origin[1])/m_dy)));
	k0 = max(0, min(m_nz, (int) floor((P[2] - m_origin[2])/m_dz)));
	
	if (P[0] >= m_xedge[i0]) 	{ ip = min(i0+1, m_nx);   }
	else                  	{ ip = max(0, i0-1);      }
	if (P[1] >= m_yedge[j0]) 	{ jp = min(j0+1, m_ny);   }
	else                  	{ jp = max(0, j0-1);      }
	if (P[2] >= m_zedge[k0]) 	{ kp = min(k0+1, m_nz);   }
	else                  	{ kp = max(0, k0-1);      }
	
	// Interpolation weights
	wx1 = max(0.0, min(1.0, abs((P[0] - m_xedge[i0])/m_dx)));     wx0 = 1.0 - wx1;
	wy1 = max(0.0, min(1.0, abs((P[1] - m_yedge[j0])/m_dy)));     wy0 = 1.0 - wy1;
	wz1 = max(0.0, min(1.0, abs((P[2] - m_zedge[k0])/m_dz)));     wz0 = 1.0 - wz1;
	
	// Interpolation
	double result = wz0 * wx0 * wy0 * (double)pointdata[accessPointData(i0,j0,k0)]
	+ wz0 * wx0 * wy1 * (double)pointdata[accessPointData(i0,jp,k0)]
	+ wz0 * wx1 * wy0 *(double) pointdata[accessPointData(ip,j0,k0)]
	+ wz0 * wx1 * wy1 *(double) pointdata[accessPointData(ip,jp,k0)]
	+ wz1 * wx0 * wy0 *(double) pointdata[accessPointData(i0,j0,kp)]
	+ wz1 * wx0 * wy1 *(double) pointdata[accessPointData(i0,jp,kp)]
	+ wz1 * wx1 * wy0 *(double) pointdata[accessPointData(ip,j0,kp)]
	+ wz1 * wx1 * wy1 *(double) pointdata[accessPointData(ip,jp,kp)];
	
	int result2 = std::floor(result+0.5);			
	return(result2);
};
/*! Interpolate value of a given data field on a target point inside the mesh
 * \param[in] point target point
 * \param[in] celldata data field defined on grid points
 * \param[out] interpResult interpolated value
 */
darray3E BASE_UStructMesh::interpolatePointData(darray3E & point, dvecarr3E & pointdata){
	
	int i0, j0, k0, ip, jp, kp;
	double wx0,wx1,wy0,wy1,wz0,wz1;
	
	darray3E P = transfToLocal(point);
	i0 = max(0, min(m_nx, (int) floor((P[0] - m_origin[0])/m_dx)));
	j0 = max(0, min(m_ny, (int) floor((P[1] - m_origin[1])/m_dy)));
	k0 = max(0, min(m_nz, (int) floor((P[2] - m_origin[2])/m_dz)));
	
	if (P[0] >= m_xedge[i0]) 	{ ip = min(i0+1, m_nx);   }
	else                  	{ ip = max(0, i0-1);      }
	if (P[1] >= m_yedge[j0]) 	{ jp = min(j0+1, m_ny);   }
	else                  	{ jp = max(0, j0-1);      }
	if (P[2] >= m_zedge[k0]) 	{ kp = min(k0+1, m_nz);   }
	else                  	{ kp = max(0, k0-1);      }
	
	// Interpolation weights
	wx1 = max(0.0, min(1.0, abs((P[0] - m_xedge[i0])/m_dx)));     wx0 = 1.0 - wx1;
	wy1 = max(0.0, min(1.0, abs((P[1] - m_yedge[j0])/m_dy)));     wy0 = 1.0 - wy1;
	wz1 = max(0.0, min(1.0, abs((P[2] - m_zedge[k0])/m_dz)));     wz0 = 1.0 - wz1;
	
	// Interpolation
	darray3E result = wz0 * wx0 * wy0 * pointdata[accessPointData(i0,j0,k0)]
	+ wz0 * wx0 * wy1 * pointdata[accessPointData(i0,jp,k0)]
	+ wz0 * wx1 * wy0 * pointdata[accessPointData(ip,j0,k0)]
	+ wz0 * wx1 * wy1 * pointdata[accessPointData(ip,jp,k0)]
	+ wz1 * wx0 * wy0 * pointdata[accessPointData(i0,j0,kp)]
	+ wz1 * wx0 * wy1 * pointdata[accessPointData(i0,jp,kp)]
	+ wz1 * wx1 * wy0 * pointdata[accessPointData(ip,j0,kp)]
	+ wz1 * wx1 * wy1 * pointdata[accessPointData(ip,jp,kp)];
	
	return(result);
};


/*! Write your grid as a point cloud in a *.vtu file. 
 * \param[in] folder output directory path
 * \param[in] outfile output namefile w/out tag
 * \param[in] counterFile integer to mark output with a counter number
 * \param[in] codexFlag boolean to distinguish between "ascii" false, and "appended" true
 * \param[in] extPoints OPTIONAL. If defined, use the current vertexList and write them as point cloud, provided that coherent dimensions with the mesh are set. 
 *			If Not, write the current mesh vertices. 
 * */
void BASE_UStructMesh::plotCloud( std::string & folder , std::string outfile, int counterFile, bool codexFlag, dvecarr3E * extPoints){
	
	std::string codex = "ascii";
	if(codexFlag){codex="appended";}
	
	ivector1D dim = getDimension();
	int sizeTot = (dim[0]+1)*(dim[1]+1)*(dim[2]+1);
	
	VTK_BASICCLOUD handle_vtk_output(folder, outfile, codex, sizeTot);
	
	dvecarr3E activeP;
	if(extPoints != NULL && extPoints->size() == sizeTot){activeP = *extPoints; }
	else{
		activeP.resize(sizeTot);
		for(int i=0; i<sizeTot; i++){
			activeP[i] = getGlobalPoint(i);
		}
	}
	
	if(counterFile>=0){handle_vtk_output.SetCounter(counterFile);}
	handle_vtk_output.SetGeomTypes("Float64","Int32", "Int32", "Int32");
	handle_vtk_output.linkData(activeP);
	handle_vtk_output.Write();
};

/*! Write your grid as a point cloud in a *.vtu file. 
 * \param[in] folder output directory path
 * \param[in] outfile output namefile w/out tag
 * \param[in] counterFile integer to mark output with a counter number
 * \param[in] codexFlag boolean to distinguish between "ascii" false, and "appended" true
 * \param[in] vertexList list of global indices of selected vertices that the user wants to write on file. 
 * \param[in] extPoints OPTIONAL. If defined, use the current vertexList and write them as point cloud, provided that coherent dimensions with the mesh are set. 
 *			If Not, write the current mesh vertices. 
 * */
void BASE_UStructMesh::plotCloud( std::string & folder , std::string outfile, int counterFile, bool codexFlag, ivector1D & vertexList, dvecarr3E * extPoints){
	
	std::string codex = "ascii";
	if(codexFlag){codex="appended";}
	
	ivector1D dim = getDimension();
	int sizeTot = (dim[0]+1)*(dim[1]+1)*(dim[2]+1);
	int sizeMap = std::min((int)vertexList.size(), sizeTot);
	
	VTK_BASICCLOUD handle_vtk_output(folder, outfile, codex, sizeMap);
	
	dvecarr3E activeP(sizeMap);
	if(extPoints != NULL && extPoints->size() == sizeTot){
		for(int i=0; i<sizeMap; i++){
			activeP[i] = (*extPoints)[vertexList[i]];
		}  
	}
	else{
		for(int i=0; i<sizeMap; i++){
			activeP[i] = getGlobalPoint(vertexList[i]);
		}
	}
	
	if(counterFile>=0){handle_vtk_output.SetCounter(counterFile);}
	handle_vtk_output.SetGeomTypes("Float64","Int32", "Int32", "Int32");
	handle_vtk_output.linkData(activeP);
	handle_vtk_output.Write();
	
};

/*! Write your grid as a hexahedrical one in a *.vtu file. 
 * \param[in] folder output directory path
 * \param[in] outfile output namefile w/out tag
 * \param[in] counterFile integer to mark output with a counter number
 * \param[in] codexFlag boolean to distinguish between "ascii" false, and "appended" true
 * \param[in] extPoints OPTIONAL. If defined, use the current vertexList and write them as point cloud, provided that coherent dimensions with the mesh are set. 
 *			If Not, write the current mesh vertices. 
 * */
void BASE_UStructMesh::plotGrid(std::string & folder, std::string outfile , int counterFile, bool codexFlag, dvecarr3E * extPoints){
	
	std::string codex = "ascii";
	if(codexFlag){codex="appended";}
	
	ivector1D dim = getDimension();
	int sizePt = (dim[0]+1)*(dim[1]+1)*(dim[2]+1);
	int sizeCl = dim[0]*dim[1]*dim[2];
	VTK_BASICMESH handle_vtk_output(folder, outfile, codex, sizePt, sizeCl, 8*sizeCl);
	
	dvecarr3E activeP(sizePt);
	ivector2D activeConn(sizeCl, ivector1D(8,0));
	
	if(extPoints != NULL && extPoints->size() == sizePt){activeP = *extPoints; }
	else{
		for(int i=0; i<sizePt; i++){
			activeP[i] = getGlobalPoint(i);
		}
	}
	
	for(int i=0; i<sizeCl; ++i){
		activeConn[i] = getCellNeighs(i); 
	}
	
	if(counterFile>=0){handle_vtk_output.SetCounter(counterFile);}
	handle_vtk_output.SetGeomTypes("Float64","Int32", "Int32", "Int32");
	handle_vtk_output.linkData(activeP,activeConn);
	handle_vtk_output.Write();
};

/*! Write your grid as a hexahedrical one in a *.vtu file. 
 * \param[in] folder output directory path
 * \param[in] outfile output namefile w/out tag
 * \param[in] counterFile integer to mark output with a counter number
 * \param[in] codexFlag boolean to distinguish between "ascii" false, and "appended" true
 * \param[in] vertexList list of global indices of selected cells that the user wants to write on file. 
 * \param[in] extPoints OPTIONAL. If defined, use the current vertexList and write them as point cloud, provided that coherent dimensions with the mesh are set. 
 *			If Not, write the current mesh vertices. 
 * */

void BASE_UStructMesh::plotGrid(std::string & folder, std::string outfile, int counterFile, bool codexFlag, ivector1D & cellList, dvecarr3E * extPoints){
	
	std::string codex = "ascii";
	if(codexFlag){codex="appended";}
	
	ivector1D dim = getDimension();
	int sizePt = (dim[0]+1)*(dim[1]+1)*(dim[2]+1);
	int sizeCl = cellList.size();
	
	if(sizeCl > dim[0]*dim[1]*dim[2]){return;}
	
	std::map<int,int> mapPoints;
	ivector2D activeConn(sizeCl, ivector1D(8,0));
	
	for(int i=0; i<sizeCl; ++i){
		activeConn[i] = getCellNeighs(i);
		for(int j=0; j<activeConn[i].size(); ++j){
			mapPoints[activeConn[i][j]] = activeConn[i][j];
		}
	}
	
	dvecarr3E activeP(mapPoints.size());
	ivector1D listP(mapPoints.size());
	int counter = 0;
	std::map<int,int>::iterator itF;
	if(extPoints != NULL && extPoints->size() == sizePt){
		for(itF = mapPoints.begin(); itF != mapPoints.end(); ++itF){
			int pos = itF->second;
			activeP[counter] = (*extPoints)[pos];
			listP[counter] = pos;
			++counter;
		}
	}
	else{
		for(itF = mapPoints.begin(); itF != mapPoints.end(); ++itF){
			int pos = itF->second;
			activeP[counter] = getGlobalPoint(pos);
			listP[counter]=  pos;
			++counter;
		}
	}
	
	//update connectivity w/ local indexing of vertex;
	
	for(int i=0; i<sizeCl; ++i){
		for(int j=0; j<activeConn[i].size(); ++j){
			int posV =posVectorFind(listP, activeConn[i][j]);
			//update local connectivity
			activeConn[i][j] = posV;
		}
	}
	
	VTK_BASICMESH handle_vtk_output(folder, outfile, codex, counter, sizeCl, 8*sizeCl);
	if(counterFile>=0){handle_vtk_output.SetCounter(counterFile);}
	handle_vtk_output.SetGeomTypes("Float64","Int32", "Int32", "Int32");
	handle_vtk_output.linkData(activeP,activeConn);
	handle_vtk_output.Write();
};




//*****************************************************************************************************************************
// UCUBICMESH IMPLEMENTATION 
/*
 *	\date			31/12/2015
 *	\authors		Gallizio Federico
 * 	\authors 		Amos Tori
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief Class for 3D uniform Cartesian Mesh 
 *
 *	Derived class from UStructMesh for uniform cartesian grids.
 */

/*! Basic Constructor */
UCubicMesh::UCubicMesh():BASE_UStructMesh() {
	m_classType = "Cartesian";
};

/*! Custom Constructor. Set mesh origin, span size, dimension and nodal data structure
 * \param[in] origin_ point origin in global reference system
 * \param[in] spanX_ width span -> must be > 0;
 * \param[in] spanY_ height span -> must be > 0;
 * \param[in] spanZ_ depth span -> must be > 0;
 * \param[in] nx_   number of cells in x-direction
 * \param[in] ny_   number of cells in y-direction
 * \param[in] nz_   number of cells in z-direction
 */
UCubicMesh::UCubicMesh(darray3E origin_, double spanX_, double spanY_, double spanZ_, int nx_, int ny_, int nz_){
	m_classType = "Cartesian";
	setMesh(origin_,spanX_, spanY_, spanZ_, nx_, ny_, nz_);
};

/*! Basic Destructor */
UCubicMesh::~UCubicMesh(){
	m_classType = "";
};

/*! Copy Constructor. 
 * \param[in] other UCubicMesh object where copy from
 */
UCubicMesh::UCubicMesh(const UCubicMesh & other){
	m_classType = other.m_classType;
	*(static_cast<BASE_UStructMesh * >(this)) = *(static_cast<const BASE_UStructMesh * >(&other));
};

/*! Copy Operator. 
 * \param[in] other UCubicMesh object where copy from
 */
UCubicMesh & UCubicMesh::operator=(const UCubicMesh & other){
	m_classType = other.m_classType;
	*(static_cast<BASE_UStructMesh * >(this)) = *(static_cast<const BASE_UStructMesh * >(&other));
	return(*this);  
};

/*! Return the class type. "Cartesian" keyword identifies the class type UCubicMesh */
std::string UCubicMesh::getClassType(){return(m_classType);};


/*! Set mesh origin, span size, dimension and nodal data structure
 * \param[in] origin_ point origin in global reference system
 * \param[in] spanX_ width span -> must be > 0;
 * \param[in] spanY_ height span -> must be > 0;
 * \param[in] spanZ_ depth span -> must be > 0;
 * \param[in] nx_   number of cells in x-direction,1 is default;
 * \param[in] ny_   number of cells in y-direction,1 is default;
 * \param[in] nz_   number of cells in z-direction,1 is default;
 */
void UCubicMesh::setMesh(darray3E origin_, double spanX_, double spanY_, double spanZ_, int nx_, int ny_, int nz_){
	
	clearMesh();
	nx_ = std::max(1, nx_); 
	ny_ = std::max(1, ny_); 
	nz_ = std::max(1, nz_); 
	//set Origin and Span
	m_origin = origin_;
	m_span[0] = std::fmax(0, spanX_); 
	m_span[1] = std::fmax(0, spanY_); 
	m_span[2] = std::fmax(0, spanZ_); 
	
	// Number of mesh cells
	m_nx = nx_; m_ny = ny_; m_nz = nz_;
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
};

/*! Scale your mesh in the 3D space 
 * \param[in] sx x coordinate translation
 * \param[in] sy y coordinate translation
 * \param[in] sz z coordinate translation
 */
void UCubicMesh::scaleMesh(double sx, double sy, double sz){
	
	if(sx <=0 || sy<=0 || sz<=0){return;}
	m_span[0]= sx*m_span[0];
	m_span[1]= sy*m_span[1];
	m_span[2]= sz*m_span[2];
	
	m_dx= sx*m_dx;
	m_dy= sy*m_dy;
	m_dz= sz*m_dz;
	
	// get point distro;
	for (int i = 0; i < m_nx+1; i++) {m_xedge[i] = ((double) i) * m_dx;} 
	for (int i = 0; i < m_ny+1; i++) {m_yedge[i] = ((double) i) * m_dy;}
	for (int i = 0; i < m_nz+1; i++) {m_zedge[i] = ((double) i) * m_dz;}
	// get cell distro
	for (int i = 0; i < m_nx; i++) {m_xnode[i] = m_xedge[i] + 0.5 * m_dx;}
	for (int i = 0; i < m_ny; i++) {m_ynode[i] = m_yedge[i] + 0.5 * m_dy;}
	for (int i = 0; i < m_nz; i++) {m_znode[i] = m_zedge[i] + 0.5 * m_dz;}
};

/*! Convert the target point from local reference system (centered on the mesh origin) to the global one
 * \param[in] point target point, in local mesh r.s.
 * \param[out] result transformed point
 */
darray3E UCubicMesh::transfToGlobal( darray3E & point){
	
	darray3E result = point + m_origin;
	return(result);  
};

/*! Convert the target point from local reference system (centered on the mesh origin) to the global one
 * \param[in] point target point, in local mesh r.s.
 * \param[out] result transformed point
 */
dvector1D UCubicMesh::transfToGlobal( dvector1D & point){
	
	darray3E result = conArray<double,3>(point) + m_origin;
	return(conVect(result));  
};

/*! Convert a target list of points from local reference system (centered on the mesh origin) to the global one
 * \param[in] point target point, in local mesh r.s.
 * \param[out] result transformed list 
 */
dvecarr3E UCubicMesh::transfToGlobal( dvecarr3E & list_points){
	
	int size = list_points.size();
	dvecarr3E result(size);
	for(int i=0; i<size; ++i){
		result[i] = transfToGlobal(list_points[i]);
	}
	return(result);
};    

/*! Convert the target point from the global reference system to the local one (centered on the mesh origin).
 * \param[in] point target point, in global r.s.
 * \param[out] result transformed point 
 */
darray3E UCubicMesh::transfToLocal( darray3E & point){
	
	darray3E result = point - m_origin;
	return(result);  
};

/*! Convert the target point from the global reference system to the local one (centered on the mesh origin).
 * \param[in] point target point, in global r.s.
 * \param[out] result transformed point
 */
dvector1D UCubicMesh::transfToLocal( dvector1D & point){
	
	darray3E result = conArray<double,3>(point) - m_origin;
	return(conVect(result));  
	
};

/*! Convert a target list of points from global reference system to the local one (centered on the mesh origin)
 * \param[in] point target point, in global r.s.
 * \param[out] result transformed list
 */
dvecarr3E UCubicMesh::transfToLocal( dvecarr3E & list_points){
	
	int size = list_points.size();
	dvecarr3E result(size);
	for(int i=0; i<size; ++i){
		result[i] = transfToLocal(list_points[i]);
	}
	return(result);
};    

/*! Get Cell Volume.
 * \param[in] index cell index
 * \param[out] result cell volume
 */
double UCubicMesh::getCellVolume(int index){
	darray3E sp = getSpacing();
	return(sp[0]*sp[1]*sp[2]);
	
};

/*! Get Scaling factors of Transformation from Mesh coordinate system to Global coordinate system 
 * \param[in] index cell global index
 * \param[out] result scaling factors
 */
darray3E UCubicMesh::getLocalScaling(int index){
	darray3E result; 
	result.fill(1.0);
	return(result);
};

/*! Get Scaling factors of Transformation from Mesh coordinate system to Global coordinate system 
 * \param[in] i_, j_, k_ cell cartesian indices
 * \param[out] result scaling factors
 */
darray3E UCubicMesh::getLocalScaling(int i_, int j_, int k_){
	darray3E result;
	result.fill(1.0);
	return(result);
};


//*****************************************************************************************************************************
// UCYLINDRICALMESH IMPLEMENTATION 
/*
 *	\date			31/12/2015
 *	\authors		Gallizio Federico
 * 	\authors 		Amos Tori
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief Class for 3D uniform Cylindrical Mesh 
 *
 *	Derived class from UStructMesh for uniform cartesian grids, in cylindrical reference frame.
 */

/*! Basic Constructor */
UCylindricalMesh::UCylindricalMesh():BASE_UStructMesh() {
	m_classType = "Cylindrical";
	m_thetaOrigin = 0.0;
};

/*! Custom Constructor. Set mesh origin, span size, dimension and nodal data structure
 * \param[in] origin_ point origin in global reference system
 * \param[in] spanR_ max base radius span: must be >0; 
 * \param[in] spanZ_ max cylinder height span: must be >0; 
 *   \param[in] thetalim lower and upper limits on the angular coordinate 
 * \param[in] nr_   number of cells in r-direction
 * \param[in] nt_   number of cells in theta-direction
 * \param[in] nz_   number of cells in z-direction
 */
UCylindricalMesh::UCylindricalMesh(darray3E origin_, double spanR_, double spanZ_, dvector1D & thetalim_, int nr_, int nt_, int nz_){
	m_classType = "Cylindrical";	      
	setMesh(origin_, spanR_, spanZ_, thetalim_, nr_, nt_, nz_);	      
};

/*! Basic Destructor */
UCylindricalMesh::~UCylindricalMesh(){
	m_classType = "";
};

/*! Copy Constructor. 
 * \param[in] other UCylindricalMesh object where copy from
 */
UCylindricalMesh::UCylindricalMesh(const UCylindricalMesh & other){
	*(static_cast<BASE_UStructMesh * >(this)) = *(static_cast<const BASE_UStructMesh * >(&other));
	m_classType = other.m_classType;
	m_thetaOrigin = other.m_thetaOrigin;
};

/*! Copy Operator. 
 * \param[in] other UCylindrical object where copy from
 */
UCylindricalMesh & UCylindricalMesh::operator=(const UCylindricalMesh & other){
	*(static_cast<BASE_UStructMesh * >(this)) = *(static_cast<const BASE_UStructMesh * >(&other));
	m_classType = other.m_classType;
	m_thetaOrigin = other.m_thetaOrigin;
	return(*this);  
};

/*! Return the class type. "Cylindrical" keyword identifies the class type UCylindicalMesh */
std::string UCylindricalMesh::getClassType(){return(m_classType);};

/*! Set mesh origin, span size, dimension and nodal data structure
 * \param[in] origin_ point origin in global reference system
 * \param[in] spanR_ max base radius span: must be >0; 
 * \param[in] spanZ_ max cylinder height span: must be >0; 
 *   \param[in] thetalim lower and upper limits on the angular coordinate 
 * \param[in] nr_   number of cells in r-direction,1 is default
 * \param[in] nt_   number of cells in theta-direction(4 is default), nt<4 check default for tangential direction
 * \param[in] nz_   number of cells in z-direction,1 is default
 */
void UCylindricalMesh::setMesh(darray3E origin_, double spanR_, double spanZ_, dvector1D & thetalim_, int nr_, int nt_, int nz_){
	
	clearMesh();
	nr_ = std::max(nr_,1);
	nz_ = std::max(nz_,1);
	nt_ = std::max(nt_,4);
	
	//check thetalim
	double hlim = 8.0*std::atan(1.0);
	dvector1D thetalim = thetalim_;
	
	if((thetalim_[1] - thetalim_[0]) < 0){	
		double dum = thetalim[1];
		thetalim[0] = thetalim[1];
		thetalim[1] = dum;
	}
	
	double spanth = thetalim[1] - thetalim[0];
	if(spanth <1.0e-12){spanth=hlim;}
	thetalim[1] = thetalim[0] + std::fmin(hlim, spanth);
	
	//set Origin and Span
	m_origin = origin_;
	m_span[0] = std::fmax(0,spanR_); 
	m_span[1] = thetalim[1] - thetalim[0]; 
	m_span[2] = std::fmax(0, spanZ_);
	
	m_thetaOrigin= thetalim[0];
	
	//preliminary check
	if(m_span[1]<=0 ){std::cout<<"Not correct span detected in setMesh"<<endl; return;};
	
	// Number of mesh cells
	m_nx = nr_; m_ny = nt_; m_nz = nz_;
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
};

/*! Scale your mesh in the 3D space 
 * \param[in] sr radial coordinate translation
 * \param[in] sz z coordinate translation
 */
void UCylindricalMesh::scaleMesh(double sr, double sz){
	
	if(sr <=0 || sz<=0){return;}
	m_span[0]= sr*m_span[0];
	m_span[2]= sz*m_span[2];
	
	m_dx= sr*m_dx;
	m_dz= sz*m_dz;
	
	// get point distro;
	for (int i = 0; i < m_nx+1; i++) {m_xedge[i] = ((double) i) * m_dx;} 
	for (int i = 0; i < m_nz+1; i++) {m_zedge[i] = ((double) i) * m_dz;}
	// get cell distro
	for (int i = 0; i < m_nx; i++) {m_xnode[i] = m_xedge[i] + 0.5 * m_dx;}
	for (int i = 0; i < m_nz; i++) {m_znode[i] = m_zedge[i] + 0.5 * m_dz;}
};

/*! Convert the target point from local reference system (centered on the mesh
 *  origin and starting from given angular origin) to the global one
 * \param[in] point target point, in local mesh r.s.
 * \param[out] result transformed point
 */
darray3E UCylindricalMesh::transfToGlobal( darray3E & point){
	
	darray3E result;
	result[0] = point[0]*std::cos(point[1] + m_thetaOrigin) + m_origin[0];
	result[1] = point[0]*std::sin(point[1] + m_thetaOrigin) + m_origin[1];
	result[2] = point[2] + m_origin[2];
	
	return(result);  
};

/*! Convert the target point from local reference system (centered on the mesh
 *  origin and starting from given angular origin) to the global one
 * \param[in] point target point, in local mesh r.s.
 * \param[out] result transformed point
 */
dvector1D UCylindricalMesh::transfToGlobal( dvector1D & point){
	
	darray3E result;
	result[0] = point[0]*std::cos(point[1] + m_thetaOrigin) + m_origin[0];
	result[1] = point[0]*std::sin(point[1] + m_thetaOrigin) + m_origin[1];
	result[2] = point[2] + m_origin[2];
	
	return(conVect(result));   
};

/*! Convert a target list of points from local reference system (centered on the mesh
 *  origin and starting from given angular origin) to the global one
 * \param[in] point target point, in local mesh r.s.
 * \param[out] result transformed list 
 */
dvecarr3E UCylindricalMesh::transfToGlobal( dvecarr3E & list_points){
	
	int size = list_points.size();
	dvecarr3E result(size);
	for(int i=0; i<size; ++i){
		result[i] = transfToGlobal(list_points[i]);
	}
	return(result);
};    

/*! Convert the target point from the global reference system to the local one (centered on the mesh
 *  origin and starting from given angular origin).
 * \param[in] point target point, in global r.s.
 * \param[out] result transformed point 
 */
darray3E UCylindricalMesh::transfToLocal( darray3E & point){
	
	darray3E result; 
	double p1, p2, pdum;
	p1 = point[0]-m_origin[0];
	p2 = point[1]-m_origin[1];
	if(p1 ==0.0 && p2 ==0.0){result[0] = 0.0; result[1] = 0.0;}
	else{
		result[0] = pow(p1*p1+p2*p2,0.5);
		pdum = std::atan2(p2,p1);
		result[1] = pdum - 4.0*(getSign(pdum)-1.0)*std::atan(1.0); 
	}
	
	//get to the correct m_thetaOrigin mark
	double param = 8*std::atan(1.0);
	result[1] = result[1] - m_thetaOrigin;
	if(result[1] < 0) 		result[1] = param + result[1];
	if(result[1] > param) 	result[1] = result[1] - param;
	
	result[2] = point[2] - m_origin[2];
	
	return(result);  
};

/*! Convert the target point from the global reference system to the local one (centered on the mesh
 *  origin and starting from given angular origin).
 * \param[in] point target point, in global r.s.
 * \param[out] result transformed point
 */
dvector1D UCylindricalMesh::transfToLocal( dvector1D & point){
	
	darray3E result; 
	double p1, p2, pdum;
	p1 = point[0]-m_origin[0];
	p2 = point[1]-m_origin[1];
	if(p1 ==0.0 && p2 ==0.0){result[0] = 0.0; result[1] = 0.0;}
	else{
		result[0] = pow(p1*p1+p2*p2,0.5);
		pdum = std::atan2(p2,p1);
		result[1] = pdum - 4.0*(getSign(pdum)-1.0)*std::atan(1.0); 
	}
	
	//get to the correct m_thetaOrigin mark
	double param = 8*std::atan(1.0);
	result[1] = result[1] - m_thetaOrigin;
	if(result[1] < 0) 		result[1] = param + result[1];
	if(result[1] > param) 	result[1] = result[1] - param;
	
	result[2] = point[2] - m_origin[2];
	return(conVect(result));  
	
};

/*! Convert a target list of points from global reference system to the local one (centered on the mesh
 *  origin and starting from given angular origin)
 * \param[in] point target point, in global r.s.
 * \param[out] result transformed list
 */
dvecarr3E UCylindricalMesh::transfToLocal( dvecarr3E & list_points){
	
	int size = list_points.size();
	dvecarr3E result(size);
	for(int i=0; i<size; ++i){
		result[i] = transfToLocal(list_points[i]);
	}
	return(result);
};    

/*! Get origin of the local azimuthal coordinate of your cylindrical reference system
 * \param[out] result actual azimuthal origin in radians
 */
double UCylindricalMesh::getAzimuthalOrigin(){
	return(m_thetaOrigin);
};

/*! Set origin of the local azimuthal coordinate of your cylindrical reference system.
 * The azimuthal span actually set is preserved.
 * \param[in] angle new actual azimuthal origin in radians
 */
void UCylindricalMesh::shiftAzimuthalOrigin(double angle){
	m_thetaOrigin=angle;
};


/*! Get Cell Volume.
 * \param[in] index cell index
 * \param[out] result cell volume
 */
double UCylindricalMesh::getCellVolume(int index){
	
	darray3E P = getGridCCell(index);
	darray3E sp = getSpacing();   
	double volume = P[0]*sp[0]*sp[1]*sp[2] ;
	
	return(volume);
	
};

/*! Get Scaling factors of Transformation from Mesh coordinate system to Global coordinate system 
 * \param[in] index cell global index
 * \param[out] result scaling factors
 */
darray3E UCylindricalMesh::getLocalScaling(int index){
	darray3E result;
	result.fill(1.0);
	darray3E P = getGridCCell(index);
	result[1] = P[0];
	return(result);
};

/*! Get Scaling factors of Transformation from Mesh coordinate system to Global coordinate system 
 * \param[in] i_, j_, k_ cell cartesian indices
 * \param[out] result scaling factors
 */
darray3E UCylindricalMesh::getLocalScaling(int i_, int j_, int k_){
	darray3E result;
	result.fill(1.0);
	darray3E P = getGridCCell(i_, j_, k_);
	result[1] = P[0];
	return(result);
	
	
	
	
};

//*****************************************************************************************************************************
// USPHERICALMESH IMPLEMENTATION 
/*
 *	\date			31/12/2015
 *	\authors		Gallizio Federico
 * 	\authors 		Amos Tori
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief Class for 3D uniform Spherical Mesh 
 *
 *	Derived class from UStructMesh for uniform cartesian grids, in spherical reference frame.
 */

/*! Basic Constructor */
USphericalMesh::USphericalMesh():BASE_UStructMesh() {
	m_classType = "Spherical";
	m_phiOrigin = 0.0;
	m_thetaOrigin=0.0;
};

/*! Custom Constructor. Set mesh origin, span size, dimension and nodal data structure
 * \param[in] origin_ point origin in global reference system
 *   \param[in] spanR_  max radius span: must be >0;
 * \param[in] thetalim_ lower and upper limits on the polar coordinate
 * \param[in] philim_ lower and upper limits on the azimuthal coordinate
 * \param[in] nr_   number of cells in r-direction
 * \param[in] nt_   number of cells in theta-direction
 * \param[in] np_   number of cells in phi-direction
 */
USphericalMesh::USphericalMesh(darray3E origin_, double spanR_, dvector1D thetalim_, dvector1D philim_, int nr_, int nt_, int np_){
	m_classType = "Spherical";	      
	setMesh(origin_, spanR_,thetalim_, philim_, nr_, nt_, np_);	      
};

/*! Basic Destructor */
USphericalMesh::~USphericalMesh(){
	m_classType = "";
};

/*! Copy Constructor. 
 * \param[in] other USphericalMesh object where copy from
 */
USphericalMesh::USphericalMesh(const USphericalMesh & other){
	*(static_cast<BASE_UStructMesh * >(this)) = *(static_cast<const BASE_UStructMesh * >(&other));
	m_classType = other.m_classType;	
};

/*! Copy Operator. 
 * \param[in] other USpherical object where copy from
 */
USphericalMesh & USphericalMesh::operator=(const USphericalMesh & other){
	*(static_cast<BASE_UStructMesh * >(this)) = *(static_cast<const BASE_UStructMesh * >(&other));
	m_classType = other.m_classType;
	return(*this);  
};

/*! Return the class type. "Spherical" keyword identifies the class type USphericalMesh */
std::string USphericalMesh::getClassType(){return(m_classType);};


/*! Set mesh origin, span size, dimension and nodal data structure
 * \param[in] origin_ point origin in global reference system
 *   \param[in] spanR_  max radius span: must be >0;
 * \param[in] thetalim_ lower and upper limits on the polar coordinate
 * \param[in] philim_ lower and upper limits on the azimuthal coordinate
 * \param[in] nr_   number of cells in r-direction, 1 is default.
 * \param[in] nt_   number of cells in theta-direction, 4 is default.
 * \param[in] np_   number of cells in phi-direction, 2 is default.
 */
void USphericalMesh::setMesh(darray3E origin_, double spanR_, dvector1D thetalim_, dvector1D philim_, int nr_, int nt_, int np_){
	
	clearMesh();
	
	nr_  = std::max(1,nr_);
	nt_  = std::max(4,nt_);
	np_  = std::max(2,np_);
	//check thetalim
	double hlim = 8.0*std::atan(1.0);
	dvector1D thetalim = thetalim_;
	dvector1D philim = philim_;
	{
		//check inversion of limits 
		if((thetalim[1] - thetalim[0])<0){
			double dum = thetalim[0];
			thetalim[0] = thetalim[1];
			thetalim[1] = dum;
		}
		
		if((philim[1] - philim[0])<0){
			double dum = philim[0];
			philim[0] = philim[1];
			philim[1] = dum;
		}
		
		
		
		//check span prerequisites;
		double spanth = (thetalim[1] - thetalim[0]); 
		if(spanth < 1.0e-12){spanth = hlim;}
		thetalim[1] = thetalim[0] + std::fmin(spanth,hlim);
		
		//phi is a fixed domain check limit minimum and maximum.
		philim[0] = std::fmax(0, philim[0]);
		philim[1] = std::fmin(0.5*hlim, philim[1]);
		
		if(abs(philim[1]-philim[0]) < 1.0e-12){
			philim[0] = 0.0;
			philim[1] = 0.5*hlim;
		}
	}
	
	//set Origin and Span
	m_origin = origin_;
	m_span[0] = std::fmax(0, spanR_); 
	m_span[1] = thetalim[1] - thetalim[0]; 
	m_span[2] = philim[1] - philim[0];
	
	m_thetaOrigin = thetalim[0];
	m_phiOrigin = philim[0];
	
	// Number of mesh cells
	m_nx = nr_; m_ny = nt_; m_nz = np_;
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
};

/*! Scale your mesh in the 3D space 
 * \param[in] sr radial coordinate translation
 */
void USphericalMesh::scaleMesh(double sr){
	
	if(sr <=0){return;}
	m_span[0]= sr*m_span[0];
	
	m_dx= sr*m_dx;
	
	// get point distro;
	for (int i = 0; i < m_nx+1; i++) {m_xedge[i] = ((double) i) * m_dx;} 
	// get cell distro
	for (int i = 0; i < m_nx; i++) {m_xnode[i] = m_xedge[i] + 0.5 * m_dx;}
};

/*! Convert the target point from local reference system (centered on the mesh origin) to the global one
 * \param[in] point target point, in local mesh r.s.
 * \param[out] result transformed point
 */
darray3E USphericalMesh::transfToGlobal( darray3E & point){
	
	darray3E result;
	result[0] = point[0]*std::cos(point[1]+m_thetaOrigin)*std::sin(point[2]+m_phiOrigin) + m_origin[0];
	result[1] = point[0]*std::sin(point[1]+m_thetaOrigin)*std::sin(point[2]+m_phiOrigin) + m_origin[1];
	result[2] = point[0]*std::cos(point[2]+m_phiOrigin) + m_origin[2];
	
	return(result);  
};

/*! Convert the target point from local reference system (centered on the mesh origin) to the global one
 * \param[in] point target point, in local mesh r.s.
 * \param[out] result transformed point
 */
dvector1D USphericalMesh::transfToGlobal( dvector1D & point){
	
	darray3E result;
	result[0] = point[0]*std::cos(point[1]+m_thetaOrigin)*std::sin(point[2]+m_phiOrigin) + m_origin[0];
	result[1] = point[0]*std::sin(point[1]+m_thetaOrigin)*std::sin(point[2]+m_phiOrigin) + m_origin[1];
	result[2] = point[0]*std::cos(point[2]+m_phiOrigin) + m_origin[2];
	
	return(conVect(result));   
};

/*! Convert a target list of points from local reference system (centered on the mesh origin) to the global one
 * \param[in] point target point, in local mesh r.s.
 * \param[out] result transformed list 
 */
dvecarr3E USphericalMesh::transfToGlobal( dvecarr3E & list_points){
	
	int size = list_points.size();
	dvecarr3E result(size);
	for(int i=0; i<size; ++i){
		result[i] = transfToGlobal(list_points[i]);
	}
	return(result);
};    

/*! Convert the target point from the global reference system to the local one (centered on the mesh origin).
 * \param[in] point target point, in global r.s.
 * \param[out] result transformed point 
 */
darray3E USphericalMesh::transfToLocal( darray3E & point){
	
	darray3E result; result.fill(0.0); 
	double pdum;
	darray3E P = point - m_origin;
	result[0] = norm2(P);
	
	if(result[0] ==0.0){return(result);}
	
	if(P[0] == 0 && P[1] ==0){ result[1] =0.0;}
	else{
		pdum = std::atan2(P[1],P[0]);
		result[1] = pdum - 4.0*(getSign(pdum)-1.0)*std::atan(1.0); 
	}
	
	//get the correct m_thetaOrigin mark
	double param = 8*std::atan(1.0);
	result[1] = result[1] - m_thetaOrigin;
	if(result[1] < 0) 		result[1] = param + result[1];
	if(result[1] > param) 	result[1] = result[1] - param;
	
	
	result[2] = std::acos(P[2]/result[0]) ;
	
	//get the correct m_phiOrigin mark
	result[2] = result[2] - m_phiOrigin;
	
	
	return(result);  
};

/*! Convert the target point from the global reference system to the local one (centered on the mesh origin).
 * \param[in] point target point, in global r.s.
 * \param[out] result transformed point
 */
dvector1D USphericalMesh::transfToLocal( dvector1D & point){
	
	darray3E result; result.fill(0.0); 
	double pdum;
	
	darray3E P = conArray<double,3>(point) - m_origin;
	result[0] = norm2(P);
	
	if(result[0] ==0.0){return(conVect(result));}
	
	if(P[0] == 0 && P[1] ==0){ result[1] =0.0;}
	else{
		pdum = std::atan2(P[1],P[0]);
		result[1] = pdum - 4.0*(getSign(pdum)-1.0)*std::atan(1.0); 
	}
	
	//get the correct m_thetaOrigin mark
	double param = 8*std::atan(1.0);
	result[1] = result[1] - m_thetaOrigin;
	if(result[1] < 0) 		result[1] = param + result[1];
	if(result[1] > param) 	result[1] = result[1] - param;
	
	
	result[2] = std::acos(P[2]/result[0]) ;
	//get the correct m_phiOrigin mark
	result[2] = result[2] - m_phiOrigin;
	
	return(conVect(result));  
	
};

/*! Convert a target list of points from global reference system to the local one (centered on the mesh origin)
 * \param[in] point target point, in global r.s.
 * \param[out] result transformed list
 */
dvecarr3E USphericalMesh::transfToLocal( dvecarr3E & list_points){
	
	int size = list_points.size();
	dvecarr3E result(size);
	for(int i=0; i<size; ++i){
		result[i] = transfToLocal(list_points[i]);
	}
	return(result);
};    


/*! Get origin of the local azimuthal coordinate of your spherical reference system
 * \param[out] result actual azimuthal origin in radians
 */
double USphericalMesh::getAzimuthalOrigin(){
	return(m_thetaOrigin);
};

/*! Set origin of the local azimuthal coordinate of your spherical reference system.
 * The azimuthal span actually set is preserved.
 * \param[in] angle new azimuthal origin in radians
 */
void USphericalMesh::shiftAzimuthalOrigin(double angle){
	m_thetaOrigin=angle;
};

/*! Get origin of the local polar coordinate of your spherical reference system
 * \param[out] result actual polar origin in radians
 */
double USphericalMesh::getPolarOrigin(){
	return(m_phiOrigin);
};

/*! Set origin of the local polar coordinate of your spherical reference system.
 * The polar span could not be preserved, because coordinate limits must be contained in [0, pi].
 * Mesh will be recalculated.
 * \param[in] angle new polar origin (radians) in [0, pi].
 */
void USphericalMesh::shiftPolarOrigin(double angle){
	double param = 4.0*std::atan(1.0);
	m_phiOrigin=std::fmin(param, std::fmax(angle, 0));
	
	//get info
	darray3E or_ = getOrigin();
	darray3E sp_ = getSpan();
	ivector1D dim_ = getDimension();
	
	dvector1D plim(2,0),tlim(2,0);
	plim[0] = m_phiOrigin;
	plim[1] = plim[0] + sp_[2];
	plim[1] = std::fmin(plim[1],param);
	
	tlim[0] = getAzimuthalOrigin();
	tlim[1] = tlim[0] + sp_[1];
	
	setMesh(or_,sp_[0], tlim,plim,dim_[0],dim_[1], dim_[2]);
};


/*! Get Cell Volume.
 * \param[in] index cell index
 * \param[out] result cell volume
 */
double USphericalMesh::getCellVolume(int index){
	
	
	darray3E P = getGridCCell(index);
	darray3E sp = getSpacing();   
	double volume = P[0]*P[0] * std::sin(P[2])*sp[0]*sp[1]*sp[2];
	return(volume);
	
};

/*! Get Scaling factors of Transformation from Mesh coordinate system to Global coordinate system 
 * \param[in] index cell global index
 * \param[out] result scaling factors
 */
darray3E USphericalMesh::getLocalScaling(int index){
	
	darray3E result;
	result.fill(1.0);
	darray3E P = getGridCCell(index);
	result[1] = P[0]*std::sin(P[2]);
	result[2] = P[0];
	return(result);
};

/*! Get Scaling factors of Transformation from Mesh coordinate system to Global coordinate system 
 * \param[in] i_, j_, k_ cell cartesian indices
 * \param[out] result scaling factors
 */
darray3E USphericalMesh::getLocalScaling(int i_, int j_, int k_){
	darray3E result;
	result.fill(1.0);
	darray3E P = getGridCCell(i_, j_, k_);
	result[1] = P[0]*std::sin(P[2]);
	result[2] = P[0];
	return(result);
	
};



















