#include "libRA_BasicMeshes2D.hpp"

//*****************************************************************************************************************************
// BASE_UStructMesh2D IMPLEMENTATION 
/*
 *	\date			31/12/2015
 *	\authors		Gallizio Federico
 * 	\authors 		Lachance Steve
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief Virtual class for 2D uniform structured mesh --> based on ucartmesh paradigm of Bitpit 
 *
 *	Base class for uniform structured grids, suitable for derivation of 2D meshes on pure cartesian or polar
 *	coordinates system.
 */

/*! Basic Constructor */
BASE_UStructMesh2D::BASE_UStructMesh2D(){
      nx = 0; ny=0;
      dx = 0; dy=0;
};

/*! Basic destructor */
BASE_UStructMesh2D::~BASE_UStructMesh2D(){
  dx = dy = 0.0;
  nx = ny = 0;
  xnode.clear(); ynode.clear();
  xedge.clear(); yedge.clear();
}

/*! Copy Constructor.
 * \param[in] other BASE_UStructMesh2D object where copy from
 */
BASE_UStructMesh2D::BASE_UStructMesh2D(const BASE_UStructMesh2D & other){
  // Number of cells
  nx = other.nx;
  ny = other.ny;
  
  // Cell Spacing
  dx = other.dx;
  dy = other.dy;
  
  // Mesh Origin & span	
  origin = other.origin;
  span = other.span;
  // Resize mesh data structure ----------------------------------------------- //
  resizeMesh();
  
  // Copy cell edges and cell centers ----------------------------------------- //
  xnode = other.xnode;
  ynode = other.ynode;
  xedge = other.xedge;
  yedge = other.yedge;
};

/*! Copy Operator.
 * \param[in] other BASE_UStructMesh2D object where copy from
 */
BASE_UStructMesh2D & BASE_UStructMesh2D::operator=(const BASE_UStructMesh2D & other){
  
  // Number of cells
  nx = other.nx;
  ny = other.ny;
  
  // Cell Spacing
  dx = other.dx;
  dy = other.dy;
  
  // Mesh Origin & span	
  origin = other.origin;
  span = other.span;
  // Resize mesh data structure ----------------------------------------------- //
  resizeMesh();
  
  // Copy cell edges and cell centers ----------------------------------------- //
  xnode = other.xnode;
  ynode = other.ynode;
  xedge = other.xedge;
  yedge = other.yedge;
  
  return(*this); 
};

/*!Clear the Mesh structure. Resize your nodal vector list to zero, but not destroy them.*/
void BASE_UStructMesh2D::clearMesh(){
  nx = 0; ny=0;;
  dx=0; dy=0; ;
  
  origin.fill(0.0);
  span.fill(0.0);
  
  resizeMesh();
};  

/*! Resize the nodal vector lists to current mesh dimensions, but don't destroy them. Please use reshape and destroy methods to handle them */
void BASE_UStructMesh2D::resizeMesh(){
  // Cell centers
  xnode.resize(nx, 0.0);
  ynode.resize(ny, 0.0);
  
  // Points
  xedge.resize(nx+1, 0.0);
  yedge.resize(ny+1, 0.0);
};

/*! Destroy the all nodal structures of the mesh. */
void BASE_UStructMesh2D::destroyNodalStructure(){
  freeContainer(xnode);
  freeContainer(ynode);
  freeContainer(xedge);
  freeContainer(yedge);  
};

/*! Destroy the all nodal structures of the mesh, and reinitialize them to current mesh dimensions. */
void BASE_UStructMesh2D::reshapeNodalStructure(){
  destroyNodalStructure();
  resizeMesh();
};


/*! Translate your mesh in the 3D space 
 * \param[in] tx x coordinate translation
 * \param[in] ty y coordinate translation
 */
void BASE_UStructMesh2D::translateMesh(double tx, double ty){
  origin[0]+= tx;
  origin[1]+= ty;
};


/*! Return current mesh origin */
darray2E BASE_UStructMesh2D::getOrigin(){return(origin);};

/*! Return current mesh span */
darray2E BASE_UStructMesh2D::getSpan(){return(span);};

/*! Return current mesh spacing */
darray2E BASE_UStructMesh2D::getSpacing(){
  darray2E res;
  res[0] = dx; res[1] =dy;
  return(res); 
};

/*! Return current dimension of the mesh (number of cells in each direction) */
ivector1D BASE_UStructMesh2D::getDimension(){
  ivector1D res(2,0);
  res[0] = nx; res[1] =ny; 
  return(res); 
};

/*! Get n-th center cell coordinates in local mesh reference frame, 
 * given its global cell index on the mesh. 
 * \param[in] index cell index in the global nodal list.
 */
darray2E BASE_UStructMesh2D::getGridCCell(int index){
  
  int i,j;
  accessCellData(index, i,j);
  darray2E res = getGridCCell(i,j);
  
  return(res);
};

/*! Get n-th center cell coordinates in local mesh reference frame, 
 * given its cartesian cell indices on the mesh. 
 * \param[in] i_ x cartesian index.
 * \param[in] j_ y cartesian index.
 */
darray2E BASE_UStructMesh2D::getGridCCell(int i_, int j_){
  
  darray2E res{0,0};
  res[0] = xnode[i_];
  res[1] = ynode[j_];
  
  return(res);
};

/*! Get n-th nodal vertex coordinates in local mesh reference frame, 
 * given its global point index on the mesh. 
 * \param[in] index point index in the global nodal list.
 */
darray2E BASE_UStructMesh2D::getGridPoint(int index){
  
  int i,j;
  accessPointData(index, i,j);
  darray2E res = getGridPoint(i,j);
  
  return(res);
};

/*! Get n-th nodal vertex coordinates in local reference frame, 
 * given its cartesian point indices on the mesh. 
 * \param[in] i_ x cartesian index.
 * \param[in] j_ y cartesian index.
 */
darray2E BASE_UStructMesh2D::getGridPoint(int i_, int j_){
  
  darray2E res{0,0};
  res[0] = xedge[i_];
  res[1] = yedge[j_];
  
  return(res);
};    

/*! Get n-th center cell coordinates in global absolute reference frame, 
 * given its global cell index on the mesh. 
 * \param[in] index cell index in the global nodal list.
 */
darray2E BASE_UStructMesh2D::getGlobalCCell(int index){
  
  darray2E res = getGridCCell(index);
  return(transfToGlobal(res));
};

/*! Get n-th center cell coordinates in global absolute reference frame, 
 * given its cartesian cell indices on the mesh. 
 * \param[in] i_ x cartesian index.
 * \param[in] j_ y cartesian index.
 */
darray2E BASE_UStructMesh2D::getGlobalCCell(int i_, int j_){
  
  darray2E res = getGridCCell(i_,j_);
  return(transfToGlobal(res));
};

/*! Get n-th nodal vertex coordinates in global absolute reference frame, 
 * given its global point index on the mesh. 
 * \param[in] index point index in the global nodal list.
 */
darray2E BASE_UStructMesh2D::getGlobalPoint(int index){
  darray2E res = getGridPoint(index);
  return(transfToGlobal(res));
};

/*! Get n-th nodal vertex coordinates in global absolute reference frame, 
 * given its cartesian point indices on the mesh. 
 * \param[in] i_ x cartesian index.
 * \param[in] j_ y cartesian index.
 */
darray2E BASE_UStructMesh2D::getGlobalPoint(int i_, int j_){
  
  darray2E res = getGridPoint(i_,j_);
  return(transfToGlobal(res));
};    

/*! Get neighbor vertices (by their global indices and in VTK-hexahedra ordered) of a given cell
 * \param[in] i_ x coordinate cell index
 * \param[in] j_ y coordinate cell index
 */ 
ivector1D BASE_UStructMesh2D::getCellNeighs(int i_, int j_){
  
  ivector1D result(4);
  
  result[0] = accessPointData(i_, j_);
  result[1] = accessPointData(i_+1, j_);
  result[2] = accessPointData(i_+1, j_+1);
  result[3] = accessPointData(i_, j_+1);
  
  return(result);  
}

/*! Get neighbor vertices (by their global indices and in VTK-hexahedra ordered) of a given cell
 * \param[in] index cell global index
 */ 
ivector1D BASE_UStructMesh2D::getCellNeighs(int index){
  
  ivector1D pp(2,0);
  accessCellData(index, pp[0],pp[1]);
  return(getCellNeighs(pp[0],pp[1]));
}


/*! Return cartesian indices of the cell containing the target point in global reference frame
 * \param[in] point 2D coordinate of target point
 * \param[out] i x cell index
 * \param[out] j y cell index 
 */ 
void BASE_UStructMesh2D::returnCellID(darray2E & point, int &i, int &j){
  
  darray2E P = transfToLocal(point);
  i = min(nx-1, max(0, (int) floor((P[0])/dx)));
  j = min(ny-1, max(0, (int) floor((P[1])/dy)));
};

/*! Return cartesian indices of the cell containing the target point in global reference framne
 * \param[in] point 2D coordinate of target point
 * \param[out] i x cell index
 * \param[out] j y cell index
 */ 
void BASE_UStructMesh2D::returnCellID(dvector1D & point, int &i, int &j){
  
  dvector1D P = transfToLocal(point);
  i = min(nx-1, max(0, (int) floor((P[0])/dx)));
  j = min(ny-1, max(0, (int) floor((P[1])/dy)));
  
};


/*! Return global index of the cell given its cartesian indices. Follows the ordering sequences y-x
 * \param[in] i x cartesian index
 *\param[in] j y cartesian index
 *\param[out] result global index 
 */
int BASE_UStructMesh2D::accessCellData(int i, int j){
  
  int index = ny * i + j ;
  return(index);
};

/*! Return cartesian indices of the cell given its global index. Follows the ordering sequences y-x
 * \param[in] N_ global index 
 *\param[out] i x cartesian index
 *\param[out] j y cartesian index
 */
void BASE_UStructMesh2D::accessCellData(int N_, int &i, int &j){
  j = N_ % ny;
  i = N_ / ny; 
};

/*! Return global index of the point given its cartesian indices. Follows the ordering sequences y-x
 * \param[in] i x cartesian index
 *\param[in] j y cartesian index
 *\param[out] result global index 
 */
int  BASE_UStructMesh2D::accessPointData(int i, int j){
  
  int index = (ny+1) * i +  j ;
  return(index);
};

/*! Return cartesian indices of the point given its global index. Follows the ordering sequences y-x
 * \param[in] N_ global index 
 *\param[out] i x cartesian index
 *\param[out] j y cartesian index
 */
void BASE_UStructMesh2D::accessPointData(int N_, int &i, int &j){
  
  j = N_ % (ny+1);
  i = N_ / (ny+1); 
};


/*! Interpolate value of a given data field on a target point inside the mesh
 * \param[in] point target point
 * \param[in] celldata data field defined on centercells
 * \param[out] interpResult interpolated value
 */
double BASE_UStructMesh2D::interpolateCellData(darray2E & point, dvector1D & celldata){
  
  int i0, j0, ip, jp;
  double wx0,wx1,wy0,wy1;
  
  returnCellID(point, i0, j0);
  darray2E P = transfToLocal(point);
  if (P[0] > xnode[i0]) 	{ ip = min(i0+1, nx-1);   }
  else                  	{ ip = max(0, i0-1);      }
  if (P[1] > ynode[j0]) 	{ jp = min(j0+1, ny-1);   }
  else                  	{ jp = max(0, j0-1);      }
  
  // Interpolation weights
  wx1 = max(0.0, min(1.0, abs((P[0] - xnode[i0])/dx)));     wx0 = 1.0 - wx1;
  wy1 = max(0.0, min(1.0, abs((P[1] - ynode[j0])/dy)));     wy0 = 1.0 - wy1;
  
  // Interpolation
  double res  = wx0 * wy0 * celldata[accessCellData(i0,j0)]
	       + wx0 * wy1 * celldata[accessCellData(i0,jp)]
	       + wx1 * wy0 * celldata[accessCellData(ip,j0)]
	       + wx1 * wy1 * celldata[accessCellData(ip,jp)];
	       
 return(res);	       
};

/*! Interpolate value of a given data field on a target point inside the mesh
 * \param[in] point target point
 * \param[in] celldata data field defined on centercells
 * \param[out] interpResult interpolated value
 */
int BASE_UStructMesh2D::interpolateCellData(darray2E & point, ivector1D & celldata){
  
   int i0, j0, ip, jp;
  double wx0,wx1,wy0,wy1;
  
  returnCellID(point, i0, j0);
  darray2E P = transfToLocal(point);
  if (P[0] > xnode[i0]) 	{ ip = min(i0+1, nx-1);   }
  else                  	{ ip = max(0, i0-1);      }
  if (P[1] > ynode[j0]) 	{ jp = min(j0+1, ny-1);   }
  else                  	{ jp = max(0, j0-1);      }
  
  // Interpolation weights
  wx1 = max(0.0, min(1.0, abs((P[0] - xnode[i0])/dx)));     wx0 = 1.0 - wx1;
  wy1 = max(0.0, min(1.0, abs((P[1] - ynode[j0])/dy)));     wy0 = 1.0 - wy1;
  
  // Interpolation
  double res1 = 	wx0 * wy0 * (double)celldata[accessCellData(i0,j0)]
		+ wx0 * wy1 *(double) celldata[accessCellData(i0,jp)]
		+ wx1 * wy0 * (double)celldata[accessCellData(ip,j0)]
		+ wx1 * wy1 *(double) celldata[accessCellData(ip,jp)];
		
	int res = std::floor(res1+0.5);
	return(res);	
};

/*! Interpolate value of a given data field on a target point inside the mesh
 * \param[in] point target point
 * \param[in] celldata data field defined on centercells
 * \param[out] interpResult interpolated value
 */
darray3E BASE_UStructMesh2D::interpolateCellData(darray2E & point, dvecarr3E & celldata){
  
    int i0, j0, ip, jp;
  double wx0,wx1,wy0,wy1;
  
  returnCellID(point, i0, j0);
  darray2E P = transfToLocal(point);
  if (P[0] > xnode[i0]) 	{ ip = min(i0+1, nx-1);   }
  else                  	{ ip = max(0, i0-1);      }
  if (P[1] > ynode[j0]) 	{ jp = min(j0+1, ny-1);   }
  else                  	{ jp = max(0, j0-1);      }
  
  // Interpolation weights
  wx1 = max(0.0, min(1.0, abs((P[0] - xnode[i0])/dx)));     wx0 = 1.0 - wx1;
  wy1 = max(0.0, min(1.0, abs((P[1] - ynode[j0])/dy)));     wy0 = 1.0 - wy1;
  
  // Interpolation
  darray3E res = wx0 * wy0 * celldata[accessCellData(i0,j0)]
	       + wx0 * wy1 * celldata[accessCellData(i0,jp)]
	       + wx1 * wy0 * celldata[accessCellData(ip,j0)]
	       + wx1 * wy1 * celldata[accessCellData(ip,jp)];
	       
	   return(res);    
};

/*! Interpolate value of a given data field on a target point inside the mesh
 * \param[in] point target point
 * \param[in] pointdata data field defined on grid points
 * \param[out] interpResult interpolated value
 */
double BASE_UStructMesh2D::interpolatePointData(darray2E & point, dvector1D & pointdata){
  
  int i0, j0, ip, jp;
  double wx0,wx1,wy0,wy1;
  
  darray2E P = transfToLocal(point);
  i0 = max(0, min(nx, (int) floor((P[0])/dx)));
  j0 = max(0, min(ny, (int) floor((P[1])/dy)));
  
  if (P[0] >= xedge[i0]) 	{ ip = min(i0+1, nx);   }
  else                  	{ ip = max(0, i0-1);      }
  if (P[1] >= yedge[j0]) 	{ jp = min(j0+1, ny);   }
  else                  	{ jp = max(0, j0-1);      }
  
  // Interpolation weights
  wx1 = max(0.0, min(1.0, abs((P[0] - xedge[i0])/dx)));     wx0 = 1.0 - wx1;
  wy1 = max(0.0, min(1.0, abs((P[1] - yedge[j0])/dy)));     wy0 = 1.0 - wy1;
  
  // Interpolation
  double res = wx0 * wy0 * pointdata[accessPointData(i0,j0)]
		+wx0 * wy1 * pointdata[accessPointData(i0,jp)]
		+wx1 * wy0 * pointdata[accessPointData(ip,j0)]
		+wx1 * wy1 * pointdata[accessPointData(ip,jp)];
		
	return(res);	
};

/*! Interpolate value of a given data field on a target point inside the mesh
 * \param[in] point target point
 * \param[in] celldata data field defined on grid points
 * \param[out] interpResult interpolated value
 */
int BASE_UStructMesh2D::interpolatePointData(darray2E & point, ivector1D & pointdata){
  
  int i0, j0, ip, jp;
  double wx0,wx1,wy0,wy1;
  
  darray2E P = transfToLocal(point);
  i0 = max(0, min(nx, (int) floor((P[0])/dx)));
  j0 = max(0, min(ny, (int) floor((P[1])/dy)));
 
  if (P[0] >= xedge[i0]) 	{ ip = min(i0+1, nx);   }
  else                  	{ ip = max(0, i0-1);      }
  if (P[1] >= yedge[j0]) 	{ jp = min(j0+1, ny);   }
  else                  	{ jp = max(0, j0-1);      }
  
  // Interpolation weights
  wx1 = max(0.0, min(1.0, abs((P[0] - xedge[i0])/dx)));     wx0 = 1.0 - wx1;
  wy1 = max(0.0, min(1.0, abs((P[1] - yedge[j0])/dy)));     wy0 = 1.0 - wy1;
  
  // Interpolation
  double res1 =wx0 * wy0 * (double)pointdata[accessPointData(i0,j0)]
		+wx0 * wy1 * (double)pointdata[accessPointData(i0,jp)]
		+wx1 * wy0 *(double) pointdata[accessPointData(ip,j0)]
		+wx1 * wy1 *(double) pointdata[accessPointData(ip,jp)];  
	int res = std::floor(res1+0.5);
	return(res);
};
/*! Interpolate value of a given data field on a target point inside the mesh
 * \param[in] point target point
 * \param[in] celldata data field defined on grid points
 * \param[out] interpResult interpolated value
 */
darray3E BASE_UStructMesh2D::interpolatePointData(darray2E & point, dvecarr3E & pointdata){
  
   int i0, j0, ip, jp;
  double wx0,wx1,wy0,wy1;
  
  darray2E P = transfToLocal(point);
  i0 = max(0, min(nx, (int) floor((P[0])/dx)));
  j0 = max(0, min(ny, (int) floor((P[1])/dy)));
 
  
  if (P[0] >= xedge[i0]) 	{ ip = min(i0+1, nx);   }
  else                  	{ ip = max(0, i0-1);      }
  if (P[1] >= yedge[j0]) 	{ jp = min(j0+1, ny);   }
  else                  	{ jp = max(0, j0-1);      }
  
  // Interpolation weights
  wx1 = max(0.0, min(1.0, abs((P[0] - xedge[i0])/dx)));     wx0 = 1.0 - wx1;
  wy1 = max(0.0, min(1.0, abs((P[1] - yedge[j0])/dy)));     wy0 = 1.0 - wy1;
  
  // Interpolation
 darray3E res = wx0 * wy0 * pointdata[accessPointData(i0,j0)]
		+wx0 * wy1 * pointdata[accessPointData(i0,jp)]
		+wx1 * wy0 * pointdata[accessPointData(ip,j0)]
		+wx1 * wy1 * pointdata[accessPointData(ip,jp)];

		return(res);
};


/*! Write your grid as a 3D point cloud in a *.vtu file. 
 * \param[in] folder output directory path
 * \param[in] outfile output namefile w/out tag
 * \param[in] counterFile integer to mark output with a counter number
 * \param[in] codexFlag boolean to distinguish between "ascii" false, and "appended" true
 * \param[in] extPoints OPTIONAL. If defined, use the current vertexList and write them as point cloud, provided that coherent dimensions with the mesh are set. 
 *			If Not, write the current mesh vertices. 
 * */
void BASE_UStructMesh2D::plotCloud( std::string & folder , std::string outfile, int counterFile, bool codexFlag, dvecarr3E * extPoints){
  
  std::string codex = "ascii";
  if(codexFlag){codex="appended";}
  
  ivector1D dim = getDimension();
  int sizeTot = (dim[0]+1)*(dim[1]+1);
  
  VTK_BASICCLOUD handle_vtk_output(folder, outfile, codex, sizeTot);
  
  dvecarr3E activeP;
  if(extPoints != NULL && extPoints->size() == sizeTot){activeP = *extPoints; }
  else{
    activeP.resize(sizeTot);
    for(int i=0; i<sizeTot; i++){
      activeP[i] = vecEnlarge(getGlobalPoint(i));
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
void BASE_UStructMesh2D::plotCloud( std::string & folder , std::string outfile, int counterFile, bool codexFlag, ivector1D & vertexList, dvecarr3E * extPoints){
  
  std::string codex = "ascii";
  if(codexFlag){codex="appended";}
  
  ivector1D dim = getDimension();
  int sizeTot = (dim[0]+1)*(dim[1]+1);
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
      activeP[i] = vecEnlarge(getGlobalPoint(vertexList[i]));
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
void BASE_UStructMesh2D::plotGrid(std::string & folder, std::string outfile , int counterFile, bool codexFlag, dvecarr3E * extPoints){
  
  std::string codex = "ascii";
  if(codexFlag){codex="appended";}
  
  ivector1D dim = getDimension();
  int sizePt = (dim[0]+1)*(dim[1]+1);
  int sizeCl = dim[0]*dim[1];

  dvecarr3E activeP;
  ivector2D activeConn(sizeCl, ivector1D(4,0));
  
if(extPoints != NULL && extPoints->size() == sizePt){
	  activeP = *extPoints; 
}else{
	activeP.resize(sizePt);  
	for(int i=0; i<sizePt; i++){
		activeP[i] = vecEnlarge(getGlobalPoint(i));
	}
}

for(int i=0; i<sizeCl; ++i){
	activeConn[i] = getCellNeighs(i); 
}

 
  VTK_BASICMESH handle_vtk_output(folder, outfile, codex, sizePt, sizeCl, 4*sizeCl);
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
void BASE_UStructMesh2D::plotGrid(std::string & folder, std::string outfile, int counterFile, bool codexFlag, ivector1D & cellList, dvecarr3E * extPoints){
  
  std::string codex = "ascii";
  if(codexFlag){codex="appended";}
  
  ivector1D dim = getDimension();
  int sizePt = (dim[0]+1)*(dim[1]+1);
  int sizeCl = cellList.size();
  
  if(sizeCl > dim[0]*dim[1]){return;}
  
  std::map<int,int> mapPoints;
  ivector2D activeConn(sizeCl, ivector1D(4,0));
  
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
  }else{
		for(itF = mapPoints.begin(); itF != mapPoints.end(); ++itF){
			int pos = itF->second;
			activeP[counter] = vecEnlarge(getGlobalPoint(pos));
			listP[counter] = pos;
			++counter;
		}
  }
  
  for(int i=0; i<sizeCl; ++i){
	  for(int j=0; j<activeConn[i].size(); ++j){
		  int posV = posVectorFind(listP, activeConn[i][j]);
		  //update local conn
		  activeConn[i][j] = posV;
	}
  }

  VTK_BASICMESH handle_vtk_output(folder, outfile, codex, counter, sizeCl, 4*sizeCl);
  if(counterFile>=0){handle_vtk_output.SetCounter(counterFile);}
  handle_vtk_output.SetGeomTypes("Float64","Int32", "Int32", "Int32");
  handle_vtk_output.linkData(activeP,activeConn);
  handle_vtk_output.Write();
};

/*! Copy your 2D point in a 3D structure with z coordinate equal to zero
 * \param[in] point 2D point
 * \param[out] result 3D point
 */
darray3E BASE_UStructMesh2D::vecEnlarge(darray2E point){
       darray3E result;
       result[0] = point[0];
       result[1] = point[1];
       result[2] = 0.0;
       return(result);
}

/*! Copy your 3D point in a 2D structure losing info on the z coordinate. 
 * \param[in] point 3D point
 * \param[out] result 2D point
 */
darray2E BASE_UStructMesh2D::vecShrink(darray3E point){
	darray2E result;
	result[0] = point[0];
	result[1] = point[1];
	return(result);
}

//*****************************************************************************************************************************
// UQUADMESH IMPLEMENTATION 
/*
 *	\date			31/12/2015
 *	\authors		Gallizio Federico
 * 	\authors 		Amos Tori
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief Class for 2D uniform Cartesian Mesh 
 *
 *	Derived class from UStructMesh2D for uniform cartesian grids.
 */

/*! Basic Constructor */
UQuadMesh::UQuadMesh() {
  classType = "Cartesian2D";
};

/*! Custom Constructor. Set mesh origin, span size, dimension and nodal data structure
 * \param[in] origin_ point origin in global reference system
 * \param[in] spanX_ width span -> must be > 0;
 * \param[in] spanY_ height span -> must be > 0;
 * \param[in] nx_   number of cells in x-direction
 * \param[in] ny_   number of cells in y-direction
 */
UQuadMesh::UQuadMesh(darray2E origin_, double spanX_, double spanY_, int nx_, int ny_){
  
  classType = "Cartesian2D";	      
  setMesh(origin_,spanX_, spanY_, nx_, ny_);
};

/*! Basic Destructor */
UQuadMesh::~UQuadMesh(){
  classType = "";
};

/*! Copy Constructor. 
 * \param[in] other UQuadMesh object where copy from
 */
UQuadMesh::UQuadMesh(const UQuadMesh & other){
	*(static_cast<BASE_UStructMesh2D * >(this)) = *(static_cast<const BASE_UStructMesh2D * >(&other));
	classType = other.classType;
};

/*! Copy Operator. 
 * \param[in] other UQuadMesh object where copy from
 */
UQuadMesh & UQuadMesh::operator=(const UQuadMesh & other){
	*(static_cast<BASE_UStructMesh2D * >(this)) = *(static_cast<const BASE_UStructMesh2D * >(&other));
	classType = other.classType;
	return(*this);  
};

/*! Return the class type. "Cartesian2D" keyword identifies the class type UQuadMesh */
std::string UQuadMesh::getClassType(){return(classType);};


/*! Set mesh origin, span size, dimension and nodal data structure
 * \param[in] origin_ point origin in global reference system
 * \param[in] spanX_ width span -> must be > 0;
 * \param[in] spanY_ height span -> must be > 0;
 * \param[in] nx_   number of cells in x-direction,1 is default
 * \param[in] ny_   number of cells in y-direction,1 is default
 */
void UQuadMesh::setMesh(darray2E origin_, double spanX_, double spanY_, int nx_, int ny_){
  
	nx_ = std::max(1, nx_);
	ny_ = std::max(1, ny_);
	
  //set Origin and Span
  origin = origin_;
  span[0] = std::fmax(0, spanX_); 
  span[1] = std::fmax(0, spanY_); 
  
  // Number of mesh cells
  nx = nx_; ny = ny_; 
  // Resize mesh data structure
  resizeMesh();
  // Create it 
  // get mesh spacing
  dx = span[0]/((double) nx);
  dy = span[1]/((double) ny);	
  // get point distro;
  for (int i = 0; i < nx+1; i++) {xedge[i] = ((double) i) * dx;} 
  for (int i = 0; i < ny+1; i++) {yedge[i] = ((double) i) * dy;}
  // get cell distro
  for (int i = 0; i < nx; i++) {xnode[i] = xedge[i] + 0.5 * dx;}
  for (int i = 0; i < ny; i++) {ynode[i] = yedge[i] + 0.5 * dy;}
 
};

/*! Scale your mesh in the 2D space 
 * \param[in] sx x coordinate translation
 * \param[in] sy y coordinate translation
 */
void UQuadMesh::scaleMesh(double sx, double sy){
  
  if(sx <=0 || sy<=0 ){return;}
  span[0]= sx*span[0];
  span[1]= sy*span[1];
  
  dx= sx*dx;
  dy= sy*dy;
  
  // get point distro;
  for (int i = 0; i < nx+1; i++) {xedge[i] = ((double) i) * dx;} 
  for (int i = 0; i < ny+1; i++) {yedge[i] = ((double) i) * dy;}
  // get cell distro
  for (int i = 0; i < nx; i++) {xnode[i] = xedge[i] + 0.5 * dx;}
  for (int i = 0; i < ny; i++) {ynode[i] = yedge[i] + 0.5 * dy;}
};

/*! Convert the target point from local reference system (centered on the mesh origin) to the global one
 * \param[in] point target point, in local mesh r.s.
 * \param[out] result transformed point
 */
darray2E UQuadMesh::transfToGlobal( darray2E & point){
  
  darray2E result = point + origin;
  return(result);  
};

/*! Convert the target point from local reference system (centered on the mesh origin) to the global one
 * \param[in] point target point, in local mesh r.s.
 * \param[out] result transformed point
 */
dvector1D UQuadMesh::transfToGlobal( dvector1D & point){
  
  darray2E result = conArray<double,2>(point) + origin;
  return(conVect(result));  
};

/*! Convert a target list of points from local reference system (centered on the mesh origin) to the global one
 * \param[in] point target point, in local mesh r.s.
 * \param[out] result transformed list 
 */
dvecarr2E UQuadMesh::transfToGlobal( dvecarr2E & list_points){
  
  int size = list_points.size();
  dvecarr2E result(size);
  for(int i=0; i<size; ++i){
    result[i] = transfToGlobal(list_points[i]);
  }
  return(result);
};    

/*! Convert the target point from the global reference system to the local one (centered on the mesh origin).
 * \param[in] point target point, in global r.s.
 * \param[out] result transformed point 
 */
darray2E UQuadMesh::transfToLocal( darray2E & point){
  
  darray2E result = point - origin;
  return(result);  
};

/*! Convert the target point from the global reference system to the local one (centered on the mesh origin).
 * \param[in] point target point, in global r.s.
 * \param[out] result transformed point
 */
dvector1D UQuadMesh::transfToLocal( dvector1D & point){
  
  darray2E result = conArray<double,2>(point) - origin;
  return(conVect(result));  
  
};

/*! Convert a target list of points from global reference system to the local one (centered on the mesh origin)
 * \param[in] point target point, in global r.s.
 * \param[out] result transformed list
 */
dvecarr2E UQuadMesh::transfToLocal( dvecarr2E & list_points){
  
  int size = list_points.size();
  dvecarr2E result(size);
  for(int i=0; i<size; ++i){
    result[i] = transfToLocal(list_points[i]);
  }
  return(result);
};    

/*! Get Cell Volume.
 * \param[in] index cell index
 * \param[out] result cell volume
 */
double UQuadMesh::getCellVolume(int index){
  darray2E sp = getSpacing();
  return(sp[0]*sp[1]);
  
};

/*! Get Scaling factors of Transformation from Mesh coordinate system to Global coordinate system 
 * \param[in] index cell global index
 * \param[out] result scaling factors
 */
darray2E UQuadMesh::getLocalScaling(int index){
  darray2E result; 
  result.fill(1.0);
  return(result);
};

/*! Get Scaling factors of Transformation from Mesh coordinate system to Global coordinate system 
 * \param[in] i_, j_, k_ cell cartesian indices
 * \param[out] result scaling factors
 */
darray2E UQuadMesh::getLocalScaling(int i_, int j_){
  darray2E result;
  result.fill(1.0);
  return(result);
};

