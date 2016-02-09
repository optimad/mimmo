#include "libRA_CorePatches.hpp"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//BASE_PATCH IMPLEMENTATION 
/*
 *	\date			31/12/2015
 *	\authors		Gallizio Federico
 * 	\authors 		Sister Rosetta
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief Virtual class for Basic Patches Representation 
 *
 *	Base class for Volumetric Core Element, suitable for interaction with Triangulation data Structure.
 */
 
/*! Basic Constructor */
BASE_PATCH::BASE_PATCH(){};
/*! Basic Destructor */
BASE_PATCH::~BASE_PATCH(){};

/*! Copy Constructor 
 * \param[in] other BASE_PATCH object where copy from
 */
BASE_PATCH::BASE_PATCH(const BASE_PATCH & other){
  patchType = other.patchType;
};

/*! Copy Operator 
 * \param[in] other BASE_PATCH object where copy from
 */
BASE_PATCH & BASE_PATCH::operator=(const BASE_PATCH & other){
   patchType = other.patchType;
};

/*! Get type of patch actually instantiated */
std::string BASE_PATCH::getPatchType(){return(patchType);};


/*! Given a triangulated tesselation, return indices of those triangles inside the patch 
 * \param[in] tri target triangulated tesselation
 * \param[out] result list-by-indices of triangles included in the volumetric patch
 */
ivector1D BASE_PATCH::includeTriangulation(Class_SurfTri * tri ){
  
  ivector1D result(tri->nSimplex); 
  int counter=0;
  
  for(int i=0; i<tri->nSimplex; ++i){
   
    if(isTriangleIncluded(tri, i)){
      result[counter] = i;
      ++counter;
    }
  }
  result.resize(counter);
  return(result);
};

/*! Given a triangulated tesselation, return indices of those triangles outside the patch 
 * \param[in] tri target triangulated tesselation
 * \param[out] result list-by-indices of triangles outside the volumetric patch
 */
ivector1D BASE_PATCH::excludeTriangulation(Class_SurfTri * tri){
  
  ivector1D result(tri->nSimplex); 
  int counter=0;
  
  for(int i=0; i<tri->nSimplex; ++i){
   
    if(!isTriangleIncluded(tri, i)){
      result[counter] = i;
      ++counter;
    }
  }
  result.resize(counter);
  return(result);
};

/*! Given a list of vertices of a point cloud, return indices of those vertices included into the patch 
 * \param[in] list list of cloud points
 * \param[out] result list-by-indices of vertices included in the volumetric patch
 */
ivector1D BASE_PATCH::includeCloudPoints(dvecarr3E & list){

  int size = list.size();
  ivector1D result(size); 
  int counter=0;
  
  for(int i=0; i<size; ++i){
   
    if(isPointIncluded(list[i])){
      result[counter] = i;
      ++counter;
    }
  }
  result.resize(counter);
  return(result);
};

/*! Given a list of vertices of a point cloud, return indices of those vertices outside the patch 
 * \param[in] list list of cloud points
 * \param[out] result list-by-indices of vertices outside the volumetric patch
 */
ivector1D BASE_PATCH::excludeCloudPoints(dvecarr3E & list){
  
  int size = list.size();
  ivector1D result(size); 
  int counter=0;
  
  for(int i=0; i<size; ++i){
   
    if(!isPointIncluded(list[i])){
      result[counter] = i;
      ++counter;
    }
  }
  result.resize(counter);
  return(result);
  
};

/*! Given a pidded part of a triangulated tesselation, return those pidded triangles inside the patch 
 * \param[in] sh target triangulated tesselation
 * \param[in] pidT target PID number
 * \param[out] result list-by-indices of pidded triangles inside the volumetric patch
 */
ivector1D BASE_PATCH::includePIDTriangulation(SHAPE * sh, int pidT ){
  
  ivector1D temp = sh->extractPidded(pidT);
  int size = temp.size();
  ivector1D result(size); 
  int counter=0;
  
  for(int i=0; i<size; ++i){
    if(isTriangleIncluded(sh, temp[i]) ){
      result[counter] = temp[i];
      ++counter;
    }
  }
  result.resize(counter);
  return(result);
};

/*! Given a pidded part of a triangulated tesselation, return those pidded triangles outside the patch 
 * \param[in] sh target triangulated tesselation
 * \param[in] pidT target PID number
 * \param[out] result list-by-indices of pidded triangles outside the volumetric patch
 */
ivector1D BASE_PATCH::excludePIDTriangulation(SHAPE * sh, int pidT){

  ivector1D temp = sh->extractPidded(pidT);
  int size = temp.size();
  ivector1D result(size); 
  int counter=0;
  
  for(int i=0; i<size; ++i){
   
    if(!isTriangleIncluded(sh, temp[i]) ){
      result[counter] = temp[i];
      ++counter;
    }
  }
  result.resize(counter);
  return(result);
};

/*! Given one or more pidded part of a triangulated tesselation, return those pidded triangles inside the patch 
 * \param[in] sh target triangulated tesselation
 * \param[in] pidT target PID number
 * \param[out] result list-by-indices of pidded triangles inside the volumetric patch
 */
ivector1D BASE_PATCH::includePIDTriangulation(SHAPE * sh, ivector1D & pidList ){
  
  int size = pidList.size();
  ivector1D result(sh->nSimplex);
  
  int counterGlobal = 0;
  for(int i=0; i<size; i++){
   
    ivector1D temp = includePIDTriangulation(sh, pidList[i]);
    
    int tempSize = temp.size();
    for(int j=0; j<tempSize; ++j){
	result[counterGlobal] = temp[j];
	++counterGlobal;
    } 
  }
    result.resize(counterGlobal);
    return(result);
};

/*! Given one or more pidded part of a triangulated tesselation, return those pidded triangles outside the patch 
 * \param[in] sh target triangulated tesselation
 * \param[in] pidT target PID number
 * \param[out] result list-by-indices of pidded triangles outside the volumetric patch
 */
ivector1D BASE_PATCH::excludePIDTriangulation(SHAPE * sh, ivector1D & pidList ){

  int size = pidList.size();
  ivector1D result(sh->nSimplex);
  
  int counterGlobal = 0;
  for(int i=0; i<size; i++){
   
    ivector1D temp = excludePIDTriangulation(sh, pidList[i]);
    
    int tempSize = temp.size();
    for(int j=0; j<tempSize; ++j){
	result[counterGlobal] = temp[j];
	++counterGlobal;
    } 
  }
    result.resize(counterGlobal);
    return(result);
};

/*! Return True if at least one vertex of a given triangle is included in the volumetric patch
 * \param[in] triangleVert 3 vertices of the given Triangle
 * \param[out] result boolean
 */   
bool BASE_PATCH::isTriangleIncluded(dvecarr3E & triangleVert){
  
  bool check = false;
  for(int i=0; i<3; ++i){
   check = check || isPointIncluded(triangleVert[i]); 
  }
  return(check);
};

/*! Return True if at least one vertex of a given triangle is included in the volumetric patch
 * \param[in] tri pointer to a Class_SurfTri tesselation 
 * \param[in] indexT triangle index of tri.
 * \param[out] result boolean
 */ 
bool BASE_PATCH::isTriangleIncluded(Class_SurfTri * tri, int indexT){

  bool check = false;
  for(int i=0; i<3; ++i){
   check = check || isPointIncluded(tri, tri->Simplex[indexT][i]); 
  }
  return(check);
};




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//PATCH_CUBE IMPLEMENTATION 
/*
 *	\date			31/12/2015
 *	\authors		Gallizio Federico
 * 	\authors 		Sister Rosetta
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief Volumetric Cubic Patch class  
 *
 *	Volumetric Core Element, of cubic shape, suitable for interaction with Triangulation data Structure.
 */

/*! Basic Constructor */
PATCH_CUBE::PATCH_CUBE(){
	patchType="patchCube";
};

/*! Custom Constructor
 * \param[in] origin_ origin of the shape
 * \param[in] w_ width of the shape
 * \param[in] h_ height of the shape
 * \param[in] d_ depth of the shape
 */
PATCH_CUBE::PATCH_CUBE(darray3E & origin_, double w_, double h_, double d_):ElementalCube(origin_, w_, h_, d_){
	patchType="patchCube";
}; 

/*! Basic Destructor */
PATCH_CUBE::~PATCH_CUBE(){};

/*! Copy Constructor 
 * \param[in] other PATCH_CUBE object where copy from
 */
PATCH_CUBE::PATCH_CUBE(const PATCH_CUBE & other){
  *(static_cast<BASE_PATCH*>(this))= *(static_cast<const BASE_PATCH*> (&other));
  *(static_cast<ElementalCube*>(this))= *(static_cast<const ElementalCube*> (&other));
};

/*! Copy Operator 
 * \param[in] other PATCH_CUBE object where copy from
 */
PATCH_CUBE & PATCH_CUBE::operator=(const PATCH_CUBE & other){
  *(static_cast<BASE_PATCH*>(this))= *(static_cast<const BASE_PATCH*> (&other));
  *(static_cast<ElementalCube*>(this))= *(static_cast<const ElementalCube*> (&other));
  return(*this);
};


/*! Return True if the given point is included in the volumetric patch
 * \param[in] point given vertex
 * \param[out] result boolean
 */
bool PATCH_CUBE::isPointIncluded(darray3E point){

  bool check = true;
  darray3E temp;
  darray3E pOr = getShapeOrigin();
  dvector1D pSp = getShapeSpan();
  
  for(int i=0; i< pSp.size(); ++i){   
    temp[i] = (point[i] -pOr[i])/pSp[i];
    
    check = check && ((temp[i] >= 0.0) && (temp[i]<=1.0));
  }
  
  return(check);  
};

/*! Return True if the given point is included in the volumetric patch
 * \param[in] tri pointer to a Class_SurfTri tesselation 
 * \param[in] indexV vertex index of tri.
 * \param[out] result boolean
 */
bool PATCH_CUBE::isPointIncluded(Class_SurfTri * tri, int indexV){

  bool check = true;
  darray3E temp;
  darray3E pOr = getShapeOrigin();
  dvector1D pSp = getShapeSpan();
  
  for(int i=0; i< pSp.size(); ++i){   
    temp[i] = (tri->Vertex[indexV][i] -pOr[i])/pSp[i];
    
    check = check && ((temp[i] >= 0.0) && (temp[i]<=1.0));
  }
  
  return(check);  
};


/*! Plot cube in vtk hexahedron. WARNING: method is deprecated, use exportSurface instead.
 * \param[in] folder destination folder
 * \param[in] name   output file name
 * \param[in] number int tag for file numbering purpose 
 */
void PATCH_CUBE::plotVTUpatch(std::string folder, std::string name, int number){
    
    int intv;
    int gdl[] = {0,0,0};
    int ii[] = {0,0,0};
    int nofHEXA=1;

    std::stringstream filenameStream;
    filenameStream.clear();  
    filenameStream << folder << "/" << name << number << ".vtu";
    std::ofstream out(filenameStream.str().c_str());
    if(!out.is_open()){
       std::cout << "The vtu can't be opened. The file " << filenameStream.str() << " can't be opened or you don't have the write permissions in the folder.";
       exit(1);
    }

    nofHEXA = 1 ;
    int nofVertices = 8;

    //extract nodes	
    dvecarr3E nodes(nofVertices);
    darray3E orig_ = getShapeOrigin();
    darray3E span_;
    {
	    dvector1D dum = getShapeSpan();
	    span_= conArray<double,3>(dum);
    }
    
    dvecarr3E eleM(3, darray3E{0,0,0});
    eleM[0][0] = eleM[1][1] = eleM[2][2]=1.0;
    
    nodes[0] = orig_;
    nodes[1] = nodes[0] + eleM[0]*span_[0];
    nodes[2] = nodes[1] + eleM[1]*span_[1];
    nodes[3] = nodes[0] + eleM[1]*span_[1];
    for(int i=0; i<4; ++i){
      nodes[i+4] = nodes[i] + eleM[2]*span_[2];
    };
    
    
    out << "<?xml version=\"1.0\"?>" << std::endl
	<< "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl
	<< "  <UnstructuredGrid>" << std::endl
	<< "    <Piece NumberOfCells=\"" << nofHEXA << "\" NumberOfPoints=\"" << nofVertices << "\">" << std::endl;
    out << "      <Points>" << std::endl
	<< "        <DataArray type=\"Float32\" Name=\"Coordinates\" NumberOfComponents=\""<< 3 <<"\" format=\"ascii\">" << std::endl
	<< "          " << std::fixed;
 
     
	
    for(int i =0; i<nofVertices; i++){ 
      for(intv=0; intv<3; intv++){
	  out << std::setprecision(6) << nodes[i][intv]<< " "; 
      }                                 
      if((i+1)%4 && i!=nofVertices-1) out << std::endl << "          ";
    }
    
    out << std::endl << "        </DataArray>" << std::endl
	<< "      </Points>" << std::endl
	<< "      <Cells>" << std::endl
	<< "        <DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl
	<< "          ";
    
    // POI CI PENSO A UN OVERLOADING PER IL 2D
    for(int i =0; i<nofVertices; i++) {
      out << i << " ";
    }
    out << std::endl << "          ";
    out << std::endl << "        </DataArray>" << std::endl
	<< "        <DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl
	<< "          ";
    
    for (int i = 0; i<nofHEXA; i++){
      out << (i+1)*8 << " ";
      if((i+1)%12==0 && i!=nofHEXA-1) out << std::endl << "          ";
    }
    out << std::endl << "        </DataArray>" << std::endl
	<< "        <DataArray type=\"UInt8\" Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl
	<< "          ";
    for(int i = 0; i < nofHEXA; i++)
    {
      int type = 12; // type segment for 2D case !! improve
      out << type << " ";
      if((i+1)%12==0 && i!=nofHEXA-1)
	out << std::endl << "          ";
    }
    out << std::endl << "        </DataArray>" << std::endl
	<< "      </Cells>" << std::endl
	<< "    </Piece>" << std::endl
	<< "  </UnstructuredGrid>" << std::endl
	<< "</VTKFile>" << std::endl;
    out.close();
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//PATCH_CYLINDER IMPLEMENTATION 
/*
 *	\date			31/12/2015
 *	\authors		Gallizio Federico
 * 	\authors 		Sister Rosetta
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief Volumetric Cylindric Patch class  
 *
 *	Volumetric Core Element, of cylindric shape, suitable for interaction with Triangulation data Structure.
 */

/*! Basic Constructor */
PATCH_CYLINDER::PATCH_CYLINDER(){
	patchType="patchCylinder";
};

/*! Custom Constructor
 * \param[in] origin_ origin of the shape
 * \param[in] r_ base radius of the shape
 * \param[in] h_ height of the shape
 */
PATCH_CYLINDER::PATCH_CYLINDER(darray3E & origin_, double r_, double h_):ElementalCylinder(origin_, r_, h_){
	patchType="patchCylinder";
}; 

/*! Basic Destructor */
PATCH_CYLINDER::~PATCH_CYLINDER(){};

/*! Copy Constructor 
 * \param[in] other PATCH_CYLINDER object where copy from
 */
PATCH_CYLINDER::PATCH_CYLINDER(const PATCH_CYLINDER & other){
  *(static_cast<BASE_PATCH*>(this))= *(static_cast<const BASE_PATCH*> (&other));
  *(static_cast<ElementalCylinder*>(this))= *(static_cast<const ElementalCylinder*> (&other));
};

/*! Copy Operator 
 * \param[in] other PATCH_CYLINDER object where copy from
 */
PATCH_CYLINDER & PATCH_CYLINDER::operator=(const PATCH_CYLINDER & other){
  *(static_cast<BASE_PATCH*>(this))= *(static_cast<const BASE_PATCH*> (&other));
  *(static_cast<ElementalCylinder*>(this))= *(static_cast<const ElementalCylinder*> (&other));
  return(*this);
};
   
/*! Return True if the given point is included in the volumetric patch
 * \param[in] point given vertex
 * \param[out] result boolean
 */
bool PATCH_CYLINDER::isPointIncluded(darray3E point){

  bool check = true;
  darray3E temp;
  darray3E pOr = getShapeOrigin();
  dvector1D pSp = getShapeSpan();
  
  temp = point - pOr;
  double radius = pow(temp[0]*temp[0] + temp[1]*temp[1], 0.5) / pSp[0];
  double hnorm  = temp[2]/pSp[1];
    
  check = check && (radius <= 1.0);
  check = check && ((hnorm >=0) && (hnorm <= 1.0));
  
  return(check);  
};

/*! Return True if the given point is included in the volumetric patch
 * \param[in] tri pointer to a Class_SurfTri tesselation 
 * \param[in] indexV vertex index of tri.
 * \param[out] result boolean
 */
bool PATCH_CYLINDER::isPointIncluded(Class_SurfTri * tri, int indexV){

  bool check = true;
  darray3E temp;
  darray3E pOr = getShapeOrigin();
  dvector1D pSp = getShapeSpan();
  
  temp = tri->Vertex[indexV] - pOr;
  double radius = pow(temp[0]*temp[0] + temp[1]*temp[1], 0.5) / pSp[0];
  double hnorm  = temp[2]/pSp[1];
    
  check = check && (radius <= 1.0);
  check = check && ((hnorm >=0) && (hnorm <= 1.0));
  
  return(check);  
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//PATCH_SPHERE IMPLEMENTATION 
/*
 *	\date			31/12/2015
 *	\authors		Gallizio Federico
 * 	\authors 		Sister Rosetta
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief Volumetric Spherical Patch class  
 *
 *	Volumetric Core Element, of spherical shape, suitable for interaction with Triangulation data Structure.
 */

/*! Basic Constructor */
PATCH_SPHERE::PATCH_SPHERE(){
	patchType="patchSphere";
};

/*! Custom Constructor
 * \param[in] origin_ origin of the shape
 * \param[in] r_ radius of the shape
 */
PATCH_SPHERE::PATCH_SPHERE(darray3E & origin_, double r_) : ElementalSphere(origin_, r_){
	patchType="patchSphere";
}; 

/*! Basic Destructor */
PATCH_SPHERE::~PATCH_SPHERE(){};

/*! Copy Constructor 
 * \param[in] other PATCH_SPHERE object where copy from
 */
PATCH_SPHERE::PATCH_SPHERE(const PATCH_SPHERE & other){
 *(static_cast<BASE_PATCH*>(this))= *(static_cast<const BASE_PATCH*> (&other));
 *(static_cast<ElementalSphere*>(this))= *(static_cast<const ElementalSphere*> (&other)); 
};

/*! Copy Operator
 * \param[in] other PATCH_SPHERE object where copy from
 */
PATCH_SPHERE & PATCH_SPHERE::operator=(const PATCH_SPHERE & other){
  *(static_cast<BASE_PATCH*>(this))= *(static_cast<const BASE_PATCH*> (&other));
  *(static_cast<ElementalSphere*>(this))= *(static_cast<const ElementalSphere*> (&other)); 
  return(*this);
};
   
   
/*! Return True if the given point is included in the volumetric patch
 * \param[in] point given vertex
 * \param[out] result boolean
 */
bool PATCH_SPHERE::isPointIncluded(darray3E point){

  bool check = false;
  darray3E temp;
  darray3E pOr = getShapeOrigin();
  double radmax = getShapeSpan();
  if(radmax ==0){return(false);}
  
  temp = point - pOr;
  double radius = norm_2(temp)/radmax;
  check = (radius <= 1.0);
  
  return(check);  
};

/*! Return True if the given point is included in the volumetric patch
 * \param[in] tri pointer to a Class_SurfTri tesselation 
 * \param[in] indexV vertex index of tri.
 * \param[out] result boolean
 */
bool PATCH_SPHERE::isPointIncluded(Class_SurfTri * tri, int indexV){

  bool check = false;
  darray3E temp;
  darray3E pOr = getShapeOrigin();
  double radmax = getShapeSpan();
  if(radmax ==0){return(false);}
  
  temp = tri->Vertex[indexV] - pOr;
  double radius = norm_2(temp)/radmax;
  
  check = (radius <= 1.0);
  return(check);
};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//PATCH_HEXA IMPLEMENTATION 
/*
 *	\date			31/12/2015
 *	\authors		Gallizio Federico
 * 	\authors 		Sister Rosetta
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief Volumetric Hexahedrical Patch class  
 *
 *	Volumetric Core Element, of generic hexahedrical shape, suitable for interaction with Triangulation data Structure.
 */

/*! Basic Constructor.  Set vertices list as those belonging to an elemental unitary cube.*/
PATCH_HEXA::PATCH_HEXA(){
	patchType="patchHexa";
};

/*! Custom Constructor. Require an external list of 8 vertices ordered according the VTK policy 
 *\param[in] vertList_ list of eight vertex defining the hexahedron
 */
PATCH_HEXA::PATCH_HEXA(dvecarr3E & vertList_):ElementalHexaHedron(vertList_){
	patchType="patchHexa";
}; 
/*! Basic Destructor */
PATCH_HEXA::~PATCH_HEXA(){};

/*! Copy Constructor 
 * \param[in] other PATCH_HEXA object where copy from
 */    
PATCH_HEXA::PATCH_HEXA(const PATCH_HEXA & other){
 *(static_cast<BASE_PATCH*>(this))= *(static_cast<const BASE_PATCH*> (&other));
 *(static_cast<ElementalHexaHedron*>(this))= *(static_cast<const ElementalHexaHedron*> (&other)); 
};

/*! Copy Operator 
 * \param[in] other PATCH_HEXA object where copy from
 */
PATCH_HEXA & PATCH_HEXA::operator=(const PATCH_HEXA & other){
 *(static_cast<BASE_PATCH*>(this))= *(static_cast<const BASE_PATCH*> (&other));
 *(static_cast<ElementalHexaHedron*>(this))= *(static_cast<const ElementalHexaHedron*> (&other)); 
  return(*this);
};

/*! Return True if the given point is included in the volumetric patch
 * \param[in] point given vertex
 * \param[out] result boolean
 */
bool PATCH_HEXA::isPointIncluded(darray3E point){
  
  bool check = true;
  dvecarr4E planes = getHexaPlanes();
  int sizeP = planes.size();
  for(int i=0; i<sizeP; ++i){
    double val = planes[i][0]*point[0]+planes[i][1]*point[1] + planes[i][2]*point[2] + planes[i][3];
    check = check && (val >=0);
  }
  
  return(check);
};

/*! Return True if the given point is included in the volumetric patch
 * \param[in] tri pointer to a Class_SurfTri tesselation 
 * \param[in] indexV vertex index of tri.
 * \param[out] result boolean
 */
bool PATCH_HEXA::isPointIncluded(Class_SurfTri * tri, int indexV){

  bool check = true;
  darray3E point = tri->Vertex[indexV];
  dvecarr4E planes = getHexaPlanes();
  int sizeP = planes.size();
  for(int i=0; i<sizeP; ++i){
    double val = planes[i][0]*point[0]+planes[i][1]*point[1] + planes[i][2]*point[2] + planes[i][3];
    check = check && (val >=0);
  }
  
  return(check);

};

/*! Plot hexahedron in vtk format . WARNING: method is deprecated, use exportSurface instead.
 * \param[in] folder destination folder
 * \param[in] name   output file name
 * \param[in] number tag for file numbering purpose 
 */
void PATCH_HEXA::plotVTUpatch(std::string folder, std::string name, int number){
    
    int intv;
    int gdl[] = {0,0,0};
    int ii[] = {0,0,0};
    int nofHEXA=1;

    std::stringstream filenameStream;
    filenameStream.clear();  
    filenameStream << folder << "/" << name << number << ".vtu";
    std::ofstream out(filenameStream.str().c_str());
    if(!out.is_open()){
       std::cout << "The vtu can't be opened. The file " << filenameStream.str() << " can't be opened or you don't have the write permissions in the folder.";
       exit(1);
    }

    nofHEXA = 1 ;
    dvecarr3E nodes = getHexaVertices();	
    int nofVertices = nodes.size();

    out << "<?xml version=\"1.0\"?>" << std::endl
	<< "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl
	<< "  <UnstructuredGrid>" << std::endl
	<< "    <Piece NumberOfCells=\"" << nofHEXA << "\" NumberOfPoints=\"" << nofVertices << "\">" << std::endl;
    out << "      <Points>" << std::endl
	<< "        <DataArray type=\"Float32\" Name=\"Coordinates\" NumberOfComponents=\""<< 3 <<"\" format=\"ascii\">" << std::endl
	<< "          " << std::fixed;
 
   	
    for(int i =0; i<nofVertices; i++){ 
      for(intv=0; intv<3; intv++){
	  out << std::setprecision(6) << nodes[i][intv]<< " "; 
      }                                 
      if((i+1)%4 && i!=nofVertices-1) out << std::endl << "          ";
    }
    
    out << std::endl << "        </DataArray>" << std::endl
	<< "      </Points>" << std::endl
	<< "      <Cells>" << std::endl
	<< "        <DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl
	<< "          ";
    
    // POI CI PENSO A UN OVERLOADING PER IL 2D
    for(int i =0; i<nofVertices; i++) {
      out << i << " ";
    }
    out << std::endl << "          ";
    out << std::endl << "        </DataArray>" << std::endl
	<< "        <DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl
	<< "          ";
    
    for (int i = 0; i<nofHEXA; i++){
      out << (i+1)*8 << " ";
      if((i+1)%12==0 && i!=nofHEXA-1) out << std::endl << "          ";
    }
    out << std::endl << "        </DataArray>" << std::endl
	<< "        <DataArray type=\"UInt8\" Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl
	<< "          ";
    for(int i = 0; i < nofHEXA; i++)
    {
      int type = 12; // type segment for 2D case !! improve
      out << type << " ";
      if((i+1)%12==0 && i!=nofHEXA-1)
	out << std::endl << "          ";
    }
    out << std::endl << "        </DataArray>" << std::endl
	<< "      </Cells>" << std::endl
	<< "    </Piece>" << std::endl
	<< "  </UnstructuredGrid>" << std::endl
	<< "</VTKFile>" << std::endl;
    out.close();
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//HULL_CUBE IMPLEMENTATION 
/*
 *	\date			31/12/2015
 *	\authors		Gallizio Federico
 * 	\authors 		Sister Rosetta
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief Volumetric Cubic HUll class  
 *
 *	Cartesian Structured Mesh, of cubic shape, suitable for interaction with Triangulation data Structure.
 */
   
/*! Basic Constructor */
HULL_CUBE::HULL_CUBE(){
	patchType="hullCube";
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
HULL_CUBE::HULL_CUBE(darray3E origin_, double spanX_, double spanY_, double spanZ_, int nx_, int ny_, int nz_):
		      UCubicMesh(origin_,spanX_,spanY_, spanZ_, nx_, ny_, nz_){
			      patchType="hullCube";
		}; 

/*! Basic Destructor */
HULL_CUBE::~HULL_CUBE(){};

/*! Copy Constructor 
 * \param[in] other HULL_CUBE object where copy from
 */
HULL_CUBE::HULL_CUBE(const HULL_CUBE & other){
  *(static_cast<BASE_PATCH*>(this))= *(static_cast<const BASE_PATCH*> (&other));
  *(static_cast<UCubicMesh*>(this))= *(static_cast<const UCubicMesh*> (&other));
};

/*! Copy Operator 
 * \param[in] other PATCH_CUBE object where copy from
 */
HULL_CUBE & HULL_CUBE::operator=(const HULL_CUBE & other){
  *(static_cast<BASE_PATCH*>(this))= *(static_cast<const BASE_PATCH*> (&other));
  *(static_cast<UCubicMesh*>(this))= *(static_cast<const UCubicMesh*> (&other));
  return(*this);
};


/*! Return True if the given point is included in the volumetric patch
 * \param[in] point given vertex
 * \param[out] result boolean
 */
bool HULL_CUBE::isPointIncluded(darray3E point){

  bool check = true;
  darray3E temp = transfToLocal(point);
  darray3E pSp = getSpan();
  
  for(int i=0; i< pSp.size(); ++i){   
    temp[i] = temp[i]/pSp[i];
    
    check = check && ((temp[i] >= 0.0) && (temp[i]<=1.0));
  }
  
  return(check);  
};

/*! Return True if the given point is included in the volumetric patch
 * \param[in] tri pointer to a Class_SurfTri tesselation 
 * \param[in] indexV vertex index of tri.
 * \param[out] result boolean
 */
bool HULL_CUBE::isPointIncluded(Class_SurfTri * tri, int indexV){

  bool check = true;
  darray3E temp = transfToLocal(tri->Vertex[indexV]);
  darray3E pSp = getSpan();
  
  for(int i=0; i< pSp.size(); ++i){   
    temp[i] = temp[i]/pSp[i];
    
    check = check && ((temp[i] >= 0.0) && (temp[i]<=1.0));
  }
  
  return(check);  
};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//HULL_CYLINDER IMPLEMENTATION 
/*
 *	\date			31/12/2015
 *	\authors		Gallizio Federico
 * 	\authors 		Sister Rosetta
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief Volumetric Cylindric HUll class  
 *
 *	Structured Mesh, of cylindrical shape and reference frame, suitable for interaction with Triangulation data Structure.
 */

/*! Basic Constructor */
HULL_CYLINDER::HULL_CYLINDER(){
	patchType="hullCylinder";
};

/*!Custom Constructor. Set mesh origin, span size, dimension and nodal data structure
 * \param[in] origin_ point origin in global reference system
 * \param[in] spanR_ max base radius span: must be >0; 
 * \param[in] spanZ_ max cylinder height span: must be >0; 
   \param[in] thetalim lower and upper limits on the angular coordinate 
 * \param[in] nr_   number of cells in r-direction
 * \param[in] nt_   number of cells in theta-direction
 * \param[in] nz_   number of cells in z-direction
 */
HULL_CYLINDER::HULL_CYLINDER(darray3E origin_, double spanR_, double spanZ_, dvector1D & thetalim_, int nr_, int nt_, int nz_):
		      UCylindricalMesh(origin_, spanR_, spanZ_, thetalim_, nr_, nt_, nz_){
			      patchType="hullCylinder";
		}; 

/*! Basic Destructor */
HULL_CYLINDER::~HULL_CYLINDER(){};

/*! Copy Constructor 
 * \param[in] other HULL_CYLINDER object where copy from
 */
HULL_CYLINDER::HULL_CYLINDER(const HULL_CYLINDER & other){
  *(static_cast<BASE_PATCH*>(this))= *(static_cast<const BASE_PATCH*> (&other));
  *(static_cast<UCylindricalMesh*>(this))= *(static_cast<const UCylindricalMesh*> (&other));
};

/*! Copy Operator 
 * \param[in] other PATCH_CYLINDER object where copy from
 */
HULL_CYLINDER & HULL_CYLINDER::operator=(const HULL_CYLINDER & other){
  *(static_cast<BASE_PATCH*>(this))= *(static_cast<const BASE_PATCH*> (&other));
  *(static_cast<UCylindricalMesh*>(this))= *(static_cast<const UCylindricalMesh*> (&other));
  return(*this);
};


/*! Return True if the given point is included in the volumetric patch
 * \param[in] point given vertex
 * \param[out] result boolean
 */
bool HULL_CYLINDER::isPointIncluded(darray3E point){

  bool check = true;
  darray3E pOr = getOrigin();
  darray3E temp2 = transfToLocal(point);
  darray3E pSp = getSpan();
  
    double radius_norm = temp2[0]/pSp[0];
    double height_norm = temp2[2]/pSp[2];
    
    check = check && (radius_norm <=1.0);
    check = check && ((height_norm >= 0.0) && (height_norm<=1.0));
  
  return(check);  
};

/*! Return True if the given point is included in the volumetric patch
 * \param[in] tri pointer to a Class_SurfTri tesselation 
 * \param[in] indexV vertex index of tri.
 * \param[out] result boolean
 */
bool HULL_CYLINDER::isPointIncluded(Class_SurfTri * tri, int indexV){

 bool check = true;
  darray3E pOr = getOrigin();
  darray3E temp2 = transfToLocal(tri->Vertex[indexV]);
  darray3E pSp = getSpan();
  
    double radius_norm = temp2[0]/pSp[0];
    double height_norm = temp2[2]/pSp[2];
    
    check = check && (radius_norm <=1.0);
    check = check && (temp2[1]>=0 && temp2[1]<=pSp[1]);
    check = check && ((height_norm >= 0.0) && (height_norm<=1.0));
  
  return(check);  
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//HULL_SPHERE IMPLEMENTATION 
/*
 *	\date			31/12/2015
 *	\authors		Gallizio Federico
 * 	\authors 		Sister Rosetta
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief Volumetric Cylindric HUll class  
 *
 *	Structured Mesh, of cylindrical shape and reference frame, suitable for interaction with Triangulation data Structure.
 */

/*! Basic Constructor */
HULL_SPHERE::HULL_SPHERE(){
	patchType="hullSphere";
};

/*! Custom Constructor. Set mesh origin, span size, dimension and nodal data structure
 * \param[in] origin_ point origin in global reference system
 * \param[in] spanR_  max radius span: must be >0;
 * \param[in] thetalim_ lower and upper limits on the polar coordinate
 * \param[in] philim_ lower and upper limits on the azimuthal coordinate
 * \param[in] nr_   number of cells in r-direction
 * \param[in] nt_   number of cells in theta-direction
 * \param[in] np_   number of cells in phi-direction
 */
HULL_SPHERE::HULL_SPHERE(darray3E origin_, double spanR_, dvector1D thetalim_, dvector1D philim_, int nr_, int nt_, int np_):
			USphericalMesh(origin_, spanR_, thetalim_, philim_, nr_, nt_, np_){
				patchType="hullSphere";
			}; 

/*! Basic Destructor */			
HULL_SPHERE::~HULL_SPHERE(){};

/*! Copy Constructor 
 * \param[in] other HULL_SPHERE object where copy from
 */
HULL_SPHERE::HULL_SPHERE(const HULL_SPHERE & other){
  *(static_cast<BASE_PATCH*>(this))= *(static_cast< const BASE_PATCH*> (&other));
  *(static_cast<USphericalMesh*>(this))= *(static_cast<const USphericalMesh*> (&other));  
};

/*! Copy Operator 
 * \param[in] other HULL_SPHERE object where copy from
 */
HULL_SPHERE & HULL_SPHERE::operator=(const HULL_SPHERE & other){
  *(static_cast<BASE_PATCH*>(this))= *(static_cast<const BASE_PATCH*> (&other));
  *(static_cast<USphericalMesh*>(this))= *(static_cast<const USphericalMesh*> (&other));
   return(*this);  
};   

/*! Return True if the given point is included in the volumetric patch
 * \param[in] point given vertex
 * \param[out] result boolean
 */
bool HULL_SPHERE::isPointIncluded(darray3E point){
   bool check = false;
   darray3E temp2 = transfToLocal(point);
   darray3E pSp = getSpan();
  
    double radius_norm = temp2[0]/pSp[0];
    
    check = (radius_norm <=1.0);
  
  return(check);  
};

/*! Return True if the given point is included in the volumetric patch
 * \param[in] tri pointer to a Class_SurfTri tesselation 
 * \param[in] indexV vertex index of tri.
 * \param[out] result boolean
 */
bool HULL_SPHERE::isPointIncluded(Class_SurfTri * tri, int indexV){
   
  bool check = false;
   darray3E temp2 = transfToLocal(tri->Vertex[indexV]);
   darray3E pSp = getSpan();
  
    double radius_norm = temp2[0]/pSp[0];
    check = (radius_norm <=1.0);
    check = check && (temp2[1]>=0 && temp2[1]<=pSp[1]);
    check = check && (temp2[2]>=0 && temp2[2]<=pSp[2]);
    
  return(check);
  
};
