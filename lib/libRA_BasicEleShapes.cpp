#include "libRA_BasicEleShapes.hpp"

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// BASE_ElementalShapes IMPLEMENTATION 
/*
*	\date			31/12/2015
*	\authors		Gallizio Federico
* 	\authors 		Ratzinger Joseph
*	\authors		Arpa Rocco
*	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
*	\par			License:\n
*	This class version is released under .
*
*	\brief Virtual class for Basic Elemental Shapes Representation 
*
*	Base class for elemental shapes definition,from cubes and hexahedra to cylinders and spheres
*/

/*! Base Constructor */	    
BASE_ElementalShapes::BASE_ElementalShapes(){
	shapeOrigin.fill(0.0);
	surfaceTessType = 0;
	surface = NULL;
	volume = NULL;
	eleShapeType="";
};

/*!Custom Constructor 
* \param[in] origin_ set the origin of the shape
*/	    
BASE_ElementalShapes::BASE_ElementalShapes(darray3E & origin_){
shapeOrigin=origin_;
surfaceTessType = 0;
surface = NULL;
volume = NULL;
eleShapeType="";
};

/*! Base Destructor */
BASE_ElementalShapes::~BASE_ElementalShapes(){
	if(surface!=NULL) {delete surface; surface=NULL;}
	if(volume != NULL) {delete volume; volume=NULL;}
};

/*! Copy Constructor. 
* \param[in] other BASE_ElementalShapes object to copyright
*/
BASE_ElementalShapes::BASE_ElementalShapes(const BASE_ElementalShapes & other){
shapeOrigin = other.shapeOrigin;
eleShapeType = other.eleShapeType;
surfaceTessType = 0;
if(other.surface != NULL){
surface = new Class_SurfTri;
*surface = *(other.surface);
surfaceTessType = other.surfaceTessType;
}
if(other.volume != NULL){
volume = new Class_VolTri;
*volume = *(other.volume);
}
};
/*! Copy Operator. 
* \param[in] other BASE_ElementalShapes object to copyright
*/
BASE_ElementalShapes & BASE_ElementalShapes::operator=(const BASE_ElementalShapes & other){

clearShape();
shapeOrigin = other.shapeOrigin;
eleShapeType = other.eleShapeType;
surfaceTessType = 0;
if(other.surface != NULL){
surface = new Class_SurfTri;
*surface = *(other.surface);
surfaceTessType = other.surfaceTessType;
}
if(other.volume != NULL){
volume = new Class_VolTri;
*volume = *(other.volume);
}

return(*this);
};

/*! Clear Class Structures */
void BASE_ElementalShapes::clearShape(){
shapeOrigin.fill(0.0);  
if(surface!=NULL) {delete surface; surface=NULL; surfaceTessType = 0;}
if(volume != NULL) {delete volume; volume=NULL;}
};

/*! Set origin of your shape
* \param[in] vx origin x coordinate
* \param[in] vy origin y coordinate
* \param[in] vz origin z coordinate
*/
void BASE_ElementalShapes::setShapeOrigin(double vx, double vy, double vz){
shapeOrigin[0] = vx;
shapeOrigin[1] = vy;
shapeOrigin[2] = vz;
};

/*! Get origin of your shape
* \param[out] result origin
*/
darray3E BASE_ElementalShapes::getShapeOrigin(){
return(shapeOrigin);
}

/*! Get your classType
* \param[out] result current class type identifying string;
*/
std::string BASE_ElementalShapes::getClassType(){return(eleShapeType);};

/*! Return your tesselated surface, by Class_SurfTri object pointer */
Class_SurfTri * BASE_ElementalShapes::getSurfaceTess(){return(surface);};

/*! Return your tesselated volume, by Class_VolTri object pointer */
Class_VolTri * BASE_ElementalShapes::getVolumeTess(){return(volume);};

/*! Check your current surface tesselation. 
* If triangulated return 1. 
* if generically tessellated return 2. 
* if none return 0 .
*/
short int BASE_ElementalShapes::checkTessellation(){return(surfaceTessType);};

/*! Export your current surface triangulation in *.stl format. If you surface member is NULL does nothing.
* If your surface is not triangulated or generically tesselated w/ class method tessellate does nothing 
* \param[in] name complete output filename w/out tag
* \param[in] codex current output codex: "ascii" or "binary"
*/
void BASE_ElementalShapes::exportSurfaceStl(std::string name, std::string codex){

bool ping = (codex=="binary");
std::string path = name + ".stl"; 
surface->Export_stl(path, ping);

};

/*! Export your current surface triangulation in *.vtu format. If you surface member is NULL does nothing.
* If your surface is not generically tesselated w/ class method tessellate does nothing 
* \param[in] name complete output filename w/out tag
*/
void BASE_ElementalShapes::exportSurface(std::string name){

std::string path = name + ".vtu"; 
surface->Export_vtu(path);

};

/*! Export your current volume tetrahedral mesh in *.vtu format. If you volume member is NULL does nothing.
* \param[in] name complete output filename w/out tag
*/
void BASE_ElementalShapes::exportVolume(std::string name){

std::string path = name + ".vtu"; 
volume->Export_vtu(path);
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// ElementalCube IMPLEMENTATION 
/*
*	\date			31/12/2015
*	\authors		Gallizio Federico
* 	\authors 		Ratzinger Joseph
*	\authors		Arpa Rocco
*	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
*	\par			License:\n
*	This class version is released under .
*
*	\brief Cube Class
*
*	Elemental shape Cube class, derived from BASE_ElementalShapes
*/
/*! Basic Constructor */
ElementalCube::ElementalCube(){
height = 0.0; 
width = 0.0;
depth = 0.0;
eleShapeType = "Cube";
};

/*! Custom Constructor. Set dimension of the Cube:
* \param[in] w_ width
* \param[in] h_ height
* \param[in] d_ depth
*/
ElementalCube::ElementalCube(darray3E origin_, double w_, double h_, double d_):BASE_ElementalShapes(origin_){
	setShapeSpan(w_,h_, d_);       
	eleShapeType = "Cube";
};

/*! Class Destructor */
ElementalCube::~ElementalCube(){eleShapeType="";};

/*! Copy Constructor.
* \param[in] other Elemental Cube object which need to be copied
*/
ElementalCube::ElementalCube(const ElementalCube & other){
	
	*(static_cast<BASE_ElementalShapes *> (this)) = *(static_cast<const BASE_ElementalShapes * >(&other));
	height = other.height;
	width = other.width;
	depth = other.depth;
};

/*! Copy Operator.
* \param[in] other Elemental Cube object which need to be copied
*/
ElementalCube & ElementalCube::operator=(const ElementalCube & other){

	*(static_cast<BASE_ElementalShapes *> (this)) = *(static_cast<const BASE_ElementalShapes * >(&other));
	height = other.height;
	width = other.width;
	depth = other.depth;
	return(*this);
};

/*! Clear Class Structures */
void ElementalCube::clearShape(){
	shapeOrigin.fill(0.0);  
	if(surface!=NULL) {delete surface; surface=NULL; surfaceTessType = 0;}
	if(volume != NULL) {delete volume; volume=NULL;}
	height = 0.0;
	width = 0.0;
	depth = 0.0;
}; 


/*! Set widht, height and depth of your cube
* \param[in] w_ width
* \param[in] h_ height
* \param[in] d_ depth
*/  
void ElementalCube::setShapeSpan(double w_, double h_, double d_){
	height = h_; 
	width = w_;	
	depth = d_;
};

/*! Return characteristic dimension of the cube.
* \param[out] result width, height, depth
*/
dvector1D ElementalCube::getShapeSpan(){
dvector1D result(3,0);
result[0] = width;
result[1] = height;
result[2] = depth;
return(result);
};

void ElementalCube::tessellate(int ){
//doing nothing for now, until i learn a smart way to tesselate
};
void ElementalCube::triangulate(int ){};

void ElementalCube::tetrahedralize(int ){};


/////////////////////////////////////////////////////////////////////////////////////////////////////////
// ElementalCylinder IMPLEMENTATION 
/*
*	\date			31/12/2015
*	\authors		Gallizio Federico
* 	\authors 		Ratzinger Joseph
*	\authors		Arpa Rocco
*	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
*	\par			License:\n
*	This class version is released under .
*
*	\brief Cylinder Class
*
*	Elemental shape Cylinder class, derived from BASE_ElementalShapes
*/

/*! Basic Constructor */
ElementalCylinder::ElementalCylinder(){
radius = 0.0;
height = 0.0;
eleShapeType = "Cylinder";
};

/*! Custom Constructor. Set dimension of the Cylinder:
* \param[in] r_ base radius
* \param[in] h_ height
*/
ElementalCylinder::ElementalCylinder(darray3E origin_, double r_, double h_):BASE_ElementalShapes(origin_){
setShapeSpan(r_, h_);
eleShapeType = "Cylinder";
};

/*! Class Destructor */
ElementalCylinder::~ElementalCylinder(){eleShapeType = "";};

/*! Copy Constructor.
* \param[in] other Elemental Cylinder object which need to be copied
*/
ElementalCylinder::ElementalCylinder(const ElementalCylinder & other){

	*(static_cast<BASE_ElementalShapes *> (this)) = *(static_cast<const BASE_ElementalShapes * >(&other));
	radius = other.radius;
	height = other.height;

};

/*! Copy Operator.
* \param[in] other Elemental Cube object which need to be copied
*/
ElementalCylinder & ElementalCylinder::operator=(const ElementalCylinder & other){
	
	*(static_cast<BASE_ElementalShapes *> (this)) = *(static_cast<const BASE_ElementalShapes * >(&other));
	radius = other.radius;
	height = other.height;
	return(*this);
};

/*! Clear Class Structures */
void ElementalCylinder::clearShape(){
	shapeOrigin.fill(0.0);  
	if(surface!=NULL) {delete surface; surface=NULL; surfaceTessType = 0;}
	if(volume != NULL) {delete volume; volume=NULL;}
	height = 0.0;
	radius = 0.0;
	
};

/*! Set base radius and height  of your cylinder
* \param[in] r_ base radius
* \param[in] h_ height
*/  
void ElementalCylinder::setShapeSpan(double r_, double h_){
radius = r_;
height = h_;
};

/*! Return characteristic dimension of cylinder.
* \param[out] result base radius, height
*/
dvector1D ElementalCylinder::getShapeSpan(){

dvector1D result(2,0);
result[0] = radius;
result[1] = height;
return(result);
};

void ElementalCylinder::tessellate(int ){};
void ElementalCylinder::triangulate(int ){};
void ElementalCylinder::tetrahedralize(int ){};

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// ElementalSphere IMPLEMENTATION 
/*
*	\date			31/12/2015
*	\authors		Gallizio Federico
* 	\authors 		Ratzinger Joseph
*	\authors		Arpa Rocco
*	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
*	\par			License:\n
*	This class version is released under .
*
*	\brief Sphere Class
*
*	Elemental shape Sphere class, derived from BASE_ElementalShapes
*/

/*! Basic Constructor */
ElementalSphere::ElementalSphere(){
radius = 0.0;
eleShapeType = "Sphere";
};

/*! Custom Constructor. Set dimension of the Cylinder:
* \param[in] r_ radius
* \param[in] h_ height
*/
ElementalSphere::ElementalSphere(darray3E origin_, double r_):BASE_ElementalShapes(origin_){
setShapeSpan(r_);
eleShapeType = "Sphere";
};

/* Class Destructor */
ElementalSphere::~ElementalSphere(){eleShapeType="";};

/*! Copy Constructor.
* \param[in] other Elemental Sphere object which need to be copied
*/
ElementalSphere::ElementalSphere(const ElementalSphere & other){
	*(static_cast<BASE_ElementalShapes *> (this)) = *(static_cast<const BASE_ElementalShapes * >(&other));
	radius = other.radius;
};

/*! Copy Operator.
* \param[in] other Elemental Sphere object which need to be copied
*/
ElementalSphere & ElementalSphere::operator=(const ElementalSphere & other){
	*(static_cast<BASE_ElementalShapes *> (this)) = *(static_cast<const BASE_ElementalShapes * >(&other));
	radius = other.radius;
	return(*this);
};

/*! Clear Class Structures */
void ElementalSphere::clearShape(){
	shapeOrigin.fill(0.0);  
	if(surface!=NULL) {delete surface; surface=NULL; surfaceTessType = 0;}
	if(volume != NULL) {delete volume; volume=NULL;}
	radius = 0.0;
};

/*! Set sphere radius 
* \param[in] r_ radius
*/
void ElementalSphere::setShapeSpan(double r_){
radius= r_;
};

/*! Return characteristic dimension of the sphere
* \param[out] result radius
*/
double ElementalSphere::getShapeSpan(){return(radius);};

void ElementalSphere::tessellate(int ){};
void ElementalSphere::triangulate(int ){};
void ElementalSphere::tetrahedralize(int ){};

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// ElementalHexaHedron IMPLEMENTATION 
/*    
*	\date			31/12/2015
*	\authors		Gallizio Federico
* 	\authors 		Ratzinger Joseph
*	\authors		Arpa Rocco
*	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
*	\par			License:\n
*	This class version is released under .
*
*	\brief Generic HexaHedron Class
*
*	Elemental shape HeaxaHedron class, derived from BASE_ElementalShapes
*/

/*Basic Constructor. Set vertices list as those belonging to an elemental unitary cube. */
ElementalHexaHedron::ElementalHexaHedron(){

vertex.resize(8, darray3E{0,0,0}); 
//fill vertex list to unitary cube;
vertex[1][0] = vertex[5][0] = 1.0;
vertex[3][1] = vertex[7][1] = 1.0;
vertex[2][0] = vertex[2][1] = vertex[6][0] = vertex[6][1] = 1.0;
vertex[4][2] = vertex[5][2] = vertex[6][2] = vertex[7][2] = 1.0;

//set planes
setFaceConnectivity();
setVertConnectivity();
setPlanes();
eleShapeType = "HexaHedron";
};

/*! Custom Constructor. Require an external list of 8 vertices ordered according the VTK policy 
*\param[in] vertList list of eight vertex defining the hexahedron
*/
ElementalHexaHedron::ElementalHexaHedron( dvecarr3E & vertList){

//fill vertex list to your custom list;
vertex = vertList;

if(vertex.size() != 8) {
	std::cout<<" Less or More than 8 vertices defined in your custom vertex List. Errors can occurr."<<endl;
}	

setFaceConnectivity();
setVertConnectivity();
setShapeOrigin(vertex[0][0], vertex[0][1], vertex[0][2]);

//setShapeOrigin already call setPlanes.
eleShapeType = "HexaHedron";
};

/*! Class Destructor */
ElementalHexaHedron::~ElementalHexaHedron(){
//destroy all structure
	if(surface!=NULL) {delete surface; surface=NULL; surfaceTessType = 0;}
	if(volume != NULL) {delete volume; volume=NULL;}
	eleShapeType = "";
	resetHexaHedron();
};

/*! Copy Constructor.
* \param[in] other Elemental Hexahedron object which need to be copied
*/
ElementalHexaHedron::ElementalHexaHedron(const ElementalHexaHedron & other){
	*(static_cast<BASE_ElementalShapes *> (this)) = *(static_cast<const BASE_ElementalShapes * >(&other));
	resetHexaHedron();
	vertex = other.vertex;
	planes = other.planes;
	face_conn = other.face_conn;
	vert_conn = other.vert_conn;
};

/*! Copy Operator.
* \param[in] other Elemental Hexahedron object which need to be copied
*/
ElementalHexaHedron & ElementalHexaHedron::operator=(const ElementalHexaHedron & other){
	*(static_cast<BASE_ElementalShapes *> (this)) = *(static_cast<const BASE_ElementalShapes * >(&other));
	resetHexaHedron();
	vertex = other.vertex;
	planes = other.planes;
	face_conn = other.face_conn;
	vert_conn = other.vert_conn;
	return(*this);
};

/*! Clear member structure of the class */  
void ElementalHexaHedron::clearShape(){
	if(surface!=NULL) {delete surface; surface=NULL; surfaceTessType = 0;}
	if(volume != NULL) {delete volume; volume=NULL;}

	resetHexaHedron();
};

/*! Reset all vertex and connectivity of the Hexahedron*/
void ElementalHexaHedron::resetHexaHedron(){
	shapeOrigin.fill(0.0);  
	freeContainer(vertex); 
	freeContainer(planes);
	freeContainer(face_conn);
	freeContainer(vert_conn);
}


/*! Set origin of your shape
* \param[in] vx origin x coordinate
* \param[in] vy origin y coordinate
* \param[in] vz origin z coordinate
*/
void ElementalHexaHedron::setShapeOrigin(double vx, double vy, double vz){
	
	if(vertex.size() == 0){return;}
	
	shapeOrigin[0] = vx;
	shapeOrigin[1] = vy;
	shapeOrigin[2] = vz;
	vertex[0] = shapeOrigin;
	setPlanes();
};

/*! Set the list of 8 vertices forming the HexaHedron, in vtk ordering format. If no vertex list is provided, 
* vertices are set to those of an elemental cube.
* \param[in] list OPTIONAL. Hexahedron vertex list.
*/ 
void ElementalHexaHedron::setBaseVertex(dvecarr3E * list){

if(list == NULL){

vertex.resize(8, darray3E{0,0,0}); 
//fill vertex list to unitary cube;
vertex[1][0] = vertex[5][0] = 1.0;
vertex[3][1] = vertex[7][1] = 1.0;
vertex[2][0] = vertex[2][1] = vertex[6][0] = vertex[6][1] = 1.0;
vertex[4][2] = vertex[5][2] = vertex[6][2] = vertex[7][2] = 1.0;

}else{
	freeContainer(vertex);
	if(list->size() != 8) {
		std::cout<<" Less or More than 8 vertices defined in your custom vertex List. Errors can occurr."<<endl;
		return;
	}
	vertex = *list;
}
if(face_conn.size() ==0) setFaceConnectivity();
if(vert_conn.size() ==0) setVertConnectivity();
//check Vertex size
setShapeOrigin(vertex[0][0], vertex[0][1], vertex[0][2]);

};

/*! Change a vertex in the prexistent vertex list. Manual check coplanarity of new point and vertices reset.
 * If index is 0, the vertex is treated as a new shape origin --> see setShapeOrigin method (coplanarity check is enforced)
* \param[in] index index point to 8 element list position
* \param[in] newPoint new point to insert
* \param[in] apply  boolean, if true activate automatic coplanarity check and reset hexahedron vertices eventually.
*/
void ElementalHexaHedron::changeVertex(int index, darray3E & newPoint, bool apply){
if(vertex.size() != 8){std::cout<<"No suitable list of vertices found. Returning"<<endl; return;}
if(index <0 || index >=8) {std::cout<<"No suitable list index provided. Returning"<<endl; return;}
if(index ==0 && !apply){std::cout<<"Changing the shape origin. Automatically apply your modifications"<<endl;}

	if(index ==0){
	setShapeOrigin(newPoint[0], newPoint[1], newPoint[2]);
	}else{
		vertex[index] = newPoint;
		if(apply){ setPlanes();}
	}
};

/*! Change a series of vertices in the prexistent vertex list,check coplanarity of new point and reset vertices.
 * If one of indices is 0, the vertex is treated as a new shape origin --> see setShapeOrigin method (coplanarity check is enforced)
 * \param[in] indexList list of index point to 8 element position list
 * \param[in] newPoints set of new points to insert
 */
void ElementalHexaHedron::changeVertices(ivector1D & indexList, dvecarr3E & newPoints){

	dvector1D newOr_;
	for(int i=0; i<indexList.size(); ++i){
		if(indexList[i] ==0){
			newOr_ = conVect(newPoints[i]);
		}else{
			changeVertex(indexList[i], newPoints[i], false);
		}
	}
	
	if(newOr_.size() !=0){
		setShapeOrigin(newOr_[0], newOr_[1], newOr_[2]);
	}else{
		setPlanes();
	}
	
	
};



/*! Get the vertex List of HexaHedron
*\param[out] result vertex list
*/
dvecarr3E ElementalHexaHedron::getHexaVertices(){return(vertex);};
/*! Get the planes list of each HexaHedron Face
*\param[out] result planes list
*/
dvecarr4E ElementalHexaHedron::getHexaPlanes(){return(planes);};
/*! Return the Face-vertex connectivity 
*\param[out] result connectivity
*/
ivector2D ElementalHexaHedron::getFaceConnectivity(){return(face_conn);};
/*! Return the Vertex-Face connectivity 
*\param[out] result connectivity
*/
ivector2D ElementalHexaHedron::getVertConnectivity(){return(vert_conn);};

/*! Get the HexaHedron vertex @ prescribed index
* \param[in] index vertex index
* \param[out] result point
*/
darray3E ElementalHexaHedron::getHexaVertex(int index){
if(index <0 || index > 7){return(darray3E{0,0,0});}
return(vertex[index]);
};

/*! Get Plane implicit equation coefficients A*x+B*y+C*z+D=0, given a certan face index. 
*  Faces are ordered according to the scheme:
*  0-West, 1-South, 2-Bottom, 3-East, 4-North, 5-Top.
* \param[in] index face index [0,5];
* \param[out] result implicit plane coefficient.
*/
darray4E ElementalHexaHedron::getPlane(int index){
if(index <0 || index > 5){return(darray4E{0,0,0,0});}
return(planes[index]);
};

/*! Get Normal to a prescribed face, given its index
*  Faces are ordered according to the scheme:
*  0-West, 1-South, 2-Bottom, 3-East, 4-North, 5-Top.
* \param[in] index face index [0,5];
* \param[out] result normal to plane/face.
*/
darray3E ElementalHexaHedron::getPlaneNormal(int index){
if(index <0 || index > 5){return(darray3E{0,0,0});}
darray3E result;
for(int i=0; i<3; ++i){result[i] = planes[index][i];}
return(result);
};

/*! Get Face connectivity of a vertex, given its index 
*\param[in] index vertex index [0,7];
* \param[out] result faces connected.
*/
ivector1D ElementalHexaHedron::getFaceConn(int index){
if(index <0 || index > 7){return(ivector1D());}
return(face_conn[index]);
};

/*! Get Vertex connectivity of a face, given its index 
*\param[in] index face index [0,5];
* \param[out] result vertex connected.
*/
ivector1D ElementalHexaHedron::getVertConn(int index){
if(index <0 || index > 5){return(ivector1D());}
return(vert_conn[index]);
};

/*! check if a point belongs to a given face of the HexaHedron 
* \param[in] point target point
* \param[in] planeIndex index fo the HexaHedra face [0,5]
* \param[out] result boolean check 
*/ 
bool ElementalHexaHedron::checkBelongPolygon(darray3E & point, int planeIndex){

// workspace
bool result = false;
double pi   = 4.0*atan(1.0);
double sum_angle = 0.0;
double dista,distb,distc, sign;
int i,j; 
darray3E normdum;
for(i=0; i<4; i++){ // counting on local vertices
	if( point == vertex[face_conn[planeIndex][i]]){
		sum_angle +=pi; 
		continue;
	} 

	j = (i+1)%4;
	
	dista = norm_2((point-vertex[face_conn[planeIndex][i]]));
	distb = norm_2((point-vertex[face_conn[planeIndex][j]]));
	distc = norm_2((vertex[face_conn[planeIndex][j]]-vertex[face_conn[planeIndex][i]]));
			
	normdum = Cross_Product((vertex[face_conn[planeIndex][i]] - point), (vertex[face_conn[planeIndex][j]] - point));
	normdum = normdum / norm_2(normdum);
	sign = Dot_Product(normdum,getPlaneNormal(planeIndex));

	sum_angle = sum_angle + sign*acos((pow(dista,2) + pow(distb,2) - pow(distc,2))/(2.0*dista*distb));                                                       
} 

if( abs(sum_angle) >= pi){result = true;}
return(result);
};

void ElementalHexaHedron::tessellate(int ){};
void ElementalHexaHedron::triangulate(int ){};
void ElementalHexaHedron::tetrahedralize(int ){};

/*! Given a set of planes describing HexaHedron faces,
*  check vertices coplanarity and reset their position eventually.*/ 
void ElementalHexaHedron::resetVertices(){

	std::array<darray3E, 3> matrix;
	dvecarr3E result(vertex.size(),darray3E{0,0,0});
	darray3E rhs;
	bool check=  true;	
	for(int i=0; i<8; i++){ //loop for intersection of epsval-translated planes to get every vertices of cornice
		for(int j=0; j<3; j++){ // loop on matrix column
			for(int k=0; k<3; k++){ // compiling a single column on element planes contributing to intersection
				matrix[k][j] = - planes[vert_conn[i][k]][j];
			}//next k element of the column  
			
			rhs[j] = planes[vert_conn[i][j]][3]; //setting translation
		}//next j column  
		
		double determinant = det(matrix);
		if(determinant == 0) {
			std::cout<< "failed to recalculate patch vertices in PATCH::resetVertices method. Vertices are untouched"<<std::endl;
			check = false;
		}else{
			darray3E dummy;
			Cramer(matrix, rhs, dummy);
			result[i] = dummy;
			check = check && true;
		}
	} // next i for  loop for intersection of vertices
	
	if(check){
		vertex = result;
	}
};

/*! Set Face-Vertex connectivity */
void ElementalHexaHedron::setFaceConnectivity(){
	freeContainer(face_conn);
	face_conn.resize(6, ivector1D(4,-1));	
face_conn[0][0] = 0;  face_conn[0][1] = 3;   face_conn[0][2] = 7; face_conn[0][3] = 4;  //west side 
face_conn[1][0] = 0;  face_conn[1][1] = 4;   face_conn[1][2] = 5; face_conn[1][3] = 1;  //south side 
face_conn[2][0] = 0;  face_conn[2][1] = 1;   face_conn[2][2] = 2; face_conn[2][3] = 3;  //front side             
face_conn[3][0] = 1;  face_conn[3][1] = 5;   face_conn[3][2] = 6; face_conn[3][3] = 2;  //east side             
face_conn[4][0] = 3;  face_conn[4][1] = 2;   face_conn[4][2] = 6; face_conn[4][3] = 7; //north side 
face_conn[5][0] = 4;  face_conn[5][1] = 7;   face_conn[5][2] = 6; face_conn[5][3] = 5; // back side 
}; 

/*! Set Vertex-Face connectivity */
void ElementalHexaHedron::setVertConnectivity(){
	freeContainer(vert_conn);
	vert_conn.resize(8, ivector1D(3,-1));
vert_conn[0][0] = 0; vert_conn[0][1] = 1; vert_conn[0][2] =2; //W-S-F plane intersection  
vert_conn[1][0] = 1; vert_conn[1][1] = 2; vert_conn[1][2] =3; //S-F-E plane intersection  
vert_conn[2][0] = 2; vert_conn[2][1] = 3; vert_conn[2][2] =4; //F-E-N plane intersection
vert_conn[3][0] = 0; vert_conn[3][1] = 2; vert_conn[3][2] =4; //W-F-N plane intersection
vert_conn[4][0] = 0; vert_conn[4][1] = 1; vert_conn[4][2] =5; //W-S-B plane intersection
vert_conn[5][0] = 1; vert_conn[5][1] = 3; vert_conn[5][2] =5; //S-E-B plane intersection
vert_conn[6][0] = 3; vert_conn[6][1] = 4; vert_conn[6][2] =5; //E-N-B plane intersection
vert_conn[7][0] = 0; vert_conn[7][1] = 4; vert_conn[7][2] =5; //W-N-B plane intersection  
}; 

/*!Set planes on Hexahedron faces once a certain vertex list is provided. It can move vertices if one or more of then
* does not satisfy coplanarity constraints on each face.
*/ 
void ElementalHexaHedron::setPlanes(){
if(vertex.size() != 8) {std::cout<<" Less or More than 8 vertices defined in your vertex List. Errors can occurr."<<endl;}

freeContainer(planes);
planes.resize(6, darray4E{0,0,0,0});	

//workspace
double toll(1.e-8);
dvecarr3E normdum(6);
darray3E vv, vaux;
double sum_up;

// define normals of patch sides - operator - between vector defined in "Operators.hpp"
		// define coefficient of ax+by+cz +d = 0 plane
		for(int i =0; i<6; i++){ 
			vv   = vertex[face_conn[i][2]] - vertex[face_conn[i][0]];
			vaux = vertex[face_conn[i][3]] - vertex[face_conn[i][1]];
			normdum[i] = Cross_Product(vv, vaux);
			normdum[i] = normdum[i] / norm_2(normdum[i]);

			// saving implicit definition of planes[i]
			sum_up = 0.0 ;
			for(int j=0; j<3; j++){
				planes[i][j] = normdum[i][j]; 
				for(int k=0; k<2; k++){
					sum_up += normdum[i][j] * vertex[face_conn[i][k]][j];
					
				}
			}
			
			planes[i][3] = - sum_up/2.0; //average regression
		}
		resetVertices(); // just in case complanarity is not satisfied, thus find the new set of vertices
};   
