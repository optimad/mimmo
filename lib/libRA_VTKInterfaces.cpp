#include "libRA_VTKInterfaces.hpp"

//#####################################################################################################
// CLASS VTK_SHAPE IMPLEMENTATION
/*
 *	\date			31/12/2015
 *	\authors		Gallizio Federico
 * 	\authors 		Allen Woody
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief VTK Class Interface for I/O managing of Unstructured Triangulated Meshes  
 *
 *	VTK_SHAPE is an I/O manager class for reading and writing Unstructured Meshes and data attached to them, 
 *	in VTK native .vtu format.
 */

/*! Default Constructor.Retains all members to null pointers */
VTK_SHAPE::VTK_SHAPE(): VTK_UnstructuredGrid<VTK_SHAPE>(){
	tri=NULL;
	ext_vertex=NULL;
	scalarField=NULL;
	vectField=NULL;
};

/*!Custom Constructor, for Output usage of the class purposes. 
 * Set link to Class_SurfTri data, set filename, directory path and writing codex. 
 * All other class members are initialized to NULL. Requires reference to:
 * \param[in] triP Pointer to Class_SurfTri data structure;
 * \param[in] dir_ string for directory path of vtu file 	
 * \param[in] name_ string for filename where vtu data will be managed	
 * \param[in] cod_  writing codex of the file (binary-1 or ascii-0);
 */
VTK_SHAPE::VTK_SHAPE(Class_SurfTri * triP, std::string dir_, std::string name_, std::string cod_) :
VTK_UnstructuredGrid<VTK_SHAPE>(dir_, name_)
{
	SetCodex(cod_);
	SetDimensions(triP->nSimplex, triP->nVertex, 3*triP->nSimplex);
	tri = triP;
	ext_vertex=NULL;
	scalarField=NULL;
	vectField=NULL;
};

/*!Custom Constructor, for Input usage of the class purposes. 
 * Set link to Class_SurfTri data, set filename, directory path and writing codex. 
 * All other class members are initialized to NULL. Requires reference to:
 * \param[in] triP Pointer to Class_SurfTri data structure;
 * \param[in] filename_ string for filename absolute path, where vtu data will be managed	
 */
VTK_SHAPE::VTK_SHAPE(Class_SurfTri * triP, std::string filename_) :VTK_UnstructuredGrid<VTK_SHAPE>(){
	tri = triP;
	ext_vertex=NULL;
	scalarField=NULL;
	vectField=NULL;
	
	std::string key1 = "/\\";
	std::string key2 = ".";
	std::size_t cut1 =filename_.find_last_of(key1);
	std::size_t cut2 =filename_.rfind(key2);
	
	std::string dir, name;
	dir  = filename_.substr(0, cut1);
	name = filename_.substr(cut1+1, cut2-cut1-1);
	
	SetNames(dir, name);
};

/*! Default Destructor of the class.*/ 
VTK_SHAPE::~VTK_SHAPE(){
	tri=NULL;
	ext_vertex=NULL;
	scalarField=NULL;
	vectField=NULL;
	freeContainer(simplTypes);
};

/*! Link Class_SurfTri data structure to current class. Requires
 * \param[in] triP Pointer to Class_SurfTri data structure;
 */
void VTK_SHAPE::linkTriMesh(Class_SurfTri * triP){
	tri = triP;
};

/*! Link current class to an external list of Vertices. Requires
 *  \param[in] evert Reference to dvecarr3E vertex structure;
 */
void VTK_SHAPE::linkExternalVertexList(dvecarr3E & evert){
	ext_vertex = &evert;
};

/*! Link current class to scalar field data. Requires
 *  \param[in] field Reference to dvector1D scalar field;
 */
void VTK_SHAPE::linkScalarField(dvector1D & field){
	scalarField = &field;
};

/*! Link current class to a vector field data. Requires
 *  \param[in] vfield Reference to dvecarr3E vector field;
 */
void VTK_SHAPE::linkVectorField(dvecarr3E & vfield){
	vectField = &vfield;
};

/*! Unlink current class to Class_SurfTri mesh data structure */
void VTK_SHAPE::unlinkTriMesh(){
	tri = NULL;
};

/*! Unlink current class to a given external vertex list*/
void VTK_SHAPE::unlinkExternalVertexList(){
	ext_vertex = NULL;
};

/*! Unlink current class to a given scalar field*/
void VTK_SHAPE::unlinkScalarField(){
	scalarField = NULL;
};
/*! Unlink current class to a given vector field*/
void VTK_SHAPE::unlinkVectorField(){
	vectField = NULL;
};

/*! CRTP method to write VTU files. Not intended for public usage purpose. 
 *  Custom implementation of VTK_UnstructuredGrid base class in BITPIT_BASE library. 
 *  Check BITPIT_BASE documentation for more information  */
void VTK_SHAPE::Flush(  fstream &str, string codex_, string name  )
{
	int n;
	
	if(tri == NULL ){std::cout<<"not linked Data Structure"<<endl; return;}
	
	string indent("         ") ;
	
	if( codex_ == "ascii"){
		if( name == "Points" ){
			for( n=0; n<tri->nVertex; n++) {
				flush_ascii( str, indent ) ;
				if(ext_vertex == NULL){flush_ascii( str, 3, tri->Vertex[n]) ;}
				else					{flush_ascii( str, 3, (*ext_vertex)[n]) ;}
				str<<endl;
			};
			
		};
	
		if( name == "connectivity" ){
			for( n=0; n<tri->nSimplex; n++) {
				flush_ascii( str, indent ) ;
				flush_ascii( str, 3, tri->Simplex[n]) ;
				str<<endl;
			};
		};
	
		if( name == "types" ){
		int type_(5) ;
			for( n=0; n<tri->nSimplex; n++) {
				flush_ascii( str, indent ) ;
				flush_ascii( str, type_  ) ;
				str<<endl;
			};
		};
		
		if( name == "offsets" ){
		int off_(0) ;
			for( n=0; n<tri->nSimplex; n++) {
				off_ += NumberOfElements( 5 ) ;
				flush_ascii( str, indent ) ;
				flush_ascii( str, off_  ) ;
				str<<endl;
			};
		};
			
	
		if (name == "scalarField" && scalarField != NULL){
		VTK::Field_C * myfield;
		bool check= GetFieldByName(name, myfield);
		std::string datatype = myfield->GetLocation();
				
			if(datatype=="Point" && scalarField->size() == tri->nVertex){
					
				for( n=0; n<scalarField->size(); n++) {
					flush_ascii( str, indent ) ;
					flush_ascii( str,(*scalarField)[n]);
					str<<endl;
				};
			};
			
			if(datatype=="Cell" && scalarField->size() == tri->nSimplex){
					
				for( n=0; n<scalarField->size(); n++) {
				flush_ascii( str, indent ) ;
				flush_ascii( str,(*scalarField)[n]);
				str<<endl;
				};
			};
		};
			
		if (name == "vectorField" && vectField != NULL){
		VTK::Field_C * myfield;
		bool check= GetFieldByName(name, myfield);
		std::string datatype = myfield->GetLocation();
				
			if(datatype=="Point" && vectField->size() == tri->nVertex){
					
				for( n=0; n<scalarField->size(); n++) {
					flush_ascii( str, indent ) ;
					flush_ascii( str, 3, (*vectField)[n]);
					str<<endl;
				};
			};
			
			if(datatype=="Cell" && vectField->size() == tri->nSimplex){
					
				for( n=0; n<vectField->size(); n++) {
					flush_ascii( str, indent ) ;
					flush_ascii( str, 3, (*vectField)[n]);
					str<<endl;
				};
			};
		};
			
	}else{		
		if( name == "Points" ){
			for( n=0; n<tri->nVertex; n++) {
				if(ext_vertex == NULL)	flush_binary( str, tri->Vertex[n]);
				else			flush_binary( str, (*ext_vertex)[n]);
			}
		};
			
		if( name == "connectivity" ){
			for( n=0; n<tri->nSimplex; n++) {flush_binary( str, tri->Simplex[n]);}
		};
				
		if( name == "types"){
		int type_(5) ;
			for( n=0; n<tri->nSimplex; n++) flush_binary( str, type_  ) ;
		};
				   
		if( name == "offsets"){
		int off_(0) ;
			for( n=0; n<tri->nSimplex; n++) {
				off_ += NumberOfElements( 5 ) ;
				flush_binary( str, off_  ) ;
			};
		};
				   
		if (name == "scalarField" && scalarField != NULL){
		VTK::Field_C * myfield;
		bool check= GetFieldByName(name, myfield);
		std::string datatype = myfield->GetLocation();
					   
			if(datatype=="Point" && scalarField->size() == tri->nVertex){
						   
				for( n=0; n<scalarField->size(); n++) {
					flush_binary( str, (*scalarField)[n]);
				};
			};
					   
			if(datatype=="Cell" && scalarField->size() == tri->nSimplex){
				for( n=0; n<scalarField->size(); n++) {
					flush_binary( str, (*scalarField)[n]);
				};
			};
		};
				   
		if (name == "vectorField" && vectField != NULL){
		VTK::Field_C * myfield;
		bool check= GetFieldByName(name, myfield);
		std::string datatype = myfield->GetLocation();
					   
			if(datatype=="Point" && vectField->size() == tri->nVertex){
						   
				myfield->SetElements(tri->nVertex);
				myfield->SetOffset(tri->nVertex);
						   
				for( n=0; n<vectField->size(); n++) {
					flush_binary( str, (*vectField)[n]);
				};
			};
			
			if(datatype=="Cell" && vectField->size() == tri->nSimplex){
						   
				myfield->SetElements(tri->nSimplex);
				myfield->SetOffset(tri->nSimplex);
						   
				for( n=0; n<vectField->size(); n++) {
					flush_binary( str, (*vectField)[n]);
				};
			};
		};
				   
	}

return ;
};

/*! CRTP method to read VTU files. Not intended for public usage purpose. 
 *  Custom implementation of VTK_UnstructuredGrid base class in BITPIT_BASE library. 
 *  Check BITPIT_BASE documentation for more information  */
void VTK_SHAPE::Absorb(  fstream &str, string codex_, string name  ){
	
	int n;
	
	if(tri == NULL){std::cout<< "Not Linked Data Structure"<<endl; return;}
	
	if( codex_ == "ascii" ){
		
		if( name == "Points"){
			tri->nVertex = nr_points;
			tri->ResizeVertex();
			
			for( n=0; n<nr_points; n++) {
				for(int jj=0; jj<3; ++jj){
					absorb_ascii( str, tri->Vertex[n][jj]) ;
				};
			};
		};

		if( name == "connectivity"){
			if(simplTypes.size() != nr_cells) {std::cout<<"Not able to reading simplicial types of your file.Exiting"<<std::endl; exit(1);}
			nconnectivity = 0;
			tri->nSimplex = nr_cells; 
			tri->Simplex.resize(tri->nSimplex);
			for(int k=0; k<simplTypes.size(); ++k) {
				nconnectivity += NumberOfElements(simplTypes[k]);
				tri->Simplex[k].resize(NumberOfElements(simplTypes[k]), 0);
			};
			
			VTK::Field_C * myfield;
			const std::string ncname = "connectivity";
			bool check= GetFieldByName(ncname, myfield);
			myfield->SetElements(nconnectivity);
			for( n=0; n<tri->nSimplex; n++) {
				for(int jj=0; jj<tri->Simplex[n].size(); ++jj){
					absorb_ascii( str, tri->Simplex[n][jj]) ;
				}
			};
		};
	
		if( name == "types"){
			if (simplTypes.size() < nr_cells){freeContainer(simplTypes);}
			simplTypes.resize(nr_cells, -1) ;
			for( n=0; n<nr_cells; n++) {
				absorb_ascii( str, simplTypes[n]) ;
			};
		};
		
		if( name == "offsets"){
			int off_ ;
			for( n=0; n<nr_cells; n++) absorb_ascii( str, off_  ) ;
			
		};
			
	}else{
		if( name == "Points"){
			tri->nVertex = nr_points;
			tri->ResizeVertex();
			
			for( n=0; n<nr_points; n++) {
				absorb_binary( str, tri->Vertex[n]) ;
			};
		};
			
		if( name == "connectivity"){
			if(simplTypes.size() != nr_cells) {std::cout<<"Not able to reading simplicial types of your file.Exiting"<<std::endl; exit(1);}
			nconnectivity = 0;
			tri->nSimplex = nr_cells; tri->Simplex.resize(tri->nSimplex);
			
			
			for(int k=0; k<simplTypes.size(); ++k) {
				
				nconnectivity += NumberOfElements(simplTypes[k]);
				tri->Simplex[k].resize(NumberOfElements(simplTypes[k]), 0);
			};
			
			VTK::Field_C * myfield;
			const std::string ncname = "connectivity";
			bool check= GetFieldByName(ncname, myfield);
			myfield->SetElements(nconnectivity);
			
			for( n=0; n<tri->nSimplex; n++) {
				for(int jj=0; jj<tri->Simplex[n].size(); ++jj){
					absorb_binary( str, tri->Simplex[n][jj]) ;
				}
			};
		};
				
		if( name == "types"){
			if (simplTypes.size() < nr_cells){freeContainer(simplTypes);}
			simplTypes.resize(nr_cells, -1) ;
			for( n=0; n<nr_cells; n++) {
				absorb_binary( str, simplTypes[n]) ;
			};
		};
					
		if( name == "offsets"){
			int off_ ;
			for( n=0; n<nr_cells; n++) {
				absorb_binary( str, off_  ) ;
			};
		};
						
	};
						
return ;
};
				   
//#####################################################################################################
// CLASS VTK_LSet IMPLEMENTATION
/*
 *	\date			31/12/2015
 *	\authors		Gallizio Federico
 * 	\authors 		Allen Woody
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief VTK Class Interface for output managing of LevelSet on 3D Cartesian Structured Meshes  
 *
 *	VTK_LSet is an Output manager class for writing Level Set scalar field and Level Set gradient field on  
 *	Cartesian Background Mesh in VTK native .vtr format.
 */

/*! Basic Constructor */
VTK_LSet::VTK_LSet():VTK_RectilinearGrid<VTK_LSet>() {
	mesh3D = NULL;
	scalarField=NULL;
	vectorField=NULL;
};

/*! Custom Constructor. Set 
 * \param[in] dir path to folder for output
 * \param[in] name file name of output, without tag 
 */
VTK_LSet::VTK_LSet(std::string dir, std::string name):VTK_RectilinearGrid<VTK_LSet>(dir, name,"ascii",0,0,0,0,0,0){
	mesh3D = NULL;
	scalarField=NULL;
	vectorField=NULL;
};
/*! Custom Constructor. Set 
 * \param[in] dir path to folder for output
 * \param[in] name file name of output, without tag 
 * \param[in] codex file codification, choose between "ascii" or "appended"
 */
VTK_LSet::VTK_LSet(std::string dir, std::string name, std::string codex):VTK_RectilinearGrid<VTK_LSet>(dir, name,codex,0,0,0,0,0,0){
	mesh3D = NULL;
	scalarField=NULL;
	vectorField=NULL;
};

/*! Destructor of the class */
VTK_LSet::~VTK_LSet(){
	mesh3D = NULL;
	scalarField=NULL;
	vectorField=NULL;
};

/*! Link Class to an external background 3D cartesian mesh.
 * \param[in] mesh reference to Class_UCartMesh3D object	
 */
void VTK_LSet::linkData(Class_UCartMesh3D & mesh){
	mesh3D = &mesh;
	SetDimensions(1, mesh.nx+1, 1, mesh.ny+1, 1, mesh.nz+1);
};

/*! Link Class to an external scalar field of doubles.
 * \param[in] mesh reference to dvector1D object	
 */
void VTK_LSet::linkScalarField(dvector1D & field){
	scalarField = &field;
};

/*! Link Class to an external 3D vector field of doubles.
 * \param[in] mesh reference to dvecarr3E object	
 */
void VTK_LSet::linkVectorField(dvecarr3E & field){
	vectorField = &field;
};

/*! Unlink Class from current background mesh*/
void VTK_LSet::unlinkData(){
	mesh3D=NULL;
};

/*! Unlink Class from current scalar field */
void VTK_LSet::unlinkScalarField(){
	scalarField=NULL;
};

/*! Unlink Class from current 3D vector field*/
void VTK_LSet::unlinkVectorField(){
	vectorField=NULL;
};

/*! CRTP method to write VTR Rectilinear Meshes files. Not intended for public usage purpose. 
 *  Custom implementation of VTK_RectilinearGrid base class in BITPIT_BASE library. 
 *  Check BITPIT_BASE documentation for more information  */
void VTK_LSet::Flush(  fstream &str, string codex_, string name  ){
	
	int n;
	if(mesh3D == NULL ){std::cout<<"not linked Data Structure"<<endl; return;}
	
	string indent("         ") ;
	
	if(codex_=="ascii") {
		if(name=="x_Coord"){
			for (n=0; n<=mesh3D->nx; ++n){
				flush_ascii(str, indent);
				flush_ascii(str, mesh3D->xedge[n]);
				str<<endl;
			};
		};

		if(name=="y_Coord"){
			for (n=0; n<=mesh3D->ny; ++n){
				flush_ascii(str, indent);
				flush_ascii(str, mesh3D->yedge[n]);
				str<<endl;
			};
		};
	
		if(name=="z_Coord"){
			for (n=0; n<=mesh3D->nz; ++n){
				flush_ascii(str, indent);
				flush_ascii(str, mesh3D->zedge[n]);
				str<<endl;
			};
		};
		
		if(name=="scalarField"){
			for (n=0; n<=scalarField->size(); ++n){
				flush_ascii(str, indent);
				flush_ascii(str, (*scalarField)[n]);
				str<<endl;
			};
		};
			
		if(name=="vectorField"){
			for (n=0; n<vectorField->size(); ++n){
				flush_ascii(str, indent);
				flush_ascii(str, 3,(*vectorField)[n]);
				str<<endl;
			};
		};
				
	}else {
		if(name=="x_Coord"){
			for (n=0; n<=mesh3D->nx; ++n){
				flush_binary(str, mesh3D->xedge[n]);
				
			};
		};
				
		if(name=="y_Coord"){
			for (n=0; n<=mesh3D->ny; ++n){
				flush_binary(str, mesh3D->yedge[n]);
			};
		};
					
		if(name=="z_Coord"){
			for (n=0; n<=mesh3D->nz; ++n){
				flush_binary(str, mesh3D->zedge[n]);
			};
		};
						
		if(name=="scalarField"){
			for (n=0; n<scalarField->size(); ++n){
				flush_binary(str, (*scalarField)[n]);
			};
		};
							
		if(name=="vectorField"){
			for (n=0; n<vectorField->size(); ++n){
				flush_binary(str,(*vectorField)[n]);
			};
		};
								
	}
return;
};

//*****************************************************************************************************************************
// VTK_BASICMESH IMPLEMENTATION 
/*
 *	\date			31/12/2015
 *	\authors		Gallizio Federico
 * 	\authors 		Amos Tori
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief Class for VTK Output purpose 
 *
 *	Class for managing output on file of uniform structured grids, based on VTK's formats (*.vtu)
 */

/*! Base Class Constructor */
VTK_BASICMESH::VTK_BASICMESH(){
	points = NULL;
	connectivity=NULL;
};
/*! Custom Constructor. Set up dimension of mesh that need to be written and its output filename and path.
 * \param[in] dir_ output directory path 
 * \param[in] name_ output filename w/out tag
 * \param[in] cod_  output codex: "ascii" or "appended"
 * \param[in] nP_   number of points in your mesh
 * \param[in] nC_   number of cells in your mesh
 * \param[in] nConn connectivity leading dimension, as sequential vector (in general 8*nC_)
 */ 
VTK_BASICMESH::VTK_BASICMESH(std::string dir_, std::string name_, std::string cod_, int nP_, int nC_, int nConn_):
VTK_UnstructuredGrid<VTK_BASICMESH>(dir_, name_) 
{
	SetCodex(cod_);
	SetDimensions(nC_, nP_, nConn_);
	points=NULL;
	connectivity=NULL;
};

/*!Basic Class Destructor */
VTK_BASICMESH::~VTK_BASICMESH(){
	points = NULL;
	connectivity = NULL;
};

/*! Link external points list and connectivity matrix to class. Data are not copied, nor stored.
 * \param[in] pp list of external points
 * \param[in] conn connectivity matrix related to pp
 */
void VTK_BASICMESH::linkData(dvecarr3E & pp, ivector2D & conn ){
	
	int typeC = 0;
	if(conn.size() >0 ){typeC = conn[0].size();}
	
	points = &pp;
	connectivity = &conn;
	//check data connection to vtk instatiated info and repair it.
	if(pp.size() != nr_points || conn.size() != nr_cells){
		SetDimensions( conn.size(), pp.size(), typeC*conn.size());
	}
};
/*! Unlink Class from points or connectivity external data*/
void VTK_BASICMESH::unlinkData(){
	points = NULL;
	connectivity = NULL;
};

/*! CRTP function to write data on file. Not suitable for a direct call. 
 * Please refer to documentation of VTK_UnstructuredGrid class in Bitpit_BASE 
 */
void VTK_BASICMESH::Flush(  fstream &str, string codex_, string name  ){
	
	int n;
	if(points == NULL || connectivity == NULL ){std::cout<<"not linked Data Structure"<<endl; return;}
	int type_ = 12;
	int lconn = 8;
	
	if((*connectivity)[0].size() == 4){
		type_ = 9; 
		lconn=4;
	}
	
	string indent("         ") ;
	
	if( codex_ == "ascii"){
		
		if( name == "Points" ){
			for( n=0; n<nr_points; n++) {
				flush_ascii( str, indent ) ;
				flush_ascii( str, 3, (*points)[n]);
				str << endl ;
			};
		};

		if( name == "connectivity" ){
			for( n=0; n<connectivity->size(); n++) {
				flush_ascii( str, indent ) ;
				flush_ascii( str, lconn, (*connectivity)[n]) ;
				str << endl ;
			};
		};
	
		if( name == "types" ){
			for( n=0; n<connectivity->size(); n++) {
				flush_ascii( str, indent ) ;
				flush_ascii( str, type_  ) ;
				str << endl ;
			};
		};
		
		if( name == "offsets" ){
			int off_(0) ;
			for( n=0; n<connectivity->size(); n++) {
				off_ += NumberOfElements( type_ ) ;
				flush_ascii( str, indent ) ;
				flush_ascii( str, off_  ) ;
				str << endl ;
			};
		};
			
	}else{
		if( name == "Points" ){
			for( n=0; n<nr_points; n++) {
				flush_binary( str, (*points)[n]);
			};
		};
			
		if( name == "connectivity" ){
			for( n=0; n<connectivity->size(); n++) {
				flush_binary( str, (*connectivity)[n]) ;
			};
		};
				
		if( name == "types" ){
			for( n=0; n<connectivity->size(); n++) {
				flush_binary( str, type_  ) ;
			};
		};
					
		if( name == "offsets" ){
			int off_(0) ;
			for( n=0; n<connectivity->size(); n++) {
				off_ += NumberOfElements( type_ ) ;
				flush_binary( str, off_  ) ;
			};
		};
						
	}
						
return ;
}; //CRTP

//*****************************************************************************************************************************
// VTK_BASICCLOUD IMPLEMENTATION 
/*
 *	\date			31/12/2015
 *	\authors		Gallizio Federico
 * 	\authors 		Amos Tori
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief Class for VTK Output purpose 
 *
 *	Class for managing output on file of point clouds, based on VTK's formats (*.vtu)
 */

/*! Basic class Constructor */
VTK_BASICCLOUD::VTK_BASICCLOUD(){
	points = NULL;
};

/*! Custom Constructor. Set up dimension of cloud points that need to be written and its output filename and path.
 * \param[in] dir_ output directory path 
 * \param[in] name_ output filename w/out tag
 * \param[in] cod_  output codex: "ascii" or "appended"
 * \param[in] nP_   number of points in your mesh
 */ 
VTK_BASICCLOUD::VTK_BASICCLOUD(std::string dir_, std::string name_, std::string cod_, int nP_):
VTK_UnstructuredGrid<VTK_BASICCLOUD>(dir_, name_)  
{
	SetCodex(cod_);
	SetDimensions(nP_, nP_, nP_);
	points=NULL;
};

/*! Basic Class Destructor */
VTK_BASICCLOUD::~VTK_BASICCLOUD(){
	points=NULL;
};

/*! Link external point lists to class
 * \param[in] pp external point cloud
 */
void VTK_BASICCLOUD::linkData(dvecarr3E & pp){
	
	points = &pp;
	//check data connection to vtk instatiated info and repair it.
	if(pp.size() != nr_points){
		SetDimensions( pp.size(), pp.size(), pp.size());
	}
};

/*! Unlink preexistent data from your class */
void VTK_BASICCLOUD::unlinkData(){
	points = NULL;
};

/*! CRTP function to write data on file. Not suitable for a direct call. 
 * Please refer to documentation of VTK_UnstructuredGrid class in Bitpit_BASE 
 */
void VTK_BASICCLOUD::Flush(  fstream &str, string codex_, string name  ){
	
	int n;
	if(points == NULL ){std::cout<<"not linked Data Structure"<<endl; return;}
	
	string indent("         ") ;
	
	if( codex_ == "ascii"){
		
		if( name == "Points" ){
			for( n=0; n<points->size(); n++) {
				flush_ascii( str, indent ) ;
				flush_ascii( str, 3, (*points)[n]);
				str << endl ;
			};
		};

		if( name == "connectivity" ){
			for( n=0; n<points->size(); n++) {
				flush_ascii( str, indent ) ;
				flush_ascii( str, n) ;
				str << endl ;
			};
		};
	
		if( name == "types" ){
			int type_(1) ;
			for( n=0; n<points->size(); n++) {
				flush_ascii( str, indent ) ;
				flush_ascii( str, type_  ) ;
				str << endl ;
			};
		};
		
		if( name == "offsets" ){
			int off_(0) ;
			for( n=0; n<points->size(); n++) {
				off_ += NumberOfElements( 1 ) ;
				flush_ascii( str, indent ) ;
				flush_ascii( str, off_  ) ;
				str << endl ;
			};
		};
			
	}else{
		if( name == "Points" ){
			for( n=0; n<points->size(); n++) {
				flush_binary( str, (*points)[n]);
			};
		};
			
		if( name == "connectivity" ){
			for( n=0; n<points->size(); n++) {
				flush_binary( str, n) ;
			};
		};
				
		if( name == "types" ){
			int type_(1) ;
			for( n=0; n<points->size(); n++) {
				flush_binary( str, type_  ) ;
			};
		};
					
		if( name == "offsets" ){
			int off_(0) ;
			for( n=0; n<points->size(); n++) {
				off_ += NumberOfElements( 1 ) ;
				flush_binary( str, off_  ) ;
			};
		};
	}
return ;
}; //CRTP


//*****************************************************************************************************************************
// VTK_3DCURVE IMPLEMENTATION 
/*
 *	\date			30/1/2016
 *	\authors		Gallizio Federico
 *	\authors 		Suzanne Vega
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief VTK Class Interface for Output managing of Boundary Contours of Unstructured Surface Meshes  
 *
 *	VTK_3DCURVE is an Output manager class for writing Cloud Points related to Unstructured Surface Meshes
 *	boundary contours, handled by a Class_SurfTri object,in VTK native .vtu format.
 */

/*! Basic Constructor */
VTK_3DCURVE::VTK_3DCURVE():VTK_UnstructuredGrid<VTK_3DCURVE>(){
	curve = NULL;
};

/*! Custom Constructor. Set Class_SurfTri Curve source, and output path, filename and codex.
 * \param[in] triP 3D curve stored as Class_Surftri object
 * \param[in] dir_ output directory
 * \param[in] name_ name of output filename
 * \param[in] cod_  output codex: "ascii" or "appended"
 */
VTK_3DCURVE::VTK_3DCURVE( Class_SurfTri & triP, std::string dir_, std::string name_, std::string cod_):
VTK_UnstructuredGrid<VTK_3DCURVE>(dir_, name_)
{
	SetCodex(cod_);
	SetDimensions(triP.nSimplex, triP.nVertex, 2*triP.nSimplex);
	curve=&triP;
};
/*! Custom Constructor. Set Class_SurfTri Curve source, and output filename as absolute path. No codex is specified.
 * \param[in] triP 3D curve stored as Class_Surftri object
 * \param[in] filename_ output filename absolute path
 */
VTK_3DCURVE::VTK_3DCURVE(Class_SurfTri & triP, std::string filename_) :
VTK_UnstructuredGrid<VTK_3DCURVE>(){
	curve = &triP;
	
	std::string key1 = "/\\";
	std::string key2 = ".";
	std::size_t cut1 =filename_.find_last_of(key1);
	std::size_t cut2 =filename_.rfind(key2);
	
	std::string dir, name;
	dir  = filename_.substr(0, cut1);
	name = filename_.substr(cut1+1, cut2-cut1-1);
	
	SetNames(dir, name);
};

/*! Basic Destructor */
VTK_3DCURVE::~VTK_3DCURVE(){
	curve=NULL;
};

/*! Link source of a 3D curve to class 
 * \param[in] triP Class_SurfTri object containing the 3D curve
 */
void VTK_3DCURVE::linkData(Class_SurfTri & triP){
	curve = &triP;
};

/*! Unlink preexistent data sources */
void VTK_3DCURVE::unlinkData(){
	curve=NULL;
};

/*! CRTP function to write data on file. Not suitable for a direct call. 
 * Please refer to documentation of VTK_UnstructuredGrid class in Bitpit_BASE 
 */
void VTK_3DCURVE::Flush(  fstream &str, string codex_, string name  )
{
	int n;
	if(curve == NULL ){std::cout<<"not linked Data Structure"<<endl; return;}
	
	string indent("         ") ;
	
	if( codex_ == "ascii"){
		
		if( name == "Points" ){
			for( n=0; n<curve->nVertex; n++) {
				flush_ascii( str, indent ) ;
				flush_ascii( str, 3, curve->Vertex[n]) ;
				str<<endl;
			};
		};

		if( name == "connectivity" ){
			for( n=0; n<curve->nSimplex; n++) {
				flush_ascii( str, indent ) ;
				flush_ascii( str, 2, curve->Simplex[n]) ;
				str<<endl;
			};
		};
	
		if( name == "types" ){
			int type_(3) ;
			for( n=0; n<curve->nSimplex; n++) {
				flush_ascii( str, indent ) ;
				flush_ascii( str, type_  ) ;
				str<<endl;
			};
		};
		
		if( name == "offsets" ){
			int off_(0) ;
			for( n=0; n<curve->nSimplex; n++) {
				off_ += NumberOfElements( 3 ) ;
				flush_ascii( str, indent ) ;
				flush_ascii( str, off_  ) ;
				str<<endl;
			};
		};
			
	}else{
		if( name == "Points" ){
			for( n=0; n<curve->nVertex; n++) {
				flush_binary( str, curve->Vertex[n]) ;
			}
		};
			
		if( name == "connectivity" ){
			for( n=0; n<curve->nSimplex; n++) {flush_binary( str, curve->Simplex[n]);}
		};
				
		if( name == "types"){
			int type_(3) ;
			for( n=0; n<curve->nSimplex; n++) flush_binary( str, type_  ) ;
		};
				   
		if( name == "offsets"){
			int off_(0) ;
			for( n=0; n<curve->nSimplex; n++) {
				off_ += NumberOfElements( 3 ) ;
				flush_binary( str, off_  ) ;
			};
		};
				   
	}
return ;
};



