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

#include "MiMMO_VTKInterfaces.hpp"

using namespace std;
using namespace mimmo;

/*
 *	\date			31/12/2015
 *	\authors		Gallizio Federico
 * 	\authors 		Amos Tori
 *	\authors		Arpa Rocco
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief Class for VTK Output of Structured Grids 
 *
 *	Class for managing output on file of uniform structured grids, based on VTK's formats (*.vtu)
 */

/*! Base Class Constructor */
VTK_BASICMESH::VTK_BASICMESH(){
	m_points = NULL;
	m_connectivity=NULL;
};
/*! Custom Constructor. Set up dimension of mesh that need to be written and its output filename and path.
 * \param[in] dir_ output directory path 
 * \param[in] name_ output filename w/out tag
 * \param[in] cod_  output codex: "ascii" or "appended"
 * \param[in] nP_   number of points in your mesh
 * \param[in] nC_   number of cells in your mesh
 * \param[in] nConn_ connectivity leading dimension, as sequential vector (in general 8*nC_)
 */ 
VTK_BASICMESH::VTK_BASICMESH(std::string dir_, std::string name_,bitpit::VTKFormat cod_, int nP_, int nC_, int nConn_):
VTKUnstructuredGrid(dir_, name_)
{
	setCodex(cod_);
	setDimensions(nC_, nP_, nConn_);
	m_points=NULL;
	m_connectivity=NULL;
};

/*!Basic Class Destructor */
VTK_BASICMESH::~VTK_BASICMESH(){
	m_points = NULL;
	m_connectivity = NULL;
};

/*! Link external points list and connectivity matrix to class. Data are not copied, nor stored.
 * \param[in] pp list of external points
 * \param[in] conn connectivity matrix related to pp
 */
void VTK_BASICMESH::linkData(dvecarr3E & pp, ivector2D & conn ){
	
	int typeC = 0;
	if(conn.size() >0 ){typeC = conn[0].size();}
	
	m_points = &pp;
	m_connectivity = &conn;
	//check data connection to vtk instatiated info and repair it.
	if(pp.size() != nr_points || conn.size() != nr_cells){
		setDimensions( conn.size(), pp.size(), typeC*conn.size());
	}
};
/*! Unlink Class from points or connectivity external data*/
void VTK_BASICMESH::unlinkData(){
	m_points = NULL;
	m_connectivity = NULL;
};

/*! CRTP function to write data on file. Not suitable for a direct call. 
 * Please refer to documentation of VTK_UnstructuredGrid class in Bitpit_BASE 
 */
void VTK_BASICMESH::flushData(  fstream &str, bitpit::VTKFormat codex_, string name  ){
	
	int n;
	if(m_points == NULL || m_connectivity == NULL ){std::cout<<"not linked Data Structure"<<endl; return;}
	bitpit::VTKElementType type_ = bitpit::VTKElementType::HEXAHEDRON;
	int lconn = 8;
	
	if((*m_connectivity)[0].size() == 4){
		type_ = bitpit::VTKElementType::QUAD;
		lconn=4;
	}
	
	string indent("         ") ;
	
	if( codex_ == bitpit::VTKFormat::ASCII){
		
		if( name == "Points" ){
			for( n=0; n<nr_points; n++) {
				bitpit::genericIO::flushASCII( str, indent ) ;
				bitpit::genericIO::flushASCII( str, 3, (*m_points)[n]);
				str << endl ;
			};
		};
		
		if( name == "connectivity" ){
			for( n=0; n<m_connectivity->size(); n++) {
				bitpit::genericIO::flushASCII( str, indent ) ;
				bitpit::genericIO::flushASCII( str, lconn, (*m_connectivity)[n]) ;
				str << endl ;
			};
		};
		
		if( name == "types" ){
			for( n=0; n<m_connectivity->size(); n++) {
				bitpit::genericIO::flushASCII( str, indent ) ;
				bitpit::genericIO::flushASCII( str, static_cast<int>(type_)  ) ;
				str << endl ;
			};
		};
		
		if( name == "offsets" ){
			int off_(0) ;
			for( n=0; n<m_connectivity->size(); n++) {
				off_ += bitpit::vtk::getNNodeInElement( type_ ) ;
				bitpit::genericIO::flushASCII( str, indent ) ;
				bitpit::genericIO::flushASCII( str, off_  ) ;
				str << endl ;
			};
		};
		
	}else{
		if( name == "Points" ){
			for( n=0; n<nr_points; n++) {
				bitpit::genericIO::flushBINARY( str, (*m_points)[n]);
			};
		};
		
		if( name == "connectivity" ){
			for( n=0; n<m_connectivity->size(); n++) {
				bitpit::genericIO::flushBINARY( str, (*m_connectivity)[n]) ;
			};
		};
		
		if( name == "types" ){
			for( n=0; n<m_connectivity->size(); n++) {
				bitpit::genericIO::flushBINARY( str, static_cast<int>(type_)  ) ;
			};
		};
		
		if( name == "offsets" ){
			int off_(0) ;
			for( n=0; n<m_connectivity->size(); n++) {
				off_ += bitpit::vtk::getNNodeInElement( type_ ) ;
				bitpit::genericIO::flushBINARY( str, off_  ) ;
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
 *	\brief Class for VTK Output of Point Clouds
 *
 *	Class for managing output on file of point clouds, based on VTK's formats (*.vtu)
 */

/*! Basic class Constructor */
VTK_BASICCLOUD::VTK_BASICCLOUD(){
	m_points = NULL;
};

/*! Custom Constructor. Set up dimension of cloud points that need to be written and its output filename and path.
 * \param[in] dir_ output directory path 
 * \param[in] name_ output filename w/out tag
 * \param[in] cod_  output codex: "ascii" or "appended"
 * \param[in] nP_   number of points in your mesh
 */ 
VTK_BASICCLOUD::VTK_BASICCLOUD(std::string dir_, std::string name_, bitpit::VTKFormat cod_, int nP_):
VTKUnstructuredGrid(dir_, name_)
{
	setCodex(cod_);
	setDimensions(nP_, nP_, nP_);
	m_points=NULL;
};

/*! Basic Class Destructor */
VTK_BASICCLOUD::~VTK_BASICCLOUD(){
	m_points=NULL;
};

/*! Link external point lists to class
 * \param[in] pp external point cloud
 */
void VTK_BASICCLOUD::linkData(dvecarr3E & pp){
	
	m_points = &pp;
	//check data connection to vtk instatiated info and repair it.
	if(pp.size() != nr_points){
		setDimensions( pp.size(), pp.size(), pp.size());
	}
};

/*! Unlink preexistent data from your class */
void VTK_BASICCLOUD::unlinkData(){
	m_points = NULL;
};

/*! CRTP function to write data on file. Not suitable for a direct call. 
 * Please refer to documentation of VTK_UnstructuredGrid class in Bitpit_BASE 
 */
void VTK_BASICCLOUD::flushData(  fstream &str,bitpit::VTKFormat codex_, string name  ){
	
	int n;
	if(m_points == NULL ){std::cout<<"not linked Data Structure"<<endl; return;}
	bitpit::VTKElementType type_ = bitpit::VTKElementType::VERTEX;
	
	string indent("         ") ;
	
	if( codex_ == bitpit::VTKFormat::ASCII){
		
		if( name == "Points" ){
			for( n=0; n<m_points->size(); n++) {
				bitpit::genericIO::flushASCII( str, indent ) ;
				bitpit::genericIO::flushASCII( str, 3, (*m_points)[n]);
				str << endl ;
			};
		};
		
		if( name == "connectivity" ){
			for( n=0; n<m_points->size(); n++) {
				bitpit::genericIO::flushASCII( str, indent ) ;
				bitpit::genericIO::flushASCII( str, n) ;
				str << endl ;
			};
		};
		
		if( name == "types" ){

			for( n=0; n<m_points->size(); n++) {
				bitpit::genericIO::flushASCII( str, indent ) ;
				bitpit::genericIO::flushASCII( str, static_cast<int>(type_)  ) ;
				str << endl ;
			};
		};
		
		if( name == "offsets" ){
			int off_(0) ;
			for( n=0; n<m_points->size(); n++) {
				off_ += bitpit::vtk::getNNodeInElement( type_ ) ;
				bitpit::genericIO::flushASCII( str, indent ) ;
				bitpit::genericIO::flushASCII( str, off_  ) ;
				str << endl ;
			};
		};
		
	}else{
		if( name == "Points" ){
			for( n=0; n<m_points->size(); n++) {
				bitpit::genericIO::flushBINARY( str, (*m_points)[n]);
			};
		};
		
		if( name == "connectivity" ){
			for( n=0; n<m_points->size(); n++) {
				bitpit::genericIO::flushBINARY( str, n) ;
			};
		};
		
		if( name == "types" ){

			for( n=0; n<m_points->size(); n++) {
				bitpit::genericIO::flushBINARY( str, static_cast<int>(type_)  ) ;
			};
		};
		
		if( name == "offsets" ){
			int off_(0) ;
			for( n=0; n<m_points->size(); n++) {
				off_ += bitpit::vtk::getNNodeInElement( type_ ) ;
				bitpit::genericIO::flushBINARY( str, off_  ) ;
			};
		};
	}
	return ;
}; //CRTP

/*! Default Constructor.Retains all members to null pointers */
VTK_DATAMESH::VTK_DATAMESH(){
	m_points=NULL;
	m_connectivity=NULL;
	m_element = bitpit::VTKElementType::UNDEFINED;
};

/*!Custom Constructor, for Output usage of the class purposes. 
 * Set link to mesh data, set filename, directory path and writing codex. 
 * Requires reference to:
 * \param[in] points pointer to list of mesh nodes coordinates;
 * \param[in] conn   pointer to mesh local connectivity;
 * \param[in] dir_ 	 string for directory path of vtu file 	
 * \param[in] name_  string for filename where vtu data will be managed	
 * \param[in] cod_   writing codex of the file (binary-1 or ascii-0);
 */
VTK_DATAMESH::VTK_DATAMESH(dvecarr3E * points, ivector2D * conn, std::string dir_, std::string name_, bitpit::VTKFormat cod_) :
VTKUnstructuredGrid(dir_, name_)
{
	if(!linkMeshData(points,conn)) return;	
	setCodex(cod_);
};

/*!Custom Constructor, for Input usage of the class purposes. 
 * Set link to mesh data, set filename, directory path and writing codex. 
 * Requires reference to:
 * \param[in] points pointer to list of mesh nodes coordinates;
 * \param[in] conn   pointer to mesh local connectivity;
 * \param[in] filename_ string for filename absolute path, where vtu data will be managed	
 */
VTK_DATAMESH::VTK_DATAMESH(dvecarr3E * points, ivector2D * conn, std::string filename_){
	
	
	if(!linkMeshData(points,conn)) return;
	
	std::string key1 = "/\\";
	std::string key2 = ".";
	std::size_t cut1 =filename_.find_last_of(key1);
	std::size_t cut2 =filename_.rfind(key2);
	
	std::string dir, name;
	dir  = filename_.substr(0, cut1);
	name = filename_.substr(cut1+1, cut2-cut1-1);
	
	setNames(dir, name);
};

/*! Default Destructor of the class.*/ 
VTK_DATAMESH::~VTK_DATAMESH(){
	m_points = NULL;
	m_connectivity = NULL;
};

/*! Link mesh data structures to current class. Requires
 * \param[in] points pointer to list of mesh nodes coordinates;
 * \param[in] conn   pointer to mesh local connectivity;
 * \return boolean, true for successfull linking
 */
bool VTK_DATAMESH::linkMeshData(dvecarr3E * points, ivector2D * conn){
	if(conn == NULL || points == NULL) return false;
	if((int)conn->size() < 1)	return false;
	
	int size = (*conn)[0].size();
	m_element = static_cast< bitpit::VTKElementType > (size);
	setDimensions((uint64_t)conn->size(), (uint64_t)points->size(), m_element);
	m_points = points;
	m_connectivity = conn;
	
	return true;
};

void  VTK_DATAMESH::unlinkMeshData(){
	m_points = NULL;
	m_connectivity = NULL;
	m_element = bitpit::VTKElementType::UNDEFINED;
}

/*! Adds a scalar field data.
 * \param[in] name	name of the field;
 * \param[in] field Reference to scalar field of doubles;
 * \param[in] loc 	choose bitpit::VTKLocation for your data on mesh(on nodes/cells)
 * \return	boolean true if field is added, false if not(try change name of your field)
 */
bool VTK_DATAMESH::addScalarField(std::string & name, dvector1D & field, bitpit::VTKLocation loc){
	std::unordered_map< std::string, dvecarr3E * >::iterator it = m_vector.find(name);
	if(it != m_vector.end()) return false;
	bool check = m_scalar.insert({name,&field}).second;
	addData(name, bitpit::VTKFieldType::SCALAR, loc, bitpit::VTKDataType::Float64);
	return check;
};

/*! Adds a vector field data.
 * \param[in] name	name of the field;
 * \param[in] field Reference to vector field of doubles;
 * \param[in] loc 	choose bitpit::VTKLocation for your data on mesh(on nodes/cells)
 * \return	boolean true if field is added, false if not(try change name of your field) 
 */
bool VTK_DATAMESH::addVectorField(std::string & name, dvecarr3E & field, bitpit::VTKLocation loc){
	std::unordered_map< std::string, dvector1D * >::iterator it = m_scalar.find(name);
	if(it != m_scalar.end()) return false;
	bool check = m_vector.insert({name,&field}).second;
	addData(name, bitpit::VTKFieldType::VECTOR, loc, bitpit::VTKDataType::Float64);
	return check;
};

/*! Remove a data field from the list */
void VTK_DATAMESH::removeField(std::string & name){
	std::unordered_map<std::string, dvector1D * >::iterator gotS = m_scalar.find(name);
	std::unordered_map<std::string, dvecarr3E * >::iterator gotV = m_vector.find(name);
	
	if(gotS != m_scalar.end())	{
		m_scalar.erase(gotS);
		removeData(name);
	}

	if(gotV != m_vector.end())	{
		m_vector.erase(gotV);
		removeData(name);
	}
};

/*! Remove all field from the list */
void VTK_DATAMESH::removeAllFields(){
	for(auto && p : m_scalar){
		removeData(p.first);
	}
	
	for(auto && v : m_vector){
		removeData(v.first);
	}
};

/*! CRTP method to write VTU files. Not intended for public usage purpose. 
 *  Custom implementation of VTK_UnstructuredGrid base class in BITPIT_BASE library. 
 *  Check BITPIT_BASE documentation for more information  */
void VTK_DATAMESH::flushData(  fstream &str, bitpit::VTKFormat codex_, string name  )
{
	int n;
	
	if(m_points == NULL || m_connectivity == NULL ){std::cout<<"not linked Mesh Data Structure"<<endl; return;}
	
	string indent("         ") ;
	
	if( codex_ == bitpit::VTKFormat::ASCII){
		if( name == "Points" ){
			for( auto && p : *m_points) {
				bitpit::genericIO::flushASCII( str, indent ) ;
				bitpit::genericIO::flushASCII( str, 3, p) ;}
				str<<endl;
		};
		
		if( name == "connectivity" ){
			for( auto && cc : *m_connectivity) {
				bitpit::genericIO::flushASCII( str, indent ) ;
				bitpit::genericIO::flushASCII( str, 3, cc) ;
				str<<endl;
			};
		};
		
		if( name == "types" ){
			for(auto && cc : *m_connectivity) {
				bitpit::genericIO::flushASCII( str, indent ) ;
				bitpit::genericIO::flushASCII( str, (int)m_element ) ;
				str<<endl;
			};
		};
		
		if( name == "offsets" ){
			int off_(0) ;
			for(auto && cc : *m_connectivity) {
				off_ += bitpit::vtk::getNNodeInElement( m_element ) ;
				bitpit::genericIO::flushASCII( str, indent ) ;
				bitpit::genericIO::flushASCII( str, off_  ) ;
				str<<endl;
			};
		};
		
	}else{		
		
		if( name == "Points" ){
			for(auto && p : *m_points) {
				bitpit::genericIO::flushBINARY( str, p);
			}
		};
		
		if( name == "connectivity" ){
			for(auto && cc : *m_connectivity) {bitpit::genericIO::flushBINARY( str,cc);}
		};
		
		if( name == "types"){
		
			for( auto && cc : *m_connectivity) bitpit::genericIO::flushBINARY( str, (int)m_element ) ;
		};
		
		if( name == "offsets"){
			int off_(0) ;
			for(auto && cc : *m_connectivity) {
				off_ += bitpit::vtk::getNNodeInElement( m_element ) ;
				bitpit::genericIO::flushBINARY( str, off_  ) ;
			};
		};
	}
	
	std::unordered_map<std::string, dvector1D * >::iterator itscalar = m_scalar.find(name);
	std::unordered_map<std::string, dvecarr3E * >::iterator itvector = m_vector.find(name);
	
	if(itscalar != m_scalar.end()){
		flushScalarField(name, codex_,str);
	}
	
	if(itvector != m_vector.end()){
		flushVectorField(name, codex_,str);
	}
	
	return ;
};


/*!
 * flush your particular scalar field data 
 */
void 	VTK_DATAMESH::flushScalarField(std::string name, bitpit::VTKFormat codex_, std::fstream & str){
	
	bitpit::VTKField **myfield;
	bool check = getFieldByName(name, myfield);
	bitpit::VTKLocation datatype = (*myfield)->getLocation();	
	
	string indent("         ") ;
	
	dvector1D * datafield = m_scalar[name]; 
	if( codex_ == bitpit::VTKFormat::ASCII){	
		
		if(datatype==bitpit::VTKLocation::POINT && datafield->size() == m_points->size()){
			
			for(auto && data : *datafield) {
				bitpit::genericIO::flushASCII( str, indent ) ;
				bitpit::genericIO::flushASCII( str, data);
				str<<endl;
			};
		};
		
		if(datatype==bitpit::VTKLocation::CELL && datafield->size() == m_connectivity->size()){
			
			for(auto && data : *datafield) {
				bitpit::genericIO::flushASCII( str, indent ) ;
				bitpit::genericIO::flushASCII( str, data);
				str<<endl;
			};
		};
	
		
	}else{
	
		if(datatype==bitpit::VTKLocation::POINT && datafield->size() == m_points->size()){
			
			for(auto && data : *datafield) {
				bitpit::genericIO::flushBINARY( str, data);
			};
		};
		
		if(datatype==bitpit::VTKLocation::CELL && datafield->size() == m_connectivity->size()){
			for(auto && data : *datafield) {
				bitpit::genericIO::flushBINARY( str, data);
			};
		};
	}
	
	
	datafield = NULL;
}

/*!
 * flush your particular vector field data 
 */
void 	VTK_DATAMESH::flushVectorField(std::string name, bitpit::VTKFormat codex_, std::fstream & str){
	
	bitpit::VTKField ** myfield;
	bool check = getFieldByName(name, myfield);
	bitpit::VTKLocation datatype = (*myfield)->getLocation();	
	
	string indent("         ") ;
	
	dvecarr3E * datafield = m_vector[name]; 
	if( codex_ == bitpit::VTKFormat::ASCII){	
		
		if(datatype==bitpit::VTKLocation::POINT && datafield->size() == m_points->size()){
			
			for(auto && data : *datafield) {
				bitpit::genericIO::flushASCII( str, indent ) ;
				bitpit::genericIO::flushASCII( str, 3, data);
				str<<endl;
			};
		};
		
		if(datatype==bitpit::VTKLocation::CELL && datafield->size() == m_connectivity->size()){
			
			for(auto && data : *datafield) {
				bitpit::genericIO::flushASCII( str, indent ) ;
				bitpit::genericIO::flushASCII( str, 3, data);
				str<<endl;
			};
		};
		
		
	}else{
		
		if(datatype==bitpit::VTKLocation::POINT && datafield->size() == m_points->size()){
			
			//myfield->setElements(m_points->size());
			//myfield->setOffset(m_points->size());
			
			for(auto && data : *datafield) {
				bitpit::genericIO::flushBINARY( str, data);
			};
		};
		
		if(datatype==bitpit::VTKLocation::CELL && datafield->size() == m_connectivity->size()){
			
			//myfield->setElements(m_connectivity->size());
			//myfield->setOffset(m_connectivity->size());
			
			for(auto && data : *datafield) {
				bitpit::genericIO::flushBINARY( str, data);
			};
		};
	};

	datafield = NULL;
}
