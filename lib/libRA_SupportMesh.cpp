#include "libRA_SupportMesh.hpp"
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// BASE_Support Implementation
/*
 *	\date			31/1/2016
 *	\authors		Federico Gallizio
 * 	\authors 		Scafroglia Ugo
 *	\authors		Rocco Arpa
 *	\copyright		Copyright 2016 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief Service Data Structure for portion of Unstructured Triangulated Meshes 
 *
 *	BASE_Support is a Service Data Structure class for Unstructured Meshes, directly derived from Class_SurfTri BitPit object of libMesh.
 *	It has no purposes of Data Structure Container, but is used to store informations of portion/patch of a preexistent Mother 
 * 	Unstructured Meshes. Retains tools for adapting and managing locally the Mother Mesh, without altering its original layout,
 *	as well as tools for extracting/handling the boundary contour of the patch.
 * 
 */

/*! Basic Constructor */ 
BASE_Support::BASE_Support(){
	name = "";
	classType = "BASE_Support";
	mother = NULL;
};

/*! Custom Constructor. Builds a patch of a Mother triangulation, given an integer flag.
 *  The integer flag identifies an ensemble of simplicies belonging to Mother Triangulation. See SHAPE::markInclusion member
 *  documentation.
 * \param[in] mother_ pointer to a SHAPE Data Structure
 * \param[in] flag   integer ID
 */ 
BASE_Support::BASE_Support(SHAPE *mother_, int flag){
	name = "";
	classType = "BASE_Support";
	setSupport(mother_, flag);
};

/*! Custom Constructor. Builds a patch of a Mother triangulation, passing a list of Simplicies belonging to it, by index.
 * \param[in] mother_ pointer to a Class_SurfTri Data Structure
 * \param[in] flag   list of Simplex indices
 */
BASE_Support::BASE_Support(Class_SurfTri *mother_, ivector1D & list){
	name = "";
	classType = "BASE_Support";
	setSupport(mother_, list);
};

/*! Basic Destructor */
BASE_Support::~BASE_Support(){
	cleanSupport();
};

/*! Copy Constructor.
 * \param[in] other BASE_Support object where copy from
 */
BASE_Support::BASE_Support(const BASE_Support &other){
	
	*(static_cast<Class_SurfTri *>(this)) = *(static_cast< const Class_SurfTri *>(&other));
	name = other.name;
	classType = other.classType;
	mother = other.mother;
	mapVertMother = other.mapVertMother;
	currentSelection = other.currentSelection;
	minInvConn = other.minInvConn;
	
};

/*! Copy Operator.
 * \param[in] other BASE_Support object where copy from
 */
BASE_Support & BASE_Support::operator=(const BASE_Support &other){
	
	*(static_cast<Class_SurfTri *>(this)) = *(static_cast< const Class_SurfTri *>(&other));
	name = other.name;
	classType = other.classType;
	mother = other.mother;
	mapVertMother = other.mapVertMother;
	currentSelection = other.currentSelection;
	minInvConn = other.minInvConn;
	return(*this);
};

/*! Clean Data stored in your BASE_Support object */
void BASE_Support::cleanSupport(){

	//free Surftri stuffs
	nVertex = 0;
	nSimplex = 0;
	freeContainer(Vertex);
	freeContainer(Simplex);
	freeContainer(Normal);
	freeContainer(Adjacency);
	
	name = "";
	mother = NULL;
	freeContainer(currentSelection);
	freeContainer(minInvConn);
	freeContainer(mapVertMother);

};

/*! Builds a patch of a Mother triangulation, given an integer flag.
 *  The integer flag identifies an ensemble of simplicies belonging to Mother Triangulation. See SHAPE::markInclusion member
 *  documentation. The method erases any class data(if any) and replace them with the new infos.
 * \param[in] mother pointer to a SHAPE Data Structure
 * \param[in] flag   integer ID
 */ 
void BASE_Support::setSupport(SHAPE * mother_, int flag){
	
	ivector1D list =mother_->extractPidded(flag);
	setSupport(mother_, list);
};

/*! Builds a patch of a Mother triangulation, given a list of its simplicies .
 *  The method erases any class data(if any) and replace them with the new infos.
 * \param[in] mother pointer to a Class_SurfTri Data Structure
 * \param[in] flag   list of indices of selected simplicies
 */ 
void BASE_Support::setSupport(Class_SurfTri *mother_, ivector1D & list){
	
	if(getMother() != NULL){cleanSupport();}
	
	int sizeList = list.size();
	//get vertex-Simplex Connectivity and Normals
	for (int i = 0; i < sizeList; ++i) {
		AddSimplex(mother_->Simplex[list[i]]);
	}	
	
	// Update normal & adjacency
	ResizeAdjacency();
	for (int i = 0; i < sizeList; ++i) {
		Adjacency[i].resize(mother_->Simplex[list[i]].size());
		for (int j = 0; j < mother_->Simplex[list[i]].size(); ++j) {
			
			if (mother_->Adjacency[list[i]][j][0] == -1) {
				Adjacency[i][j].push_back(-1);
			}else {
				for ( int k=0; k < mother_->Adjacency[list[i]][j].size(); ++k){
					
					int neigh = posVectorFind(list, mother_->Adjacency[list[i]][j][k]); 
					if(neigh != -1){	
						Adjacency[i][j].push_back(neigh);
					}
				}
				if(Adjacency[i][j].size() ==0){
					Adjacency[i][j].push_back(-1);
				}	
			}
		}
	}	
	
	//set class own members
	mother = mother_;
	currentSelection = list;
	compMinInvConn();
	//done
};

/*! Set an optional name to recognize your current support
 * \param[in] name 
 */
void BASE_Support::setName(std::string name_){name = name_;};

/*! Return current name of the support  
 * \param[out] result optional name of the object 
 */
std::string BASE_Support::getName(){return(name);};

/*! Return current class Type of the object  
 * \param[out] result get class type string 
 */
std::string BASE_Support::getClassType(){return(classType);};

/*! Return pointer to current Mother triangulation  
 * \param[out] result Class_SurfTri Mother 
 */
Class_SurfTri * BASE_Support::getMother(){return(mother);};

/*! Return a list of neighbor vertices in the 1-Ring of a given target vertex. I/O indexing 
 * is referred to the Mother tessellation. 
 * \param[in] iVert  target vertex in Mother tessellation numbering
 * \param[out] result list of 1-Ring Vertices in Mother tessellation numbering
 */ 
ivector1D BASE_Support::getVNeighs( int iVert)
{
	ivector1D result;
	ivector1D listS = getSNeighs(iVert);
	
	std::map<int,int> mapV;
	std::map<int,int>::iterator itM;
	
	for(int i=0; i<listS.size(); ++i){
		for(int j=0; j<Simplex[listS[i]].size(); ++j){
				mapV[Simplex[listS[i]][j]] = Simplex[listS[i]][j];
		}
	}
	mapV.erase(iVert);
	
	result.resize(mapV.size());
	int counter=0;
	for(itM=mapV.begin(); itM!=mapV.end(); ++itM){
		result[counter] = itM->second;
		++counter;
	}
	return(result);	
};

/*! Return a list of simplex in the 1-Ring of a given target vertex. Input vertex indexing 
 * is referred to the Mother tessellation, output simplex indexing is local, referred to the current support. 
 * \param[in] iVert  target vertex in Mother tessellation numbering
 * \param[out] result list of 1-Ring Simplex in local numbering
 */ 
ivector1D BASE_Support::getSNeighs( int iVert){
	
	int vPos = posVectorFind(mapVertMother,iVert);
	if(vPos == -1){return(ivector1D());}
	
	//get belonging local triangle of iVert
	int iT = minInvConn[vPos];
	//get local indexing
	int vT = posVectorFind(Simplex[iT], iVert);
	
	bool check;
	ivector1D result = Ring_1(iT,vT, check);
	return(result);
};

/*! Return list of simplicies composing support
 * \param[out] result list of simplicies -> mother indexing
 */
ivector1D BASE_Support::getSimplexMap(){return(currentSelection);};

/*! Return list of vertices composing support
 * \param[out] result list of vertices -> mother indexing
 */
ivector1D BASE_Support::getVertexMap(){return(mapVertMother);};

/*! Return Normals associated to each support simplex
 * \param[out] result normals
 */
dvecarr3E BASE_Support::getNormals(){
	
	dvecarr3E result(nSimplex);
	for(int i=0; i<nSimplex; ++i){
		result[i] = mother->Normal[currentSelection[i]];
	}
	return(result);
};

/*! Return Normals associated to each support vertex. 
 * If Vertex Normals are not defined in your Mother tesselation, generate them
 * \param[out] result vertex normals
 */
dvecarr3E BASE_Support::getVNormals(){
	
	if(mother->VNormal.size() == 0){mother->GenerateVNormals();}
	
	dvecarr3E result(mapVertMother.size());
	for(int i=0; i<mapVertMother.size(); ++i){
		result[i] = mother->VNormal[mapVertMother[i]];
	}
	return(result);
};


/*! Return Normal associated to a given support simplex in local support simplex-indexing
 * \param[out] result normals
 */
darray3E BASE_Support::getNormals(int sIndex){
	
	darray3E zero; zero.fill(0.0);
	if(sIndex<0 || sIndex >=nSimplex){return zero;}
	return(mother->Normal[currentSelection[sIndex]]);
};

/*! Return Normal associated to a given support vertex in mother tesselation indexing. 
 * If Vertex Normals are not defined in your Mother tesselation, generate them
 * \param[out] result vertex normals
 */
darray3E BASE_Support::getVNormals(int iVert){
	
	if(mother->VNormal.size() == 0){mother->GenerateVNormals();}
	darray3E zero; zero.fill(0.0);
	if(!checkVectorFind(mapVertMother, iVert)){return zero;}
	return(mother->VNormal[iVert]);
};

/*! Return list of vertex composing the whole boundary contour of your support
 * \param[out] result list of boundary vertex indices
 */
ivector1D BASE_Support::getWholeBoundaryIndex(){

	ivector1D result;
	std::map<int,int> map;
	std::map<int,int>::iterator it;
	
	for(int T=0; T<nSimplex; ++T){
		int n=Adjacency[T].size();
		for(int i=0; i<n; ++i){
			if(Adjacency[T][i][0]< 0){
				map[Simplex[T][i]] = Simplex[T][i];
			}
		}
		
	}
	
	int counter=0;
	result.resize(map.size());
	for(it=map.begin(); it!= map.end(); ++it){
		result[counter]= it->second;
		++counter;
	}

	return(result);
};

/*! Get index of those support boundary vertices, internal to the mother tessellation
 * \param[out] result list of constrained boundary vertices, by index
 */  
ivector1D BASE_Support::getConstrBoundaryIndex(){
	
	ivector1D mBV = getMotherBoundaryVIndex();
	ivector1D lBV = getWholeBoundaryIndex();
	
	ivector1D result(lBV.size());
	int counter=0;
	for(int i=0; i<lBV.size(); ++i){
		if(!checkVectorFind(mBV, lBV[i])){
			result[counter]= lBV[i];
			++counter;
		}
	}
	result.resize(counter);
	return(result);
};

/*! Get those support boundary vertices belonging also to the mother's boundary contour
 * \param[out] result list of free boundary vertices, by index
 */  
ivector1D BASE_Support::getFreeBoundaryIndex(){
	ivector1D mBV = getMotherBoundaryVIndex();
	ivector1D lBV = getWholeBoundaryIndex();
	
	ivector1D result(lBV.size());
	int counter=0;
	for(int i=0; i<lBV.size(); ++i){
		if(checkVectorFind(mBV, lBV[i])){
			result[counter]= lBV[i];
			++counter;
		}
	}
	result.resize(counter);
	return(result);
};

/*! Return list of vertex composing the whole boundary contour of your support
 * \param[out] result list of boundary vertices
 */
dvecarr3E BASE_Support::getWholeBoundary(){
	
	ivector1D list = getWholeBoundaryIndex();
	return(indexToPoints(list));
};

/*! Get index of those support boundary vertices, internal to the mother tessellation
 * \param[out] result list of constrained boundary vertices
 */  
dvecarr3E BASE_Support::getConstrBoundary(){
	
	ivector1D list = getConstrBoundaryIndex();
	return(indexToPoints(list));
};

/*! Get those support boundary vertices belonging also to the mother's boundary contour
 * \param[out] result list of free boundary vertices
 */  
dvecarr3E BASE_Support::getFreeBoundary(){
	
	ivector1D list = getFreeBoundaryIndex();
	return(indexToPoints(list));
};

/*! Extract Points from Mother Triangulation, given a selected list of vertex indices 
 * \param[in] map list of selected vertex indices
 * \param[out] result return list of vertex coordinates 
 */
dvecarr3E BASE_Support::indexToPoints(ivector1D & map){
	
	//preliminary check
	int sMap = map.size();
	if(sMap ==0){
		return(dvecarr3E());
	}
	int minMap, maxMap;
	minval(map,minMap);
	maxval(map,maxMap);
	
	if(minMap < 0 && maxMap >= mother->nVertex){
		std::cout<<"Invalid argument passed in BASE_Support::indexToPoints"<<endl;
		return(dvecarr3E());
		
	}
	
	
	dvecarr3E result(sMap);
	
	for(int i=0; i<sMap; ++i){
		result[i] = mother->Vertex[map[i]];
	}
	return(result);
};

/*! Get index of simplex in Mother tessellation, passing its local support index
 *\param[in] loc_target local simplex index
 *\param[out] result simplex index in Mother  
 */ 
int BASE_Support::simp_LocToMother(int loc_target){
	if(loc_target < 0 && loc_target >= nSimplex){
		std::cout<<"Invalid argument passed in BASE_Support::simp_LocToMother"<<endl;
		return(-1);
	}
	
	return(currentSelection[loc_target]);
};

/*! Get index of simplex in local support, passing its global index in Mother tessellation.
 * If simplex is not available in current support, returns -1.
 * \param[in] glb_target simplex index in Mother
 *\param[out] result   local simplex index
 */ 
int BASE_Support::simp_MotherToLoc(int glb_target){
	return(posVectorFind(currentSelection, glb_target));
};

/*! Get indices of a simplex list in Mother tessellation, passing its local support index list.
 * Irregular input values are skipped in the output list.
 * \param[in] loc_list local simplex index list
 *\param[out] result simplex index list in Mother  
 */ 
ivector1D BASE_Support::simp_LocToMother(ivector1D &loc_list){
	
	ivector1D result(loc_list.size());
	int counter =0;
	
	for(int i=0; i<loc_list.size(); ++i){
		result[counter] = simp_LocToMother(loc_list[i]);
		++counter;
	}
	result.resize(counter);
	return(result);
};

/*! Get index of simplex in local support, passing its global index in Mother tessellation.
 * rregular input values are skipped in the output list.
 * \param[in] glb_target simplex index in Mother
 *\param[out] result   local simplex index
 */ 
ivector1D BASE_Support::simp_MotherToLoc(ivector1D &glb_list){
	
	ivector1D result(glb_list.size());
	int counter =0;
	
	for(int i=0; i<glb_list.size(); ++i){
		result[counter] = simp_MotherToLoc(glb_list[i]);
		++counter;
	}
	result.resize(counter);
	return(result);
};

/*! check a List of Vertex indices and cut off all those vertices than do not belong to Support tessellation
 * \param[in] list  original list;
 * \param[out] result cleaned list;
 */
ivector1D BASE_Support::cleanVertexList(ivector1D & list){
	
	ivector1D result(list.size());
	int counter=0;
	for(int i=0; i<list.size(); ++i){
		if(checkVectorFind(mapVertMother, list[i])){
			result[counter]=list[i];
			++counter;
		}
	}
	
	result.resize(counter);
	return(result);
};

/*! check a List of Simplex indices (local indexing) and cut off all those simplicies than do not belong to Support tessellation
 * \param[in] list  original list;
 * \param[out] result cleaned list;
 */
ivector1D BASE_Support::cleanSimplexList(ivector1D & list){
	
	ivector1D result(list.size());
	int counter=0;
	for(int i=0; i<list.size(); ++i){
		if(checkVectorFind(currentSelection, list[i])){
			result[counter]=list[i];
			++counter;
		}
	}
	
	result.resize(counter);
	return(result);
};

/*! Create a boundary 3D curve on your current support, passing a list of Support vertices.
 * \param[in] list indices of support vertices
 * \param[out] result 3DCurve stored in a SurfTri object
 */
Class_SurfTri BASE_Support::createBoundary3DCurve(ivector1D & vList){
	
	
	//check input list
	ivector1D vListClean = vList;
	int mapT = vListClean.size();
	if( mapT == 0){return(Class_SurfTri());}
	
	// Local variables
	Class_SurfTri boundaryCurve;
	
	for (int iV = 0; iV < mapT; ++iV) {
		//extract neighbours
		int V = vListClean[iV];
		ivector1D Vneighs = getVNeighs(V);
		int n = Vneighs.size();
		
		ivector1D trueNeighs(n,-1);
		int counter=0;
		for(int i=0; i<n; ++i){
			if(checkVectorFind(vListClean, Vneighs[i])){
				trueNeighs[counter]=Vneighs[i];
				++counter;
			}
			
		}
		trueNeighs.resize(counter);
		
		//refine trueNeighs search. From all possible connection 
		// 1) get the boundary edges, 
		if(counter>2){
			
			ivector1D Sneighs = getSNeighs(V);
			int nT=Sneighs.size();
			
			std::map<int,int> map;
			std::map<int,int>::iterator it;
			
			for(int iK=0; iK<counter;++iK){
				for(int iT=0; iT<nT; ++iT){
					if(checkVectorFind(Simplex[Sneighs[iT]], trueNeighs[iK])){
						int target = Sneighs[iT];
						int locEdge = edge(target,V,trueNeighs[iK]);
						if(Adjacency[target][locEdge][0] == -1){
							map[trueNeighs[iK]] = trueNeighs[iK];
							//dummy_vec[countNew] = trueNeighs[iK];
							//++countNew;
						}
					}
				}
			}
			
			int countNew=0;
			trueNeighs.resize(map.size());
			for(it=map.begin(); it !=map.end(); it++){
				trueNeighs[countNew] = it->second;
				++countNew;
			}
			counter = countNew;
		}
		
		//adding vertex
		boundaryCurve.AddVertex(mother->Vertex[V]);
		//create simplicies
		ivector1D simplexD(2, -1);
		simplexD[0]=iV;
		for(int i=0; i<counter; ++i){
			simplexD[1]=posVectorFind(vListClean,trueNeighs[i]);
			boundaryCurve.AddSimplex(simplexD);
		}
	} //next iV
	
  	boundaryCurve.RemoveTrueDoubleSimplex();
	boundaryCurve.RemoveIsolatedVertex();
 	boundaryCurve.ResizeVertex();
 	boundaryCurve.ResizeSimplex();
	return(boundaryCurve);
};

/*! Plot current Support in VTK Unstructured Grid format
 * \param[in] dir folder path
 * \param[in] name file name
 * \param[in] flag codex boolean 0-"ascii", 1-"appended"
 * \param[in] vertices OPTIONAL, external vertex list to be plotted. If NULL, plot original vertex list
 */ 
void BASE_Support::plotVTU(std::string dir, std::string name, bool flag, dvecarr3E * vertices){
	
	Class_SurfTri supportDum = extractWholeTess(vertices);
	std::string codex = "ascii";
	if(flag){codex="appended";}
	VTK_SHAPE handle_vtk_output(&supportDum, dir, name, codex);
	handle_vtk_output.SetGeomTypes("Float64", "Int32", "Int32", "Int32");
	handle_vtk_output.Write();
};

/*! Plot scalar fields on current Support in VTK Unstructured Grid format
 * \param[in] dir folder path
 * \param[in] name file name
 * \param[in] flag codex boolean 0-"ascii", 1-"appended"
 * \param[in] counter integer counter to tag filename
 * \param[in] field scalar field of doubles
 * \param[in] vertices OPTIONAL, external vertex list to be plotted. If NULL, plot original vertex list
 */ 
void BASE_Support::plotScalarVTU(std::string dir, std::string name, bool flag, int counter, dvector1D & field, dvecarr3E * vertices){
	
	Class_SurfTri supportDum = extractWholeTess(vertices);
	std::string codex = "ascii";
	if(flag){codex="appended";}
	VTK_SHAPE handle_vtk_output(&supportDum, dir, name, codex);
	
	handle_vtk_output.SetGeomTypes("Float64", "Int32", "Int32", "Int32");
	
	if(counter>=0){handle_vtk_output.SetCounter(counter);}
	
	std::string loc_;
	if(field.size() == supportDum.nVertex){loc_="Point";}
	if(field.size() == supportDum.nSimplex){loc_="Cell";}
	
	if(loc_.empty()){std::cout<<"scalar field not suitable to be written on mesh, nor by points nor by cells"<<endl; return;}
	
	handle_vtk_output.AddData("scalarField",1, "Float64", loc_, codex);
	handle_vtk_output.linkScalarField(field);
	handle_vtk_output.Write();
};

/*! Plot vector fields on current Support in VTK Unstructured Grid format
 * \param[in] dir folder path
 * \param[in] name file name
 * \param[in] flag codex boolean 0-"ascii", 1-"appended"
 * \param[in] counter integer counter to tag filename
 * \param[in] field vector field of 3-elements arrays of doubles
 * \param[in] vertices OPTIONAL, external vertex list to be plotted. If NULL, plot original vertex list
 */ 
void BASE_Support::plotVectorVTU(std::string dir, std::string name, bool flag, int counter, dvecarr3E & field, dvecarr3E * vertices){

	Class_SurfTri supportDum = extractWholeTess(vertices);
	std::string codex = "ascii";
	if(flag){codex="appended";}
	VTK_SHAPE handle_vtk_output(&supportDum, dir, name, codex);
	
	handle_vtk_output.SetGeomTypes("Float64", "Int32", "Int32", "Int32");
	
	if(counter>=0){handle_vtk_output.SetCounter(counter);}
	
	std::string loc_;
	if(field.size() == supportDum.nVertex){loc_="Point";}
	if(field.size() == supportDum.nSimplex){loc_="Cell";}
	
	if(loc_.empty()){std::cout<<"scalar field not suitable to be written on mesh, nor by points nor by cells"<<endl; return;}
	
	handle_vtk_output.AddData("vectorField",3, "Float64", loc_, codex);
	handle_vtk_output.linkVectorField(field);
	handle_vtk_output.Write();
};

/*! Plot current Support in STL format
 * \param[in] filename output filename -> absolute
 * \param[in] vertices vertex list to be plotted
 * \param[in] flag codex boolean 0-"ascii", 1-"binary"
 */ 
void BASE_Support::plotSTL(string filename,  dvecarr3E & vertices, bool flag){
	
	STL_obj STL_solid(trim(filename), flag);
	Class_SurfTri supportDum = extractWholeTess(&vertices);
	supportDum.GenerateNormals();
	STL_solid.save("", supportDum.nVertex, supportDum.nSimplex, supportDum.Vertex, supportDum.Normal,
		       supportDum.Simplex);
};

/*! Plot current Support in Dune grid format format
 * \param[in] filename output filename -> absolute
 * \param[in] vertices vertex list to be plotted
 */ 
void BASE_Support::plotDGF(string filename,  dvecarr3E & vertices){

	DGF_obj DGF(trim(filename));
	Class_SurfTri supportDum = extractWholeTess(&vertices);
	DGF.save(supportDum.nVertex, supportDum.nSimplex, supportDum.Vertex,supportDum.Simplex);
};

/*! Plot the whole boundary in VTK cloud points or Unstructured 3D curve 
 * \param[in] dir output directory
 * \param[in] name output filename
 * \param[in] flag boolean for codex 0-"ascii",1-"appended"
 * \param[in] type  plot in "3DCurve" style or "Cloud". Default is "Cloud". Other specified types get the default one. 
 */
void BASE_Support::plotWholeBoundary(std::string dir, std::string name, bool flag, std::string type){

	//extractBoundary
	ivector1D bound = getWholeBoundaryIndex();
	plotBoundary(dir, name, flag, type, bound);
};

/*! Plot the constrained boundary in VTK cloud points or Unstructured 3D curve 
 * \param[in] dir output directory
 * \param[in] name output filename
 * \param[in] flag boolean for codex 0-"ascii",1-"appended"
 * \param[in] type  plot in "3DCurve" style or "Cloud". Default is "Cloud". Other specified types get the default one. 
 */
void BASE_Support::plotConstrBoundary(std::string dir, std::string name, bool flag, std::string type){
	
	//extractBoundary
	ivector1D bound = getConstrBoundaryIndex();
	plotBoundary(dir, name, flag, type, bound);
};

/*! Plot the free boundary in VTK cloud points or Unstructured 3D curve 
 * \param[in] dir output directory
 * \param[in] name output filename
 * \param[in] flag boolean for codex 0-"ascii",1-"appended"
 * \param[in] type  plot in "3DCurve" style or "Cloud". Default is "Cloud". Other specified types get the default one. 
 */
void BASE_Support::plotFreeBoundary(std::string dir, std::string name, bool flag, std::string type){
	
	//extractBoundary
	ivector1D bound = getFreeBoundaryIndex();
	plotBoundary(dir, name, flag, type, bound);
};

/*! Plot a generic portion of the boundary, passed by an external vertex list, in VTK cloud points or Unstructured 3D curve 
 * \param[in] dir output directory
 * \param[in] name output filename
 * \param[in] flag boolean for codex 0-"ascii",1-"appended"
 * \param[in] type  plot in "3DCurve" style or "Cloud". Default is "Cloud". Other specified types get the default one. 
 * \param[in] list list of vertex of boundary subset, passed by their indices(mother indexing)
 */
void BASE_Support::plotBoundary(std::string dir, std::string name, bool flag, std::string type, ivector1D & list){
	
	ivector1D list2 = cleanVertexList(list);
	if(list2.size() == 0){return;}
	
	std::string codex = "ascii";
	if(flag){codex="appended";}
	
	//check type
	bool typ_int = (type =="3DCurve");
	if(typ_int){
		Class_SurfTri boundaryCurve = createBoundary3DCurve(list2); 
		VTK_3DCURVE handle_vtk_output(boundaryCurve, dir, name, codex);
		handle_vtk_output.SetGeomTypes("Float64", "Int32", "Int32", "Int32");
		handle_vtk_output.Write();
		
	}else{
		dvecarr3E pt = indexToPoints(list2);
		VTK_BASICCLOUD handle_vtk_output(dir, name, codex, pt.size());
		handle_vtk_output.SetGeomTypes("Float64", "Int32", "Int32", "Int32");
		handle_vtk_output.linkData(pt);
		handle_vtk_output.Write();
	}
	
};

/*! Compile the minimum inverse Connectivity map (local Vertex->local Simplex) of the current support and
 *  extract the list of Vertex indices currently present in the Support.  
 *  The inverse connectivity map is minimal since only one index of the possible simplex the vertex belongs to is stored.
 *  Results are stored in the class members minInvConn and mapVertMother */
void BASE_Support::compMinInvConn(){

	freeContainer(minInvConn);
	std::map<int,int> tempMap;
	std::map<int,int>::iterator itMap;
	
	for(int i=0; i<nSimplex; ++i){
		for(int j=0; j<Simplex[i].size(); ++j){
			tempMap[Simplex[i][j]] = i;
		}	
	}
	minInvConn.resize(tempMap.size(), -1);
	mapVertMother.resize(tempMap.size(),-1);
	
	int counter=0;
	for(itMap=tempMap.begin(); itMap !=tempMap.end(); ++itMap){
		minInvConn[counter] = itMap->second;
		mapVertMother[counter]=itMap->first;
		++counter;
	}
};

/*! Set classType string 
 * \param[in] name_ string naming class;
 */
void BASE_Support::setClassType(std::string name_){
	classType = name_;
}


/*! Extract list of vertices on boundary countour of your Mother tessellation */ 
ivector1D BASE_Support::getMotherBoundaryVIndex(){
	
	ivector1D result;
	std::map<int,int> map;
	std::map<int,int>::iterator it;
	
	for(int T=0; T<mother->nSimplex; ++T){
		int n=mother->Simplex[T].size();
		for(int i=0; i<n; ++i){
			if(mother->Adjacency[T][i][0]<0){
				map[mother->Simplex[T][i]]=mother->Simplex[T][i];
			}
		}
		
	}
	
	int counter = 0;
	result.resize(map.size());
	for(it=map.begin(); it != map.end(); ++it){
		result[counter] = it->second;
		++counter;
	}
	
	return(result);
};

/*! Utility to get an essential copy of a BASE_Support object in full indipendent SurfTri object
 * \param[in] optional external list of vertex, equal in dimension to mother->Vertex structure
 * \param[out] result Class_SurfTri copy
 */
Class_SurfTri BASE_Support::extractWholeTess(dvecarr3E * vertices){
	
	Class_SurfTri result;
	result.nVertex = mapVertMother.size();
	result.nSimplex = nSimplex;
	result.ResizeVertex();
	result.ResizeSimplex();
	
	if(vertices != NULL && vertices->size() != mother->nVertex){
		for(int i=0; i<result.nVertex; ++i){
			result.Vertex[i] = (*vertices)[mapVertMother[i]];
		}	
	}else{	
		for(int i=0; i<result.nVertex; ++i){
			result.Vertex[i] = mother->Vertex[mapVertMother[i]];
		}	
	}
	
	for(int i=0; i<nSimplex; ++i){
		ivector1D dum;
		for(int j=0; j<Simplex[i].size(); ++j){
			int pos = posVectorFind(mapVertMother, Simplex[i][j]);
			dum.push_back(pos);
		}
		result.Simplex[i]=dum;
	}
	return(result);
};


//TODO capire se serve
// void  SupportMesh::cleanCornerNodes()
// { //clean automatically corner nodes and corner simplicies associated to it
// 
// 	ivector1D bmap = Daughter.FindFreeSimplex();
// 	ivector1D cornerList(bmap.size(),0);
// 	int counter = 0;
// 	int bmap_size= bmap.size();
// 	for (int i=0; i<bmap_size; i++)
// 	{
// 		int count = 0;
// 		for(int j=0; j<3; j++) {if(Daughter.Adjacency[bmap[i]][j][0] == -1) count++;}
// 		if(count > 1) {cornerList[counter] = bmap[i]; counter++;}
// 	}//next i
// 
// 	cornerList.resize(counter);
// 
// 	if(cornerList.size() > 0){
// 
// 	ivector1D cleanSim;
// 	Daughter.RemoveSimplex(cornerList,cleanSim);
// 	Daughter.RemapCellData(cleanSim,currentSelection);
// 
// 	ivector1D cleanMap;
// 	Daughter.RemoveIsolatedVertex(cleanMap);
// 	Daughter.RemapPointData(cleanMap, mapMotherV);
// 
// 	Daughter.ResizeVertex();
// 	Daughter.ResizeSimplex();
// 	Daughter.ResizeNormal();
// 	Daughter.ResizeAdjacency();
// 
// 	//refresh inverseconnectivity
// 	//evaluate the inverse connectivity map
// 	freeContainer(inverseConn);
// 	inverseConn.resize(Daughter.nVertex);
// 
// 	for (int iT = 0; iT < Daughter.nSimplex; ++iT) {
// 		for (int locV = 0; locV < 3; ++locV) {
// 			inverseConn[Daughter.Simplex[iT][locV]].push_back(iT);
// 		} //next local indices locV;
// 	} //next triangle iT;
// 
// 	}
// }//end of cleanCornerNodes
// 

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SupportSplitB Implementation
/*
 *	\date			31/1/2016
 *	\authors		Federico Gallizio
 * 	\authors 		Scafroglia Ugo
 *	\authors		Rocco Arpa
 *	\copyright		Copyright 2016 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief Service Data Structure for portion of Unstructured Triangulated Meshes, with boundary splitting methods 
 *
 *	SupportSplitB is a class derived from BASE_Support. It adds further features to the base class, for splitting the 
 * 	boundary contour in four sub-branches, whenever is possible.
 * 
 */

/*! Basic Constructor */
SupportSplitB::SupportSplitB(){
	setClassType("SupportSplitB");
	typeSegmentation=2;
	closedLoop=true;
	switchBoundary=true;
	seed.fill(1e18);
};
/*! Custom Constructor. Builds a patch of a Mother triangulation, given an integer flag.
 *  The integer flag identifies an ensemble of simplicies belonging to Mother Triangulation. See SHAPE::markInclusion member
 *  documentation.
 * \param[in] mother pointer to a SHAPE Data Structure
 * \param[in] flag   integer ID
 */ 
SupportSplitB::SupportSplitB(SHAPE *mother_, int flag){
	setClassType("SupportSplitB");
	typeSegmentation=2;
	switchBoundary=true;
	setSupport(mother_, flag);
	seed.fill(1e18);
};

/*! Custom Constructor. Builds a patch of a Mother triangulation, passing a list of Simplicies belonging to it, by index.
 * \param[in] mother pointer to a Class_SurfTri Data Structure
 * \param[in] flag   list of Simplex indices
 */
SupportSplitB::SupportSplitB(Class_SurfTri *mother_, ivector1D & list){
	setClassType("SupportSplitB");
	typeSegmentation=2;
	switchBoundary=true;
	setSupport(mother_, list);
	seed.fill(1e18);
};

/*! Basic Destructor */
SupportSplitB::~SupportSplitB(){
	cleanSupport();
};

/*! Copy Constructor.
 * \param[in] other SupportSplitB object where copy from
 */
SupportSplitB::SupportSplitB(const SupportSplitB &other){
	
	*(static_cast<BASE_Support *>(this)) = *(static_cast< const BASE_Support *>(&other));
	boundSegmentedMap = other.boundSegmentedMap;
	boundCurve = other.boundCurve;
	mapOrderBoundary = other.mapOrderBoundary;
	typeSegmentation = other.typeSegmentation;
	seed = other.seed;
	list_geo = other.list_geo;
	splitPoints = other.splitPoints;
	switchBoundary= other.switchBoundary;
	closedLoop= other.closedLoop;  
};

/*! Copy Operator.
 * \param[in] other SupportSplitB object where copy from
 */
SupportSplitB & SupportSplitB::operator=(const SupportSplitB &other){
	
	*(static_cast<BASE_Support *>(this)) = *(static_cast< const BASE_Support *>(&other));
	boundSegmentedMap = other.boundSegmentedMap;
	boundCurve = other.boundCurve;
	mapOrderBoundary = other.mapOrderBoundary;
	typeSegmentation = other.typeSegmentation;
	seed = other.seed;
	list_geo = other.list_geo;
	splitPoints = other.splitPoints;
	switchBoundary= other.switchBoundary;
	closedLoop= other.closedLoop;
	return(*this);
};

/*! Clean Data stored in your SupportSplitB object */
void SupportSplitB::cleanSupport(){
	
	//free Surftri stuffs
	nVertex = 0;
	nSimplex = 0;
	freeContainer(Vertex);
	freeContainer(Simplex);
	freeContainer(Normal);
	freeContainer(Adjacency);
	
	
	setName("");
	mother = NULL;
	freeContainer(currentSelection);
	freeContainer(minInvConn);
	freeContainer(mapVertMother);
	
	boundCurve.nVertex = 0;
	boundCurve.nSimplex= 0;
	freeContainer(boundCurve.Vertex);
	freeContainer(boundCurve.Simplex);
	freeContainer(boundCurve.Normal);
	freeContainer(boundCurve.Adjacency);
	
	cleanSplitBoundaries();
	freeContainer(mapOrderBoundary);
	freeContainer(boundSegmentedMap);
	closedLoop = true;
	switchBoundary=true;
	
};

/*! Builds a patch of a Mother triangulation, given an integer flag.
 *  The integer flag identifies an ensemble of simplicies belonging to Mother Triangulation. See SHAPE::markInclusion member
 *  documentation. The method erases any class data(if any) and replace them with the new infos. Create the boundary structure 
 *  boundCurve and sort it. If Boundary does not exist an internal flag will inhibit all boundary manipulation methods 
 * \param[in] mother pointer to a SHAPE Data Structure
 * \param[in] flag   integer ID
 */ 
void SupportSplitB::setSupport(SHAPE *mother_, int flag){
	
	ivector1D list =mother_->extractPidded(flag);
	setSupport(mother_, list);
};

/*! Builds a patch of a Mother triangulation, given a list of its simplicies .
 *  The method erases any class data(if any) and replace them with the new infos.Create the boundary structure 
 *  boundCurve and sort it. If Boundary does not exist an internal flag will inhibit all boundary manipulation methods
 * \param[in] mother pointer to a Class_SurfTri Data Structure
 * \param[in] flag   list of indices of selected simplicies
 */ 
void SupportSplitB::setSupport(Class_SurfTri *mother_, ivector1D & list){
	
	if(getMother() != NULL){cleanSupport();}
	int sizeList = list.size();
	//get vertex-Simplex Connectivity and Normals
	for (int i = 0; i < sizeList; ++i) {
		AddSimplex(mother_->Simplex[list[i]]);
	}
	
	// Update normal & adjacency
	ResizeAdjacency();
	for (int i = 0; i < sizeList; ++i) {
		Adjacency[i].resize(mother_->Simplex[list[i]].size());
		for (int j=0; j < mother_->Simplex[list[i]].size(); ++j) {
						
			if (mother_->Adjacency[list[i]][j][0] == -1) {
				Adjacency[i][j].push_back(-1);
			}else {
				for ( int k=0; k < mother_->Adjacency[list[i]][j].size(); ++k){
								
					int neigh = posVectorFind(list, mother_->Adjacency[list[i]][j][k]); 
					if(neigh != -1){	
						Adjacency[i][j].push_back(neigh);
					}
				}
				if(Adjacency[i][j].size() ==0) Adjacency[i][j].push_back(-1);
			}
		}
	}	

	//set class own members
	mother = mother_;
	currentSelection = list;
	compMinInvConn();
	//done
	
	//now istantiate your boundary curve
	createBoundary();
	//check for closed loop or ambiguos boundaries
	if(boundCurve.nVertex != 0){closedLoop=true;}
	if(!CG_PLCurve::IsClosed(boundCurve.Adjacency) ) {
		std::cout<<"the curve selected is not closed "<<endl;
		closedLoop=false;}
	
	//then sort it
	if(closedLoop ){
		sortBoundary(mapOrderBoundary[0]);
		boundSegmentedMap.resize(mapOrderBoundary.size(), -1);
	}
};

/*!Clean boundary management structure only*/
void SupportSplitB::cleanSplitBoundaries(){
	
	freeContainer(boundSegmentedMap);
	boundSegmentedMap.resize(mapOrderBoundary.size(), -1);
	seed.fill(1e18);
	freeContainer(list_geo);
	freeContainer(splitPoints);
	typeSegmentation=2;
};

/*! Segment your boundary according to the current list of split points stored in member SplitPoints
 * \param[out] result  result boolean true/false for "clean exit w/ errors?" monitoring  
 */
bool SupportSplitB::splitBoundary(){
	
	if(!closedLoop || !switchBoundary){return(false);}
		
	typeSegmentation=2;
	int sP = splitPoints.size();
	if(sP == 0){return(false);}
	
	freeContainer(boundSegmentedMap);
	boundSegmentedMap.resize(mapOrderBoundary.size(),-1);
	
	//found local position of splitPoints on the ordered Vertex boundary map
	ivector1D locIndex(sP);
	int iMark = 0;
	for(int i=0; i<sP; ++i){
		locIndex[i] = posVectorFind(mapOrderBoundary, splitPoints[i] );
		if(boundSegmentedMap[locIndex[i]] == 0){
			iMark = i;
		}
	}
	
	for(int i=0; i<sP; ++i){
		int j = i;
		int k = (i+1)%sP;
		int cVal = (sP + i -iMark)%sP;
		fillVectorSubset(locIndex[j], locIndex[k], boundSegmentedMap, cVal);
	}
		
	return(true);
};

/*! Segment your boundary getting common boundaries of adjacent triangulated meshes
 * \param[in] geoList list of external superficial mesh, adjacent to current support
 * \param[out] result boolean true/false for "clean exit w/ errors?" monitoring
 */
bool SupportSplitB::splitBoundary_geoList(svector1D & geoList){
	
	if(!closedLoop || !switchBoundary){return(false);}
	
	cleanSplitBoundaries();
	typeSegmentation=1;
	
	bool check = true;
	for(int ll=0; ll<geoList.size(); ++ll){
		
		Class_SurfTri local_geo = extractBoundaryFromGeoFile(geoList[ll]);
		
		ivector1D match = compareCurves(boundCurve,local_geo);
		int size_match = match.size();
		
		for(int k=0; k<size_match; ++k){
			
			boundSegmentedMap[match[k]]=ll;
		}
	}
	
	//extract splitPoints
	int sz = boundSegmentedMap.size();
	int valSx = boundSegmentedMap[sz-1]; 
	for(int i=0; i<sz; i++){
		if(boundSegmentedMap[i] !=valSx){
			splitPoints.push_back(mapOrderBoundary[i]);
			valSx = boundSegmentedMap[i];
		}
	}
	
	//resplit boundary for non perfect geometry adjacency 
	return((check && splitBoundary()));	
};

/*! Segment your boundary automatically in an indicative number N_ of branches.
 * It tries to find N_ significant split points in your contour, according to a criterium of 
 * maximum angle formed by two consecutive edges. If a number of split points less then N_ is found,
 * store those found split points, but does not segment the boundary and return w/ a boolean error flag.
 * \param[in] seed starting seed of your segmentation
 * \param[in] N_   number of desired branch
 * \param[out] n_eff number of effective split points found
 * \param[out] result boolean true/false for "clean exit w/ errors?" monitoring
 */
bool SupportSplitB::splitBoundary_auto(darray3E seed_, int N_, int n_eff){
	
	if(!closedLoop){return(false);}
	
	cleanSplitBoundaries();
	typeSegmentation=0;
	
	sortBoundary(seed_);
	
	dvector1D angles;
	ivector1D mapV = regularizeBoundary(angles) ;
	cout<<angles<<endl;	
	ivector1D extremaLoc = findSPoints(mapV, angles, N_);
	cout<<extremaLoc<<endl;
	n_eff = extremaLoc.size();
	cout<<n_eff<<endl;
	if(N_ != n_eff){return(false);}	
	
	//get extremaLoc vertex in mother support indexing & reorder them 
	std::map<int,int> mExtr;
	std::map<int,int>::iterator itExtr;
	for(int i=0; i<n_eff; ++i){mExtr[extremaLoc[i]] = extremaLoc[i];}
	
	splitPoints.resize(mExtr.size());
	int counter=0;
	for(itExtr=mExtr.begin(); itExtr != mExtr.end(); ++itExtr){
		int loc = itExtr->second;
		splitPoints[counter] = mapOrderBoundary[loc];
		++counter;
	}
	
	
	//split your boundary
	bool check = splitBoundary();	

	return(check);
};

/*! Add a split point to your boundary and segment it.
 *\param[in] point  point coordinates. If point not belong to boundary, take the nearest boundary vertex to it 
 *\param[out] result boolean true/false for "clean exit w/ errors" monitoring 
 */
bool SupportSplitB::addSplitPoint(darray3E point){
	if(!closedLoop || !switchBoundary){return(false);}
	return(addSplitPoint(point,true));
};

/*! Add a split point to your boundary.
 * \param[in] point  point coordinates. If point not belong to boundary, take the nearest boundary vertex to it 
 *\param[in] force  boolean flag, if true apply also segmentation of your boundary, if false store only your point. 
 *\param[out] result boolean true/false for "clean exit w/ errors" monitoring 
 */
bool SupportSplitB::addSplitPoint(darray3E point, bool force){
	if(!closedLoop || !switchBoundary){return(false);}
	int iV = getNearestBPointIndex(point);
	return(addSplitPoint(iV, force));
};

/*! Add a split point to your boundary and segment it.
 * \param[in] index boundary vertex index ->mother indexing. If point not belong to boundary,does nothing 
 * \param[out] result boolean true/false for "clean exit w/ errors?" monitoring
 */
bool SupportSplitB::addSplitPoint(int index){
	if(!closedLoop || !switchBoundary){return(false);}
	return(addSplitPoint(index, true));
};

/*! Add a split point to your boundary.
 * \param[in] index boundary vertex index ->mother indexing. If point not belong to boundary,does nothing 
 * \param[in] force  boolean flag, if true apply also segmentation of your boundary, if false store only your point.  
 * \param[out] result boolean true/false for "clean exit w/ errors?" monitoring
 */
bool SupportSplitB::addSplitPoint(int index, bool force ){
	if(!closedLoop || !switchBoundary){return(false);}
	
	typeSegmentation =2;
	int posIndex  = posVectorFind(mapOrderBoundary, index);
	if(posIndex == -1){return(false);}
	
	std::map<int, int> mapSplit;
	std::map<int, int>::iterator itSplit;
	mapSplit[posIndex] = index;
	for(int i=0; i<splitPoints.size(); ++i){
		int loc = posVectorFind(mapOrderBoundary, splitPoints[i]);
		mapSplit[loc] = splitPoints[i];
	}
	
	freeContainer(splitPoints);
	splitPoints.resize(mapSplit.size());
	int counter=0;
	for(itSplit=mapSplit.begin(); itSplit !=mapSplit.end(); ++itSplit){
		splitPoints[counter] = itSplit->second;
		++counter;
	}

	bool check = true;
	if(force){check = check && splitBoundary();}
	return(check);
};

/*! Remove a split point from your boundary segmentation and stitch related branches.
 * \param[in] point split vertex coordinate. If point not exactly belong to split point list, does nothing 
 * \param[out] result boolean true/false for "clean exit w/ errors?" monitoring
 */
bool SupportSplitB::removeSplitPoint(darray3E point){
	if(!closedLoop || !switchBoundary){return(false);}
	return(removeSplitPoint(point,true));
};

/*! Remove a split point from your boundary segmentation and stitch related branches.
 * \param[in] index split vertex index ->mother indexing. If point not belong to split point list, does nothing 
 * \param[out] result boolean true/false for "clean exit w/ errors?" monitoring
 */
bool SupportSplitB::removeSplitPoint(int index){
	if(!closedLoop || !switchBoundary){return(false);}
	return(removeSplitPoint(index,true));
};

/*! Remove a split point from your boundary segmentation.
 * \param[in] point split vertex coordinate. If point does not belong to split point list, does nothing 
 * \param[in] force  boolean flag, if true recalculate also segmentation of your boundary, if false throw away only your point.  
 * \param[out] result boolean true/false for "clean exit w/ errors?" monitoring 
 */
bool SupportSplitB::removeSplitPoint(darray3E point, bool force){
	if(!closedLoop || !switchBoundary){return(false);}
	
	int candidate = getNearestBPointIndex(point);
	int pos = posVectorFind(splitPoints, candidate);
	if(pos == -1) return(false);
	return(removeSplitPoint(splitPoints[pos], force));
};


/*! Remove a split point from your boundary segmentation.
 * \param[in] index split vertex index ->mother indexing. If point not belong to split point list, does nothing 
 * \param[in] force  boolean flag, if true recalculate also segmentation of your boundary, if false throw away only your point.  
 * \param[out] result boolean true/false for "clean exit w/ errors?" monitoring 
 */
bool SupportSplitB::removeSplitPoint(int index, bool force){
	if(!closedLoop || !switchBoundary){return(false);}
	
	typeSegmentation = 2;
	int posIndex  = posVectorFind(splitPoints, index);
	if(posIndex == -1){return(false);}
	
	ivector1D::iterator mSp = splitPoints.begin();
	for(int i=0; i<posIndex; ++i){++mSp;}
	splitPoints.erase(mSp);

	std::map<int, int> mapSplit;
	std::map<int, int>::iterator itSplit;
	for(int i=0; i<splitPoints.size(); ++i){
		int loc = posVectorFind(mapOrderBoundary, splitPoints[i]);
		mapSplit[loc] = splitPoints[i];
	}
	
	freeContainer(splitPoints);
	splitPoints.resize(mapSplit.size());
	int counter=0;
	for(itSplit=mapSplit.begin(); itSplit !=mapSplit.end(); ++itSplit){
		splitPoints[counter] = itSplit->second;
		++counter;
	}
	
	bool check = true;
	if(force){check = check && splitBoundary();}
	return(check);
};

/*! Add a list of split points to your boundary and segment it. Given points will be ordered automatically according to
 * the actual sorting of your boundary vertices.
 * \param[in] pointlist list point coordinates. If a point not belong to boundary, take the nearest boundary vertex to it 
 * \param[out] result boolean true/false for "clean exit w/ errors?" monitoring 
 */
bool SupportSplitB::addSplitPoints(dvecarr3E & pointlist){
	if(!closedLoop || !switchBoundary){return(false);}
	bool check = true;
	for(int i=0; i<pointlist.size(); ++i){
		check = check && addSplitPoint(pointlist[i], false);
	}
	check = check && splitBoundary();
	return(check);
};

/*! Add a list of split point to your boundary and segment it. Given points will be ordered automatically according to
 * the actual sorting of your boundary vertices.
 * \param[in] pointlist boundary vertex index list ->mother indexing. If point not belong to boundary,does nothing 
 * \param[out] result boolean true/false for "clean exit w/ errors?" monitoring
 */
bool SupportSplitB::addSplitPoints(ivector1D & pointlist){
	if(!closedLoop || !switchBoundary){return(false);}
	
	bool check = true;
	for(int i=0; i<pointlist.size(); ++i){
		check = check && addSplitPoint(pointlist[i], false);
	}
	check = check && splitBoundary();
	return(check);
};

/*! Remove a list of split points from your boundary segmentation and stitch related branches.
 * \param[in] pointlist split vertex coordinate list. If a point not belong to split point list, does nothing 
 * \param[out] result boolean true/false for "clean exit w/ errors?" monitoring
 */
bool SupportSplitB::removeSplitPoints(dvecarr3E & pointlist){
	if(!closedLoop || !switchBoundary){return(false);}
	
	bool check = true;
	for(int i=0; i<pointlist.size(); ++i){
		check = check && removeSplitPoint(pointlist[i], false);
	}
	check = check && splitBoundary();
	return(check);
};


/*! Remove a list of split points from your boundary segmentation and stitch related branches.
 * \param[in] pointlist split vertex index list ->mother indexing. If a point not belong to split point list, does nothing 
 * \param[out] result boolean true/false for "clean exit w/ errors?" monitoring
 */
bool SupportSplitB::removeSplitPoints(ivector1D & pointlist){
	if(!closedLoop || !switchBoundary){return(false);}
	
	bool check = true;
	for(int i=0; i<pointlist.size(); ++i){
		check = check && removeSplitPoint(pointlist[i], false);
	}
	check = check && splitBoundary();
	return(check);
};

/*! Return branch which a given point is nearest to/belongs to
 * \param[in] point point coordinates
 * \param[out] result branch index. If index <0, errors occurred 
 */
int SupportSplitB::whichBranch(darray3E point){
	if(!closedLoop || !switchBoundary){return(-1);}
	int iV = getNearestBPointIndex(point);
	return(whichBranch(iV));
};

/*! Return branch which a given boundary vertex belongs to
 *\param[in] index vertex index
 *\param[out] result branch index. If index <0, errors occurred
 */
int SupportSplitB::whichBranch(int index ){
	if(!closedLoop || !switchBoundary){return(-1);}
	int pos = posVectorFind(mapOrderBoundary, index);
	if(pos < 0){return(-1);}
	return(boundSegmentedMap[pos]);
};

/*! Change numeration of consecutive branches of your boundary segmentation starting 
 * from a target branch (counter-clockwise orientation) 
 * \param[in] branchIndex target branch index
 * \param[out] result error flag: if 0-successfull exit, if <0, errors occurred
 */
int SupportSplitB::orientSegmentation(int branchIndex){
	if(!closedLoop || !switchBoundary){return(-1);}
	
	int size = getBoundaryBranchSize();
	if(branchIndex<0 || branchIndex >= size){return(-1);}
	
	//update your boundSegmentedMap
	for(int i=0; i<boundSegmentedMap.size(); ++i){
		boundSegmentedMap[i] = (size + i - branchIndex)%size;
	}
	return(0);
};

/*! Sort your boundary curve, in counter-clockwise direction w.r.t. tessellation normal distribution, starting
 * from a given point. Member boundCurve is manipulated, new ordering is stored in mapOrderBoundary member.
 *\param[in] startVert starting point 
 *\param[out] result return 0-for successfull exit, -1 - errors occurred 
 */
int SupportSplitB::sortBoundary(darray3E &startVert ){
	
	if(!closedLoop || !switchBoundary){return(-1);}
	
	int index = getNearestBPointIndex(startVert);
	return(sortBoundary(index));
};

/*!Sort your boundary curve, in counter-clockwise direction w.r.t. tessellation normal distribution, starting
 * from a given point. Member boundCurve is manipulated, new ordering is stored in mapOrderBoundary member.
 *\param[in] startVert starting point index -> mother indexing
 *\param[out] result return 0-for successfull exit, -1 - errors occurred  
 */
int SupportSplitB::sortBoundary(int startVert){
	if(!closedLoop || !switchBoundary){return(-1);}
	
	int posVert = posVectorFind(mapOrderBoundary, startVert);
	
	if(	(boundCurve.nSimplex != boundCurve.nVertex) || boundCurve.nSimplex == 0 || 
		boundCurve.nVertex == 0 || posVert <0){
		return(-1);
	}
	
	if(boundCurve.Adjacency.size() == 0){boundCurve.BuildAdjacency();}
	if(boundSegmentedMap.size() ==0){boundSegmentedMap.resize(boundCurve.nSimplex, -1);}
	
	//get the new map of vertex
	ivector1D newMap(boundCurve.nVertex);
	newMap[0] = posVert;
	
	int tSim = 0;
	int counter = 0;
	{	//find Simplex containing the starting vertex
		ivector1D dumS(boundCurve.nSimplex);
		for(int i=0; i<boundCurve.nSimplex; ++i){dumS[i] = boundCurve.Simplex[i][0];}
		tSim = posVectorFind(dumS, posVert);
	}
	//compile map
	while(tSim != -1 && counter<(boundCurve.nSimplex-1) ){
		newMap[counter+1]   = boundCurve.Simplex[tSim][1];
		tSim = boundCurve.Adjacency[tSim][1][0];
		++counter;
	}
	//check problems
	if(counter != boundCurve.nSimplex-1){return(-1);}
	
	
	//reset BoundCurve member.
	//VertexList;
	//update mapOrderBoundary and boundSegmentedMap
	{ 
		dvecarr3E tempV = boundCurve.Vertex;
		ivector1D mapSeg = boundSegmentedMap;
		
		ivector1D mapOrdNew(boundCurve.nSimplex);
		
		for(int i=0; i<boundCurve.nSimplex; ++i){
			boundCurve.Vertex[i]=  tempV[newMap[i]];
			boundSegmentedMap[i] = mapSeg[newMap[i]];
			mapOrdNew[i] = mapOrderBoundary[newMap[i]];
		}
		mapOrderBoundary = mapOrdNew;
	}
	
	
	//recreate Connectivity and adjacency of an ordered curve loop
	for(int i=0; i<boundCurve.nSimplex; ++i){
		boundCurve.Simplex[i][0] = i;
		boundCurve.Simplex[i][1] = (i+1)%boundCurve.nSimplex;
		boundCurve.Adjacency[i][0][0] = (boundCurve.nSimplex + i -1)%boundCurve.nSimplex;
		boundCurve.Adjacency[i][1][0] = (i+1)%boundCurve.nSimplex;
	}
	
	return(0);
};

/*! Return the nearest boundary vertex to a given seed point
 *\param[in] target point
 *\param[out] result nearest boundary vertex index -> mother indexing
 */ 
int SupportSplitB::getNearestBPointIndex(darray3E target){
	if(!closedLoop || !switchBoundary){return(-1);}
	
	dvector1D dist(boundCurve.nVertex);
	for(int i=0; i<boundCurve.nVertex; ++i){
		dist[i] = norm_2(boundCurve.Vertex[i] - target);
	}
	
	double minValue;
	minval(dist, minValue);
	int pos = posVectorFind(dist, minValue);
	if(pos <0){return(pos);}
	return(mapOrderBoundary[pos]);
};

/*! Return type of boundary splitting tech currently used. 0-autoSplitting, 1-external adjacent Geometries Splitting,
 * 2-Manual Splitting. Default option is 2.
 * \param[out] result type of splitting
 */
int SupportSplitB::getSplittingType(){return(typeSegmentation);};

/*! Return list of ordered boundary points defining your splitted contour
 *\param[out] result split points of boundary contour
 */ 
dvecarr3E SupportSplitB::getSplitPoints(){
	
	dvecarr3E result(splitPoints.size());
	Class_SurfTri * mum = getMother();
	for(int i=0; i<splitPoints.size(); ++i){
		result[i] = mum->Vertex[splitPoints[i]];
	}
	mum=NULL;
	return(result);
};

/*! Return list of ordered boundary points defining your splitted contour
 * \param[out] result split points of boundary contour, by their mother index
 */ 
ivector1D SupportSplitB::getSplitPointsIndex(){return(splitPoints);};

/*! Return the initial seed point set for automatic Splitting. If type segmentation is set to manual,
 * or w/ external adjacent geometries, return {0,0,0}.
 * \param[out] result seed point
 */
darray3E SupportSplitB::getSeed(){return(seed);};

/*! Return the list of external adjacent geoemtries used for Splitting. If type segmentation is manual or
 * automatic, return a void list.
 * \param[out] result list of external geometries
 */
svector1D SupportSplitB::getGeoList(){return(list_geo);};

/*! Return current number of splitted branch in your boundary
 * \param[out] result current number of boundary branches
 */
int SupportSplitB::getBoundaryBranchSize(){return(splitPoints.size());};

/*! Return list of vertex contained in the "flag" branch.
 * \param[in] branch index
 * \param[out] result list of vertices, by index --> mother indexing
 */
ivector1D SupportSplitB::getBoundaryBranchIndex(int flag){

	if(flag <0 || flag >= getBoundaryBranchSize() ){return(ivector1D());}
	
	ivector1D result(boundSegmentedMap.size());
	int counter=0;
	for(int i=0; i<boundSegmentedMap.size(); ++i){
		if(boundSegmentedMap[i]==flag){
			result[counter] = mapOrderBoundary[i];
			++counter;
		}
	}
	result.resize(counter);
	return(result);
};

/*! Return list of vertex contained in the "flag" branch.
 * \param[in] branch index
 * \param[out] result list of vertex coordinates 
 */
dvecarr3E SupportSplitB::getBoundaryBranch(int flag){
	
	Class_SurfTri * mum = getMother();
	ivector1D list = getBoundaryBranchIndex(flag);
	dvecarr3E result(list.size());
	for(int i=0; i<list.size(); ++i){
		result[i] = mum->Vertex[list[i]];
	}
	mum=NULL;
	return(result);
};

/*! Return the whole branch map of your boundaries.
  * \param[out] result branch boundary map of vertex --> mother indexing 
 */
ivector2D SupportSplitB::getBoundaryBranchMap(){
	int size = getBoundaryBranchSize();
	
	ivector2D result(size);
	for(int i=0; i<size; i++){
		result[i] = getBoundaryBranchIndex(i);
	}
	return(result);
};

/*! Return a copy of your boundary structure
 * \param[out] result 3D curve stored in a surfTri object
 */
Class_SurfTri * SupportSplitB::getBoundaryCurve(){return(&boundCurve);};


/*! Activate/Deactivate Boundary manipulation methods. If Not active, the class behaves like a BASE_Support class essentially 
 * \param[in] OnOff boolean true-active, false- not active
 */
void SupportSplitB::switchBoundaryManipulation(bool OnOff){
	switchBoundary = OnOff;
};

/*! Return Boundary Manipulation status */
bool SupportSplitB::getBoundaryManipulationStatus(){return(switchBoundary);};

/*! Plot all branches of segmentation in different files, in VTK cloud points or Unstructured 3D curves 
 * \param[in] dir output directory
 * \param[in] name output root filename
 * \param[in] flag boolean for codex 0-"ascii",1-"appended"
 * \param[in] type  plot in "3DCurve" style or "Cloud". Default is "Cloud". Other specified types get the default one. 
 */
void SupportSplitB::plotBoundaryBranches(std::string dir, std::string name, bool flag, std::string type){
	
	int size = getBoundaryBranchSize();
	
	for(int i=0; i<size; ++i){
		//manipulate new loc name for i-th branch
		std::stringstream ss;
		ss<<name<<"_branch_"<<i;
		
		ivector1D list = getBoundaryBranchIndex(i);
		plotBoundary(dir, ss.str(), flag, type, list);
	}
	return;
};

/*! Decimate the 3D curve boundary boundCurve, so that angles between adjacent edges will be all convex(>=0).
 * w.r.t current boundary edge ordering.
 *\param[out] angles  list of angles between i-th and i-th +1 edges
 *\param[out] result map of decimated boundary vertices, in local boundCurve indexing  
 */
ivector1D SupportSplitB::regularizeBoundary(dvector1D & angles){
	
	ivector1D result (boundCurve.nVertex);
	for(int i=0; i<boundCurve.nVertex; ++i){result[i]=i;}
	
	angles.resize(boundCurve.nVertex, -1.0);

	//evaluate an average normal of the current support surface which the boundary belongs to
	darray3E SurfaceNormal;
	SurfaceNormal.fill(0.0);
		
	dvecarr3E locVNorm = getVNormals();

	for(int k=0; k<mapVertMother.size(); ++k){
		SurfaceNormal = SurfaceNormal + locVNorm[k]/((double) mapVertMother.size());
	}
		
	int counter = result.size();
	int result_size;
	//start decimating
	do{
		//calculate angles
		for(int i=0; i<counter; ++i){
				
			int iV = result[i];
			int iL = result[(i-1+counter)%counter];
			int iR = result[(i+1)%counter];
				
			darray3E n1, n2, wNorm, locNormal;
			n1.fill(0);
			n2.fill(0);
			wNorm.fill(0);
			locNormal = SurfaceNormal; 
				
			n1 = boundCurve.Vertex[iV] - boundCurve.Vertex[iL];
			n2 = boundCurve.Vertex[iR] - boundCurve.Vertex[iV];
				
			n1 = n1 - Dot_Product(n1,locNormal)*locNormal;
			n2 = n2 - Dot_Product(n2,locNormal)*locNormal;
				
			n1 = n1/norm_2(n1);
			n2 = n2/norm_2(n2);
				
			double proj = Dot_Product(n1,n2);
			wNorm = Cross_Product(n1,n2);
			if(abs(proj)>0.99) {wNorm = locNormal; proj=1.0;} 
				
			double sign = Dot_Product(wNorm,locNormal);
			if(abs(sign)< 0.01){
				sign=1.0;
				
			}else{
				sign = sign/abs(sign);
			}

			angles[i] = sign * acos(proj);
		}
		
 		cout<<angles<<endl;
 		cout<<"-----------------------"<<endl;		
		result_size = counter;
		
		counter = 0;
		//decimate results vector when angles value is negative
		for(int k=0; k<result_size; ++k){
			if(angles[k] >= 0){
				result[counter]=result[k];
				angles[counter]=angles[k];
				++counter;
			}
		}
		
		result.resize(counter);
		angles.resize(counter);
// 		cout<<result<<endl;
	}while(result_size != counter);

	return(result);
};

/*! Find the first guessN boundary vertices who share the greater edge-edge angle. 
 * Available points are always returned, even if they are lesser than guessN.
 * \param[in] map list of boundCurve vertex indices --> local indexing
 * \param[in] angles list of edge-edge angle associated to map.
 * \param[in] guessN number of desidered points.
 * \param[in] result avalaible points list, descending order w/ angles --> local boundCurve indexing
 */
ivector1D SupportSplitB::findSPoints(ivector1D & map, dvector1D & angles, int guessN){

	ivector1D result;
	int size = angles.size();
	bvector1D gradientAngles(size, true);
	dvector1D dummy(size);
	
	for (int i=0; i<size; ++i){
		int iL = (size+i-1)%size;
		int iR = (i+1)%size;
		dummy[i] =angles[iR] - 2*angles[i] + angles[iL];
		if((angles[i] - angles[iL])< 0.0){gradientAngles[i] = false;}
	}
	cout<<dummy<<endl;
	ofstream out;
	out.open("angles.dat");
	for(int i=0; i<size; i++){
	out<<i<<'\t'<<angles[i]<<'\t'<<dummy[i]<<endl;
	}
	out.close();
	ivector1D candidates;
	for (int i=0; i<size; ++i){
		int iR = (i+1)%size;
		if(gradientAngles[i] && !gradientAngles[iR]){
			candidates.push_back(i);
		}
	}
	
	std::map<double, int> mapCand;
	for (int i=0; i<candidates.size(); ++i){
		mapCand[angles[candidates[i]]] = candidates[i];
	}
	
	int checkSize = std::min((int)mapCand.size(), guessN);
	freeContainer(candidates);
	candidates.resize(checkSize);
	
	std::map<double,int>::iterator itFind = mapCand.end();
	int counter = 0;
	
	while(counter< checkSize){
		--itFind;
		candidates.push_back((*itFind).second);
		++counter;
	}
		
	//get the list of candidates;
	result.resize(checkSize);
	for(int k=0; k<checkSize; ++k){
		result[k] = map[candidates[k]];
	}

	return(result);
};

/*! Extract boundary contour from Surface Tessellation given from file. 
 *  Formats admissible are those reported in class SHAPE.
 * \param[in] filename abs path of your geometry file
 * \param[out] result boundary curve stored in SurfTri object
 */
Class_SurfTri SupportSplitB::extractBoundaryFromGeoFile(std::string & filename){
	
	Class_SurfTri curve;
	SHAPE geo;
	geo.init(filename);
	geo.cleaning();
	
	geo.Boundaries(curve);
	if(curve.Vertex.size()==0){
		std::cout<<"WARNING! your geometry "<<filename<< " has no free boundaries"<<endl;
	}
	
	return(curve);
};

/*! Extract common vertices between two 3D curves, mother and daughter. 
 * \param[in] mother first 3D curve 
 * \param[in] daughter second 3D curve 
 * \param[out] result list of common vertices, in mother indexing
 */
ivector1D SupportSplitB::compareCurves(Class_SurfTri & mother_,Class_SurfTri & daughter_){
	
	ivector1D result(mother_.nVertex);
	double tol = 1.0E-6;
	int counter=0;
	for(int i=0; i<mother_.nVertex; ++i){
		int j = 0;
		bool check = false;
		while(j<daughter_.nVertex && !check){
			
			double norm = norm_2(mother_.Vertex[i] - daughter_.Vertex[j]);
			if(norm <= tol){
				check=true;
				result[counter]=i;
				counter++;
			}
			++j;
		}
	}
	
	result.resize(counter);
	return(result);
};
	
/*! Adjust orientation of boundCurve member connectivity
 * \param[in] seed  index of a seed simplex
 */ 
void SupportSplitB::adjust3DCurveOrientation( int seed_)
{
	// check if your seed is available;
	if(boundCurve.nSimplex == 0 || seed_<0 ||seed_>= boundCurve.nSimplex){
		std::cout<<"seed not found in your curve. proceed with default value"<<endl; 
		seed_=0;
	}

	//check for built adjacency
	if(boundCurve.Adjacency.size()!= boundCurve.nSimplex){boundCurve.BuildAdjacency();}

	adjustSeedOrientation(seed_);
	
	//makes curve orientation homogeneous
	CG_PLCurve::AdjustOrder(boundCurve.Simplex, boundCurve.Adjacency, seed_);
	boundCurve.GenerateNormals();	
	
	return;
}//end method

/*! Adjust orientation of a seed simplex belonging to boundCurve
 * \param[in] seed  index of a seed simplex
 */ 
void SupportSplitB::adjustSeedOrientation( int seed_)
{
	
	 int V1 = mapOrderBoundary[boundCurve.Simplex[seed_][0]];	
	 int V2 = mapOrderBoundary[boundCurve.Simplex[seed_][1]];
	 
	 ivector1D S_V1 = getSNeighs(V1);
	 ivector1D S_V2 = getSNeighs(V2);
	 
	 int pos=-1;
	 int count=0;
	 while(pos==-1 && count<S_V1.size()){
		pos =posVectorFind(S_V2, S_V1[count]);
		count++;
	}
	 
	 if(pos==-1){return;}
	 int targetS = S_V2[pos];
	 
	 
	 int checkV1 = posVectorFind(Simplex[targetS], V1);
	 int checkV2 = (checkV1+1)%Simplex[targetS].size();
	 
	 if(V2 != Simplex[targetS][checkV2]){
		 int conn_temp = boundCurve.Simplex[seed_][0];
		 ivector1D adj_temp  = boundCurve.Adjacency[seed_][0];

		 boundCurve.Simplex[seed_][0] = boundCurve.Simplex[seed_][1];
		 freeContainer(boundCurve.Adjacency[seed_][0]);
		 boundCurve.Adjacency[seed_][0] = boundCurve.Adjacency[seed_][1];
		 
		 boundCurve.Simplex[seed_][1] = conn_temp;
		 freeContainer(boundCurve.Adjacency[seed_][1]);
		 boundCurve.Adjacency[seed_][1] = adj_temp;
		 
	}

	
}//end method

/*! Fill the boundCurve member of the class, with the boundary of the current support */
void SupportSplitB::createBoundary(){
	
	//check input list
	ivector1D vListClean= getWholeBoundaryIndex();
	int mapT = vListClean.size();
	if( mapT == 0){return;}
	
	// Local variables
	for (int iV = 0; iV < mapT; ++iV) {
		//extract neighbours
		int V = vListClean[iV];
		ivector1D Vneighs = getVNeighs(V);
		int n = Vneighs.size();
		
		ivector1D trueNeighs(n,-1);
		int counter=0;
		for(int i=0; i<n; ++i){
			if(checkVectorFind(vListClean, Vneighs[i])){
				trueNeighs[counter]=Vneighs[i];
				++counter;
			}
			
		}
		trueNeighs.resize(counter);
		
		//refine trueNeighs search. From all possible connection 
		// 1) get the boundary edges, 
		if(counter>2){
			ivector1D Sneighs = getSNeighs(V);
			int nT=Sneighs.size();
			std::map<int,int> map;
			std::map<int,int>::iterator it;
			
			for(int iK=0; iK<counter;++iK){
				for(int iT=0; iT<nT; ++iT){
					if(checkVectorFind(Simplex[Sneighs[iT]], trueNeighs[iK])){
						int target = Sneighs[iT];
						int locEdge = edge(target,V,trueNeighs[iK]);
						if(Adjacency[target][locEdge][0] == -1){
							map[trueNeighs[iK]] = trueNeighs[iK];
							//dummy_vec[countNew] = trueNeighs[iK];
							//++countNew;
						}
					}
				}
			}
			
			
			int countNew=0;
			freeContainer(trueNeighs);
			trueNeighs.resize(map.size());
			for(it=map.begin(); it !=map.end(); it++){
				trueNeighs[countNew] = it->second;
				++countNew;
			}
			counter = countNew;
		}
		//adding vertex
		boundCurve.AddVertex(mother->Vertex[V]);
		
		//create simplicies
		ivector1D simplexD(2, -1);
		simplexD[0]=iV;
		for(int i=0; i<counter; ++i){
			simplexD[1]=posVectorFind(vListClean,trueNeighs[i]);
			boundCurve.AddSimplex(simplexD);
		}
	} //next iV
	
	boundCurve.RemoveTrueDoubleSimplex();
	boundCurve.RemoveIsolatedVertex();
	boundCurve.ResizeVertex();
	boundCurve.ResizeSimplex();
	boundCurve.ResizeNormal();
	boundCurve.GenerateNormals();
	
 	mapOrderBoundary = vListClean;
 	if(boundCurve.nVertex != 0){
 		adjust3DCurveOrientation(0);
 	}	
}

