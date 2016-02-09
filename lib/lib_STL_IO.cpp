#include "lib_STL_IO.hpp"

//#####################################################################################################
// CLASS SHAPE IMPLEMENTATION
/*
 *	\date			31/12/2015
 *	\authors		Federico Gallizio
 * 	\authors		Allen Woody
 *	\authors		Rocco Arpa
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief Data Structure Manager Class for Unstructured Triangulated Meshes 
 *
 *	SHAPE is a Data Structure class for Unstructured Meshes, directly derived from Class_SurfTri BitPit object of libMesh.
 *	I/O can extend to Standard Tasselation Language,  *.stl(ASCII or bynary), BDF NASTRAN native format *.nas, *.dgf and VTK 
 *     	*.vtu formats. Retains tool for checking/cleaning the meshes and handling refinement.
 */	

/*! Default Constructor. Doing Nothing */
SHAPE::SHAPE() {
	filetype = 0;
	cleanData = false;
};

/*! Copy Constructor of the Class. Retains info from
 * \param[in] other SHAPE class
 */ 
SHAPE::SHAPE(const SHAPE & other){
	 filename = other.filename;
         filetype = other.filetype;  
              
         *(static_cast<Class_SurfTri * > (this) ) = *(static_cast<const Class_SurfTri * > (&other));       
         
	 InverseSimplexMap = other.InverseSimplexMap;
         deform_OR = other.deform_OR;
	 markRefinement = other.markRefinement;
	 markInclusion = other.markInclusion;
	 pid_surface = other.pid_surface;
	 moll_filter = other.moll_filter;
	 pidType = other.pidType;
	 cleanData = other.cleanData;
};

/*! Default Destructor. Free all structure of the class */
SHAPE::~SHAPE(){
    freeSHAPE();
}

/*! Copy Operator of the Class. Retains info from
 * \param[in] other SHAPE class
 */ 
SHAPE & SHAPE::operator=(const SHAPE & other){
  
	 freeSHAPE(); 
  	 filename = other.filename;
         filetype = other.filetype;  
              
         *(static_cast<Class_SurfTri * > (this) ) = *(static_cast<const Class_SurfTri * > (&other));       
         
	 InverseSimplexMap = other.InverseSimplexMap;
         deform_OR = other.deform_OR;
	 markRefinement = other.markRefinement;
	 markInclusion = other.markInclusion;
	 pid_surface = other.pid_surface;
	 moll_filter = other.moll_filter;
	 pidType = other.pidType;
	 cleanData = other.cleanData;
  
	 return(*this);
};

/*! Cleaner of the class. Retain all data to void structures */
void SHAPE::freeSHAPE() {

	filename = "";
	freeContainer(InverseSimplexMap);
	freeContainer(deform_OR);
	freeContainer(markRefinement);
	freeContainer(markInclusion);
	freeContainer(pid_surface);
	freeContainer(moll_filter);
	freeContainer(pidType);
	cleanData = false;
};

/*! Initialize Shape Class reading triangulated mesh from file and 
 *  store it in data structure. Formats allowed are *.stl, *.nas, *.vtu, *dgf.  Requires:
 * \param[in] input_file absolute path to geometry file
 */ 
void SHAPE::init(std::string input_file) {

	//setting input_file;
	filename = input_file;
	VTK_SHAPE handle_VTK_input((static_cast<Class_SurfTri * > (this)), input_file);
	
	std::string key2 = ".";
	std::size_t cut2 =input_file.find_last_of(key2);
	std::string ftype = input_file.substr(cut2);
	bool STLstatus;
	
	filetype = -1;
	if (ftype == ".stl" || ftype == ".STL") {
		filetype = 0;
		STLstatus = checkSTLfile(input_file);
	}
	if (ftype == ".nas" || ftype == ".NAS") {
		filetype = 1;
	}
	if (ftype == ".dgf" || ftype == ".DGF") {
		filetype = 2;
	}
	if (ftype == ".vtu" || ftype == ".VTU") {
		filetype = 3;
	}

	bool dummy;
	switch (filetype) {
	case 0: // reading the STL ASCII/BINARY................
		Import_stl(filename, STLstatus);
		break;
	case 1: // reading the nas format......................
		Import_nas(filename);
		break;
	case 2: //reading the dgf  format
		Import_dgf(filename);
		break;
	case 3: //reading the vtu format
		handle_VTK_input.linkTriMesh((static_cast<Class_SurfTri *> (this)));
		handle_VTK_input.Read();
		break;

	default: //read nothing
		std::cout << " No useful file format retrieved for reading..." << '\n';
		break;
	} //end switch

	if (nSimplex == 0) {
		std::cout << " No geometry available in file, or non existing file."
				<< '\n';
		std::cout << " Please check the file " << filename
				<< " or check your main control dictionary" << '\n';
		std::cout << "Now exiting..." << '\n';
		exit(1);
	}


	std::cout << "number of triangles read : " << nSimplex << endl;
	std::cout << "number of nodes read     : " << nVertex << endl;


	// building object triangulation(conn tri/tri, ricalculate normal)
	GenerateNormals();

	// Initializing data structure of the class appended to SurfTri
	deform_OR.resize(nVertex, darray3E{0.0,0.0,0.0});

	markRefinement.resize(nSimplex, 0); //initialize all simplex to status of blank/no refinement
	markInclusion.resize(nSimplex, -1); //initialize all simplex to status of not belong to lattices/selections

	//only in case pid_surface is not instantiated
	if (pid_surface.size() == 0)
		pid_surface.resize(nSimplex, -1);

		// check Cleaning troubles

		int nDV = CountDoubleVertex();

		std::cout << "Cleaning Alert--> found in the mother data structure:"<< '\n';
		std::cout << "Repeated  vertices  : " << nDV << '\n';
		if (nDV != 0)
			cleanData = true;
		if (nDV == 0) {
			BuildAdjacency();

			ivector1D List0Area = Find0AreaSimplex();
			int nIV = List0Area.size();
			int nIS = CountIsolatedSimplex();
			int nDS = CountTrueDoubleSimplex();

			std::cout << "Zero-Area triangles : " << nIV << '\n';
			std::cout << "Repeated  triangles : " << nDS << '\n';
			std::cout << "Isolated  triangles : " << nIS << '\n';

			if ((nIV + nIS + nDS) > 0) {
				cleanData = true;
			} else {
				GenerateVNormals();
				setInverseConnectivityMap();
			}
		}

}; // end of initializer

/*! Cleaning the data structure and create Normals and Triangle-Triangle connectivity */

void SHAPE::cleaning() { // cleaning triangulation (double vertex, double triangles, isolated vertex and isolated simplex.)
// plus build adjacency and recalculate Normals

	std::cout << "cleaning up geometry..." << '\n';

	ivector1D cleaningMap;
	bool check_normal, check_adjacency;

	 	 	 	 SetTolerance();
	 	 	 	 // cleaning up

	    	                 RemoveDoubleVertex();
				 BuildAdjacency();

				 { ivector1D list, cleaningMap;
				   list = Find0AreaSimplex();

				   if(list.size()>0){
	   	                     RemoveSimplex(list, cleaningMap);
	   	                     RemapCellData(cleaningMap, pid_surface);
			             }
				 }

				 { ivector1D list, cleaningMap;
				   list = FindIsolatedSimplex();
	    	                   RemoveSimplex(list, cleaningMap);
	   	                   RemapCellData(cleaningMap, pid_surface);
				 }

				    RemoveIsolatedVertex();

				 { ivector1D list, cleaningMap;
				   list = FindTrueDoubleSimplex();
				   RemoveSimplex(list, cleaningMap);
				   RemapCellData(cleaningMap, pid_surface);
				 }

	    ResizeVertex();
	    ResizeSimplex();
	    ResizeNormal();
	    ResizeAdjacency();

	    FixNodeNumb();
	    GenerateVNormals();


	std::cout << "number of triangles after cleaning : " << nSimplex << endl;
	std::cout << "number of nodes after cleaning:      " << nVertex  << endl;
	std::cout << "Normals & Adjacency regeneration done." << endl;
	std::cout << "Generated the Inverse Connectivity Map." << endl;

	// Resizing external data structure appended to SurfTri
	deform_OR.resize(nVertex);

	markRefinement.resize(nSimplex);
	markInclusion.resize(nSimplex);
	cleanData = false;
	// create inverse connectivity map
	setInverseConnectivityMap();
};
//end cleaning

 /*! Utility function for reading BDF Nastran *.nas file.
  * \param[in] fileinput string for input file;
  */
void SHAPE::Import_nas(std::string fileinput) {
	//====================================================================================================
	//READ THE NAS FILE
	//===================================================================================================
	// workspace
	std::string type, line, start, endF;
	std::string work2, work1;
	int n1, sorting_point, n2;
	int maxvalPID = 0;

	ivector1D pid_node, pid_triangle;

	int number_vertex = 0;
	int number_simplex = 0;

	std::ifstream indata(fileinput.c_str());
	// read first time to allocate structures;
	if (indata.is_open()) {
		do {
			std::getline(indata, start);
			endF = trim(start);
			std::stringstream ss;
			ss.str(start);
			ss >> type;
		} while (trim(type) != "GRID");

		// starting reading GRID and CTRIA3 info
		std::stringstream ss1;
		ss1.str(start);
		ss1 >> type;

		while (trim(type) == "GRID") {
			number_vertex++;
			std::stringstream ss2;
			getline(indata, line);
			ss2.str(line);
			ss2 >> type;
		}

		while (trim(type) == "CTRIA3") {
			number_simplex++;
			std::stringstream ss2;
			getline(indata, line);
			ss2.str(line);
			ss2 >> type;
		}

		// further reading for complete PID list.
		std::cout << "Recognized PID:             " << '\n';
		while (trim(line) != "ENDDATA") {
			if (trim(type) == "PSHELL") {
				std::stringstream ss;
				ss.str(line);
				std::string dummy;
				ss >> dummy >> n1;
				maxvalPID = max(maxvalPID, n1);
				std::cout << n1 << '\t';
				pidType.push_back(n1);
			}

			std::stringstream ss2;
			getline(indata, line);
			ss2.str(line);
			ss2 >> type;
		}
		std::cout << '\n';
	}

	indata.close();

	//Resizing class_surfTri.

	nVertex = number_vertex;
	nSimplex = number_simplex;
	ResizeVertex();
	ResizeSimplex(3);


	pid_surface.resize(number_simplex);
	std::map<int, int> nodeNames;

	//back to the file and reading vertices and connectivity;
	std::ifstream indata2(fileinput.c_str());
	if (indata2.is_open()) {

		do {
			std::stringstream ss3;
			getline(indata2, line);
			ss3.str(line);
			ss3 >> type;
		} while (trim(type) != "GRID");

		for (int iV = 0; iV < number_vertex; iV++) {
			std::string dummy;
			std::stringstream ss3;
			ss3.str(line);
			ss3 >> dummy >> work1;
			int etiquette = atoi(work1.c_str());
			nodeNames[etiquette] = iV;

			line = trim(line);
			work2 = line.substr(line.size() - 24, 24);

			darray3E vecDum = NASconvertVertex(work2);
			Vertex[iV] = vecDum;
			getline(indata2, line);
		}

		std::stringstream ss4;
		ss4.str(line);
		ss4 >> type;

		while (trim(type) != "CTRIA3") {
			std::stringstream ss3;
			getline(indata2, line);
			ss3.str(line);
			ss3 >> type;
		}

		for (int iT = 0; iT < number_simplex; iT++) {
			std::string dummy;
			std::stringstream ss3;
			ss3.str(line);
			ss3 >> dummy;
			ss3 >> n1 >> n2;
			pid_surface[iT] = (short int) n2;
			
			for (int k = 0; k < 3; k++) {
				ss3 >> sorting_point;
				Simplex[iT][k] = nodeNames[sorting_point];
			}
			getline(indata2, line);
		}

	}	   //end of indata2 reading
	indata2.close();

};

/*! Utility function for writing triangular mesh in VTU - Paraview format. Optional external vertex list can be specified. 
 *  If not, original vertices Vertex of the will be accounted for.
 *  \param[in] dir   directory path where the file will be written
 *  \param[in] name  name of the file
 *  \param[in] flag  boolean to activate ascii-writing (false) or binary-writing (true) codex
 *  \param[in] vertices	 external list of vertices, alternative to original vertices of the data structure (OPTIONAL)	  
 */
void SHAPE::plotVTU(std::string dir, std::string name, bool flag, dvecarr3E * vertices){
	std::string codex = "ascii";
	if(flag){codex="appended";}
	VTK_SHAPE handle_vtk_output((static_cast<Class_SurfTri *> (this)), dir, name, codex);
	if(vertices != NULL){handle_vtk_output.linkExternalVertexList(*vertices);}
	handle_vtk_output.SetGeomTypes("Float64", "Int32", "Int32", "Int32");
	handle_vtk_output.Write();
};
/*! Utility function for plotting a scalar field on triangular mesh in VTU - Paraview format. 
 *   Optional external vertex list for mesh can be specified. 
 *  If not, original vertices Vertex of the will be accounted for.
 *  \param[in] dir   directory path where the file will be written
 *  \param[in] name  name of the file
 *  \param[in] flag  boolean to activate ascii-writing (false) or binary-writing (true) codex
 *  \param[in] counter integer flag to mark name of your files.
 *  \param[in] field	scalar field to be plotted (can be given both by Cells or Points)
 *  \param[in] vertices	 external list of vertices, alternative to original vertices of the data structure (OPTIONAL)	  
 */
void SHAPE::plotScalarVTU(std::string dir, std::string name, bool flag, int counter, dvector1D& field, dvecarr3E * vertices)
{
	std::string codex = "ascii";
	if(flag){codex="appended";}
	VTK_SHAPE handle_vtk_output((static_cast<Class_SurfTri *> (this)), dir, name, codex);
	if(vertices != NULL){handle_vtk_output.linkExternalVertexList(*vertices);}

	if(counter>=0){handle_vtk_output.SetCounter(counter);}
	handle_vtk_output.SetGeomTypes("Float64", "Int32", "Int32", "Int32");

	std::string loc_;
	if(field.size() == nVertex){loc_="Point";}
	if(field.size() == nSimplex){loc_="Cell";}

	if(loc_.empty()){std::cout<<"scalar field not suitable to be written on mesh, nor by points nor by cells"<<endl; return;}
	handle_vtk_output.AddData("scalarField",1, "Float64", loc_, codex);
	handle_vtk_output.linkScalarField(field);
	handle_vtk_output.Write();
};

/*! Utility function for plotting a vector field on triangular mesh in VTU - Paraview format. 
 *   Optional external vertex list for mesh can be specified. 
 *  If not, original vertices Vertex of the will be accounted for.
 *  \param[in] dir   directory path where the file will be written
 *  \param[in] name  name of the file
 *  \param[in] flag  boolean to activate ascii-writing (false) or binary-writing (true) codex
 *  \param[in] counter integer flag to mark name of your files.
 *  \param[in] field	vector field to be plotted (can be given both by Cells or Points)
 *  \param[in] vertices	 external list of vertices, alternative to original vertices of the data structure (OPTIONAL)	
 * */
void SHAPE::plotVectorVTU(std::string dir, std::string name, bool flag, int counter, dvecarr3E & field, dvecarr3E * vertices)
{
	std::string codex = "ascii";
	if(flag){codex="appended";}
	VTK_SHAPE handle_vtk_output((static_cast<Class_SurfTri *> (this)), dir, name, codex);
	if(vertices != NULL){handle_vtk_output.linkExternalVertexList(*vertices);}
	handle_vtk_output.SetGeomTypes("Float64", "Int32", "Int32", "Int32");

	if(counter>=0){handle_vtk_output.SetCounter(counter);}

	std::string loc_;
	if(field.size() == nVertex){loc_="Point";}
	if(field.size() == nSimplex){loc_="Cell";}

	if(loc_.empty()){std::cout<<"scalar field not suitable to be written on mesh, nor by points nor by cells"<<endl; return;}

	handle_vtk_output.AddData("vectorField",3, "Float64", loc_, codex);
	handle_vtk_output.linkVectorField(field);
	handle_vtk_output.Write();
};

/*! Utility function for plotting triangular mesh in STL StereoLithography format. 
 *  \param[in] fileinput  destination file of your output 
 *  \param[in] vertices	 list of triangulation vertices	
 *  \param[in] flag  boolean to activate ascii-writing (false) or binary-writing (true) codex
 * */
void SHAPE::plotSTL(string fileinput, dvecarr3E & vertices, bool flag) {

	STL_obj STL_solid(trim(fileinput), flag);
	STL_solid.save("", nVertex, nSimplex, vertices, Normal,	Simplex);
};

/*! Utility function for plotting triangular mesh in dgf Dune Grid format. 
 *  \param[in] fileinput  destination file of your output 
 *  \param[in] vertices	 list of triangulation vertices	
 * */
void SHAPE::plotDGF(string fileinput,  dvecarr3E & vertices){

	DGF_obj DGF(trim(fileinput));
	DGF.save(nVertex,nSimplex, vertices, Simplex);
}
/*! Utility function to modify vertices list and connectivity in a pre-read structure of BDF Nastran *.nas.
 * The method is not suitable to transfer grid read in other formats in a *.nas one. 
 *  \param[in] fileinput  source file or your original input mesh. (*.nas format is mandatory). 
 *  \param[in] outfilename  destination file of your output 
 *  \param[in] vertices	 list of triangulation vertices	
*/ 
void SHAPE::modifyNAS(std::string fileinput, std::string outfilename,dvecarr3E & vertex) {
	
	
	std::size_t found = fileinput.find_last_of(".");
	if(fileinput.substr(found) !=".nas" && fileinput.substr(found) !=".NAS") {return;}
	
	std::string type, line, end;
	std::stringstream linedef;
	int kcount;
	
	std::ifstream indata(fileinput.c_str());
	std::ofstream outdata(outfilename.c_str());

	std::vector < std::string > number(3);
	std::size_t fd_e;
	std::string work1, work2, result, result2;
	int size_result, size_exp;

	int nvertices = nVertex;

	if (fileinput == outfilename) {
		std::cout
				<< " WARNING/ERROR: input *.NAS and output *.NAS files coincide. Overwriting is not possible!!"
				<< '\n';
		return;
	}

	//Open source file
	if (indata.is_open()) {

		do {
			std::getline(indata, line);
			end = trim(line);
			outdata << end << '\n'; //roughly copying infos on target outfile
		} while (trim(line) != "BEGIN BULK");

		// take pointer of indata to the first line after the last element of connectivity CTRIA3
		do {
			std::getline(indata, line);
			stringstream ss;
			ss.str(trim(line));
			ss >> type;
		} while (trim(type) == "GRID" || trim(type) == "CTRIA3");

		// Meanwhile, start writing nodes and connectivity on outdata

		// writing grid points
		for (int j = 0; j < nvertices; j++) {

			// writing 3 vertices kcount in a unique record of 24 lenght.
			for (int i = 0; i < 3; i++) {
				result2 = "";
				if ((abs(vertex[j][i]) > 1.0e4)
						|| (abs(vertex[j][i]) < 1.0e-2)) {
					std::stringstream ss_num;
					ss_num << std::setprecision(4) << std::scientific
							<< vertex[j][i];
					result = ss_num.str();
					size_result = result.size();

					fd_e = result.find("e");
					work1 = result.substr(fd_e + 1, size_result - fd_e);

					if (work1.substr(1, 1) == "0")
						work1.erase(1, 1);
					size_exp = work1.size();

					work2 = result.substr(0, 8 - size_exp);

					work2 += work1;

					if (work2.size() < 8)
						for (int j = 0; j < (8 - work2.size()); j++)
							result2 += " ";

					result2 += work2;
					number[i] = result2;
					result2 = "";
				} else {
					std::stringstream ss_num;
					ss_num << std::setprecision(8) << vertex[j][i];
					result = ss_num.str().substr(0, 8);

					if (result.size() < 8)
						for (int j = 0; j < (8 - result.size()); j++)
							result2 += " ";

					result2 += result;
					number[i] = result2;
					result2 = "";
				}
			} // end of loop on coords of vertex

			outdata << "GRID    " << resize_to_record8(j + 1) << "        "
					<< number[0] << number[1] << number[2] << '\n';

		} // end of loop on grid points

		for (int j = 0; j < nSimplex; j++) {
			outdata << "CTRIA3  " << resize_to_record8(j + 1)
					<< resize_to_record8(pid_surface[j])
					<< resize_to_record8(Simplex[j][0] + 1)
					<< resize_to_record8(Simplex[j][1] + 1)
					<< resize_to_record8(Simplex[j][2] + 1) << '\n';
		}

		// writing personal data finished. Come back copying the last part of indata file
		do {
			end = trim(line);
			outdata << end << '\n';
			std::getline(indata, line);
		} while (!line.empty());

	}
	indata.close();
	outdata.close();
};
// end of modifyNAS ....................................................................


/*! Utility function for writing portion of triangular mesh in VTU - Paraview format. Optional external vertex list can be specified. 
 *  If not, original vertices Vertex of the will be accounted for. 
 *  \param[in] dir   directory path where the file will be written
 *  \param[in] name  name of the file
 *  \param[in] flag  boolean to activate ascii-writing (false) or binary-writing (true) codex
 *  \param[in] counter counter of output files
 *  \param[in] Tlist   list of subportion triangle indices.
 *  \param[in] vertices	 external list of vertices, alternative to original vertices of the data structure (OPTIONAL)	  
 */
void SHAPE::plotExtraction(std::string dir, std::string name, bool flag, int counter, ivector1D & TList2, dvecarr3E * vertices){
	
	std::string codex = "ascii";
	if(flag){
		codex="appended";
	}

	ivector1D TList = checkTriangleList(TList2);
	//count and reorder vertices and simplicies
	std::map<int,int> mapV;
	std::map<int,int>::iterator itV;
	
	for(int i=0; i<TList.size(); ++i){
		for(int j=0; j<Simplex[TList[i]].size(); ++j){
			mapV[Simplex[TList[i]][j]] = Simplex[TList[i]][j]; 
		}
	}
	
	ivector1D listP(mapV.size());
	dvecarr3E activeP(mapV.size());
	int ctt = 0;
	if(vertices != NULL && vertices->size() == nVertex){
		for(itV = mapV.begin(); itV != mapV.end(); ++itV){
			int pos = itV->second;
			activeP[ctt] = (*vertices)[pos];
			listP[ctt] = pos;
			++ctt;
		}
	}else{
		for(itV = mapV.begin(); itV != mapV.end(); ++itV){
			int pos = itV->second;
			activeP[ctt] = Vertex[pos];
			listP[ctt]=  pos;
			++ctt;
		}
	}
	
	Class_SurfTri temp;
	temp.nVertex = ctt;
	temp.nSimplex = TList.size();
	temp.ResizeVertex();
	temp.Simplex.resize(temp.nSimplex);
	
	//push vertices;
	temp.Vertex = activeP;
	//push Simplicies	
	for(int i=0; i<TList.size(); ++i){
		
		ivector1D dummy(Simplex[TList[i]].size());
		
		for(int j=0; j<Simplex[TList[i]].size(); ++j){
		
			int posV =posVectorFind(listP, Simplex[TList[i]][j]);
			//update local connectivity
			dummy[j] = posV;
		}
		
		temp.Simplex[i] = dummy;
	}
	
	VTK_SHAPE handle_vtk_output(&temp, dir, name, codex);
	if(counter>=0){handle_vtk_output.SetCounter(counter);}
	handle_vtk_output.SetGeomTypes("Float64", "Int32", "Int32", "Int32");
	handle_vtk_output.Write();
};

/*! Utility function for writing portion of triangular mesh in a Point Cloud of VTU - Paraview format. Optional external vertex list can be specified. 
 *  If not, original vertices Vertex of the will be accounted for.
 *  \param[in] dir   directory path where the file will be written
 *  \param[in] name  name of the file
 *  \param[in] flag  boolean to activate ascii-writing (false) or binary-writing (true) codex
 *  \param[in] counter counter of output files
 *  \param[in] Vlist   list of subportion vertex indices.
 *  \param[in] vertices	 external list of vertices, alternative to original vertices of the data structure (OPTIONAL)	  
 */
void SHAPE::plotCloudExtraction(std::string dir, std::string name, bool flag, int counter, ivector1D & VList2, dvecarr3E * vertices){
	
	std::string codex = "ascii";
	if(flag){codex="appended";}
	
	ivector1D VList = checkVertexList(VList2);
	
	VTK_BASICCLOUD handle_vtk_output(dir, name, codex, VList.size());
	
	dvecarr3E activeP(VList.size());
	if(vertices != NULL && vertices->size() == nVertex){
		for(int i=0; i<VList.size(); i++){
			activeP[i] = (*vertices)[VList[i]];
		}  
	}
	else{
		for(int i=0; i<VList.size(); i++){
			activeP[i] = Vertex[VList[i]];
		}
	}
	
	if(counter>=0){handle_vtk_output.SetCounter(counter);}
	handle_vtk_output.SetGeomTypes("Float64","Int32", "Int32", "Int32");
	handle_vtk_output.linkData(activeP);
	handle_vtk_output.Write();
};


/*!Given a selection of Simplex indices, identify sub-patches (continous clusters of neighbor simplicies)
  * and returns them in ordered map of indices
  * \param[in] selection list of selected simplicies of mesh triangulation
  * \param[in] lumen     scalar field to guide the choice of best point candidate for each sub-patch.
  * \param[out] iCheck   boolean to handle successfull (true) sub-patch decomposition.	
  * \param[out] result   return sub-patches in ivector2D Matrix structure. Each sub-patches i return in first position Matrix[i][0]
  * 			 the best candidate simplex index.		
  */
ivector2D SHAPE::createPatchesGraph(ivector1D & selection,fvector1D & lumen, bool & iCheck) {

	ivector2D result(selection.size());

	bvector1D unChecked(selection.size() + 1, true); //+1 is due to cut an if statement in the graph search
	unChecked[selection.size()] = false;

	iCheck = true;

	//initialize the graph search
	int counterPatch = 0;

	while (summing_Vector(unChecked)) {
		int indexT = 0;

		ivector1D ivTarget; //((actualVertex-nVisited));

		while (!unChecked[indexT])
			indexT++;

		//ivTarget[counter] = selection[indexT];
		ivTarget.push_back(selection[indexT]);
		unChecked[indexT] = false;
		float lumenRef = lumen[ivTarget[0]];
		int pickRayCandidate = ivTarget[0];

		//checking a whole patch

		for (int iV = 0; iV < ivTarget.size(); iV++) {
			bool dummyFlag, checkFlag;
			int ITriangle = InverseSimplexMap[ivTarget[iV]][0];
			int locV = vertex(ITriangle, ivTarget[iV]);

			//TODO check c'e' forse un metodo uguale in surftri. Usa quello
			ivector1D ListRing = localVRing_1(ITriangle, locV, dummyFlag,
					checkFlag);
			if (!checkFlag) {
				iCheck = false;
				return (result);
			}

			for (int iN = 0; iN < ListRing.size(); iN++) {
				int iTarget = ListRing[iN];

				// check if iTarget belongs to the selection and if already visited;
				ivector1D::iterator itCheck = std::find(selection.begin(),
						selection.end(), iTarget);
				int pos = std::distance(selection.begin(), itCheck);

				if (itCheck != selection.end() && unChecked[pos]) {
					ivTarget.push_back(iTarget); // = iTarget;
					unChecked[pos] = false;
					if (lumen[pos] < lumenRef) {
						pickRayCandidate = iTarget;
						lumenRef = lumen[pos];
					}

				} // end if
			} //next iN
		} //loop on ivTarget map

		// calculate the pickrayCandidate point and rearrange the graph patch.

		result[counterPatch].resize(ivTarget.size() + 1);
		result[counterPatch][0] = pickRayCandidate;
		for (int j = 0; j < ivTarget.size(); j++)
			result[counterPatch][j + 1] = ivTarget[j];

		counterPatch++;
	}

	result.resize(counterPatch);
	return (result);
} //end of createPatchGraph

/*! Evaluate the Inverse connectivity map of the mesh triangulation (Simplex-Vertex) and store it
 * in the InverseSimplexMap member. It needs to be recalculated each time the data structure of the 
 * mesh changes. So Be careful when using it. Automatically clean the triangulation if not.
 */
void SHAPE::setInverseConnectivityMap() { 
  
	if (nVertex == 0 && nSimplex == 0) {
		std::cout
				<< "Warning! attempt a generation of An Inverse connectivity Matrix"
						"without a suitable data structure! program exiting....";
		exit(1);
	}
	if (cleanData)
		cleaning(); //make sure of warking with cleaned stuffs;

	if(InverseSimplexMap.size()!= 0){ freeContainer(InverseSimplexMap);};
	InverseSimplexMap.resize(nVertex);

	for (int iT = 0; iT < nSimplex; iT++) {
		for (int locV = 0; locV < 3; locV++) {
			InverseSimplexMap[Simplex[iT][locV]].push_back(iT);
		} //next local indices locV;
	} //next triangle iT;
} //end of setInverseConnectivityMap

/*!Extract those vertex-indices belonging to an ensemble of free edges in the mesh 
 *  triangulation. If empty vector is returned, triangulation is closed.*/
ivector1D SHAPE::extractFrontier() { 

	bvector1D temp(Vertex.size(), false);
	ivector1D result(Vertex.size(), -1);
	bool flag;
	int counter=0;

	//almighty loop on Simplex of Father tasselation
	for (int nT = 0; nT < nSimplex; nT++) {
		for (int ilocV = 0; ilocV < 3; ilocV++) {
			int iVert = Simplex[nT][ilocV];
			if (!temp[iVert] && (Adjacency[nT][ilocV][0] == -1)) {

					temp[iVert] = true;
					result[counter] = iVert; //store node of frontier
					++counter;
			}
		} //next iLocV
	} //next nT

	result.resize(counter);
	return result;
} // end of markVertFrontier

/*!Extract all available included patch ID's, present in the markInclusion member of the class 
 * \param[out] list of available ID's
 */
ivector1D SHAPE::getPatchIDIncluded(){

	ivector1D result;
	std::map<int, short int> IDmap;
	std::map<int, short int>::iterator IDit;

	for (int i=0; i<nSimplex; ++i)	IDmap[markInclusion[i]]=markInclusion[i];

	result.resize(IDmap.size());
	
	int counter=0;
	for (IDit=IDmap.begin(); IDit != IDmap.end(); ++IDit){
		result[counter]= IDit->second;
		++counter;
	}
	return result;
};

/*! Extract indices of triangles marked by a target PID Flag(if any)
 * \param[in] PIDflag target PIDflag
 * \param[out] result  list of extracted triangle indices 
 */
ivector1D SHAPE::extractPidded(int PIDflag){
  
	if(!checkVectorFind(pidType,PIDflag)){
		return ivector1D();
	}
   
	int size = pid_surface.size();
	ivector1D result(size);
	int counter=0;
	for(int i=0; i<size; i++){
		if(pid_surface[i] == (short int)PIDflag){
		result[counter] = i;
		++counter;
		}
	}

	result.resize(counter);
	return(result);
};

/*! Mark with a flag all triangle indices provided, and store them in the class member markInclusion
 * \param[in] map list of triangle indices that need to be marked
 * \param[in] flag integer number marking the current inclusion
 */ 
void SHAPE::markIncluded(ivector1D & map, int flag){
  
  int size = map.size();
  for(int i=0; i<size; ++i){
    markInclusion[map[i]] = (short int) flag;
  }
};

/*! Extract indices of triangles marked by a target flag (if any) in class member markInclusion
 * \param[in] flag target flag
 * \param[out] result  list of extracted triangle indices  
 */
ivector1D SHAPE::extractIncluded(int flag){
  
   short int target = (short int) flag;
   
   int check = posVectorFind(markInclusion,target); 
   if(check <0){return ivector1D();}
   
   ivector1D result(nSimplex);
   int counter=0;
   
   for(int i=0; i<nSimplex; i++){
      if(markInclusion[i] == target){
	result[counter] = i;
        ++counter;
      }
   }

   result.resize(counter);
   return(result); 
};


/*!TODO moving this to Refinement library 
 * Utility to return edges in common between a target Green Simplex 
 * having a single RED neighbour in its adjacency.
 * \param[in] Nt index on triangulation of the Green Simplex
 * \param[out] result returning number of edges in common 
 */
int SHAPE::marking_RED_edge(int Nt) { 
	ivector1D nN(3, 0);
	std::vector<bool> work(3, false);
	int sum;

	for (int j = 0; j < 3; j++) {
		nN[j] = Adjacency[Nt][j][0];
		if (markRefinement[nN[j]] == 2)
			work[j] = true;
	}

	sum = work[0] * 0 + work[1] * 1 + work[2] * 2;
	if ((work[0] + work[1] + work[2]) != 1) {
		std::cout
		<< "Warning on marking_RED_edge of class SHAPE, 2 marked red simplex at least neighbouring the green simplex:"
		<< Nt << '\n';
		exit(1);
	}

	return (sum);
}

//TODO is already present something similar in surftri. Please check and erase if so. 
ivector1D SHAPE::localVRing_1(int T, int j, bool & flag, bool & check) {
	ivector1D TList, VList;

	TList = Ring_1(T, j, flag, check);
	if(!check) return(VList);
	//TODO list need to write a fail safe Ring_1 in class surf tri

	std::map<int, int> appendVertex;

	for (int i = 0; i < TList.size(); i++) {
		int iTN = TList[i];
		for (int j = 0; j < 3; j++)
			appendVertex[Simplex[iTN][j]] =	Simplex[iTN][j];
	}

	VList.resize(appendVertex.size());
	// convert map to ivector list
	int counter = 0;
	for (std::map<int, int>::iterator iTapp = appendVertex.begin();
			iTapp != appendVertex.end(); iTapp++) {
		VList[counter] = (*iTapp).second;
		counter++;
	}
	return (VList);
}		// end of localVRing_1;


//******PRIVATE STUFFS***************************************************************

 /*! Utility for converting an integer number in a string of fixed 8 character size. 
    The method is referenced by SHAPE::modifyNas */ 
std::string SHAPE::resize_to_record8(int number) { 
	std::string work, readnum;
	std::stringstream ss;
	int size;

	work = "";
	ss << number;
	readnum = ss.str();
	trim(readnum);

	size = readnum.size();

	for (int i = 0; i < (8 - size); i++)
		work += " ";
	work += readnum;
	return (work);
};

/*! Convert BDF string of 24 characters of BDF Nastran format in an array of 3 doubles*/
darray3E SHAPE::NASconvertVertex(std::string & work) {
	
	darray3E point{0.0,0.0,0.0};
	std::string number, work3;

	std::string plus("+");
	std::string minus("-");
	std::size_t found_plus, found_minus;

	double value, base, esponente;
	int itype;

	for (int i = 0; i < 3; i++) {
		number = work.substr(i * 8, 8);
		number = trim(number);
		work3 = number.substr(1, number.size() - 1);

		found_minus = work3.find(minus);
		found_plus = work3.find(plus);

		std::stringstream ss1, ss2;

		if (found_minus != std::string::npos)
			itype = 0;
		if (found_plus != std::string::npos)
			itype = 1;
		if ((found_plus == std::string::npos)
				&& (found_minus == std::string::npos))
			itype = 2;
		switch (itype) {
		case 0:

			ss1.str(number.substr(0, found_minus));
			ss2.str(number.substr(found_minus + 1, 8 - found_minus - 1));
			ss1 >> base;
			ss2 >> esponente;
			value = base * pow(10.0, esponente);
			break;
		case 1:

			ss1.str(number.substr(0, found_plus));
			ss2.str(number.substr(found_plus + 1, 8 - found_plus - 1));
			ss1 >> base;
			ss2 >> esponente;
			value = base * pow(10.0, esponente);
			break;
		case 2:
			ss1.str(number);
			ss1 >> value;
			break;
		default: //do nothing
			break;
		} //end switch;

		point[i] = value;

	} //next i;
	return (point);
} //end of NASconvertVertex

/*! Utility to check if a *.stl file is written in ASCII or BINARY formats.
 * \param[in] input string containg absolute path to your *.stl file;
 * \param[out] result boolean FALSE for ASCII, TRUE for binary;
 * Unpredictable behaviours are not handled.  
*/
bool SHAPE::checkSTLfile(std::string & input){
  
  std::string line, out;
  std::ifstream read(input.c_str());
  if(read.is_open()){
	std::getline(read, line);
	stringstream ss;
	ss.str(trim(line));
	ss >> out;
  }
  else {std::cout<<"CANNOT OPEN FILE "<<input<<endl;}  
  read.close();
  
  return (out != "solid"); 
}

/*! check input list of triangle indices of the mother tessellation and epurate those indices 
 * not corresponding to a real simplex.
 * \param[in] Tlist list of triangle indices
 * \param[out] result cleaned list
 */
ivector1D SHAPE::checkTriangleList(ivector1D & TList){
	ivector1D result(TList.size());
	ivector1D::iterator it;
	int counter= 0 ;
	for(it=TList.begin(); it != TList.end(); ++it){
		if(*it >=0 && *it<nSimplex){
			result[counter] = *it;
			counter++;
		}
	}
	result.resize(counter);
	return(result);
};

/*! check input list of vertex indices of the mother tessellation and epurate those indices 
 * not corresponding to a real vertex.
 * \param[in] Vlist list of vertex indices
 * \param[out] result cleaned list
 */
ivector1D SHAPE::checkVertexList(ivector1D & VList){
	ivector1D result(VList.size());
	ivector1D::iterator it;
	int counter= 0 ;
	for(it=VList.begin(); it != VList.end(); ++it){
		if(*it >=0 && *it<nVertex){
			result[counter] = *it;
			counter++;
		}
	}
	result.resize(counter);
	return(result);
};

//#####################################################################################################
// NAMESPACE LS_UTILS_SHAPE IMPLEMENTATION
/*
 *	\date			31/12/2015
 *	\authors		Federico Gallizio
 * 	\authors 		Allen Woody
 *	\authors		Rocco Arpa
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This namespace version is released under .
 *
 *	\brief Utilities for Level Set calculation of Unstructured Triangulated Meshes on background Cartesian Mesh 
 *
 *	Collection of methods to generate Level Set informations on a Cartesian 3D Background Mesh of a 3D object
 *	passed as SHAPE container  .
 */

/*! Basic Constructor */
LS_UTILS_SHAPE::VoronoiBox::VoronoiBox(){
	demispan={0,0,0};
	origin={0,0,0};
}

/*! Put your cartesian mesh in a vector of cells of VoronoiBox structure.
 * \param[in] mesh pointer to Class_UcartMesh3D object
 * \param[out] result return a vector of VoronoiBox structs.
 */
std::vector<LS_UTILS_SHAPE::VoronoiBox> LS_UTILS_SHAPE::evaluateCartesianVoronoi(Class_UCartMesh3D * mesh)
{
  // i' m going to specialize finite volume mesh in a VoronoiBox structure.
    int dim = (mesh->nx)*(mesh->ny)*(mesh->nz);
    std::vector<VoronoiBox> result;
    result.resize(dim);

    darray3E dspan;
    dspan[0] = mesh->dx/2.0;
    dspan[1] = mesh->dy/2.0;
    dspan[2] = mesh->dz/2.0;

  for(int i=0; i<dim; i++)
  {
	 int l1,l2,l3;
	 l1 = i/(mesh->ny)/(mesh->nz);
	 l2 = (i - l1*(mesh->ny)*(mesh->nz)) / (mesh->nz);
	 l3 =  i - l1*(mesh->ny)*(mesh->nz) - l2*(mesh->nz);
	 result[i].demispan = dspan;
	 result[i].origin[0] = mesh->xnode[l1];
	 result[i].origin[1] = mesh->ynode[l2];
	 result[i].origin[2] = mesh->znode[l3];
  }

  return result;
} //end of evaluateCartesianVoronoi

/*! Generate a map of Mesh Node indices, reordering them according to Camilo's numeration.
 * \param[in] mesh pointer to Class_UcartMesh3D object
 * \param[out] result map[i] = j where i follow Camilo's indexing and j Class_UcartMesh3D's one.
 */
ivector1D LS_UTILS_SHAPE::generateNodeList(Class_UCartMesh3D * mesh)
{
  ivector1D totList;
  int i,j,k;
  int dim = (mesh->nx+1)*(mesh->ny+1)*(mesh->nz+1);

  totList.resize(dim);

  for (int ll=0; ll<dim; ll++)
  {
	i = ll / (mesh->nz+1)/(mesh->ny+1);
	j = (ll - i*(mesh->nz+1)*(mesh->ny+1))/(mesh->nz+1) ;
	k =  ll - i*(mesh->nz+1)*(mesh->ny+1) - j*(mesh->nz+1);

	totList[ll] = k*(mesh->nx+1)*(mesh->ny+1) + j*(mesh->nx+1) + i;
  }
	return totList;
} //end of generateNodeList

/*! Check if a point belongs to a VoronoiBox Cell (true). Return false if not.
 * \param[in] point 3D point-vertex 
 * \param[in] vbox  pointer to VoronoiBox cell
 * \param[out] result boolean flag 
 */
bool LS_UTILS_SHAPE::checkPointVBox(darray3E & point, VoronoiBox * vbox)
{
	bool check;
	int checkIntegral;
	darray3E dumP{0,0,0};
	dumP = point - vbox->origin;
	double tolerance = 0.05; // 5% canonical on lattice stuff bbox

	checkIntegral =  (int)(abs(dumP[0])< (1.0+tolerance) * vbox->demispan[0]) *
					 (int)(abs(dumP[1])< (1.0+tolerance) * vbox->demispan[1]) *
					 (int)(abs(dumP[2])< (1.0+tolerance) * vbox->demispan[2]) ;


	check = (bool)checkIntegral;
	return check;
}//end of checkPointVBox

/*! Generate Signed Level Set field of a given 3D surface on a given Background mesh(Internal-negative value of level set).
 *  Mesh limits are automatically calculated by using the bounding box of the original geometry. The method is intended for  
 *  closed surfaces only. Please refer to GenerateLevelSetOPEN for open surface version.
 * \param[in] tri 3D surface pointer
 * \param[in] mesh background cartesian mesh pointer
 * \param[in] scale scaling factor >=1 for surface bounding box
 * \param[in] n_tot total number of cells in the background mesh
 * \param[out] result signed level set field 
 */
dvector1D LS_UTILS_SHAPE::generateLevelSet(SHAPE* tri, Class_UCartMesh3D * mesh,double scale, int n_tot )
{
	std::cout<<"entered in LS closed surface"<<endl;
	dvector1D result;
	dvector1D x_ext(2,0),y_ext(2,0),z_ext(2,0);
	int nx,ny,nz;

	if(tri->Adjacency.size()==0){tri->BuildAdjacency();}
	if(CG_PLSurf::IsOpen(tri->Adjacency))
		{
		std::cout<<"WARNING! Open geometry/Superficial flaws on tessellation detected."<<endl;
		std::cout<<"Error may occurs in level set computation..."<<endl;
		return result;
		}

	//calculate bounding box of the 3D-object
	tri->BoundingBox(x_ext,y_ext,z_ext);
	// rescale your grid limits
	{
	 double span, delta;

         span = fmax(1, scale) * (x_ext[1] - x_ext[0]);
	 delta = 0.5*(span - (x_ext[1] - x_ext[0]));
	 x_ext[0] = x_ext[0] - delta;
	 x_ext[1] = x_ext[1] + delta;

	 span = fmax(1, scale) * (y_ext[1] - y_ext[0]);
	 delta = 0.5*(span - (y_ext[1] - y_ext[0]));
	 y_ext[0] = y_ext[0] - delta;
	 y_ext[1] = y_ext[1] + delta;

	 span = fmax(1, scale) * (z_ext[1] - z_ext[0]);
	 delta = 0.5*(span - (z_ext[1] - z_ext[0]));
	 
	 z_ext[0] = z_ext[0] - delta;
	 z_ext[1] = z_ext[1] + delta;

	}
	//compute the indicative spacing lenght for a uniform cartesian mesh to be created w/ the n_tot point available
	double delta = (x_ext[1] - x_ext[0]) * (y_ext[1] - y_ext[0]) * 	(z_ext[1] - z_ext[0]) / n_tot;
	delta = pow(delta,1.0/3.0);

	std::cout<<"Grid stepsize :   "<<delta<<std::endl;
	//set the grid dimension
	nx = (int)((x_ext[1] - x_ext[0])/delta + 0.5);
	ny = (int)((y_ext[1] - y_ext[0])/delta + 0.5);
	nz = (int)((z_ext[1] - z_ext[0])/delta + 0.5);

	std::cout<<"Background Grid real dimensions : "<<nx<<'\t'<<ny<<'\t'<<nz<<'\t'<<std::endl;

	//initialize the Cartesian mesh:
	mesh->SetMesh(x_ext,y_ext,z_ext, nx, ny,nz);
	
	Class_LevelSet_Stl<Class_UCartMesh3D> levSet(mesh, tri);
	levSet.setSdfPropagation(true);
	
	levSet.sdfNarrowBand();
	result = levSet.getSdf();
	
	for(int k=0; k<result.size(); ++k){
		if(std::isnan(result[k])){std::cout<<"find nan in "<<k<<"  of level set"<<endl;}
	}
return(result);
};

/*! Generate Level Set field of a given 3D surface on a given Background mesh. Mesh limits
 *  are automatically calculated by using the bounding box of the original geometry. The method is intended for  
 *  opened surfaces only. If a closed surface is selected internal/external distances are not distinguished by sign.
 * \param[in] tri 3D surface pointer
 * \param[in] mesh background cartesian mesh pointer
 * \param[in] scale scaling factor >=1  to enlarge surface bounding box
 * \param[in] n_tot total number of cells in the background mesh
 * \param[out] result level set field 
 */
dvector1D LS_UTILS_SHAPE::generateLevelSet_OPEN(SHAPE* tri, Class_UCartMesh3D * mesh, double distance, int n_tot )
{
	std::cout<<"entered in LS open surface"<<endl;
	dvector1D result;
	dvector1D x_ext(2,0),y_ext(2,0),z_ext(2,0);
	int nx,ny,nz;

	//calculate bounding box of the 3D-object
	tri->BoundingBox(x_ext,y_ext,z_ext);
	// rescale your grid limits
	{
	 x_ext[0] = x_ext[0] - 2*distance;
	 x_ext[1] = x_ext[1] + 2*distance;

	 y_ext[0] = y_ext[0] - 2*distance;
	 y_ext[1] = y_ext[1] + 2*distance;

	 z_ext[0] = z_ext[0] - 2*distance;
	 z_ext[1] = z_ext[1] + 2*distance;
	}
// 	{
// 	 double span;
// 
//          span = fmax(1, (scale - 1.0)) * (x_ext[1] - x_ext[0]);
// 	 x_ext[0] = x_ext[0] - 0.5*span;
// 	 x_ext[1] = x_ext[1] + 0.5*span;
// 
//          span = fmax(1, (scale - 1.0)) * (y_ext[1] - y_ext[0]);
// 	 y_ext[0] = y_ext[0] - 0.5*span;
// 	 y_ext[1] = y_ext[1] + 0.5*span;
// 
//          span = fmax(1, (scale - 1.0)) * (z_ext[1] - z_ext[0]);
// 	 z_ext[0] = z_ext[0] - 0.5*span;
// 	 z_ext[1] = z_ext[1] + 0.5*span;
// 
// 	}

	//compute the indicative spacing lenght for a uniform cartesian mesh to be created w/ the n_tot point available
	double delta = (x_ext[1] - x_ext[0]) * (y_ext[1] - y_ext[0]) * 	(z_ext[1] - z_ext[0]) / n_tot;
	delta = pow(delta,1.0/3.0);

	std::cout<<"Grid stepsize :   "<<delta<<std::endl;

	//if(delta>distance) std::cout<<"WARNING! Grid stepsized is greater than the control distance: this may impar the quality of result." << endl;

	//set the grid dimension
	nx = (int)((x_ext[1] - x_ext[0])/delta + 0.5);
	ny = (int)((y_ext[1] - y_ext[0])/delta + 0.5);
	nz = (int)((z_ext[1] - z_ext[0])/delta + 0.5);


	std::cout<<"Background Grid real dimensions : "<<nx<<'\t'<<ny<<'\t'<<nz<<'\t'<<std::endl;

	//initialize the Cartesian mesh:
	mesh->SetMesh(x_ext,y_ext,z_ext, nx,ny,nz);


	Class_LevelSet_Stl<Class_UCartMesh3D> levSet(mesh, tri);
	levSet.setSdfPropagation(true);
	levSet.dfNarrowBand();

	result = levSet.getSdf();

	for(int k=0; k<result.size(); ++k)
		{
			if(std::isnan(result[k])){std::cout<<"find nan in "<<k<<"  of level set"<<endl;}
		}

	return(result);
};

/*! Plot A Cartesian Grid on VTK-VTR file w/ an optional scalar field attached. 
 * \param[in] dir path to output directory
 * \param[in] name output file name, without tag
 * \param[in] cod  codex for writing, "ascii" or "appended" (binary)
 * \param[in] openSurface flag true if Level Set is referred to an open surface, false otherwise
 * \param[in] mesh cartesian mesh
 * \param[in] field scalar field to be plotted, OPTIONAL
 */
	
void LS_UTILS_SHAPE::plotScalarVTR(std::string dir, std::string name,
							 bool cod,
							 bool openSurface, 
							 Class_UCartMesh3D & mesh,
							 dvector1D * field)
{
	std::string codex="ascii";
	if(cod){codex="appended";}
	VTK_LSet handle_out(dir,name,codex);
	handle_out.linkData(mesh);
	dvector1D result;
	if(field != NULL){
		std::string loc;
		if(field->size() == (mesh.nx+1)*(mesh.ny+1)*(mesh.nz+1)) {loc="Point";}
		if(field->size() == (mesh.nx)*(mesh.ny)*(mesh.nz)) {loc="Cell";}

		if(loc.empty()){std::cout<<"ScalarField not compatible w/ current mesh"<<endl; return;}
	
 		if(openSurface){
		//invert the ordering of the mesh
		int sizeR = field->size();
		result.resize(sizeR,0.0);
			for(int i=0; i<sizeR; ++i){
				ivector1D coord(3,0);
				int I_;
				mesh.AccessPointData(i, coord[0],coord[1],coord[2]);
				I_ = (mesh.nx+1)*(mesh.ny+1)*coord[2] + (mesh.nx+1)*coord[1] + coord[0];
				result[I_] = (*field)[i];
			}
			handle_out.AddData("scalarField",1,"Float64", loc, codex);
			handle_out.linkScalarField(result);
		}
		else{
		handle_out.AddData("scalarField",1,"Float64", loc, codex);
		handle_out.linkScalarField(*field);
 		}
	}

	handle_out.Write();
	return;
};

/*! Plot A Cartesian Grid on VTK-VTR file w/ an optional scalar field attached. 
 * \param[in] dir path to output directory
 * \param[in] name output file name, without tag
 * \param[in] cod  codex for writing, "ascii" or "appended" (binary)
 * \param[in] openSurface flag true if Level Set is referred to an open surface, false otherwise
 * \param[in] mesh cartesian mesh
 * \param[in] field v field to be plotted, OPTIONAL
 */
void LS_UTILS_SHAPE::plotVectorVTR(std::string dir, std::string name,
							 bool cod,
							 bool openSurface, 
							 Class_UCartMesh3D & mesh,
							 dvecarr3E * field)
{
	std::string codex="ascii";
	if(cod){codex="appended";}
	VTK_LSet handle_out(dir,name,codex);
	handle_out.linkData(mesh);
	dvecarr3E result;
	if(field != NULL){
		std::string loc;
		if(field->size() == (mesh.nx+1)*(mesh.ny+1)*(mesh.nz+1)) {loc="Point";}
		if(field->size() == (mesh.nx)*(mesh.ny)*(mesh.nz)) {loc="Cell";}

		if(loc.empty()){std::cout<<"VectorField not compatible w/ current mesh"<<endl; return;}
 		if(openSurface){
		//invert the ordering of the mesh
		int sizeR = (*field).size();
		result.resize(sizeR,darray3E{0,0,0});

			for(int i=0; i<sizeR; ++i){
				ivector1D coord(3,0);
				int I_;
				mesh.AccessPointData(i, coord[0],coord[1],coord[2]);

				I_ = (mesh.nx+1)*(mesh.ny+1)*coord[2] + (mesh.nx+1)*coord[1] + coord[0];

				result[I_] = (*field)[i];
			}
			handle_out.AddData("vectorField",3,"Float64", loc, codex);
			handle_out.linkVectorField(result);
		}
		else{
			handle_out.AddData("vectorField",3,"Float64", loc, codex);
			handle_out.linkVectorField(*field);
 		}
	}

	handle_out.Write();
	return;
};

//#####################################################################################################
// NAMESPACE PLOT_DISTANCE_SHAPE IMPLEMENTATION
/*
 *	\date			31/12/2015
 *	\authors		Federico Gallizio
 * 	\authors 		Allen Woody
 *	\authors		Rocco Arpa
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This namespace version is released under .
 *
 *	\brief Utilities for managing Distance evaluation on vertices of an unstructured mesh 
 *
 *	Collection of methods to post-process level set information and plot them as a "Distance" scalar
 *	field wrt to vertices of a SHAPE Data Structure  .
 */

/*! Generate a list of triangulated mesh points belonging/relative to a 3D surface, given an ordered list of their indices.
 * \param[in] tri 3D surface structure
 * \param[in] map list of point indices
 * \param[in] vertex external list of vertices, OPTIONAL
 * \param[out] result selected points
 */
dvecarr3E PLOT_DISTANCE_SHAPE::generateListPoints(SHAPE & tri, ivector1D & map, dvecarr3E * vertex)
{
	 dvecarr3E result;
	 result.resize(map.size(), darray3E{0,0,0});
	 int itype = 1;

	 if(vertex != NULL && vertex->size() == tri.nVertex) {itype = 0;}

	 switch(itype) {

	 case 0: for(int k=0; k<map.size(); ++k){result[k] = (*vertex)[map[k]];}
		 break;
	 default:
		 for(int k=0; k<map.size(); ++k){result[k] = tri.Vertex[map[k]];}
		 break;
	}
	 return(result);
};

/*! Interpolate values of a given scalar field on a cartesian mesh, for each point provided by an external list
 * \param[in] mesh reference cartesian mesh
 * \param[in] field scalar field defined on cartesian mesh 
 * \param[in] vertex list of points where field need to be interpolated
 * \param[out] result return interpolated values
 */
dvector1D PLOT_DISTANCE_SHAPE::interpolateCartMeshField(Class_UCartMesh3D & mesh, dvector1D & field, dvecarr3E & vertex)
{
 dvector1D result(vertex.size(), 0);

 for(int k=0; k<result.size(); ++k){
	dvector1D target = conVect(vertex[k]);
	mesh.interpolatePointData(target,field,result[k]);
 }//next k;

 return(result);
};

/*! Plot list of vertices and their associated scalar field in xyz format 
 * \param[in] dir output directory path
 * \param[in] name output filename 
 * \param[in] vertex list of vertices
 * \param[out] field scalar field associated to the cloud points vertex
 */
void  PLOT_DISTANCE_SHAPE::plotOnFile(std::string dir, std::string name, dvecarr3E & vertex, dvector1D & field )
{
  std::string path = dir+"/"+name;
  std::ofstream out(path.c_str());
  if(!out.is_open()){
   std::cout << "The file " << path << " can't be opened or you don't have the write permissions in the folder.";
   exit(1);
  }
  for(int K=0; K<vertex.size(); ++K){
	for(int j=0;j<vertex[K].size(); ++j){
		out<<scientific<<vertex[K][j]<<'\t';
		}
	out<<scientific<<field[K]<<std::endl;
 }//next K
 out.close();
};

	
	
















