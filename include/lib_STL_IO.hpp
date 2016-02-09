#ifndef EL_MC_STL_IO
#define EL_MC_STL_IO

// local libs
#include "Operators.hpp"
#include "customOperators.hpp"
#include "STL_IOFunct.hpp"
#include "DGF_IOFunct.hpp"
#include "libRA_VTKInterfaces.hpp"
#include "Class_UCartMesh.hpp"
#include "Class_LevelSet_Stl.hpp"


/*!
 *	\date			31/12/2015
 *	\authors		Federico Gallizio
 * 	\authors 		Allen Woody
 *	\authors		Rocco Arpa
 *	\copyright		Copyright 2015 Optimad engineering srl. All rights reserved.
 *	\par			License:\n
 *	This class version is released under .
 *
 *	\brief Data Structure Manager Class for Unstructured Triangulated Meshes 
 *
 *	SHAPE is a Data Structure class for Unstructured Meshes, directly derived from Class_SurfTri BitPit object of libMesh.
 *	I/O can extend to Standard Tesselation Language,  *.stl(ASCII or bynary), BDF NASTRAN native format *.nas, *.dgf and VTK 
 *     	*.vtu formats. Retains tool for checking/cleaning the meshes and handling refinement.
 */
class SHAPE:public Class_SurfTri{

		public:
                std::string filename;		/**< std::string. File Name of the original geometry*/
                short int filetype;		/**< Integer. Flag to identify the type of geometry file.
                							stl 0 - NAS 1 - dgf 2 - vtu 3 */
                ivector2D InverseSimplexMap; 	/**< Vertex To Simplex Connectivity */

                dvecarr3E deform_OR;		/**< Accessory mesh data structure, interfacing with Deformation class/functions.
						  Coordinates of each deformed/displaced triangulation vertex expressed 
								  in the unique CAD reference system. */ 
              
                shivector1D   markRefinement;      /**< Accessory mesh data structure, interfacing with Refinement class/methods.
							     It is a list of all simplicies in the tasselation. Such simplicies are marked 
							     with an integer flag necessary to guide the local Refinement. Admissible flags are:

		                                                         -    0: no refinement on simplex

									 -    1: Green refinement on simplex

									 -    2: Red refinement on simplex.

								see class REFINEMENT for further details. */
               shivector1D   markInclusion;	/**< Accessory mesh data structure, interfacing with Selection class/methods, Lattice  
							    and Patch class. It is a list of all simplicies in the tasselation, marked with 
							    an integer flag that corresponds to the "name" of a Lattice. See LATTICE<*,*>.
							    Each sub-group of simplicies marked with the name of the same lattice constitues 
							    the deformable support of that lattice. Each defined Lattice has its own 
							    deformable support contained in markInclusion and called with its own name.	*/
		
		shivector1D  pid_surface;	/**< Accessory meturnesh data structure, interfacing with Selection class/methods.
						   It is a list of all simplicies in the tasselation, marked with the PID integer flag 
						   readed from the geometry file in *.nas format. It is void for all other formats.*/ 

		dvector1D moll_filter;		/**< Structure containing the deformability field of the original geometry. It resized
						    and filled in the method FFD_Ext::deformation*/

		ivector1D pidType; /**< save the total number of different PID found in the structure */
		bool cleanData; /**< cleaning flag for tracking cleaning stuffs*/

      SHAPE();
      SHAPE(const SHAPE &);
      ~SHAPE();
      void freeSHAPE();
      
     SHAPE & operator=(const SHAPE &);
     void init(std::string input_file);
     void cleaning();
     
     // structural methods
     void Import_nas(std::string file);
     void setInverseConnectivityMap(); 
     
    //plotting
    void plotVTU(std::string dir, std::string name, bool flag, dvecarr3E * vertices=NULL);
    void plotScalarVTU(std::string dir, std::string name, bool flag, int counter, dvector1D &, dvecarr3E * vertices=NULL);
    void plotVectorVTU(std::string dir, std::string name, bool flag, int counter, dvecarr3E &, dvecarr3E * vertices=NULL);
    void plotSTL(string filename,  dvecarr3E & vertices, bool flag); 
    void plotDGF(string filename,  dvecarr3E & vertices);
    void modifyNAS (std::string filename, std::string outfilename, dvecarr3E & vertex);
    
    
    void plotExtraction(std::string dir, std::string name, bool flag, int counter, ivector1D & Tlist, dvecarr3E * vertices=NULL);
    void plotCloudExtraction(std::string dir, std::string name, bool flag, int counter, ivector1D & Vlist, dvecarr3E * vertices=NULL);
    
    //utilities	
    ivector1D getPatchIDIncluded();
    ivector1D extractPidded(int PIDflag);
    void markIncluded(ivector1D & map, int flag);
    ivector1D extractIncluded(int flag);
    
    //service utilities
    ivector2D createPatchesGraph(ivector1D &, fvector1D &, bool &);
    ivector1D extractFrontier();
     
    //TODO moving to Refinement library
    int marking_RED_edge(int Nt);
    //TODO non so se mi serve davvero, con la derivazione in public di surfTri LO USO SOLO IN createPatchesGraph.
    ivector1D localVRing_1(int , int, bool &, bool &);

     private:
        std::string resize_to_record8(int number);
	darray3E NASconvertVertex(std::string &) ;
	bool checkSTLfile(std::string &);
	ivector1D checkVertexList(ivector1D &);
	ivector1D checkTriangleList(ivector1D &);

}; // end of Class STL_I/O


/*!
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
namespace LS_UTILS_SHAPE{
  
  struct VoronoiBox {
	darray3E origin;
	darray3E demispan;
	//constructor
	VoronoiBox();
  };	

	std::vector<VoronoiBox> evaluateCartesianVoronoi(Class_UCartMesh3D * );
	ivector1D    generateNodeList(Class_UCartMesh3D *);
	
	bool checkPointVBox(darray3E &, VoronoiBox *);

	dvector1D generateLevelSet(SHAPE*, Class_UCartMesh3D *,double, int);
	dvector1D generateLevelSet_OPEN(SHAPE*, Class_UCartMesh3D *,double, int);

	void plotScalarVTR(std::string dir, std::string name, bool cod, bool openSurface, Class_UCartMesh3D & mesh, dvector1D * field=NULL);
	void plotVectorVTR(std::string dir, std::string name, bool cod, bool openSurface,Class_UCartMesh3D & mesh, dvecarr3E * field=NULL);
};

/*!
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
namespace PLOT_DISTANCE_SHAPE{

	dvecarr3E generateListPoints(SHAPE &, ivector1D &, dvecarr3E * );
	dvector1D interpolateCartMeshField(Class_UCartMesh3D &, dvector1D &, dvecarr3E &);
	void 	  plotOnFile(std::string, std::string, dvecarr3E &, dvector1D &);
};

#endif
