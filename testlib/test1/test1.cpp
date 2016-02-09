/* test 001 -> class SHAPE testing
 * Testing all methods of the General Unstructured Mesh manager
 */

#include "lib_STL_IO.hpp"
#include <iostream>
#include <fstream>


svector1D decodeMainInput(int n_info, char* info[]){
	
	svector1D readS(n_info -1);
	svector1D result;
	for(int i=1; i<n_info; ++i ){
		readS[i-1] = info[i];
		readS[i-1] = trim(readS[i-1]);
	}
	
	//std::cout<<readS<<std::endl;
	std::string keypath, keytype, keyhelp;
	keypath ="-p";
	keytype ="-t";
	keyhelp ="-h";
	
	int posP = posVectorFind(readS, keypath);
	int posT = posVectorFind(readS, keytype);
	int posH = posVectorFind(readS, keyhelp);
	
	bool check = ( (2*(posP != -1)+2*(posT != -1) + (posH !=-1) ) == readS.size());
	
	//std::cout<<posP<<'\t'<<posT<<'\t'<<posH<<'\t'<<check<<std::endl;
	
	check = check && ((posP+1 != posH) && (posP+1 != posT) && 
			(posT+1 != posH) && (posT+1 != posP)); 
	
	//std::cout<<check<<std::endl;
	
	if(posH == -1 && posP != -1 && posT != -1 && check){
		result.push_back(readS[posP+1]);
		result.push_back(readS[posT+1]);
	}else{
		std::cout<<"=================================================="<<std::endl;
		std::cout<<"test1 help				  "<<std::endl;
		std::cout<<"==================================================="<<std::endl;
		std::cout<<""<<std::endl;
		std::cout<<"-h:		print this help		   "<<std::endl;
		std::cout<<""<<std::endl;
		std::cout<<"-p <path>:	define abs path to your geometry "<<std::endl;
		std::cout<<""<<std::endl;
		std::cout<<"-t <int>:	define your current test, 0-Shape, 1-LevelSetUtils, 2-PlotDistance "<<std::endl;
		std::cout<<""<<std::endl;
		std::cout<<"=================================================="<<std::endl;
		std::cout<<"WARNING				  "<<std::endl;
		std::cout<<"==================================================="<<std::endl;
		std::cout<<"To run tests you need to specify a valid path and valid type of test				  "<<std::endl;
		std::cout<<"==================================================="<<std::endl;
		exit(0);
	}
	return(result);
}

void testShape(std::string &filename, std::ofstream & out ){
	
	//CLASS SHAPE testing I/O
	SHAPE geo;
	geo.init(filename);
	out<<"Loading Geometry is ...OK"<<endl;
		
	geo.cleaning();
	out<<"Cleaning Geometry is ...OK"<<endl;
	
	{	
		//testing markIncluded/extractPidded/extractIncluded
		ivector1D map;
		if(geo.filetype != 1 && geo.pidType.size() == 0){
			int target = (int)geo.nSimplex/10;
			map.resize(target,1);
			for(int i=1; i<target; ++i){map[i] = map[i-1]+1;}
		}else{
			map = geo.extractPidded(geo.pidType[0]);
		}
		
		bool check = true; 
		geo.markIncluded(map, 8);
		
		ivector1D inclPD = geo.getPatchIDIncluded(); 
		
		int mark = 8;
		check = check && (checkVectorFind(inclPD,mark));
		ivector1D result = geo.extractIncluded(8);
		check = check && (map == result);
		if(!check){
			out<<"Inclusion/Extraction methods...FAILED"<<endl;
		}else{
			out<<"Inclusion/Extraction methods...OK"<<endl;
		}
	}

	{
		//plotting
		geo.plotVTU(".", "VTUconversion", true, &geo.Vertex);
		geo.plotSTL("./asciiSTLconversion.stl", geo.Vertex, false);
		geo.plotSTL("./binarySTLconversion.stl", geo.Vertex, true);
		geo.plotDGF("./DGFconversion.dgf", geo.Vertex);
		if(geo.filetype ==1){
			geo.modifyNAS(geo.filename, "./NASreplot.nas", geo.Vertex);
		}
		//create a fake scalar and vector fields
		dvector1D sfield(geo.nSimplex, 4);
		dvecarr3E vfield(geo.nVertex, darray3E{5,5,5});
		
		geo.plotScalarVTU(".", "VTUScalarField", true, 4, sfield, &geo.Vertex);
		geo.plotVectorVTU(".", "VTUVectorField", true, 5, vfield, &geo.Vertex);
		
		out<<"Plotting DONE...check Output"<<endl;
	}
	
	return;
}

void testLS(std::string &filename, std::ofstream & out){
	
	Class_UCartMesh3D * mesh = new Class_UCartMesh3D;
	SHAPE geo;
	geo.init(filename);
	geo.cleaning();
	
	dvector1D LS_Closed = LS_UTILS_SHAPE::generateLevelSet(&geo, mesh, 1.2, 100000);
	
	if(LS_Closed.size() !=0){
		LS_UTILS_SHAPE::plotScalarVTR(".", "signedLS", true, false, *mesh, &LS_Closed);
		out<<"signed level set...OK"<<endl;
	}else{
		out<<"signed level set...NOT AVAILABLE FOR THIS GEOMETRY"<<endl;
	}
	dvector1D LS_Open = LS_UTILS_SHAPE::generateLevelSet_OPEN(&geo, mesh, 0.1, 100000);
	LS_UTILS_SHAPE::plotScalarVTR(".", "unsignedLS", true, true, *mesh, &LS_Open);	
	
	out<<"unsigned level set...OK"<<endl;
	out<<"Level Set test DONE...check Output"<<endl;
	delete mesh; mesh==NULL;
}

void testPlotDist(std::string &filename, std::ofstream & out){
	
	Class_UCartMesh3D * mesh = new Class_UCartMesh3D;
	SHAPE geo; 
	geo.init(filename);
	geo.cleaning();
	ivector1D indexList(geo.nVertex, 1);
	
	for(int i=1; i<geo.nVertex; ++i){
		indexList[i] = indexList[i-1]+1;
	}
	
	dvector1D LS_Closed = LS_UTILS_SHAPE::generateLevelSet(&geo, mesh, 1.2, 100000);
	
	if(LS_Closed.size() !=0){
		LS_UTILS_SHAPE::plotScalarVTR(".", "signedLS", true, false, *mesh, &LS_Closed);
		out<<"signed level set...OK"<<endl;
		dvecarr3E points = PLOT_DISTANCE_SHAPE::generateListPoints(geo, indexList, &geo.Vertex);
		dvector1D distance = PLOT_DISTANCE_SHAPE::interpolateCartMeshField(*mesh, LS_Closed, points);
		PLOT_DISTANCE_SHAPE::plotOnFile(".", "plotDistanceTest.dat", points, distance);
		out<<"Plot Distance test DONE...check Output"<<endl;
	}else{
		
		out<<"Plot Distance test FAILED...error occurred or open surface Detected"<<endl;
		out<<"Please Note PLOT_DISTANCE_SHAPE works only with closed surfaces"<<endl;
	}
	return;
	
	delete mesh; mesh==NULL;
}


int main(int argc, char*argv[])
{
	svector1D infoShape =  decodeMainInput(argc, argv);

	int valType;
	{
		std::stringstream ss;
		ss<<infoShape[1];
		ss>>valType;
	}

	std::ofstream out;
	out.open("testingSHAPE.log");
	if(out.is_open()){

		switch(valType){
			case 1 :testLS(infoShape[0], out);
				break;
				
			case 2 :testPlotDist(infoShape[0], out);
				break;
				
			default:testShape(infoShape[0], out);
				break;
		}
			
	}else{
		cout<<"Not able to open log file"<<endl;
	}
	
	out.close();
	return(0);
} //end program 
