/* test 006 -> Support Mesh Classes methods of libRA_SupportMesh
 * Testing all methods of all classes defining a class SupportMesh
 */

 using namespace std;
#include <iostream>
#include <fstream> 
#include "customOperators.hpp" 
#include "libRA_CorePatches.hpp"
#include "libRA_SupportMesh.hpp" 

int main(int argc, char*argv[]){
	
std::ofstream out;
out.open("testingSupportMesh.log");
if(out.is_open()){
	
	
	SHAPE sphere;
	sphere.init("./sphere.nas");
	sphere.cleaning();
	int pidsphere;
	if(sphere.pidType.size() != 0){
		pidsphere = sphere.pidType[0];
	}
	
	SHAPE plate;
	plate.init("./plate.stl");
	plate.cleaning();
	dvector1D xlim(2,0),ylim(2,0), zlim(2,0);
	plate.BoundingBox(xlim,ylim,zlim);
	
	std::string folder = ".";
	ivector1D tList;
	ivector1D vList;
	dvecarr3E checkNormals;
	dvecarr3E checkVNormals;
	
	{
		darray3E or_;
		or_[0] = (xlim[0] + xlim[1])*0.5;
		or_[1] = (ylim[0] + ylim[1])*0.5;
		or_[2] = (zlim[0] + zlim[1])*0.5;
		
		double radius = std::fmax((xlim[1] - xlim[0]), (ylim[1] - ylim[0]));
		radius = 0.25 *std::fmax(radius, (zlim[1] - zlim[0])); 
		
		PATCH_SPHERE * select = new PATCH_SPHERE(or_, radius); 
		tList = select->includeTriangulation(&plate);
		delete select; select=NULL;
	}
	
	{
		std::map<int,int> map;
		std::map<int,int>::iterator it;
		
		for(int i=0; i<tList.size(); ++i){
			for(int j=0; j<plate.Simplex[tList[i]].size(); ++j){
				map[plate.Simplex[tList[i]][j]] = plate.Simplex[tList[i]][j];
			}
		}
		
		int counter=0;
		vList.resize(map.size());
		for(it=map.begin(); it !=map.end(); ++it){
			vList[counter] = it->second;
			++counter;
		}
	}
	
	checkNormals.resize(tList.size());	
	for(int i=0; i<tList.size(); ++i){
		checkNormals[i] = plate.Normal[tList[i]];	
	}
	
	checkVNormals.resize(vList.size());	
	for(int i=0; i<vList.size(); ++i){
		checkVNormals[i] = plate.VNormal[vList[i]];	
	}
	

// 	//testing BASE_Support methods.
// 	{
// 		out<<"Testing BASE_Support methods..."<<endl;
// 		
// 		BASE_Support * support = new BASE_Support();
// 		{
// 			support->setSupport(&sphere, pidsphere);
// 			support->plotVTU(folder, "supportWshape&PID", true);
// 			out<<"...test setSupport method w/ SHAPE and PID...see output file"<<endl;
// 		}
// 		
// 		{
// 			support->cleanSupport();
// 			bool check = (support->getMother()== NULL);
// 			check = check &&((support->getSimplexMap()).size() == 0);
// 			check = check &&((support->getVertexMap()).size() == 0);
// 			out<<"...test cleanSupport method ...";
// 			if(check){
// 				out<<"OK"<<endl;
// 			}else{
// 				out<<"FAILED"<<endl;		
// 			}
// 		}
// 		{
// 			support->setSupport(&plate, tList);
// 			support->plotVTU(folder, "supportWithExtList", true);
// 			out<<"...test setSupport method w/ ext List...see output file"<<endl;
// 		}
// 		{
// 			support->setName("trialSupport");
// 			bool check = (support->getName() == "trialSupport");
// 			out<<"...test set/getName methods ...";
// 			if(check){
// 				out<<"OK"<<endl;
// 			}else{
// 				out<<"FAILED"<<endl;		
// 			}
// 		}
// 		{
// 			std::string prova = "BASE_Support";
// 			bool check = (support->getClassType()==prova );
// 			out<<"...test getClassType method ...";
// 			if(check){
// 				out<<"OK"<<endl;
// 			}else{
// 				out<<"FAILED"<<endl;		
// 			}
// 			
// 		}
// 		{
// 			bool check = (support->getMother()== &plate);
// 			out<<"...test getMother method ...";
// 			if(check){
// 				out<<"OK"<<endl;
// 			}else{
// 				out<<"FAILED"<<endl;		
// 			}
// 		}
// 		{
// 			bool check = (support->getSimplexMap()== tList);
// 			out<<"...test getSimplexMap method ...";
// 			if(check){
// 				out<<"OK"<<endl;
// 			}else{
// 				out<<"FAILED"<<endl;		
// 			}
// 		}
// 		{
// 			bool check = (support->getVertexMap()== vList);
// 			out<<"...test getVertexMap method ...";
// 			if(check){
// 				out<<"OK"<<endl;
// 			}else{
// 				out<<"FAILED"<<endl;		
// 			}
// 		}
// 		{
// 			dvecarr3E cN = support->getNormals();
// 			bool check = true;
// 			for(int i=0; i<cN.size(); ++i){
// 				check = check &&(cN[i]== checkNormals[i]);
// 			}	
// 			out<<"...test getNormals method ...";
// 			if(check){
// 				out<<"OK"<<endl;
// 			}else{
// 				out<<"FAILED"<<endl;		
// 			}
// 		}
// 		{
// 			dvecarr3E cVN = support->getVNormals();
// 			bool check = true;
// 			for(int i=0; i<cVN.size(); ++i){
// 				check = check &&(cVN[i]== checkVNormals[i]);
// 			}	
// 			out<<"...test getVNormals method ...";
// 			if(check){
// 				out<<"OK"<<endl;
// 			}else{
// 				out<<"FAILED"<<endl;		
// 			}
// 		}
// 		{
// 			support->plotWholeBoundary(folder, "plotWholeBoundary",true, "Cloud");
// 			out<<"...test plotWholeBoundary...see output file"<<endl;;
// 		}
// 		{
// 			support->plotFreeBoundary(folder, "plotFreeBoundary",true, "3DCurve");
// 			out<<"...test plotFreeBoundary...see output file"<<endl;
// 		}
// 		{
// 			support->plotConstrBoundary(folder, "plotConstrBoundary",true, "3DCurve");
// 			out<<"...test plotConstrBoundary...see output file"<<endl;
// 		}
// 		
// 		delete support; support = NULL;
// 	}	
// 	
	
	//testing SupportsplitB methods.
	{
 		out<<"Testing SupportSplitB methods..."<<endl;
 		
 		SupportSplitB * support = new SupportSplitB();
// 		{
// 			support->setSupport(&sphere, pidsphere);
// 			support->plotVTU(folder, "supportSplitB_Wshape&PID", true);
// 			out<<"...test setSupport method w/ SHAPE and PID...see output file"<<endl;
// 		}
// 		{
// 			support->cleanSupport();
// 			bool check = (support->getMother()== NULL);
// 			check = check &&((support->getSimplexMap()).size() == 0);
// 			check = check &&((support->getVertexMap()).size() == 0);
// 			out<<"...test cleanSupport method ...";
// 			if(check){
// 				out<<"OK"<<endl;
// 			}else{
// 				out<<"FAILED"<<endl;		
// 			}
// 		}
// 		{
// 			support->setSupport(&plate, tList);
// 			support->plotVTU(folder, "supportSplitBWithExtList", true);
// 			out<<"...test setSupport method w/ ext List...see output file"<<endl;
// 		}
// 		{
// 			std::string prova = "SupportSplitB";
// 			bool check = (support->getClassType()==prova );
// 			out<<"...test getClassType method ...";
// 			if(check){
// 				out<<"OK"<<endl;
// 			}else{
// 				out<<"FAILED"<<endl;		
// 			}
// 			
// 		}
		
// 		support->cleanSupport();
		SHAPE plate2;
		plate2.init("./centralPlate.vtu");
		plate2.cleaning();
		ivector1D listPlate2(plate2.nSimplex);
		for(int i=0; i<plate2.nSimplex; ++i)	listPlate2[i]=i;
		
		support->setSupport(&plate2, listPlate2);
			
// 		{
// 			out<<"...TESTING MANUAL SPLITTING..."<<endl;
// 			
// 		
// 			{
// 				bool check = support->getBoundaryManipulationStatus();
// 				out<<"......test getBoundaryManipulationStatus method ...";
// 				if(check){
// 					out<<"OK"<<endl;
// 				}else{
// 					out<<"FAILED"<<endl;		
// 				}
// 			}
// 			
// 			{
// 				support->switchBoundaryManipulation(false);
// 				support->addSplitPoint(4, true);
// 				dvecarr3E prova = support->getSplitPoints();
// 				bool check = (prova.size() == 0);
// 				out<<"......test switchBoundaryManipulation method ...";
// 				if(check){
// 					out<<"OK"<<endl;
// 				}else{
// 					out<<"FAILED"<<endl;		
// 				}
// 			}
// 			
// 			{
// 				support->switchBoundaryManipulation(true);
// 				
// 				darray3E p1,p2,p3,p4,p5;
// 				p1[0] = 2.6; p1[1] = 1.0; p1[2] = 0.0;
// 				p2[0] = 2.8; p2[1] = 0.0; p2[2] = 0.0;				
// 				p3[0] = 5.4; p3[1] = 1.0; p3[2] = 0.0;
// 				p4[0] = 3.8; p4[1] = 0.0; p4[2] = 0.0;
// 				p5[0] = 5.05; p5[1] = 0.0; p5[2] = 0.0;
// 				
// 				dvecarr3E pp(3);
// 				pp[0] = p1; pp[1] = p2; pp[2] = p3;
// 				bool check1 = support->addSplitPoint(p4, false);
//  				check1 = check1 && support->addSplitPoint(p5, false);
//  				check1 = check1 && support->addSplitPoints(pp);
// 				check1 = check1 && support->removeSplitPoint(p4, false);
// 				check1 = check1 && support->splitBoundary();
// 				dvecarr3E prova = support->getSplitPoints();
// 				bool check = (prova.size() == 4);
// 				
// 				out<<"......test add/remove/get Split Points methods ...";
// 				if(check){
// 					out<<"OK"<<endl;
// 				}else{
// 					out<<"FAILED"<<endl;		
// 				}
// 				
// 				support->plotBoundaryBranches(folder, "manualBsplitting", true, "3DCurve");
// 				out<<"......test manual boundary splitting...";
// 				if(check1){
// 					out<<"OK->see output result"<<endl;
// 				}else{
// 					out<<"FAILED"<<endl;		
// 				}
// 				
// 			}
// 		}
// 		{
// 			support->cleanSplitBoundaries();
// 			dvecarr3E prova = support->getSplitPoints();
// 			bool check = (prova.size() == 0);
// 			out<<"......test cleanSplitBoundaries method ...";
// 			if(check){
// 				out<<"OK"<<endl;
// 			}else{
// 				out<<"FAILED"<<endl;		
// 			}
// 		}
		{
			out<<"...TESTING AUTO SPLITTING..."<<endl;
			
			darray3E p1;
			p1[0] = 2.6; p1[1] = 1.0; p1[2] = 0.0;
			int N_=4;
			int n_eff;
			
			support->plotVTU(folder, "autoTrueSupport", true);
			
			bool check1 = support->splitBoundary_auto(p1,N_, n_eff);
			
//			bool check = (norm_2(support->getSeed() -  p1)<1.0e-12);
//			out<<"......test getSeed method ...";
//			if(check){
//				out<<"OK"<<endl;
//			}else{
//				out<<"FAILED"<<endl;		
//			}
			
			
//			int branch = support->whichBranch(p1);
//			check = (branch != -1);
//			out<<"......test whichBranch and getNearestBPoint methods ...";
//			if(check){
//				out<<"OK"<<endl;
//			}else{
//				out<<"FAILED"<<endl;		
//			}
			
//			int new_branch = (branch+2)%support->getBoundaryBranchSize();
//			int cc = support->orientSegmentation(new_branch);
//			check = (cc != -1);
//			out<<"......test orientSegmentation method ...";
//			if(check){
//				out<<"OK"<<endl;
//			}else{
//				out<<"FAILED"<<endl;		
//			}
			
			support->plotBoundaryBranches(folder, "autoBsplitting", true, "3DCurve");
			out<<"......test auto boundary splitting...";
			if(check1){
				out<<"OK->see output result"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		exit(1);
		support->cleanSplitBoundaries();	
		{
			out<<"...TESTING EXT GEOMETRIES SPLITTING..."<<endl;
			darray3E p1;
			p1[0] = 2.6; p1[1] = 1.0; p1[2] = 0.0;
			int cc = support->sortBoundary(p1);
			bool check = (cc != -1);
			out<<"......test sortBoundary methods ...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
			
			svector1D geoList(2,"");
			geoList[1] = "./leftNeighs.vtu";
			geoList[2] = "./rightNeighs.vtu";
			bool check1 = support->splitBoundary_geoList(geoList);
			
			check = (support->getSplittingType() ==1);
			out<<"......test getSplittingType method ...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
			
			check = (support->getGeoList() ==geoList);
			out<<"......test getGeoList method ...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
			
			
			support->plotBoundaryBranches(folder, "geoListBsplitting", true, "3DCurve");
			out<<"......test by ext geometries boundary splitting...";
			if(check1){
				out<<"OK->see output result"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}	
				
		out<<"...test getBoundaryBranch methods...OK"<<endl;
		delete support; support = NULL;
	}	
}else{
	cout<<"Not able to open log file"<<endl;
}

out.close();	
return(0);

} //end program 
