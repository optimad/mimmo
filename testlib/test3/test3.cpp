/* test 003 -> Basic Structured Mesh of libRA_BasicMeshes
 * Testing all methods of all classes defining a basic structured mesh.
 */

 using namespace std;
#include <iostream>
#include <fstream> 
#include "customOperators.hpp" 
#include "libRA_BasicMeshes.hpp"

int main(int argc, char*argv[]){
	
std::ofstream out;
out.open("testingBasicMeshes.log");
if(out.is_open()){
		
	//testing BASE methods/UCubicMesh methods.
	{
		out<<"Testing BASE_UstructMesh/UCubicMesh methods..."<<endl;
		
		darray3E origin; origin.fill(1.0);
		double spanX = 3;
		double spanY = 4;
		double spanZ = 5;
		int nx = 3;
		int ny = 4;
		int nz = 5;
		
		BASE_UStructMesh * mesh = new UCubicMesh(origin, spanX,spanY,spanZ, nx,ny,nz);
		
		bool check;
		
		{
			check = (mesh->getOrigin()==origin);
			out<<"...test getOrigin method...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			darray3E sp;
			sp[0] = spanX; sp[1]=spanY; sp[2] = spanZ;
			check = (mesh->getSpan()==sp);
			out<<"...test getSpan method...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			darray3E sp;
			sp[0] = spanX/nx; sp[1]=spanY/ny; sp[2] = spanZ/nz;
			check = (mesh->getSpacing()==sp);
			out<<"...test getSpacing method...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			ivector1D sp(3);
			sp[0] = nx; sp[1]=ny; sp[2] = nz;
			check = (mesh->getDimension()== sp);
			out<<"...test getDimension method...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			int p1 = mesh->accessCellData(1,2,3);
			
			int i,j,k;
			mesh->accessCellData(p1,i,j,k);
					
			check = (i== 1 && j==2 && k==3);
			out<<"...test accessCellData methods...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			int p1 = mesh->accessPointData(3,2,4);
			int i,j,k;
			mesh->accessPointData(p1,i,j,k);
			check = (i== 3 && j==2 && k==4);
			out<<"...test accessPointData methods...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			darray3E point1; 
			point1[0]=1.5; point1[1]=2.5; point1[2]=3.5;
			dvector1D point2 = conVect(point1);
			int i1,j1,k1, i2,j2,k2;
			mesh->returnCellID(point1,i1,j1,k1);
			mesh->returnCellID(point2,i2,j2,k2);
			
			check = (i1== 0 && j1==1 && k1==2);
			check = check && (i2== 0 && j2==1 && k2==2);
			out<<"...test returnCellID methods...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			darray3E ccell{1.5,2.5,3.5};
		        int val = mesh->accessCellData(1,2,3);
			check = (mesh->getGridCCell(1,2,3)== ccell);
			check = check &&(mesh->getGridCCell(val)== ccell);
			out<<"...test getGridCCell methods...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			darray3E point{1.0,2.0,3.0};
			int val = mesh->accessPointData(1,2,3);
			check = (mesh->getGridPoint(1,2,3)== point);
			check = check &&(mesh->getGridPoint(val)== point);
			out<<"...test getGridPoint methods...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			darray3E ccell{1.5,2.5,3.5};
			ccell = ccell + mesh->getOrigin();
			int val = mesh->accessCellData(1,2,3);
			check = (mesh->getGlobalCCell(1,2,3)== ccell);
			check = check &&(mesh->getGlobalCCell(val)== ccell);
			out<<"...test getGlobalCCell methods...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			darray3E point{1.0,2.0,3.0};
			point = point + mesh->getOrigin();
			int val = mesh->accessPointData(1,2,3);
			check = (mesh->getGlobalPoint(1,2,3)== point);
			check = check &&(mesh->getGlobalPoint(val)== point);
			out<<"...test getGlobalPoint methods...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			ivector1D  test;
			test.push_back(mesh->accessPointData(0,0,0));
			test.push_back(mesh->accessPointData(1,0,0));
			test.push_back(mesh->accessPointData(1,1,0));
			test.push_back(mesh->accessPointData(0,1,0));
			test.push_back(mesh->accessPointData(0,0,1));
			test.push_back(mesh->accessPointData(1,0,1));
			test.push_back(mesh->accessPointData(1,1,1));
			test.push_back(mesh->accessPointData(0,1,1));
			
			check = (mesh->getCellNeighs(0,0,0)== test);
			check = check &&(mesh->getCellNeighs(0)== test);
			out<<"...test getCellNeighs methods...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			ivector1D dim = mesh->getDimension();
			int size = dim[0]*dim[1]*dim[2];
			dvector1D celldata(size, 3.42);
			ivector1D celldata2(size, -2);
			darray3E ptarg; ptarg.fill(2.12345);
		
			double result1 = mesh->interpolateCellData(ptarg, celldata);
			int result2 = mesh->interpolateCellData(ptarg, celldata2);
			
			check = (abs(result1-3.42)<1.0e-12) && (result2 == -2);
			out<<"...test interpolateCellData methods...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			ivector1D dim = mesh->getDimension();
			int size = (dim[0]+1)*(dim[1]+1)*(dim[2]+1);
			dvector1D ptdata(size, 3.42);
			ivector1D ptdata2(size, -2);
			darray3E ptarg; ptarg.fill(1.12345);
		
			double result1 = mesh->interpolatePointData(ptarg, ptdata);
			int result2 = mesh->interpolatePointData(ptarg, ptdata2);
			
			
			check = (abs(result1-3.42)<1.0e-12) && (result2 == -2);
			out<<"...test interpolatePointData methods...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			std::string folder = ".";
			ivector1D dim = mesh->getDimension();
			int sizeskC = std::floor(0.3*dim[0]*dim[1]*dim[2]+0.5);
			int sizeskP = std::floor(0.6*(dim[0]+1)*(dim[1]+1)*(dim[2]+1)+0.5);
		
			ivector1D skC(sizeskC), skP(sizeskP);
			
			for(int i=0; i<sizeskC; ++i){
				skC[i] = i;
			}
			for(int i=0; i<sizeskP; ++i){
				skP[i] = i;
			}
			
			mesh->plotCloud(folder,"meshCloud", 0, true);
			mesh->plotCloud(folder,"meshCloudPart", 1, false, skP);
			mesh->plotGrid(folder,"meshGrid", 2, true);
			mesh->plotGrid(folder,"meshGridPart", 3, false, skC);

			out<<"...test interpolatePointData methods...";
			out<<"see Output results"<<endl;
		}
		{
			BASE_UStructMesh *mesh2 = new UCubicMesh();
			*mesh2 = *mesh;
			int i1,i2,i3;
			ivector1D dim = mesh2->getDimension();
			
			i1 = dim[0];
			i2 = dim[2];
			i3 = dim[0]+dim[1]+dim[2];
			
			check = (mesh2->getOrigin() == mesh->getOrigin());
			check = check && (mesh2->getSpan() == mesh->getSpan());
			check = check && (mesh2->getDimension() == mesh->getDimension());
			check = check && (mesh2->getSpacing() == mesh->getSpacing());
			
			//sparse check on nodes and point
			check = check && (mesh2->getGridCCell(i1) == mesh->getGridCCell(i1));
			check = check && (mesh2->getGridCCell(i2) == mesh->getGridCCell(i2));
			check = check && (mesh2->getGridPoint(i3) == mesh->getGridPoint(i3));
			
			out<<"...test BASE copy Operator...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
			
			delete mesh2; mesh2=NULL;
		}
		{
			mesh->clearMesh();
			darray3E zero{0,0,0};
			ivector1D szero(3,0);
			
			check = (mesh->getOrigin() == zero);
			check = check && (mesh->getSpan() == zero);
			check = check && (mesh->getDimension() == szero);
			check = check && (mesh->getSpacing() == zero);

			out<<"...test clear Mesh...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		delete mesh; mesh = NULL;
	}	
	/////////////////////////////////////////////////////////////
	//testing UCubicMesh methods.
	{
		out<<"Testing UCubicMesh methods..."<<endl;
		
		darray3E origin; origin.fill(1.0);
		double spanX = 3;
		double spanY = 4;
		double spanZ = 5;
		int nx = 3;
		int ny = 4;
		int nz = 5;
		
		UCubicMesh * mesh = new UCubicMesh();
		
		bool check;
		
		{
			mesh->setMesh(origin, spanX, spanY, spanZ, nx,ny,nz);
			darray3E spanC; 
			spanC[0] = spanX; spanC[1] = spanY; spanC[2] = spanZ;
			ivector1D dimC(3,nx);
			dimC[1] = ny; dimC[2]=nz;
			darray3E spaceC;
			spaceC[0] = spanX/nx; spaceC[1] = spanY/ny; spaceC[2] = spanZ/nz;
			
			check = (mesh->getOrigin()==origin);
			check = check && (mesh->getSpan()== spanC);
			check = check && (mesh->getDimension()== dimC);
			check = check && (mesh->getSpacing()== spaceC);
			out<<"...test setMesh method...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			check = (mesh->getClassType()=="Cartesian");
			out<<"...test getClassType method...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			darray3E spanC{6,8,10};
			mesh->scaleMesh(2,2,2);
			check = check && (mesh->getSpan() == spanC);
			out<<"...test scaleMesh method...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
			
			mesh->scaleMesh(0.5,0.5,0.5);
		}
		{
			dvecarr3E test(2, darray3E{0,0,0});
			test[1][0] = 4.0;
			dvecarr3E result(2,darray3E{1,1,1});
			result[1][0] = 5.0;
			
			dvecarr3E target = mesh->transfToGlobal(test);
			
			check = (target[0] == result[0]);
			check = check && (target[1] == result[1]);
			out<<"...test transfToGlobal method...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
			
			dvecarr3E target2 = mesh->transfToLocal(result);
			
			check = (target2[0] == test[0]);
			check = check && (target2[1] == test[1]);
			out<<"...test transfToLocal method...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			darray3E space = mesh->getSpacing();
			double volume = space[0]*space[1]*space[2];
			
			check = (mesh->getCellVolume(0) == volume);
			out<<"...test getCellVolume method...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			darray3E result{1,1,1};
			check = (mesh->getLocalScaling(0) == result);
			out<<"...test getLocalScaling method...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			UCubicMesh *mesh2 = new UCubicMesh();
			*mesh2 = *mesh;
			int i1,i2,i3;
			ivector1D dim = mesh2->getDimension();
			
			i1 = dim[0];
			i2 = dim[2];
			i3 = dim[0]+dim[1]+dim[2];
			
			check = (mesh2->getOrigin() == mesh->getOrigin());
			check = check && (mesh2->getSpan() == mesh->getSpan());
			check = check && (mesh2->getDimension() == mesh->getDimension());
			check = check && (mesh2->getSpacing() == mesh->getSpacing());
			check = check && (mesh2->getClassType() == "Cartesian");
			
			//sparse check on nodes and point
			check = check && (mesh2->getGridCCell(i1) == mesh->getGridCCell(i1));
			check = check && (mesh2->getGridCCell(i2) == mesh->getGridCCell(i2));
			check = check && (mesh2->getGridPoint(i3) == mesh->getGridPoint(i3));
			
			out<<"...test UCubicMesh copy Operator...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
			
			delete mesh2; mesh2=NULL;
		}

		delete mesh; mesh = NULL;
	}	
	/////////////////////////////////////////////////////////////
	//testing UCylindricalMesh methods.
	{
		out<<"Testing UCylindricalMesh methods..."<<endl;
		
		darray3E origin; origin.fill(1.0);
		double spanR = 2;
		double spanZ = 6;
		dvector1D thLim(2,0);
		thLim[1] = 8*std::atan(1.0);
		int nr = 2;
		int nt = 15;
		int nz = 4;
		
		UCylindricalMesh * mesh = new UCylindricalMesh();
		
		bool check;
		
		{
			mesh->setMesh(origin, spanR, spanZ, thLim, nr,nt,nz);
			darray3E spanC; 
			spanC[0] = spanR; spanC[1] = thLim[1]; spanC[2] = spanZ;
			ivector1D dimC(3,nr);
			dimC[1] = nt; dimC[2]=nz;
			darray3E spaceC;
			spaceC[0] = spanR/nr; spaceC[1] = thLim[1]/nt; spaceC[2] = spanZ/nz;
			
			check = (mesh->getOrigin()==origin);
			check = check && (mesh->getSpan()== spanC);
			check = check && (mesh->getDimension()== dimC);
			check = check && (mesh->getSpacing()== spaceC);
			out<<"...test setMesh method...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			check = (mesh->getClassType()=="Cylindrical");
			out<<"...test getClassType method...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			darray3E spanC;
			spanC[0] = 4; spanC[1] = thLim[1]; spanC[2] = 12;
			mesh->scaleMesh(2,2);
			check = check && (mesh->getSpan() == spanC);
			out<<"...test scaleMesh method...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
			
			mesh->scaleMesh(0.5,0.5);
		}
		{
			dvecarr3E test(2, darray3E{0,0,0});
			test[0][0] = 1; test[0][2] = 3;
			test[1][0] = 2; test[1][1] = 2*std::atan(1.0); test[1][2] = 5;
			dvecarr3E result(2,darray3E{1,1,1});
			result[0][0] = 2; result[0][2] = 4;
			result[1][1] = 3; result[1][2] = 6;
			
			dvecarr3E target = mesh->transfToGlobal(test);
			check = (norm_2(target[0]-result[0])<1.e-12);
			check = check && ((norm_2(target[1]-result[1])<1.e-12));

			out<<"...test transfToGlobal method...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
			
			dvecarr3E target2 = mesh->transfToLocal(result);
			check = ((norm_2(target2[0]-test[0])<1.e-12));
			check = check && ((norm_2(target2[1]-test[1])<1.e-12));

			out<<"...test transfToLocal method...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			darray3E space = mesh->getSpacing();
			darray3E point = mesh->getGridCCell(0);
			double volume = space[0]*space[1]*space[2]*point[0];
			
			check = (mesh->getCellVolume(0) == volume);
			out<<"...test getCellVolume method...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			darray3E result{1,1,1};
			darray3E point = mesh->getGridCCell(0);
			result[1] = point[0];
			check = (mesh->getLocalScaling(0) == result);
			out<<"...test getLocalScaling method...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			UCylindricalMesh *mesh2 = new UCylindricalMesh();
			*mesh2 = *mesh;
			int i1,i2,i3;
			ivector1D dim = mesh2->getDimension();
			
			i1 = dim[0];
			i2 = dim[2];
			i3 = dim[0]+dim[1]+dim[2];
			
			check = (mesh2->getOrigin() == mesh->getOrigin());
			check = check && (mesh2->getSpan() == mesh->getSpan());
			check = check && (mesh2->getDimension() == mesh->getDimension());
			check = check && (mesh2->getSpacing() == mesh->getSpacing());
			check = check && (mesh2->getClassType() == "Cylindrical");
			
			//sparse check on nodes and point
			check = check && (mesh2->getGridCCell(i1) == mesh->getGridCCell(i1));
			check = check && (mesh2->getGridCCell(i2) == mesh->getGridCCell(i2));
			check = check && (mesh2->getGridPoint(i3) == mesh->getGridPoint(i3));
			
			out<<"...test UCylindricalMesh copy Operator...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
			
			delete mesh2; mesh2=NULL;
		}
		
		
		{
			std::string folder = ".";
			ivector1D dim = mesh->getDimension();
			int sizeskC = std::floor(0.3*dim[0]*dim[1]*dim[2]+0.5);
			int sizeskP = std::floor(0.6*(dim[0]+1)*(dim[1]+1)*(dim[2]+1)+0.5);
			
			ivector1D skC(sizeskC), skP(sizeskP);
			
			for(int i=0; i<sizeskC; ++i){
				skC[i] = i;
			}
			for(int i=0; i<sizeskP; ++i){
				skP[i] = i;
			}
			
			mesh->plotGrid(folder,"cylinderGrid", 2, true);
			mesh->plotGrid(folder,"cylinderGridPart", 3, false, skC);
			
			out<<"...test plotting Cylindrical mesh ...";
			out<<"see Output results"<<endl;	
		}
		
		delete mesh; mesh = NULL;
	}	
	/////////////////////////////////////////////////////////////	
	//testing UCylindricalMesh methods.
	{
		out<<"Testing USphericalMesh methods..."<<endl;
		
		darray3E origin; origin.fill(1.0);
		double spanR = 2;
		dvector1D thLim(2,0);
		dvector1D phLim(2,0);
		thLim[1] = 8*std::atan(1.0);
		phLim[1] = 4*std::atan(1.0);
		int nr = 2;
		int nt = 15;
		int np = 8;
		
		USphericalMesh * mesh = new USphericalMesh();
		
		bool check;
		
		{
			mesh->setMesh(origin, spanR, thLim, phLim, nr,nt,np);
			darray3E spanC; 
			spanC[0] = spanR; spanC[1] = thLim[1]; spanC[2] = phLim[1];
			ivector1D dimC(3,nr);
			dimC[1] = nt; dimC[2]=np;
			darray3E spaceC;
			spaceC[0] = spanR/nr; spaceC[1] = thLim[1]/nt; spaceC[2] = phLim[1]/np;
			
			check = (mesh->getOrigin()==origin);
			check = check && (mesh->getSpan()== spanC);
			check = check && (mesh->getDimension()== dimC);
			check = check && (mesh->getSpacing()== spaceC);
			out<<"...test setMesh method...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			check = (mesh->getClassType()=="Spherical");
			out<<"...test getClassType method...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			darray3E spanC;
			spanC[0] = 4; spanC[1] = thLim[1]; spanC[2] = phLim[1];
			mesh->scaleMesh(2);
			check = check && (mesh->getSpan() == spanC);
			out<<"...test scaleMesh method...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
			
			mesh->scaleMesh(0.5);
		}
		{
			dvecarr3E test(2, darray3E{0,0,0});
			test[0][0] = 1; 
			test[1][0] = 2; test[1][1] = 2*std::atan(1.0); test[1][2] = 2*std::atan(1.0);
			dvecarr3E result(2,darray3E{1,1,1});
			result[0][2] = 2;
			result[1][1] = 3; 
			
			dvecarr3E target = mesh->transfToGlobal(test);
			check = (norm_2(target[0]-result[0])<1.e-12);
			check = check && ((norm_2(target[1]-result[1])<1.e-12));
			
			out<<"...test transfToGlobal method...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
			
			dvecarr3E target2 = mesh->transfToLocal(result);
			check = ((norm_2(target2[0]-test[0])<1.e-12));
			check = check && ((norm_2(target2[1]-test[1])<1.e-12));
			
			out<<"...test transfToLocal method...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			darray3E space = mesh->getSpacing();
			darray3E point = mesh->getGridCCell(0);
			double volume = space[0]*space[1]*space[2]*point[0]*point[0]*std::sin(point[2]);
			
			check = (mesh->getCellVolume(0) == volume);
			out<<"...test getCellVolume method...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			darray3E result{1,1,1};
			darray3E point = mesh->getGridCCell(0);
			result[1] = point[0]*std::sin(point[2]);
			result[2] = point[0];
			check = (mesh->getLocalScaling(0) == result);
			out<<"...test getLocalScaling method...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			USphericalMesh *mesh2 = new USphericalMesh();
			*mesh2 = *mesh;
			int i1,i2,i3;
			ivector1D dim = mesh2->getDimension();
			
			i1 = dim[0];
			i2 = dim[2];
			i3 = dim[0]+dim[1]+dim[2];
			
			check = (mesh2->getOrigin() == mesh->getOrigin());
			check = check && (mesh2->getSpan() == mesh->getSpan());
			check = check && (mesh2->getDimension() == mesh->getDimension());
			check = check && (mesh2->getSpacing() == mesh->getSpacing());
			check = check && (mesh2->getClassType() == "Spherical");
			
			//sparse check on nodes and point
			check = check && (mesh2->getGridCCell(i1) == mesh->getGridCCell(i1));
			check = check && (mesh2->getGridCCell(i2) == mesh->getGridCCell(i2));
			check = check && (mesh2->getGridPoint(i3) == mesh->getGridPoint(i3));
			
			out<<"...test USpericalMesh copy Operator...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
			
			delete mesh2; mesh2=NULL;
		}
		
		
		{
			std::string folder = ".";
			ivector1D dim = mesh->getDimension();
			int sizeskC = std::floor(0.3*dim[0]*dim[1]*dim[2]+0.5);
			int sizeskP = std::floor(0.6*(dim[0]+1)*(dim[1]+1)*(dim[2]+1)+0.5);
			
			ivector1D skC(sizeskC), skP(sizeskP);
			
			for(int i=0; i<sizeskC; ++i){
				skC[i] = i;
			}
			for(int i=0; i<sizeskP; ++i){
				skP[i] = i;
			}
			
			mesh->plotGrid(folder,"sphereGrid", 2, true);
			mesh->plotGrid(folder,"sphereGridPart", 3, false, skC);
			
			out<<"...test plotting Spherical mesh ...";
			out<<"see Output results"<<endl;	
		}
		
		delete mesh; mesh = NULL;
	}	
	/////////////////////////////////////////////////////////////	
	
}else{
	cout<<"Not able to open log file"<<endl;
}

out.close();	
return(0);

} //end program 
