/* test 004 -> Basic Structured Mesh of libRA_BasicMeshes2D
 * Testing all methods of all classes defining a basic structured mesh 2D.
 */

 using namespace std;
#include <iostream>
#include <fstream> 
#include "customOperators.hpp" 
#include "libRA_BasicMeshes2D.hpp"

int main(int argc, char*argv[]){
	
std::ofstream out;
out.open("testingBasicMeshes2D.log");
if(out.is_open()){
		
	//testing BASE methods/UCubicMesh methods.
	{
		out<<"Testing BASE_UstructMesh2D/UQuadMesh2D methods..."<<endl;
		
		darray2E origin; origin.fill(1.0);
		double spanX = 3;
		double spanY = 4;
		int nx = 3;
		int ny = 4;
		
		BASE_UStructMesh2D * mesh = new UQuadMesh(origin, spanX,spanY, nx,ny);
		
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
			darray2E sp;
			sp[0] = spanX; sp[1]=spanY; 
			check = (mesh->getSpan()==sp);
			out<<"...test getSpan method...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			darray2E sp;
			sp[0] = spanX/nx; sp[1]=spanY/ny;
			check = (mesh->getSpacing()==sp);
			out<<"...test getSpacing method...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			ivector1D sp(2);
			sp[0] = nx; sp[1]=ny;
			check = (mesh->getDimension()== sp);
			out<<"...test getDimension method...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			int p1 = mesh->accessCellData(1,2);
			
			int i,j;
			mesh->accessCellData(p1,i,j);
					
			check = (i== 1 && j==2);
			out<<"...test accessCellData methods...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			int p1 = mesh->accessPointData(3,2);
			int i,j;
			mesh->accessPointData(p1,i,j);
			check = (i== 3 && j==2);
			out<<"...test accessPointData methods...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			darray2E point1; 
			point1[0]=1.5; point1[1]=2.5;
			dvector1D point2 = conVect(point1);
			int i1,j1,i2,j2;
			mesh->returnCellID(point1,i1,j1);
			mesh->returnCellID(point2,i2,j2);
			
			check = (i1== 0 && j1==1);
			check = check && (i2== 0 && j2==1 );
			out<<"...test returnCellID methods...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			darray2E ccell{1.5,2.5};
		        int val = mesh->accessCellData(1,2);
			check = (mesh->getGridCCell(1,2)== ccell);
			check = check &&(mesh->getGridCCell(val)== ccell);
			out<<"...test getGridCCell methods...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			darray2E point{1.0,2.0};
			int val = mesh->accessPointData(1,2);
			check = (mesh->getGridPoint(1,2)== point);
			check = check &&(mesh->getGridPoint(val)== point);
			out<<"...test getGridPoint methods...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			darray2E ccell{1.5,2.5};
			ccell = ccell + mesh->getOrigin();
			int val = mesh->accessCellData(1,2);
			check = (mesh->getGlobalCCell(1,2)== ccell);
			check = check &&(mesh->getGlobalCCell(val)== ccell);
			out<<"...test getGlobalCCell methods...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			darray2E point{1.0,2.0};
			point = point + mesh->getOrigin();
			int val = mesh->accessPointData(1,2);
			check = (mesh->getGlobalPoint(1,2)== point);
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
			test.push_back(mesh->accessPointData(0,0));
			test.push_back(mesh->accessPointData(1,0));
			test.push_back(mesh->accessPointData(1,1));
			test.push_back(mesh->accessPointData(0,1));
			
			check = (mesh->getCellNeighs(0,0)== test);
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
			int size = dim[0]*dim[1];
			dvector1D celldata(size, 3.42);
			ivector1D celldata2(size, -2);
			darray2E ptarg; ptarg.fill(2.12345);
		
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
			int size = (dim[0]+1)*(dim[1]+1);
			dvector1D ptdata(size, 3.42);
			ivector1D ptdata2(size, -2);
			darray2E ptarg; ptarg.fill(1.12345);
		
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
			int sizeskC = std::floor(0.3*dim[0]*dim[1]+0.5);
			int sizeskP = std::floor(0.6*(dim[0]+1)*(dim[1]+1)+0.5);
		
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

			out<<"...test Plot Grid methods...";
			out<<"see Output results"<<endl;
		}
		{
			BASE_UStructMesh2D *mesh2 = new UQuadMesh();
			*mesh2 = *mesh;
			int i1,i2,i3;
			ivector1D dim = mesh2->getDimension();
			
			i1 = dim[0];
			i2 = dim[1];
			i3 = dim[0]+dim[1];
			
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
			darray2E zero{0,0};
			ivector1D szero(2,0);
			
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
		out<<"Testing UQuadMesh methods..."<<endl;
		
		darray2E origin; origin.fill(1.0);
		double spanX = 3;
		double spanY = 4;
		int nx = 3;
		int ny = 4;
		
		UQuadMesh * mesh = new UQuadMesh();
		
		bool check;
		
		{
			mesh->setMesh(origin, spanX, spanY,nx,ny);
			darray2E spanC; 
			spanC[0] = spanX; spanC[1] = spanY;
			ivector1D dimC(2,nx);
			dimC[1] = ny; 
			darray2E spaceC;
			spaceC[0] = spanX/nx; spaceC[1] = spanY/ny; 
			
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
			check = (mesh->getClassType()=="Cartesian2D");
			out<<"...test getClassType method...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			darray2E spanC{6,8};
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
			dvecarr2E test(2, darray2E{0,0});
			test[1][0] = 4.0;
			dvecarr2E result(2,darray2E{1,1});
			result[1][0] = 5.0;
			
			dvecarr2E target = mesh->transfToGlobal(test);
			
			check = (target[0] == result[0]);
			check = check && (target[1] == result[1]);
			out<<"...test transfToGlobal method...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
			
			dvecarr2E target2 = mesh->transfToLocal(result);
			
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
			darray2E space = mesh->getSpacing();
			double volume = space[0]*space[1];
			
			check = (mesh->getCellVolume(0) == volume);
			out<<"...test getCellVolume method...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			darray2E result{1,1};
			check = (mesh->getLocalScaling(0) == result);
			out<<"...test getLocalScaling method...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			UQuadMesh *mesh2 = new UQuadMesh();
			*mesh2 = *mesh;
			int i1,i2,i3;
			ivector1D dim = mesh2->getDimension();
			
			i1 = dim[0];
			i2 = dim[1];
			i3 = dim[0]+dim[1];
			
			check = (mesh2->getOrigin() == mesh->getOrigin());
			check = check && (mesh2->getSpan() == mesh->getSpan());
			check = check && (mesh2->getDimension() == mesh->getDimension());
			check = check && (mesh2->getSpacing() == mesh->getSpacing());
			check = check && (mesh2->getClassType() == "Cartesian2D");
			
			//sparse check on nodes and point
			check = check && (mesh2->getGridCCell(i1) == mesh->getGridCCell(i1));
			check = check && (mesh2->getGridCCell(i2) == mesh->getGridCCell(i2));
			check = check && (mesh2->getGridPoint(i3) == mesh->getGridPoint(i3));
			
			out<<"...test UQuadMesh copy Operator...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
			
			delete mesh2; mesh2=NULL;
		}

		delete mesh; mesh = NULL;
	}	
}else{
	cout<<"Not able to open log file"<<endl;
}

out.close();	
return(0);

} //end program 
