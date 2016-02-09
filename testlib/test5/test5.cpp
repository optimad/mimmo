/* test 005 -> Core Patches of libRA_CorePatches
 * Testing all methods of all classes defining a an elementary PATCH.
 */

 using namespace std;
#include <iostream>
#include <fstream> 
#include "customOperators.hpp" 
#include "libRA_CorePatches.hpp"

int main(int argc, char*argv[]){
	
std::ofstream out;
out.open("testingCorePatches.log");
if(out.is_open()){
	
	
	SHAPE geo;
	geo.init("./sphere.nas");
	geo.cleaning();
	dvector1D xlim(2,0),ylim(2,0), zlim(2,0);
	geo.BoundingBox(xlim,ylim,zlim);
	std::string folder = ".";
	int pidtarg;
	if(geo.pidType.size() != 0){
		pidtarg = geo.pidType[0];
	}
	
	
	//testing BASE methods/PatchCube methods.
	{
		out<<"Testing BASE_PATCH/PATCH_CUBE methods..."<<endl;
		
		darray3E origin; 
		origin[0] = xlim[0] - 0.3*(xlim[1] - xlim[0]);
		origin[1] = ylim[0] - 0.3*(ylim[1] - ylim[0]);
		origin[2] = zlim[0] + 0.3*(zlim[1] - zlim[0]);
		
		double spanX = (xlim[1] - xlim[0]);
		double spanY = (ylim[1] - ylim[0]);
		double spanZ = (zlim[1] - zlim[0]);	
		
		BASE_PATCH * patch = new PATCH_CUBE(origin, spanX,spanY,spanZ);
		{
			bool check = (patch->getPatchType()=="patchCube");
			out<<"...test getPatchType method...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			ivector1D list = patch->includeTriangulation(&geo);
			geo.plotExtraction(folder, "includeCube", true, 0, list);
			out<<"...test includeTriangulation method...see output file"<<endl;
		}
		{
			ivector1D list = patch->excludeTriangulation(&geo);
			geo.plotExtraction(folder, "excludeCube", true, 1, list);
			out<<"...test excludeTriangulation method...see output file"<<endl;
		}
		{
			ivector1D list = patch->includeCloudPoints(geo.Vertex);
			geo.plotCloudExtraction(folder, "includeCloudCube", true, 2, list);
			out<<"...test includeCloudPoints method...see output file"<<endl;
				
		}
		
		{
			ivector1D list = patch->excludeCloudPoints(geo.Vertex);
			geo.plotCloudExtraction(folder, "excludeCloudCube", true, 3, list);
			out<<"...test excludeCloudPoints method...see output file"<<endl;
		}
		
		{
			ivector1D list = patch->includePIDTriangulation(&geo, pidtarg);
			geo.plotExtraction(folder, "includePIDCube", true, 4, list);
			out<<"...test includePIDTriangulation method...see output file"<<endl;
			
		}
		{
			ivector1D list = patch->excludePIDTriangulation(&geo, pidtarg);
			geo.plotExtraction(folder, "excludePIDCube", true, 5, list);
			out<<"...test excludePIDTriangulation method...see output file"<<endl;
		}
		{
			(dynamic_cast<PATCH_CUBE * >(patch))->plotVTUpatch(".", "patchCubeExample",6);
			out<<"...test plotVTUpatch method...see output file"<<endl;
		}
		
		delete patch; patch = NULL;
	}	
	///////////////////////////////////////////
	//testing PatchCylinder methods.
	{
		out<<"Testing PATCH_CYLINDER methods..."<<endl;
		
		darray3E origin; 
		origin[0] = xlim[0] +0.5*(xlim[1] - xlim[0]);
		origin[1] = ylim[0] +0.5*(ylim[1] - ylim[0]);
		origin[2] = zlim[0] +0.5*(zlim[1] - zlim[0]);
		
		double spanX = 0.32*(xlim[1] - xlim[0]);
		double spanZ = 0.6*(zlim[1] - zlim[0]);	
		
		PATCH_CYLINDER * patch = new PATCH_CYLINDER(origin, spanX,spanZ);
		{
			bool check = (patch->getPatchType()=="patchCylinder");
			out<<"...test getPatchType method...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			ivector1D list = patch->includeTriangulation(&geo);
			geo.plotExtraction(folder, "includeCylinder", true, 7, list);
			out<<"...test includeTriangulation method...see output file"<<endl;
		}
		
		delete patch; patch = NULL;
	}	
	///////////////////////////////////////////
	//testing PatchSphere methods.
	{
		out<<"Testing PATCH_SPHERE methods..."<<endl;
		
		darray3E origin; 
		origin[0] = xlim[0] -0.15*(xlim[1] - xlim[0]);
		origin[1] = ylim[0] -0.21*(ylim[1] - ylim[0]);
		origin[2] = zlim[0] +0.25*(zlim[1] - zlim[0]);
		
		double spanR = ((xlim[1] - xlim[0])+(ylim[1] - ylim[0])+ (zlim[1] - zlim[0]))/3.0 ;
		
		
		PATCH_SPHERE * patch = new PATCH_SPHERE(origin, spanR);
		{
			bool check = (patch->getPatchType()=="patchSphere");
			out<<"...test getPatchType method...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			ivector1D list = patch->includeTriangulation(&geo);
			geo.plotExtraction(folder, "includeSphere", true, 8, list);
			out<<"...test includeTriangulation method...see output file"<<endl;
		}
		delete patch; patch = NULL;
	}	
	///////////////////////////////////////////
	//testing PatchHexa methods.
	{
		out<<"Testing PATCH_HEXAHEDRON methods..."<<endl;
		
		darray3E origin; 
		double sx,sy,sz;
		darray3E xv{1,0,0}, yv{0,1,0}, zv{0,0,1};
		
		sx = (xlim[1] - xlim[0]);
		sy = (ylim[1] - ylim[0]);
		sz = (zlim[1] - zlim[0]);
		
		origin[0] = xlim[0] +0.5*sx;
		origin[1] = ylim[0] +0.5*sy;
		origin[2] = zlim[0] +0.5*sz;
		
		dvecarr3E vertList(8);
		vertList[0] = origin - 0.2*sx*xv - 0.2*sy*yv;
		vertList[1] = origin + 0.2*sx*xv - 0.2*sy*yv;
		vertList[2] = origin + 0.2*sx*xv + 0.2*sy*yv;
		vertList[3] = origin - 0.2*sx*xv + 0.2*sy*yv;
		vertList[4] = origin - 0.6*sx*xv - 0.6*sy*yv + 0.6*sz*zv;
		vertList[5] = origin + 0.6*sx*xv - 0.6*sy*yv + 0.6*sz*zv;
		vertList[6] = origin + 0.6*sx*xv + 0.6*sy*yv + 0.6*sz*zv;
		vertList[7] = origin - 0.6*sx*xv + 0.6*sy*yv + 0.6*sz*zv;
		
		
		PATCH_HEXA * patch = new PATCH_HEXA(vertList);
		{
			bool check = (patch->getPatchType()=="patchHexa");
			out<<"...test getPatchType method...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			ivector1D list = patch->includeTriangulation(&geo);
			geo.plotExtraction(folder, "includeHexaHedron", true, 9, list);
			out<<"...test includeTriangulation method...see output file"<<endl;
		}
		
		{
			patch->plotVTUpatch(".", "patchHexaExample",10);
			out<<"...test plotVTUpatch method...see output file"<<endl;
		}
		delete patch; patch = NULL;
	}	

	///////////////////////////////////////////
	//testing HullCube methods.
	{
		out<<"Testing HULL_CUBE methods..."<<endl;
		
		darray3E origin; 
		origin[0] = xlim[0] +0.5*(xlim[1] - xlim[0]);
		origin[1] = ylim[0] +0.5*(ylim[1] - ylim[0]);
		origin[2] = zlim[0] +0.5*(zlim[1] - zlim[0]);
		
		double spanX = 0.6*(xlim[1] - xlim[0]);
		double spanY = 0.4*(xlim[1] - xlim[0]);
		double spanZ = 0.6*(zlim[1] - zlim[0]);	
		
		HULL_CUBE * patch = new HULL_CUBE(origin, spanX, spanY,spanZ, 10,6,8);
		{
			bool check = (patch->getPatchType()=="hullCube");
			out<<"...test getPatchType method...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			ivector1D list = patch->includeTriangulation(&geo);
			geo.plotExtraction(folder, "includeCubicMesh", true, 11, list);
			out<<"...test includeTriangulation method...see output file"<<endl;
			std::string folder = ".";
			patch->plotGrid(folder,"cubicMesh",12,true);
		}
		
		delete patch; patch = NULL;
	}	
	//testing HullCylinder methods.
	{
		out<<"Testing HULL_CYLINDER methods..."<<endl;
		
		darray3E origin; 
		origin[0] = xlim[0] +0.5*(xlim[1] - xlim[0]);
		origin[1] = ylim[0] +0.5*(ylim[1] - ylim[0]);
		origin[2] = zlim[0] +0.5*(zlim[1] - zlim[0]);
		
		double spanR = 0.23*(xlim[1] - xlim[0]);
		double spanZ = 0.6*(zlim[1] - zlim[0]);	
		dvector1D philim(2,0);
		philim[1] = 4.5*std::atan(1.0);
		
		HULL_CYLINDER * patch = new HULL_CYLINDER(origin, spanR,spanZ, philim, 3,20,10);
		{
			bool check = (patch->getPatchType()=="hullCylinder");
			out<<"...test getPatchType method...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			ivector1D list = patch->includeTriangulation(&geo);
			geo.plotExtraction(folder, "includeCylindricalMesh", true, 12, list);
			out<<"...test includeTriangulation method...see output file"<<endl;
			std::string folder = ".";
			patch->plotGrid(folder,"cylindricalMesh",13,true);
		}
		
		{
			patch->shiftAzimuthalOrigin(atan(1.0));
			ivector1D list = patch->includeTriangulation(&geo);
			geo.plotExtraction(folder, "includeCylindricalMeshANGULARSHIFTED", true, 14, list);
			out<<"...test includeTriangulation method...see output file"<<endl;
			std::string folder = ".";
			patch->plotGrid(folder,"cylindricalMeshANGULARSHIFTED",15,true);
		}
		
		delete patch; patch = NULL;
	}	

	//testing HullSphere methods.
	{
		out<<"Testing HULL_SPHERE methods..."<<endl;
		
		darray3E origin; 
		origin[0] = xlim[0] +0.5*(xlim[1] - xlim[0]);
		origin[1] = ylim[0] +0.5*(ylim[1] - ylim[0]);
		origin[2] = zlim[0] +0.75*(zlim[1] - zlim[0]);
		
		double spanR = 0.42*(zlim[1] - zlim[0]);
		dvector1D thetalim(2,0), philim(2,0);
		philim[0] = 0;//2.0*std::atan(1.0);
		philim[1] = 2.0*std::atan(1.0);//4.0*std::atan(1.0);
		thetalim[1] = 6*std::atan(1.0);
		
		HULL_SPHERE * patch = new HULL_SPHERE(origin, spanR,thetalim, philim, 3,25,25);
		{
			bool check = (patch->getPatchType()=="hullSphere");
			out<<"...test getPatchType method...";
			if(check){
				out<<"OK"<<endl;
			}else{
				out<<"FAILED"<<endl;		
			}
		}
		{
			ivector1D list = patch->includeTriangulation(&geo);
			geo.plotExtraction(folder, "includeSphericalMesh", true, 16, list);
			out<<"...test includeTriangulation method...see output file"<<endl;
			std::string folder = ".";
			patch->plotGrid(folder,"sphericalMesh",17,true);
		}
		
		{
			patch->shiftPolarOrigin(atan(1.0));
			ivector1D list = patch->includeTriangulation(&geo);
			geo.plotExtraction(folder, "includeSphericalMeshANGULARSHIFTED", true, 18, list);
			out<<"...test includeTriangulation method...see output file"<<endl;
			std::string folder = ".";
			patch->plotGrid(folder,"sphericalMeshANGULARSHIFTED",19,true);
		}
		
		delete patch; patch = NULL;
	}	
}else{
	cout<<"Not able to open log file"<<endl;
}

out.close();	
return(0);

} //end program 
