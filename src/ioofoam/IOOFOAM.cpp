/*---------------------------------------------------------------------------*\
 *
 *  mimmo
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of mimmo.
 *
 *  mimmo is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  mimmo is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with mimmo. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/
#include "IOOFOAM.hpp"
#include <vtkSmartPointer.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkTriangleFilter.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkPolyDataWriter.h>

using namespace std;
using namespace bitpit;
namespace mimmo{

/*!Default constructor of IOOFOAM.
 */
IOOFOAM::IOOFOAM(){
    m_name          = "mimmo.IOOFOAM";
    setDefaults();
}

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
IOOFOAM::IOOFOAM(const bitpit::Config::Section & rootXML){

    m_name = "mimmo.IOOFOAM";
    setDefaults();

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::trim(input);
    if(input == "mimmo.IOOFOAM"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*!Default destructor of IOOFOAM.
 */
IOOFOAM::~IOOFOAM(){};

/*!Copy constructor of IOOFOAM.
 */
IOOFOAM::IOOFOAM(const IOOFOAM & other):BaseManipulation(){
    *this = other;
};

/*!Assignement operator of IOOFOAM.
 */
IOOFOAM & IOOFOAM::operator=(const IOOFOAM & other){
    m_read = other.m_read;
    m_rfilenameS= other.m_rfilenameS;
    m_write = other.m_write;
    m_wfilenameV = other.m_wfilenameV;
    m_rdirS = other.m_rdirS;
    m_wdirV = other.m_wdirV;
    m_surfmesh_ext = other.m_surfmesh_ext;
    return *this;
};

/*!Default values for IOOFOAM.
 */
void
IOOFOAM::setDefaults(){

    m_read       = false;
    m_rfilenameS.clear();
    m_rdirS.clear();
    m_write      = false;
    m_wfilenameV = "mimmoPoints";
    m_wdirV      = "./";
    m_surfmesh_ext = NULL;
    m_maxf      = 0.0;
    m_scaling   = 1.0;
    m_normalize = true;
}

/*! It builds the input/output ports of the object
 */
void
IOOFOAM::buildPorts(){
    bool built = true;
    built = (built && createPortIn<MimmoObject*, IOOFOAM>(this, &IOOFOAM::setGeometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
    built = (built && createPortIn<MimmoObject*, IOOFOAM>(this, &IOOFOAM::setSurfaceBoundary, PortType::M_GEOM2, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
    built = (built && createPortIn<dvector1D, IOOFOAM>(this, &IOOFOAM::setField, PortType::M_SCALARFIELD, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::FLOAT));

    built = (built && createPortOut<MimmoObject*, IOOFOAM>(this, &IOOFOAM::getGeometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
    built = (built && createPortOut<MimmoObject*, IOOFOAM>(this, &IOOFOAM::getSurfaceBoundary, PortType::M_GEOM2, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
    built = (built && createPortOut<dvector1D, IOOFOAM>(this, &IOOFOAM::getField, PortType::M_SCALARFIELD, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::FLOAT));
    m_arePortsBuilt = built;
};

/*!It sets the condition to read the geometries (both surface and volume patches) on file during the execution.
 * \param[in] read Does it read the geometry in execution?
 */
void
IOOFOAM::setRead(bool read){
    m_read = read;
}

/*!It sets the name of directory to read the geometry (surface patches).
 * \param[in] dir Vector of names of input directories.
 */
void
IOOFOAM::setVTKReadDir(vector<string> dir){
    m_rdirS = dir;
}

/*!It sets the name of file to read the geometry (surface patches).
 * \param[in] filename Vector of names of input files.
 */
void
IOOFOAM::setVTKReadFilename(vector<string> filename){
    m_rfilenameS = filename;
}

/*!It adds a name of directory to read a geometry (surface patch).
 * \param[in] dir Name of input directory.
 */
void
IOOFOAM::addVTKReadDir(string dir){
    m_rdirS.push_back(dir);
}

/*!It adds a name of file to read a geometry (surface patch).
 * \param[in] filename Name of input file.
 */
void
IOOFOAM::addVTKReadFilename(string filename){
    m_rfilenameS.push_back(filename);
}

/*!It sets the name of directory to read the points cloud (volume patch).
 * \param[in] dir Name of input directory.
 */
void
IOOFOAM::setPointsReadDir(string dir){
    m_rdirV = dir;
}

/*!It sets the name of file to read the points cloud (volume patch).
 * \param[in] filename Name of input file.
 */
void
IOOFOAM::setPointsReadFilename(string filename){
    m_rfilenameV = filename;
}


/*!It sets the condition to write the geometry (volume points cloud) on file during the execution.
 * \param[in] write Does it write the geometry in execution?
 */
void
IOOFOAM::setWrite(bool write){
    m_write = write;
}

/*!It sets the name of directory to write the geometry (surface patch).
 * \param[in] dir Name of directory.
 */
void
IOOFOAM::setVTKWriteDir(string dir){
    m_wdirS = dir;
}

/*!It sets the name of file to write the geometry (surface patch).
 * \param[in] filename Name of output file.
 */
void
IOOFOAM::setVTKWriteFilename(string filename){
    m_wfilenameS = filename;
}

/*!It sets the name of directory to write the geometry (volume points cloud).
 * \param[in] dir Name of directory.
 */
void
IOOFOAM::setPointsWriteDir(string dir){
    m_wdirV = dir;
}

/*!It sets the name of file to write the geometry (volume points cloud).
 * \param[in] filename Name of output file.
 */
void
IOOFOAM::setPointsWriteFilename(string filename){
    m_wfilenameV = filename;
}


/*!
 * Set current geometry to an external points cloud mesh.
 */
void
IOOFOAM::setGeometry(MimmoObject * geom){
    if(geom == NULL || m_read)   return;
    if(geom->isEmpty())  return;
    if(geom->getType() != 3) return;

    BaseManipulation::setGeometry(geom);
    //m_volmesh.reset(nullptr); //TODO clean local geometry if present?
}

/*!
 * Set boundary surface relative to the volume mesh.Option active only in writing mode.
 * Pre-existent boundary surfaces read from file, when class is set in read mode,  will be erased.
 * \param[in] geosurf pointer to surface boundary mesh
 */
void
IOOFOAM::setSurfaceBoundary(MimmoObject* geosurf){
    if(geosurf == NULL || m_read)   return;
    if(geosurf->isEmpty())  return;
    if(geosurf->getType() != 1) return;

    m_surfmesh_ext = geosurf;
    //m_surfmesh.reset(nullptr); //TODO clean local geometry if present?
};

/*!It sets the scaling factor to be applied to the scalar field
 * related to surface boundary mesh.
 * \param[in] scaling Scaling factor.
 */
void
IOOFOAM::setScaling(double scaling){
    m_scaling = scaling;
}

/*!It sets the scalar field
 * related to surface boundary mesh.
 * \param[in] field Scalar field.
 */
void
IOOFOAM::setField(dvector1D field){
    m_field = field;
}

/*!It sets if the scalar field has to be normalized with its maximum absolute value.
 * \param[in] normalize Has the scalar field to be normalized?
 */
void
IOOFOAM::setNormalize(bool normalize){
    m_normalize = normalize;
}

/*!
 * Return all the surface bounding the current volume mesh.
 * If reading mode is active returns info contained in vtk, marking with internal PID all the possible boundary patches,
 * otherwise refers to the actual object pointed by the User.
 * \return pointer to surface boundary mesh
 */
MimmoObject*
IOOFOAM::getSurfaceBoundary(){
    if(m_read)  return m_surfmesh.get();
    else    return m_surfmesh_ext;
}


/*!
 * Return current pointer to geometry.If read mode is active return local read volumetric mesh, else
 * otherwise return pointed externally geometry
 * \return pointer to volume points cloud mesh
 */
MimmoObject*
IOOFOAM::getGeometry(){
    if(m_read) return m_volmesh.get();
    else return BaseManipulation::getGeometry();
}

/*!It gets the scalar field
 * related to surface boundary mesh.
 * \return Scalar field.
 */
dvector1D
IOOFOAM::getField(){
    return (m_field);
}

/*!It reads the mesh geometries from input file.
 * It reads even the scalar fields trelated to surface patches if they are present.
 * \return False if files don't exist or are not a polydata (surface) or OpenFOAM points format (volume).
 */
bool
IOOFOAM::read(){


    //Read OpenFOAM Points
    {

        std::ifstream infile(m_rdirV+"/"+m_rfilenameV);
        bool check = infile.good();
        if (!check){
            m_stopat = SHRT_MAX;
            return false;
        }
        infile.close();

        dvecarr3E   Ipoints;
        readOFP(m_rdirV, m_rfilenameV, Ipoints);

        //Reverse info in your grids.
        std::unique_ptr<MimmoObject> patchVol(new MimmoObject(3));

        int sizeV = Ipoints.size();
        patchVol->getPatch()->reserveVertices(sizeV);

        //  Stock vertices in volume grid
        for(auto & vv : Ipoints)    patchVol->addVertex(vv);

        //release the points cloud mesh
        m_volmesh = std::move(patchVol);

    }

    //compute kdtree for points cloud
    m_volmesh->buildKdTree();


    //read vtk patches and assign ordered PID starting from 0
    //TODO input PID for each patch has to be implemented.

    std::unique_ptr<MimmoObject> patchBnd(new MimmoObject(1));

    short nPID = m_rfilenameS.size();
    bool check;
    for (short iPID = 0; iPID < nPID; iPID++){
        check = readVTK(m_rdirS[iPID], m_rfilenameS[iPID], iPID, patchBnd.get());
        if (!check) return check;
    }

    if (m_normalize && m_maxf > 0.0){
        for (auto vertex : patchBnd->getVertices()){
            m_field[patchBnd->getMapDataInv(vertex.getId())] /= m_maxf;
            m_field[patchBnd->getMapDataInv(vertex.getId())] *= m_scaling;
        }
    }
    else{
        for (auto vertex : patchBnd->getVertices()){
            m_field[patchBnd->getMapDataInv(vertex.getId())] *= m_scaling;
        }
    }

    //release the surface mesh
    m_surfmesh = std::move(patchBnd);

    return true;

}

/*!It writes the mesh geometries on output .vtk (surface patche) and points format (points cloud) files.
 * It unifies the polydata input files with a unique surface patch.
 *\return False if one geometry is empty.
 */
bool
IOOFOAM::write(){

    if (getGeometry() == NULL && getSurfaceBoundary() == NULL){
        m_stopat = 2;
        return false;
    }
    if (getSurfaceBoundary() != NULL){
        //write surface mesh
        string outputFilename = m_wdirS+"/"+m_wfilenameS;
        getSurfaceBoundary()->getPatch()->write(outputFilename);
    }

    if (getGeometry() != NULL){
        //write point cloud mesh
        writeOFP(m_wdirV, m_wfilenameV, getGeometry()->getVertices());
    }
    return true;
}


//===============================//
//====== OFOAM INTERFACE ========//
//===============================//

/*!
 *  Read openFoam format geometry file and absorb it as a point cloud ONLY.
 *\param[in]    inputDir    folder of file
 *\param[in]    pointsName  name of file
 *\param[out]   points      list of points in the cloud
 *
 */
void IOOFOAM::readOFP(string& inputDir, string& pointsName, dvecarr3E& points){

    ifstream is(inputDir +"/"+pointsName);

    points.clear();
    int ip = 0;
    int np;
    darray3E point;
    string sread;
    char par;

    for (int i=0; i<18; i++){
        getline(is,sread);
    }
    is >> np;
    getline(is,sread);
    getline(is,sread);

    points.resize(np);
    while(!is.eof() && ip<np){
        is.get(par);
        for (int i=0; i<3; i++) is >> point[i];
        is.get(par);
        getline(is,sread);
        points[ip] = point;
        ip++;
    }
    is.close();
    return;

}

/*!
 *  Write geometry file in openFoam format as a point cloud ONLY.
 *\param[in]    outputDir    folder of file
 *\param[in]    pointsName  name of file
 *\param[out]   vertices    list of points in the cloud
 *
 */
void IOOFOAM::writeOFP(string& outputDir, string& pointsName, PiercedVector<Vertex>& vertices){

    ofstream os(outputDir +"/"+pointsName);
    char nl = '\n';

    string separator(" ");
    string parl("("), parr(")");
    string hline;

    hline = "/*--------------------------------*- C++ -*----------------------------------*\\" ;
    os << hline << nl;
    hline = "| =========                 |                                                 |";
    os << hline << nl;
    hline = "| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |";
    os << hline << nl;
    hline = "|  \\\\    /   O peration     | Version:  2.4.x                                 |";
    os << hline << nl;
    hline = "|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |";
    os << hline << nl;
    hline = "|    \\\\/     M anipulation  |                                                 |";
    os << hline << nl;
    hline = "\\*---------------------------------------------------------------------------*/";
    os << hline << nl;
    hline = "FoamFile";
    os << hline << nl;
    hline = "{";
    os << hline << nl;
    hline = "    version     2.0;";
    os << hline << nl;
    hline = "    format      ascii;";
    os << hline << nl;
    hline = "    class       vectorField;";
    os << hline << nl;
    hline = "    location    \"constant/polyMesh\";";
    os << hline << nl;
    hline = "    object      points;";
    os << hline << nl;
    hline = "}";
    os << hline << nl;
    hline = "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //";
    os << hline << nl;
    os << nl;
    os << nl;
    int np = vertices.size();
    os << np << nl;
    os << parl << nl;
    for (auto vertex : vertices){
        os << parl;
        for (int j=0; j<2; j++){
            os << setprecision(16) << vertex[j] << separator;
        }
        os << setprecision(16) << vertex[2] << parr << nl;
    }
    os << parr << nl;
    os << nl;
    os << nl;
    hline = "// ************************************************************************* //";
    os << hline << nl;

    os.close();
}

//===============================//
//===============================//


/*!It reads a vtk surface patch from file. The structure is stored in a MimmoObject
 * by matching the IDs of the vertices wth the vertices stored in the local volume
 * points cloud patch (we are in reading modality).
 *\param[in]    inputDir    folder of file
 *\param[in]    surfaceName name of file
 *\param[in]    PID         PID of the patch
 *\param[in]    patchBnd    Actual boundary surface patch to be filled
 */
bool IOOFOAM::readVTK(string& inputDir, string& surfaceName, short PID, MimmoObject* patchBnd){


    std::ifstream infile(inputDir+"/"+surfaceName+".vtk");
    bool check = infile.good();
    if (!check) return false;
    infile.close();


    string inputFilename = inputDir+"/"+surfaceName+".vtk";
    //int np = 0;
    int nt = 0;
    darray3E point;
    bitpit::Vertex vertex;

    //mapper for connectivity
    map<vtkIdType, long> mapID;
    map<vtkIdType, long> mapid;

    // Get all data from the file
    vtkSmartPointer<vtkGenericDataObjectReader> reader =
            vtkSmartPointer<vtkGenericDataObjectReader>::New();
    reader->SetFileName(inputFilename.c_str());
    reader->Update();

    // All of the standard data types can be checked and obtained like this:
    if(reader->IsFilePolyData())
    {
        vtkPolyData* output = reader->GetPolyDataOutput();

        vtkSmartPointer<vtkTriangleFilter> tri= vtkSmartPointer<vtkTriangleFilter>::New();
        tri->SetInputData(output);
        tri->Update();

        vtkPolyData* output2 = tri->GetOutput();

        vtkSmartPointer<vtkPoints> points = output2->GetPoints();
        vtkCellArray *cells = output2->GetPolys();
        //vtkCellData *cdata = output2->GetCellData();
        vtkPointData *pdata = output2->GetPointData();


        bitpit::KdTree<3, bitpit::Vertex, long> * kdtree = m_volmesh->getKdTree();
        long ID;

        double point_[3], tol = 1.0e-08;
        livector1D ids;
        livector1D noids;
        bool check;

        for (vtkIdType id=0; id<points->GetNumberOfPoints(); id++ ){
            points->GetPoint(id, point_);
            for (int i=0; i<3; i++) point[i] = point_[i];
            vertex.setCoords(point);
            check = false;
            while(!check){
                ids.clear();
                kdtree->hNeighbors(&vertex, tol, &ids, &noids);
                if (ids.size() == 0){
                    tol *= 1.5;
                }
                else if (ids.size() > 1){
                    tol /= 1.5;
                }
                else{
                    ID = ids[0];
                    check = true;
                }
            }

            if (!(patchBnd->getVertices().exists(ID))){
                patchBnd->addVertex(point, ID);
            }
            mapID[id] = ID;
            mapid[id] = patchBnd->getMapDataInv(ID);
        }

        vtkIdType* conn_;
        vtkIdType npts;
        nt = cells->GetNumberOfCells();
        livector1D connectivity(3);
        bitpit::ElementInfo::Type eltype = bitpit::ElementInfo::TRIANGLE;
        cells->InitTraversal();
        for (vtkIdType id=0; id<nt; id++ ){
            cells->GetNextCell(npts, conn_);
            for (int i=0; i<npts; i++) connectivity[i] = mapID[conn_[i]];
            patchBnd->addConnectedCell(connectivity, eltype, PID);
        }

        patchBnd->cleanGeometry();


        //IOOFOAM reads only a scalar field on each patch if present (if not add 0 values)
        //If m_normalize is active the data are normalized with the maximum value on the patches
        vtkDataArray* data = pdata->GetArray(0);
        m_field.resize(patchBnd->getNVertex());
        if (data != NULL){
            for (vtkIdType id=0; id<points->GetNumberOfPoints(); id++){
                m_field[mapid[id]] = (data->GetComponent(id,0));
                m_maxf = max(m_maxf, abs(m_field[mapid[id]]));
            }
        }
        else{
            for (vtkIdType id=0; id<points->GetNumberOfPoints(); id++){
                m_field[mapid[id]] = 0.0;
            }
        }
    }
    else{
        (*m_log) << m_name << " error: polydata not found in : "<< inputFilename << std::endl;
        m_stopat = PID;
        return false;
    }

    return true;
}



/*!Execution command.
 * It reads the geometry if the condition m_read is true.
 * It writes the geometry if the condition m_write is true.
 */
void
IOOFOAM::execute(){
    bool check = true;
    if (m_read) check = read();
    if (!check){
        if (m_stopat == SHRT_MAX){
            (*m_log) << m_name << " error: file not found : "<< m_rfilenameV << std::endl;
            (*m_log) << " " << std::endl;
            throw std::runtime_error (m_name + ": file not found : " + m_rfilenameV);
        }
        (*m_log) << m_name << " error: file not found : "<< m_rfilenameS[m_stopat] << std::endl;
        (*m_log) << " " << std::endl;
        throw std::runtime_error (m_name + ": file not found : " + m_rfilenameS[m_stopat]);
    }
    if (m_write) check = write();
    if (!check){
        if (m_stopat == 2){
            (*m_log) << m_name << " error: write not done : surface and volume geometry not linked " << std::endl;
            (*m_log) << " " << std::endl;
            throw std::runtime_error (m_name + ": write not done : surface and volume geometry not linked ");
        }
    }
}


/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
IOOFOAM::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    std::string input;

    if(slotXML.hasOption("Priority")){
        input = slotXML.get("Priority");
        int value =0;
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::trim(input));
            ss>>value;
        }
        setPriority(value);
    };

    if(slotXML.hasOption("ReadFlag")){
        input = slotXML.get("ReadFlag");
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::trim(input));
            ss >> value;
        }
        setRead(value);
    };

    std::vector<std::string> temp_dirs;
    if(slotXML.hasSection("VTKReadDirs")){

        const bitpit::Config::Section & filesXML = slotXML.getSection("VTKReadDirs");

        for(auto & subfile : filesXML.getSections()){
            std::string dir;

            if(subfile.second->hasOption("dir"))   {
                dir = subfile.second->get("dir");
            }

            if(!dir.empty()){
                temp_dirs.push_back(dir);
            }
        }
        setVTKReadDir(temp_dirs);
    }


    std::vector<std::string> temp_files;
    if(slotXML.hasSection("VTKReadFilenames")){

        const bitpit::Config::Section & filesXML = slotXML.getSection("VTKReadFilenames");

        for(auto & subfile : filesXML.getSections()){
            std::string file;

            if(subfile.second->hasOption("filename"))   {
                file = subfile.second->get("filename");
            }

            if(!file.empty()){
                temp_files.push_back(file);
            }
        }
        setVTKReadFilename(temp_files);
    }

    if(slotXML.hasOption("PointsReadDir")){
        input = slotXML.get("PointsReadDir");
        input = bitpit::utils::trim(input);
        if(input.empty())   input = "./";
        setPointsReadDir(input);
    };


    if(slotXML.hasOption("PointsReadFilename")){
        input = slotXML.get("PointsReadFilename");
        input = bitpit::utils::trim(input);
        if(input.empty())   input = "./";
        setPointsReadFilename(input);
    };


    if(slotXML.hasOption("WriteFlag")){
        input = slotXML.get("WriteFlag");
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::trim(input));
            ss >> value;
        }
        setWrite(value);
    };


    if(slotXML.hasOption("PointsWriteDir")){
        input = slotXML.get("PointsWriteDir");
        input = bitpit::utils::trim(input);
        if(input.empty())   input = "./";
        setPointsWriteDir(input);
    };


    if(slotXML.hasOption("PointsWriteFilename")){
        input = slotXML.get("PointsWriteFilename");
        input = bitpit::utils::trim(input);
        if(input.empty())   input = "./";
        setPointsWriteFilename(input);
    };


    if(slotXML.hasOption("VTKWriteDir")){
        input = slotXML.get("VTKWriteDir");
        input = bitpit::utils::trim(input);
        if(input.empty())   input = "./";
        setVTKWriteDir(input);
    };


    if(slotXML.hasOption("VTKWriteFilename")){
        input = slotXML.get("VTKWriteFilename");
        input = bitpit::utils::trim(input);
        if(input.empty())   input = "./";
        setVTKWriteFilename(input);
    };


    if(slotXML.hasOption("Normalize")){
        input = slotXML.get("Normalize");
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::trim(input));
            ss >> value;
        }
        setNormalize(value);
    };


    if(slotXML.hasOption("Scaling")){
        input = slotXML.get("Scaling");
        double value = 1.0;
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::trim(input));
            ss >> value;
        }
        setScaling(value);
    };

};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
IOOFOAM::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    slotXML.set("ClassName", m_name);
    slotXML.set("Priority", std::to_string(getPriority()));


    std::string output;

    output = std::to_string(m_read);
    slotXML.set("ReadFlag", output);

    {
        bitpit::Config::Section & filesXML = slotXML.addSection("VTKReadDirs");

        int counter = 0;
        for(auto & dir : m_rdirS){
            std::string name = "dir"+std::to_string(counter);
            bitpit::Config::Section & local = filesXML.addSection(name);
            local.set("dir", dir);
            ++counter;
        }
    }

    {
        bitpit::Config::Section & filesXML = slotXML.addSection("VTKReadFilenames");

        int counter = 0;
        for(auto & file : m_rfilenameS){
            std::string name = "file"+std::to_string(counter);
            bitpit::Config::Section & local = filesXML.addSection(name);
            local.set("filename", file);
            ++counter;
        }
    }

    slotXML.set("PointsReadDir", m_rdirV);
    slotXML.set("PointsReadFilename", m_rfilenameV);

    output = std::to_string(m_write);
    slotXML.set("WriteFlag", output);

    slotXML.set("VTKWriteDir", m_wdirS);
    slotXML.set("VTKWriteFilename", m_wfilenameS);
    slotXML.set("PointsWriteDir", m_wdirV);
    slotXML.set("PointsWriteFilename", m_wfilenameV);

    output = std::to_string(m_normalize);
    slotXML.set("Normalize", output);

    if (m_scaling != 1.0){
        std::stringstream ss;
        ss<<std::scientific<<m_scaling;
        slotXML.set("Scaling", ss.str());
    }

};


}
