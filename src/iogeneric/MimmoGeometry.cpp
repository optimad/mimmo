/*---------------------------------------------------------------------------*\
 *
 *  mimmo
 *
 *  Copyright (C) 2015-2021 OPTIMAD engineering Srl
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
#include "MimmoGeometry.hpp"
#include "VTUGridReader.hpp"
#include "VTUGridWriterASCII.hpp"
#include <iostream>

namespace mimmo {

/*!
Default constructor of MimmoGeometry.
Require to specify mode of the block as reader, writer or converter.
NOTE: in case of converter the directory and filename will be the same
for input & output
(the I/O extensions have to be set by setReadFileType & setWriteFileType).

\param[in] mode Modality of the block -READ/WRITE/CONVERT
*/
MimmoGeometry::MimmoGeometry(MimmoGeometry::IOMode mode){
    m_name         = "mimmo.Geometry";
    setDefaults();
    setIOMode(mode);
}

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
MimmoGeometry::MimmoGeometry(const bitpit::Config::Section & rootXML){

    m_name = "mimmo.Geometry";
    setDefaults();

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::string::trim(input);

    std::string fallback_name2 = "READ";
    std::string input2 = rootXML.get("IOMode", fallback_name2);
    input2 = bitpit::utils::string::trim(input2);

    setIOMode(IOMode::READ);
    if(input2 == "WRITE") setIOMode(IOMode::WRITE);
    if(input2 == "CONVERT") setIOMode(IOMode::CONVERT);

    if(input == "mimmo.Geometry"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*!Default destructor of MimmoGeometry.
 */
MimmoGeometry::~MimmoGeometry(){
    clear();
};

/*!Copy constructor of MimmoGeometry.If to-be-copied object has an
 * internal MimmoObject instantiated, it will be soft linked
 * in the current class (its pointer only will be copied in m_geometry member);
 */
MimmoGeometry::MimmoGeometry(const MimmoGeometry & other):BaseManipulation(other){
    m_rinfo = other.m_rinfo;
    m_winfo = other.m_winfo;
    m_read = other.m_read;
    m_write = other.m_write;
    m_wformat = other.m_wformat;
    m_codex = other.m_codex;
    m_buildSkdTree = other.m_buildSkdTree;
    m_buildKdTree = other.m_buildKdTree;
    m_refPID = other.m_refPID;
    m_multiSolidSTL = other.m_multiSolidSTL;
    m_tolerance = other.m_tolerance;
    m_clean = other.m_clean;
};

/*!
 * Assignement operator of MimmoGeometry.If to-be-copied object has an
 * internal MimmoObject instantiated, it will be soft linked
 * in the current class (its pointer only will be copied in m_geometry member);
 */
MimmoGeometry & MimmoGeometry::operator=(MimmoGeometry other){
    swap(other);
    return *this;
};

/*!
 * Swap method.
 * \param[in] x object to be swapped
 */
void MimmoGeometry::swap(MimmoGeometry & x) noexcept
{
    std::swap(m_rinfo, x.m_rinfo);
    std::swap(m_winfo, x.m_winfo);
    std::swap(m_read, x.m_read);
    std::swap(m_write, x.m_write);
    std::swap(m_wformat, x.m_wformat);
    std::swap(m_codex, x.m_codex);
    std::swap(m_buildSkdTree, x.m_buildSkdTree);
    std::swap(m_buildKdTree, x.m_buildKdTree);
    std::swap(m_refPID, x.m_refPID);
    std::swap(m_multiSolidSTL, x.m_multiSolidSTL);
    std::swap(m_tolerance, x.m_tolerance);
    std::swap(m_clean, x.m_clean);
    BaseManipulation::swap(x);
}

/*!
 * Building the ports available in the class
 */
void
MimmoGeometry::buildPorts(){
    bool built = true;
    bool mandatory_input = (m_write && !m_read);
    built = (built && createPortIn<MimmoSharedPointer<MimmoObject>, MimmoGeometry>(this, &mimmo::BaseManipulation::setGeometry, M_GEOM, mandatory_input));
    built = (built && createPortOut<MimmoSharedPointer<MimmoObject>, MimmoGeometry>(this, &mimmo::BaseManipulation::getGeometry, M_GEOM));
    m_arePortsBuilt = built;
}

/*!
 * Return a const pointer to your current class.
 */
const MimmoGeometry * MimmoGeometry::getCopy(){
    return this;
}

/*!
 * Set proper member of the class to defaults
 */
void
MimmoGeometry::setDefaults(){

    m_read            = false;
    m_rinfo.fname    = "mimmoGeometry";
    m_rinfo.ftype    = static_cast<int>(FileType::STL);
    m_rinfo.fdir    = "./";

    m_write            = false;
    m_winfo.fname    = "mimmoGeometry";
    m_winfo.ftype    = static_cast<int>(FileType::STL);
    m_winfo.fdir    = "./";

    m_wformat        = Short;
    m_codex            = true;
    m_buildSkdTree    = false;
    m_buildKdTree    = false;
    m_refPID = 0;
    m_multiSolidSTL = false;
    m_tolerance = 1.0e-06;
    m_clean = true;
    m_parallelRestore = MIMMO_ENABLE_MPI;
}


/*!It sets the condition to read the geometry on file during the execution.
 * \param[in] read Does it read the geometry in execution?
 */
void
MimmoGeometry::_setRead(bool read){
    m_read = read;
}

/*!It sets the type of file to read the geometry during the execution.
 * \param[in] type Extension of file.
 * NOTE: use it only for converter mode, otherwise setFileType is recommended
 */
void
MimmoGeometry::setReadFileType(FileType type){
    m_rinfo.ftype = type._to_integral();
}

/*!It sets the type of file to read the geometry during the execution.
 * \param[in] type Label of file type (0 = bynary STL).
 * NOTE: use it only for converter mode, otherwise setFileType is recommended
 */
void
MimmoGeometry::setReadFileType(int type){
    auto maybe_type = FileType::_from_integral_nothrow(type);
    if(!maybe_type) type = 0;
    m_rinfo.ftype = type;
}

/*!It sets the name of directory to read the geometry.
 * \param[in] dir Name of directory.
 * NOTE: use it only for converter mode, otherwise setDir is recommended
 */
void
MimmoGeometry::setReadDir(std::string dir){
    m_rinfo.fdir = dir;
}

/*!It sets the name of file to read the geometry.
 * \param[in] filename Name of input file.
 * NOTE: use it only for converter mode, otherwise setFilename is recommended
 */
void
MimmoGeometry::setReadFilename(std::string filename){
    m_rinfo.fname = filename;
}

/*!It sets the type of file to write the geometry during the execution.
 * \param[in] type Extension of file.
 * NOTE: use it only for converter mode, otherwise setFileType is recommended
 */
void
MimmoGeometry::setWriteFileType(FileType type){
    m_winfo.ftype = type._to_integral();
}

/*!It sets the type of file to write the geometry during the execution.
 * \param[in] type Label of file type (0 = bynary STL).
 * NOTE: use it only for converter mode, otherwise setFileType is recommended
 */
void
MimmoGeometry::setWriteFileType(int type){
    auto maybe_type = FileType::_from_integral_nothrow(type);
    if(!maybe_type) type = 0;
    m_winfo.ftype = type;
}


/*!It sets the condition to write the geometry on file during the execution.
 * \param[in] write Does it write the geometry in execution?
 */
void
MimmoGeometry::_setWrite(bool write){
    m_write = write;
}

/*!It sets the name of directory to write the geometry.
 * \param[in] dir Name of directory.
 * NOTE: use it only for converter mode, otherwise setDir is recommended
 */
void
MimmoGeometry::setWriteDir(std::string dir){
    m_winfo.fdir = dir;
}

/*!It sets the name of file to write the geometry.
 * \param[in] filename Name of output file.
 * NOTE: use it only for converter mode, otherwise setFilename is recommended
 */
void
MimmoGeometry::setWriteFilename(std::string filename){
    m_winfo.fname = filename;

}

/*!Set the mode of the block as reader, writer or converter.
 * \param[in] mode Modality of the block
 * NOTE: in case of converter the directory and filename will be the same
 * for input & output
 * (the I/O extensions have to be set by setReadFileType & setWriteFileType).
 */
void
MimmoGeometry::setIOMode(MimmoGeometry::IOMode mode){
    switch(mode){
        case MimmoGeometry::IOMode::READ :
            _setRead();
            _setWrite(false);
            break;
        case MimmoGeometry::IOMode::WRITE :
            _setRead(false);
            _setWrite();
            break;
        case MimmoGeometry::IOMode::CONVERT :
            _setRead();
            _setWrite();
            break;
        default:
            //never been reached
            break;
        }
}

/*!Set the mode of the block as reader, writer or converter.
 * \param[in] mode Modality of the block
 * NOTE: in case of converter the directory and filename will be the same
 * for input & output
 * (the I/O extensions have to be set by setReadFileType & setWriteFileType).
 */
void
MimmoGeometry::setIOMode(int mode){
    setIOMode(static_cast<MimmoGeometry::IOMode>(mode));
}

/*!It sets the name of directory to read/write the geometry.
 * \param[in] dir Name of directory.
 */
void
MimmoGeometry::setDir(std::string dir){
    m_rinfo.fdir = dir;
    m_winfo.fdir = dir;
}

/*!It sets the name of file to read/write the geometry.
 * \param[in] filename Name of output file.
 */
void
MimmoGeometry::setFilename(std::string filename){
    m_winfo.fname = filename;
    m_rinfo.fname = filename;
}

/*!It sets the type of file to read/write the geometry during the execution.
 * \param[in] type Extension of file.
 * NOTE: in case of converter mode use separate setReadFileType and setWriteFileType
 */
void
MimmoGeometry::setFileType(FileType type){
    m_winfo.ftype = type._to_integral();
    m_rinfo.ftype = type._to_integral();
}

/*!It sets the type of file to read/write the geometry during the execution.
 * \param[in] type Label of file type (0 = bynary STL).
 * NOTE: in case of converter mode use separate setReadFileType and setWriteFileType
 */
void
MimmoGeometry::setFileType(int type){
    auto maybe_type = FileType::_from_integral_nothrow(type);
    if(!maybe_type) type = 0;
    m_winfo.ftype = type;
    m_rinfo.ftype = type;
}

/*!
 * set codex ASCII false, BINARY true for writing sessions ONLY.
 * Default is Binary/Appended. Pay attention, binary writing is effective
 * only those file formats which support it.(ex STL, STVTU, SQVTU, VTVTU, VHVTU)
 * \param[in] binary codex flag.
 */
void MimmoGeometry::setCodex(bool binary){
    m_codex = binary;
}

/*!
 * Set ASCII Multi Solid STL writing. The method has effect only while writing
 * standard surface triangulations in STL format type.
 * \param[in] multi boolean, true activate Multi Solid STL writing.
 */
void MimmoGeometry::setMultiSolidSTL(bool multi){
    m_multiSolidSTL = multi;
}

/*!
 * Set geometric tolerance used to perform geometric operations on the mimmo object.
 * param[in] tol input geometric tolerance
 */
void
MimmoGeometry::setTolerance(double tol){
	m_tolerance = std::max(1.0e-15, tol);
}

/*!
 * Set if the geometry has to cleaned after reading.
 * param[in] clean cleaning flag
 */
void
MimmoGeometry::setClean(bool clean){
	m_clean = clean;
}

/*!
 * Set if the geometry to read is parallel or not.
 * In case of reading a geometry dumped in mimmo (*.geomimmo) format
 * this function has to be called in order to restore correctly a parallel patch.
 * param[in] parallelRestore is parallel flag
 */
void
MimmoGeometry::setParallelRestore(bool parallelRestore){
    m_parallelRestore = parallelRestore;
}

/*!
 * Force your class to allocate an internal MimmoObject of type 1-Superficial mesh
 * 2-Volume Mesh,3-Point Cloud, 4-3DCurve. Other internal object allocated or externally linked geometries
 * will be destroyed/unlinked.
 * \param[in] type 1-Surface MimmoObject, 2-Volume MimmoObject, 3-Point Cloud, 4-3DCurve. Default is 1, no other type are supported
 */
void
MimmoGeometry::setGeometry(int type){
    if(type > 4)    type = 1;
    int type_ = std::max(type,1);
    m_geometry.reset(new MimmoObject(type_));
};

/*!
 * Return a pointer to the Vertex structure of the MimmoObject geometry actually pointed or allocated by
 * the class. If no geometry is actually available return a nullptr
 * \return pointer to the vertices structure
 */
bitpit::PiercedVector<bitpit::Vertex> *
MimmoGeometry::getVertices(){
    return &(getGeometry()->getVertices());
};

/*!
 * Return a pointer to the Cell structure of the MimmoObject geometry actually pointed or allocated by
 * the class. If no geometry is actually available return a nullptr
 * \return pointer to the cells structure
 */
bitpit::PiercedVector<bitpit::Cell> * MimmoGeometry::getCells(){
    return    &(getGeometry()->getCells());

};

/*!It sets the PIDs of all the cells of the geometry Patch.
 * \param[in] pids PIDs of the cells of geometry mesh, in compact sequential order. If pids size does not match number of current cell does nothing
 */
void
MimmoGeometry::setPID(livector1D pids){
    getGeometry()->setPID(pids);
};

/*!It sets the PIDs of part of/all the cells of the geometry Patch.
 * \param[in] pidsMap PIDs map list w/ id of the cell as first value, and pid as second value.
 */
void
MimmoGeometry::setPID(std::unordered_map<long, long> pidsMap){
    getGeometry()->setPID(pidsMap);
};

/*!
 * Set a reference PID for the whole geometry cell. If the geometry is already pidded, all existent pid's
 * will be translated w.r.t. the value of the reference pid. Default value is 0. Maximum value allowed is 30000.
 * The effect of reference pidding will be available only after the execution of the class. After execution, reference PID will
 * be reset ot zero.
 * \param[in] pid reference pid to assign to the geometry cells
 */
void
MimmoGeometry::setReferencePID(long pid){
    m_refPID = std::max(long(0),pid);
}


/*!It sets if the SkdTree of the patch has to be built during execution.
 * \param[in] build If true the SkdTree is built in execution and stored in
 * the related MimmoObject member.
 */
void
MimmoGeometry::setBuildSkdTree(bool build){
    m_buildSkdTree = build;
}

/*!It sets if the KdTree of the patch has to be built during execution.
 * \param[in] build If true the KdTree is built in execution and stored in
 * the related MimmoObject member.
 */
void
MimmoGeometry::setBuildKdTree(bool build){
    m_buildKdTree = build;
}

/*!
 * Check if geometry is not linked or not locally instantiated in your class.
 * True - no geometry present, False otherwise.
 * \return is the geometry empty?
 */
bool MimmoGeometry::isEmpty(){
    if(getGeometry() == nullptr)   return true;
    return getGeometry()->isEmpty();
}

/*!
 * Clear all stuffs in your class
 */
void
MimmoGeometry::clear(){
    setDefaults();
    BaseManipulation::clear();
};

/*!It sets the format to export .nas file.
 * \param[in] wform Format of .nas file (Short/Long).
 */
void
MimmoGeometry::setFormatNAS(WFORMAT wform){
    m_wformat = wform;
}

/*!It writes the mesh geometry on output .vtu file.
 *\return False if geometry is not linked.
 */
bool
MimmoGeometry::write(){

    //adjusting with reference pid;
    {
        auto locpids = getGeometry()->getCompactPID();
        auto pidsMap = getGeometry()->getPIDTypeListWNames();
        for( auto & val: locpids) val+=m_refPID;
        setPID(locpids);
        auto pids = getGeometry()->getPIDTypeList();
        for( auto & val: pids) getGeometry()->setPIDName(val, pidsMap[val - m_refPID]);
        //action completed, reset m_refPID to zero.
        m_refPID = 0;
    }

    //check writing filetype
    if (!FileType::_from_integral_nothrow(m_winfo.ftype)){
        return false;
    }

    switch(FileType::_from_integral(m_winfo.ftype)){

    case FileType::STL :
        //Export STL
    {
        auto pidsMap = getGeometry()->getPIDTypeListWNames();
        std::unordered_map<int, std::string> mpp;
        for(const auto & touple : pidsMap){
            mpp[touple.first] = touple.second;
        }
        std::string name = (m_winfo.fdir+"/"+m_winfo.fname+".stl");
        dynamic_cast<bitpit::SurfUnstructured*>(getGeometry()->getPatch())->exportSTL(name, m_codex, m_multiSolidSTL, &mpp);
        return true;
    }
    break;


    case FileType::SURFVTU :
    case FileType::VOLVTU :
    case FileType::CURVEVTU:
    case FileType::PCVTU:
        //Export Surface/Volume/3DCurve VTU
    {
        if(!m_codex){
            VTUFlushStreamerASCII streamer;
            VTUGridWriterASCII vtkascii(streamer, *(getGeometry()->getPatch()) );
            vtkascii.write(m_winfo.fdir+"/", m_winfo.fname);
        }
        else{
            getGeometry()->getPatch()->getVTK().setCodex(bitpit::VTKFormat::APPENDED);
            getGeometry()->getPatch()->getVTK().setDirectory(m_winfo.fdir+"/");
            getGeometry()->getPatch()->getVTK().setName(m_winfo.fname);
            getGeometry()->getPatch()->write();
        }
        return true;
    }
    break;

    case FileType::NAS :
        //Export Nastran file - Id of vertex and cell cannot be 0 or negative. So all
        // nas vertex/cell id etiquette are incremented by 1;
        // same for PID: if 0 exists, pids are incremented by 1; if -1 exists are incremented by 2;
        // Beware if id >= 10^8 nas format does not support it.
    {
        lilimap mapDataInv;
        dvecarr3E    points = getGeometry()->getVerticesCoords();
        livector1D   pointsID;
        pointsID.reserve(points.size());
        bitpit::PiercedVector<bitpit::Vertex> & vertices = getGeometry()->getVertices();
        long id;
        for(bitpit::Vertex & vv : vertices){
            id = vv.getId();
            mapDataInv[id] = id+1;
            pointsID.push_back(id+1);
        }
        //return compactconnectivity using vertex map
        livector2D   connectivity = getGeometry()->getCompactConnectivity(mapDataInv);
        livector1D   elementsID;
        elementsID.reserve(connectivity.size());
        for(bitpit::Cell & cell : getGeometry()->getCells()){
            elementsID.push_back(cell.getId()+1);
        }
        NastranInterface nastran;
        nastran.setWFormat(m_wformat);
        livector1D pids = getGeometry()->getCompactPID();
        std::unordered_set<long>  pidsset = getGeometry()->getPIDTypeList();
        // pid cannot be negative or 0 in NAS.
        if(pidsset.count(-1) > 0 || pidsset.count(0) > 0){
            std::unordered_set<long> temp;
            long offset = 1+ long(pidsset.count(-1)>0);
            for(auto & val : pidsset){
                temp.insert(val+offset);
            }
            std::swap(temp, pidsset);
            for(long & val :pids){
                val+=offset;
            }
        }
        std::string namefile = m_winfo.fname;
#if MIMMO_ENABLE_MPI
        // Only master rank 0 writes on file
        if (getRank() == 0)
#endif
        {
            if (pids.size() == connectivity.size()){
                nastran.write(m_winfo.fdir,namefile,points, pointsID, connectivity,elementsID, &pids, &pidsset);
            }else{
                nastran.write(m_winfo.fdir,namefile,points, pointsID, connectivity, elementsID);
            }
        }
        return true;
    }
    break;
//
//    case FileType::PCVTU :
//        //Export Point Cloud VTU
//    {
//        dvecarr3E    points = getGeometry()->getVerticesCoords();
//        ivector2D    connectivity(points.size(), ivector1D(1,0));
//        int counter = 0;
//        for(auto & c:connectivity){
//            c[0] = counter;
//            counter++;
//        }
//        bitpit::VTKUnstructuredGrid  vtk(m_winfo.fdir, m_winfo.fname, bitpit::VTKElementType::VERTEX);
//        vtk.setGeomData( bitpit::VTKUnstructuredField::POINTS, points) ;
//        vtk.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, connectivity) ;
//        vtk.setDimensions(connectivity.size(), points.size());
//        if(!m_codex)    vtk.setCodex(bitpit::VTKFormat::ASCII);
//        else            vtk.setCodex(bitpit::VTKFormat::APPENDED);
//        vtk.write() ;
//        return true;
//    }
//    break;

    case FileType::MIMMO :
    	//Export in mimmo (bitpit) dump format
    {
    	int archiveVersion = 1;
    	std::string header = m_name;
    	std::string filename = (m_winfo.fdir+"/"+m_winfo.fname);
#if MIMMO_ENABLE_MPI
    	bitpit::OBinaryArchive binaryWriter(filename, "geomimmo", archiveVersion, header, getRank());
#else
    	bitpit::OBinaryArchive binaryWriter(filename, "geomimmo", archiveVersion, header);
#endif
    	getGeometry()->dump(binaryWriter.getStream());
    	binaryWriter.close();
    	return true;
    }
    break;

    default: //never been reached
        break;
    }

    return false;
};

/*!It reads the mesh geometry from an input file and reverse it in the internal
 * MimmoObject container. If an external container is linked skip reading and doing nothing.
 * \return False if file doesn't exists or not found geometry container address.
 */
bool
MimmoGeometry::read(){

    //check reading filetype
    if (!FileType::_from_integral_nothrow(m_rinfo.ftype)){
        return false;
    }

    switch(FileType::_from_integral(m_rinfo.ftype)){

    //Import STL
    case FileType::STL :
    {
        setGeometry(1);
        std::string name;
#if MIMMO_ENABLE_MPI
        if (getRank() == 0) {
#endif
            {
                name = m_rinfo.fdir+"/"+m_rinfo.fname+".stl";
        		bool check = fileExist(name);
        		if (!check){
                    name = m_rinfo.fdir+"/"+m_rinfo.fname+".STL";
        			check = fileExist(name);
        			if (!check) return false;
        		}
        	}

        	std::ifstream in(name);
        	std::string    ss, sstype;
        	std::getline(in,ss);
        	std::stringstream ins;
        	ins << ss;
        	ins >> sstype;
        	bool binary = true;
        	if (sstype == "solid" || sstype == "SOLID") binary = false;
        	in.close();

        	std::unordered_map<int,std::string> mapPIDSolid;
        	dynamic_cast<bitpit::SurfUnstructured*>(getGeometry()->getPatch())->importSTL(name, binary, 0, false, &mapPIDSolid);

        	//count PID if multi-solid
        	auto & mapset = getGeometry()->getPIDTypeList();
        	for(const auto & cell : getGeometry()->getCells() ){
        		mapset.insert(cell.getPID());
        	}
        	std::unordered_map<long, std::string> & mapWNames = getGeometry()->getPIDTypeListWNames();
        	for(const auto & val :mapset ){
        		mapWNames[val] = mapPIDSolid[val];
        	}
#if MIMMO_ENABLE_MPI
        }
#endif
    }
    break;

    case FileType::SURFVTU :
        //Import Surface VTU meshes
    {
        setGeometry(1);
        bool masterRankOnly = true;
        std::string name = m_rinfo.fdir+"/"+m_rinfo.fname;
        std::string extension = ".vtu";
#if MIMMO_ENABLE_MPI
        if (getGeometry()->getProcessorCount() > 1){
            //check for pvtu
            extension = ".pvtu";
            if(fileExist(name+extension)){
                masterRankOnly = false;
            }else{
                //fallback to normal vtu.
                extension = ".vtu";
            }
        }
#endif

        if (!fileExist(name+extension)) return false;

        VTUGridStreamer vtustreamer;
        VTUGridReader  input(m_rinfo.fdir, m_rinfo.fname, vtustreamer, *(getGeometry()->getPatch()), masterRankOnly);

        input.read() ;

        getGeometry()->resyncPID();
    }
    break;

    case FileType::VOLVTU :
        //Import Volume VTU meshes
    {
        setGeometry(2);
        bool masterRankOnly = true;
        std::string name = m_rinfo.fdir+"/"+m_rinfo.fname;
        std::string extension = ".vtu";
#if MIMMO_ENABLE_MPI
        if (getGeometry()->getProcessorCount() > 1){
            //check for pvtu
            extension = ".pvtu";
            if(fileExist(name+extension)){
                masterRankOnly = false;
            }else{
                //fallback to normal vtu.
                extension = ".vtu";
            }
        }
#endif

        if (!fileExist(name+extension)) return false;

        VTUGridStreamer vtustreamer;
        VTUGridReader  input(m_rinfo.fdir, m_rinfo.fname, vtustreamer, *(getGeometry()->getPatch()), masterRankOnly);

        input.read() ;

        getGeometry()->resyncPID();
    }
    break;

    case FileType::NAS :
        //Import Surface NAS
    {
        setGeometry(1);
#if MIMMO_ENABLE_MPI
    	if (getRank() == 0) {
#endif
        bool check = fileExist(m_rinfo.fdir+"/"+m_rinfo.fname+".nas");
        if (!check) return false;

        dvecarr3E    Ipoints ;
        livector2D    Iconnectivity ;
        livector1D    pointsID;
        livector1D    cellsID;

        NastranInterface nastran;
        nastran.setWFormat(m_wformat);

        livector1D pids;
        nastran.read(m_rinfo.fdir, m_rinfo.fname, Ipoints, pointsID, Iconnectivity, cellsID, pids );

        bitpit::ElementType eltype;

        int sizeV, sizeC;
        sizeV = Ipoints.size();
        sizeC = Iconnectivity.size();
        getGeometry()->getPatch()->reserveVertices(sizeV);
        getGeometry()->getPatch()->reserveCells(sizeC);

        int counter = 0;
        for(const auto & vv : Ipoints){
            getGeometry()->addVertex(vv, pointsID[counter]);
            ++counter;
        }
        counter = 0;
        for(const auto & cc : Iconnectivity)    {
            eltype = bitpit::ElementType::UNDEFINED;
            std::size_t ccsize = cc.size();
            if(ccsize == 1)  eltype = bitpit::ElementType::VERTEX;
            if(ccsize == 2)  eltype = bitpit::ElementType::LINE;
            if(ccsize == 3)  eltype = bitpit::ElementType::TRIANGLE;
            if(ccsize == 4)  eltype = bitpit::ElementType::QUAD;
            if(ccsize > 4)   eltype = bitpit::ElementType::POLYGON;

            getGeometry()->addConnectedCell(cc, eltype, cellsID[counter]);
            ++counter;
        }
        getGeometry()->setPID(pids);
#if MIMMO_ENABLE_MPI
    	}
#endif
    }
    break;

    //Import point cloud from vtu
    case FileType::PCVTU :
    {
        setGeometry(3);
        bool masterRankOnly = true;
        std::string name = m_rinfo.fdir+"/"+m_rinfo.fname;
        std::string extension = ".vtu";
#if MIMMO_ENABLE_MPI
        if (getGeometry()->getProcessorCount() > 1){
            //check for pvtu
            extension = ".pvtu";
            if(fileExist(name+extension)){
                masterRankOnly = false;
            }else{
                //fallback to normal vtu.
                extension = ".vtu";
            }
        }
#endif

        if (!fileExist(name+extension)) return false;

        VTUGridStreamer vtustreamer;
        VTUGridReader  input(m_rinfo.fdir, m_rinfo.fname, vtustreamer, *(getGeometry()->getPatch()), masterRankOnly);

        input.read() ;
        getGeometry()->resyncPID();
    }
    break;

    case FileType::CURVEVTU :
        //Import 3D Curve VTU
    {
        setGeometry(4);
        bool masterRankOnly = true;
        std::string name = m_rinfo.fdir+"/"+m_rinfo.fname;
        std::string extension = ".vtu";
#if MIMMO_ENABLE_MPI
        if (getGeometry()->getProcessorCount() > 1){
            //check for pvtu
            extension = ".pvtu";
            if(fileExist(name+extension)){
                masterRankOnly = false;
            }else{
                //fallback to normal vtu.
                extension = ".vtu";
            }
        }
#endif

        if (!fileExist(name+extension)) return false;

        VTUGridStreamer vtustreamer;
        VTUGridReader  input(m_rinfo.fdir, m_rinfo.fname, vtustreamer, *(getGeometry()->getPatch()), masterRankOnly);

        input.read() ;

        getGeometry()->resyncPID();
    }
    break;

    case FileType::MIMMO :
    	//Import in mimmo (bitpit) restore format
    {
    	std::string filename = (m_rinfo.fdir+"/"+m_rinfo.fname);
#if MIMMO_ENABLE_MPI
    	bitpit::IBinaryArchive binaryReader(filename,"geomimmo", getRank());
#else
    	bitpit::IBinaryArchive binaryReader(filename, "geomimmo");
#endif

    	// Reset to a generic geometry with the correct parallel propriety
        getGeometryReference().reset(new MimmoObject(0, m_parallelRestore));
        m_geometry->restore(binaryReader.getStream());
    	binaryReader.close();
    }
    break;

    default: //never been reached
        break;

    }

    // set Tolerance and clean geometry
	getGeometry()->setTolerance(m_tolerance);
	if (m_clean){
		getGeometry()->cleanGeometry();
	}

    //adjusting with reference pid;
    auto locpids = getGeometry()->getCompactPID();
    auto pidsMap = getGeometry()->getPIDTypeListWNames();
    for( auto & val: locpids) val+=m_refPID;
    setPID(locpids);
    auto pids = getGeometry()->getPIDTypeList();
    for( auto & val: pids) getGeometry()->setPIDName(val, pidsMap[val - m_refPID]);

    //action completed, reset m_refPID to zero.
    m_refPID = 0;

    // Set all the sync status variables to unsync
    getGeometry()->setUnsyncAll();

    return true;
};

/*!
    Return true if a file exists on the filesystem and it is readable
    (check successfull opening with ifstream)
    \param[in] filename full path to file to be checked
    \ return true/false if the file does/does-not exist.
*/
bool
MimmoGeometry::fileExist(const std::string & filename){
    std::ifstream tryfile(filename);
    bool check = tryfile.good();
    tryfile.close();
    return  check;
}

/*!Execution command.
 * It reads the geometry if the condition m_read is true.
 * It writes the geometry if the condition m_write is true.
 */
void
MimmoGeometry::execute(){
    bool check = true;
    // Read geometry
    if (m_read){
        check = read();
        if (!check){
            (*m_log) << m_name << " error: file not found : "<< m_rinfo.fname << std::endl;
            (*m_log) << " " << std::endl;
            throw std::runtime_error (m_name + " : file not found : " + m_rinfo.fname);
        }

        // Update geometry
        getGeometry()->update();

        if (m_buildSkdTree) getGeometry()->buildSkdTree();
        if (m_buildKdTree) getGeometry()->buildKdTree();
    }
    check = true;
    // Write geometry
    if (m_write){
        check = write();
        if (!check){
            (*m_log) << m_name << " error: write not done : geometry not linked " << std::endl;
            (*m_log) << " " << std::endl;
            return;
        }
    }
}

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
MimmoGeometry::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    std::string input;

    BaseManipulation::absorbSectionXML(slotXML, name);

    if(slotXML.hasOption("FileType")){
        input = slotXML.get("FileType");
        input = bitpit::utils::string::trim(input);
        if(!input.empty()){
            for(auto c: FileType::_values()){
                if(input == c._to_string()){
                    setFileType(c);
                }
            }
        }else{
            setFileType(0);
        }
    };

    if(slotXML.hasOption("ReadFileType")){
        input = slotXML.get("ReadFileType");
        input = bitpit::utils::string::trim(input);
        if(!input.empty()){
            for(auto c: FileType::_values()){
                if(input == c._to_string()){
                    setReadFileType(c);
                }
            }
        }
    };

    if(slotXML.hasOption("WriteFileType")){
        input = slotXML.get("WriteFileType");
        input = bitpit::utils::string::trim(input);
        if(!input.empty()){
            for(auto c: FileType::_values()){
                if(input == c._to_string()){
                    setWriteFileType(c);
                }
            }
        }
    };

    if(slotXML.hasOption("Dir")){
        input = slotXML.get("Dir");
        input = bitpit::utils::string::trim(input);
        if(input.empty())    input = "./";
        setDir(input);
    };


    if(slotXML.hasOption("Filename")){
        input = slotXML.get("Filename");
        input = bitpit::utils::string::trim(input);
        if(input.empty())    input = "mimmoGeometry";
        setFilename(input);
    };

    if(slotXML.hasOption("ReadDir")){
        input = slotXML.get("ReadDir");
        input = bitpit::utils::string::trim(input);
        if(!input.empty())
            setReadDir(input);
    };


    if(slotXML.hasOption("ReadFilename")){
        input = slotXML.get("ReadFilename");
        input = bitpit::utils::string::trim(input);
        if(!input.empty())
            setReadFilename(input);
    };

    if(slotXML.hasOption("WriteDir")){
        input = slotXML.get("WriteDir");
        input = bitpit::utils::string::trim(input);
        if(!input.empty())
            setWriteDir(input);
    };


    if(slotXML.hasOption("WriteFilename")){
        input = slotXML.get("WriteFilename");
        input = bitpit::utils::string::trim(input);
        if(!input.empty())
            setWriteFilename(input);
    };

    if(slotXML.hasOption("Codex")){
        input = slotXML.get("Codex");
        bool value = true;
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::string::trim(input));
            ss >> value;
        }
        setCodex(value);
    };

    //to be deprecated BvTree field option. Substituted by SkdTree.
    if(slotXML.hasOption("BvTree")){
        input = slotXML.get("BvTree");
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::string::trim(input));
            ss >> value;
        }
        setBuildSkdTree(value);
    };

    if(slotXML.hasOption("SkdTree")){
        input = slotXML.get("SkdTree");
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::string::trim(input));
            ss >> value;
        }
        setBuildSkdTree(value);
    };


    if(slotXML.hasOption("KdTree")){
        input = slotXML.get("KdTree");
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::string::trim(input));
            ss >> value;
        }
        setBuildKdTree(value);
    };

    if(slotXML.hasOption("AssignRefPID")){
        input = slotXML.get("AssignRefPID");
        long value = 0;
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::string::trim(input));
            ss >> value;
        }
        setReferencePID(value);
    };

    if(slotXML.hasOption("WriteMultiSolidSTL")){
        input = slotXML.get("WriteMultiSolidSTL");
        bool value = true;
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::string::trim(input));
            ss >> value;
        }
        setMultiSolidSTL(value);
    };

	if(slotXML.hasOption("Tolerance")){
		input = slotXML.get("Tolerance");
		double value = 1.0e-06;
		if(!input.empty()){
			std::stringstream ss(bitpit::utils::string::trim(input));
			ss >> value;
			setTolerance(value);
		}
	};

    if(slotXML.hasOption("Clean")){
        input = slotXML.get("Clean");
        bool value = true;
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::string::trim(input));
            ss >> value;
        }
        setClean(value);
    };

    if(slotXML.hasOption("FormatNAS")){
        input = slotXML.get("FormatNAS");
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::string::trim(input));
            ss >> value;
        }
        if(value){
            setFormatNAS(WFORMAT::Long);
        }else{
            setFormatNAS(WFORMAT::Short);
        }
    };

    if(slotXML.hasOption("ParallelRestore")){
        input = slotXML.get("ParallelRestore");
        bool value = true;
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::string::trim(input));
            ss >> value;
        }
        setParallelRestore(value);
    };

};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
MimmoGeometry::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    BaseManipulation::flushSectionXML(slotXML, name);

    std::string output;

    {
        std::string temp = "READ";
        if(m_read && m_write)   temp = "CONVERT";
        if(!m_read && m_write)  temp = "WRITE";
        slotXML.set("IOMode", temp);
    }

    if (m_read && m_write)
    {
        slotXML.set("ReadDir", m_rinfo.fdir);
        slotXML.set("ReadFilename", m_rinfo.fname);
        slotXML.set("WriteDir", m_winfo.fdir);
        slotXML.set("WriteFilename", m_winfo.fname);
    }
    else{
        slotXML.set("Dir", m_rinfo.fdir);
        slotXML.set("Filename", m_rinfo.fname);
    }

    if (m_read && m_write)
    {
        std::string temp = (FileType::_from_integral(m_rinfo.ftype))._to_string();
        slotXML.set("ReadFileType", temp);
        temp = (FileType::_from_integral(m_winfo.ftype))._to_string();
        slotXML.set("WriteFileType", temp);
    }
    else
    {
        std::string temp = (FileType::_from_integral(m_rinfo.ftype))._to_string();
        slotXML.set("FileType", temp);
    }

    output = std::to_string(m_codex);
    slotXML.set("Codex", output);

    output = std::to_string(m_buildSkdTree);
    slotXML.set("SkdTree", output);

    output = std::to_string(m_buildKdTree);
    slotXML.set("KdTree", output);
    slotXML.set("AssignRefPID", std::to_string(m_refPID));
    slotXML.set("WriteMultiSolidSTL", std::to_string(m_multiSolidSTL));

	std::stringstream ss;
	ss<<std::scientific<<m_tolerance;
	slotXML.set("Tolerance", ss.str());

    output = std::to_string(m_clean);
    slotXML.set("Clean", output);
    slotXML.set("FormatNAS", std::to_string(bool(m_wformat == WFORMAT::Long)));

    output = std::to_string(m_parallelRestore);
    slotXML.set("ParallelRestore", output);

};



//===============================//
//====== NASTRAN INTERFACE ======//
//===============================//

/*!
 * Default constructor
 */
NastranInterface::NastranInterface(){
    // Enable grid points and surface elements by default
    enable(NastranElementType::GRID);
    enable(NastranElementType::CTRIA);
    enable(NastranElementType::CQUAD);
    disable(NastranElementType::CBAR);
    disable(NastranElementType::RBE2);
    disable(NastranElementType::RBE3);
}

/*!
 * Set the format Short/Long of data in your nastran file
 * \param[in] wf WFORMAT Short/Long of your nastran file
 */
void NastranInterface::setWFormat(WFORMAT wf){
    m_wformat = wf;
}

/*!
 * Write a keyword string according to WFORMAT chosen by the User
 * \param[in]     key    keyword
 * \param[in,out] os    ofstream where the keyword is written
 */
void NastranInterface::writeKeyword(std::string key, std::ofstream& os){
    os.setf(std::ios_base::left);
    switch (m_wformat)
    {
    case Short:
    {
        os << std::setw(8) << key;
        break;
    }
    case Long:
    {
        os << std::setw(8) << std::string(key + '*');
        break;
    }
    }
    os.unsetf(std::ios_base::left);
}

/*!
 * Write a 3D point to an ofstream
 * \param[in]     p         point
 * \param[in]       pointI    point label
 * \param[in,out] os        ofstream where the point is written
 */
void NastranInterface::writeCoord(const darray3E& p, const long& pointI, std::ofstream& os){
    // Fixed short/long formats:
    // 1 GRID
    // 2 ID : point ID - requires starting index of 1
    // 3 CP : co-ordinate system ID (blank)
    // 4 X1 : point x cp-ordinate
    // 5 X2 : point x cp-ordinate
    // 6 X3 : point x cp-ordinate
    // 7 CD : co-ordinate system for displacements (blank)
    // 8 PS : single point constraints (blank)
    // 9 SEID : super-element ID

    std::string separator_("");
    writeKeyword(std::string("GRID"), os);
    os << separator_;
    os.setf(std::ios_base::right);
    writeValue(pointI, os);
    os << separator_;
    writeValue("", os);
    os << separator_;
    writeValue(p[0], os);
    os << separator_;
    writeValue(p[1], os);
    os << separator_;
    switch (m_wformat)
    {
    case Short:
    {
        writeValue(p[2], os);
        os << nl;
        //        os << setw(8) << p[2]
        //                           << nl;
        os.unsetf(std::ios_base::right);
        break;
    }
    case Long:
    {
        os << nl;
        os.unsetf(std::ios_base::right);
        writeKeyword("", os);
        os.setf(std::ios_base::right);
        writeValue(p[2], os);
        os << nl;
        break;
    }
    default:
    {
        throw std::runtime_error ("NastranInterface : Unknown writeFormat enumeration");
    }
    }
    os.unsetf(std::ios_base::right);
}

/*!
 * Write a face to an ofstream
 * \param[in]     faceType     string reporting label fo the type of element
 * \param[in]       facePts    Element Points connectivity
 * \param[in]       nFace        element label
 * \param[in,out] os        ofstream where the face is written
 * \param[in]     PID       part identifier associated to the element
 */
void NastranInterface::writeFace(std::string faceType, const livector1D& facePts, const long& nFace, std::ofstream& os,long PID){
    // Only valid surface elements are CTRIA3 and CQUAD4

    // Fixed short/long formats:
    // 1 CQUAD4
    // 2 EID : element ID
    // 3 PID : property element ID; default = EID (blank)
    // 4 G1 : grid point index - requires starting index of 1
    // 5 G2 : grid point index
    // 6 G3 : grid point index
    // 7 G4 : grid point index
    // 8 onwards - not used

    // For CTRIA3 elements, cols 7 onwards are not used

    WFORMAT wformat_ = m_wformat;
    m_wformat = Short;
    std::string separator_("");
    writeKeyword(faceType, os);
    os << separator_;
    os.setf(std::ios_base::right);
    writeValue(nFace, os);
    os << separator_;
    writeValue(PID, os);
    int fp1;
    switch (m_wformat)
    {
    case Short:
    {
        for (int i=0; i< (int)facePts.size(); i++)
        {
            fp1 = facePts[i];
            writeValue(fp1, os);
        }

        break;
    }
    case Long:
    {
        for (int i=0; i< (int)facePts.size(); i++)
        {
            fp1 = facePts[i];
            writeValue(fp1, os);
            //            if (i == 1)
            //            {
            //                os << nl;
            //                os.unsetf(std::ios_base::right);
            //                writeKeyword("", os);
            //                os.setf(std::ios_base::right);
            //            }
        }
        break;
    }
    default:
    {
        std::cout << "Unknown writeFormat enumeration" << std::endl;
        throw std::runtime_error ("Unknown writeFormat enumeration");
    }
    }
    os << nl;
    os.unsetf(std::ios_base::right);
    m_wformat = wformat_;
}

/*!
 * Write a face to an std::ofstream
 * \param[in]     points     list of vertices
 * \param[in]     pointsID   list of vertices labels
 * \param[in]     faces        element points connectivity
 * \param[in]     facesID   list of element labels
 * \param[in,out] os        ofstream where the geometry is written
 * \param[in]     PIDS      list of Part Identifiers for each cells (optional)
 * \return true if the geometry is correctly written.
 */
bool NastranInterface::writeGeometry(dvecarr3E& points, livector1D& pointsID, livector2D& faces, livector1D & facesID,
                                     std::ofstream& os, livector1D* PIDS){

    if(points.size() != pointsID.size())    return false;
    if(faces.size() != facesID.size())    return false;
    bool flagpid = (PIDS != nullptr);
    if(flagpid && (PIDS->size() != faces.size())) return false;
    int counter =  0;
    for (const auto & p: points){
        writeCoord(p, pointsID[counter], os);
        ++counter;
    }

    counter = 0;
    for (const auto & f: faces){
        long PID = 1;
        if (flagpid) PID = (*PIDS)[counter];
        if (f.size() == 3){
            writeFace("CTRIA3", f, facesID[counter], os, PID);
            ++counter;
        }
        else if (f.size() == 4){
            writeFace("CQUAD4", f, facesID[counter], os, PID);
            ++counter;
        }else{
            std::cerr << "Warning NastranInterface::writeGeometry : unknown face format" << std::endl;
        }
    }

    return true;
}

/*!
 * Write nastran footer to an ofstream
 * \param[in,out] os        ofstream where the footer is written
 * \param[in]     PIDSSET   list of overall available part identifiers
 */
void NastranInterface::writeFooter(std::ofstream& os, std::unordered_set<long>* PIDSSET){
    std::string separator_("");
    if (PIDSSET == nullptr){
        int PID = 1;
        writeKeyword("PSHELL", os);
        os << separator_;
        writeValue(PID, os);
        for (int i = 0; i < 6; i++)
        {
            // Dummy values
            long uno = 1;
            os << separator_;
            writeValue(uno, os);
        }
        os << nl;
    }
    else{
        for (std::unordered_set<long>::iterator it = PIDSSET->begin(); it != PIDSSET->end(); it++)
        {
            writeKeyword("PSHELL", os);
            int PID = (*it);
            os << separator_;
            writeValue(PID, os);
            for (int i = 0; i < 6; i++)
            {
                // Dummy values
                long uno = 1;
                os << separator_;
                writeValue(uno, os);
            }
            os << nl;

        }
    }

    writeKeyword("MAT1", os);
    os << separator_;
    int MID = 1;
    writeValue(MID, os);
    for (int i = 0; i < 7; i++)
    {
        // Dummy values
        os << separator_;
        writeValue("", os);
    }
    os << nl;
}

/*!
 * Write a bdf nastran file
 * \param[in] outputDir    output directory
 * \param[in] surfaceName    output filename
 * \param[in] points    reference of a point container that has to be written
 * \param[in] pointsID  reference of a label point container
 * \param[in] faces     reference to element-point connectivity that has to be written
 * \param[in] facesID   reference to a label element container
 * \param[in] PIDS        reference to long int vector for Part Identifiers that need to be associated to elements
 * \param[in] PIDSSET    map of overall available Part Identifiers
 */
void NastranInterface::write(std::string& outputDir, std::string& surfaceName, dvecarr3E& points, livector1D & pointsID,
                             livector2D& faces, livector1D & facesID, livector1D* PIDS, std::unordered_set<long>* PIDSSET){

    std::ofstream os(outputDir +"/"+surfaceName + ".nas");
    os << "TITLE=mimmo " << surfaceName << " mesh" << nl
            << "$" << nl
            << "BEGIN BULK" << nl;
    bool check = writeGeometry(points, pointsID,  faces, facesID, os, PIDS);
    if(!check){
        throw std::runtime_error ("Error in NastranInterface::write : geometry mesh cannot be written");
    }
    writeFooter(os, PIDSSET);
    os << "ENDDATA" << std::endl;
}

//========READ====//
/*!
 * Read a bdf nastran file
 * \param[in] inputDir    input directory
 * \param[in] surfaceName    input filename
 * \param[out] points    reference of a point container that has to be filled
 * \param[out] pointsID  reference to a label point container to be filled
 * \param[out] faces     reference to element-point connectivity that has to be filled
 * \param[out] facesID  reference to a label element container to be filled
 * \param[out] PIDS        reference to long int vector for Part Identifier storage
 */
void NastranInterface::read(std::string& inputDir, std::string& surfaceName, dvecarr3E& points, livector1D & pointsID,
                            livector2D& faces, livector1D & facesID, livector1D& PIDS){

    std::ifstream is(inputDir +"/"+surfaceName + ".nas");

    points.clear();
    pointsID.clear();
    faces.clear();
    facesID.clear();
    PIDS.clear();

    long ipoint = 0;
    long iface  = 0;
    darray3E point;
    livector1D face;
    long pid;
    std::string sread;
    std::string ssub = trim(sread.substr(0,8));

    while(!is.eof()){
        while(((ssub != "GRID" && ssub != "GRID*") || !isEnabled(NastranElementType::GRID)) &&
                ((ssub != "CTRIA3" && ssub != "CTRIA3*") || !isEnabled(NastranElementType::CTRIA)) &&
                ((ssub != "CQUAD4" && ssub != "CQUAD4*") || !isEnabled(NastranElementType::CQUAD)) &&
                (ssub != "RBE2" || !isEnabled(NastranElementType::RBE2)) &&
                (ssub != "RBE3" || !isEnabled(NastranElementType::RBE3)) &&
                (ssub != "CBAR" || !isEnabled(NastranElementType::CBAR)) &&
                !is.eof()){
            std::getline(is,sread);
            ssub = trim(sread.substr(0,8));
        }
        if(ssub == "GRID"){
            ipoint = stoi(sread.substr(8,8));
            point[0] = stod(convertVertex(trim(sread.substr(24,8))));
            point[1] = stod(convertVertex(trim(sread.substr(32,8))));
            point[2] = stod(convertVertex(trim(sread.substr(40,8))));
            points.push_back(point);
            pointsID.push_back(ipoint);
            std::getline(is,sread);
            ssub = trim(sread.substr(0,8));
        }
        else if(ssub == "GRID*"){
            ipoint = stoi(sread.substr(16,16));
            point[0] = stod(convertVertex(trim(sread.substr(48,16))));
            point[1] = stod(convertVertex(trim(sread.substr(64,16))));
            point[2] = stod(convertVertex(trim(sread.substr(80,16))));
            points.push_back(point);
            pointsID.push_back(ipoint);
            std::getline(is,sread);
            ssub = trim(sread.substr(0,8));
        }
        else if(ssub == "CTRIA3"){
            face.resize(3);
            iface = stoi(sread.substr(8,8));
            pid   = stoi(sread.substr(16,8));
            face[0] = stoi(sread.substr(24,8));
            face[1] = stoi(sread.substr(32,8));
            face[2] = stoi(sread.substr(40,8));
            faces.push_back(face);
            facesID.push_back(iface);
            PIDS.push_back(pid);
            std::getline(is,sread);
            ssub = trim(sread.substr(0,8));
        }
        else if(ssub == "CTRIA3*"){
            face.resize(3);
            iface = stoi(sread.substr(16,16));
            pid   = stoi(sread.substr(32,16));
            face[0] = stoi(sread.substr(48,16));
            face[1] = stoi(sread.substr(64,16));
            face[2] = stoi(sread.substr(80,16));
            faces.push_back(face);
            facesID.push_back(iface);
            PIDS.push_back(pid);
            std::getline(is,sread);
            ssub = trim(sread.substr(0,8));
        }
        else if(ssub == "CQUAD4"){
            face.resize(4);
            iface = stoi(sread.substr(8,8));
            pid = stoi(sread.substr(16,8));
            face[0] = stoi(sread.substr(24,8));
            face[1] = stoi(sread.substr(32,8));
            face[2] = stoi(sread.substr(40,8));
            face[3] = stoi(sread.substr(48,8));
            faces.push_back(face);
            facesID.push_back(iface);
            PIDS.push_back(pid);
            std::getline(is,sread);
            ssub = trim(sread.substr(0,8));
        }
        else if(ssub == "CQUAD4*"){
            face.resize(4);
            iface = stoi(sread.substr(16,16));
            pid = stoi(sread.substr(32,16));
            face[0] = stoi(sread.substr(48,16));
            face[1] = stoi(sread.substr(64,16));
            face[2] = stoi(sread.substr(80,16));
            face[3] = stoi(sread.substr(96,16));
            faces.push_back(face);
            facesID.push_back(iface);
            PIDS.push_back(pid);
            std::getline(is,sread);
            ssub = trim(sread.substr(0,8));
        }
        else if(ssub == "CBAR"){
            face.resize(2);
            iface = stoi(sread.substr(8,8));
            pid = stoi(sread.substr(16,8));
            face[0] = stoi(sread.substr(24,8));
            face[1] = stoi(sread.substr(32,8));
            faces.push_back(face);
            facesID.push_back(iface);
            PIDS.push_back(pid);
            std::getline(is,sread);
            ssub = trim(sread.substr(0,8));
        }
        else if(ssub == "RBE3"){

            long ID;
            std::vector<long> conn;
            long PID;
            std::string components;

            std::vector<std::string> records;
            bool toexit = false;
            while (!toexit){
                appendLineRecords(sread, 8, records);
                getline(is,sread);
                ssub = sread.substr(0,8);
                if ((ssub != "        " && ssub[0] != '+'))
                    toexit = true;
            }

            absorbRBE3(records, ID, PID, conn, components);

            //If surface elements build polygon
            std::size_t nv = conn.size();

            int ispolygon = int(nv > 4);
            face.resize(nv+ispolygon);
            face[0] = nv;
            for (std::size_t i=0; i<nv; i++){
                face[i+ispolygon] = conn[i];
            }
            faces.push_back(face);
            facesID.push_back(ID);
            PIDS.push_back(PID);

        }
        else if(ssub == "RBE2"){

            long ID;
            std::vector<long> conn;
            long PID;
            std::string components;

            std::vector<std::string> records;
            bool toexit = false;
            while (!toexit){
                appendLineRecords(sread, 8, records);
                getline(is,sread);
                ssub = sread.substr(0,8);
                if ((ssub != "        " && ssub[0] != '+')){
                    toexit = true;
                }
            }

            absorbRBE2(records, ID, PID, conn, components);

            // build polygon
            std::size_t nv = conn.size();

            int ispolygon = int(nv > 4);
            face.resize(nv+ispolygon);
            face[0] = nv;
            for (std::size_t i=0; i<nv; i++){
                face[i+ispolygon] = conn[i];
            }
            faces.push_back(face);
            facesID.push_back(ID);
            PIDS.push_back(PID);

        }

    }
    is.close();

}

/*!
 * Custom trimming of nas string
 * \return trimmed string
 */
std::string
NastranInterface::trim(std::string in){

    return bitpit::utils::string::trim(in);

/*
 * Old implementation for CentOS 5 system

    std::stringstream out;
    out << in;
    out >> in;
    return in;
*/

}

/*!
 * Convert nas string of numerical floats
 * \return string of regular floats
 */
std::string
NastranInterface::convertVertex(std::string in){
    int pos = in.find_last_of("-");
    if (pos<(int)in.size() && pos > 0){
        in.insert(pos, "E");
    }
    else{
        pos = in.find_last_of("+");
        if (pos<(int)in.size() && pos > 0){
            in.insert(pos, "E");
        }
    }
    return in;
}

/*!
 * Append records of given length of an input line in a records structure
 * \param[in] sread input line given as string
 * \param[in] recordlength length of single record string
 * \param[out] records vector of string records to be filled
 */
void
NastranInterface::appendLineRecords(std::string & sread, std::size_t recordlength, std::vector<std::string> & records)
{
    std::size_t readlength = sread.size();
    std::size_t readpos = 0;
    while (readpos < readlength){
        std::string record = sread.substr(readpos, recordlength);
        record = bitpit::utils::string::trim(record);
        // Do not append empty records or blank spaces or new line tag (ANSA format, starting with character '+')
        if (!record.empty() && record != "        " && record[0] != '+'){
            records.push_back(record);
        }
        readpos += recordlength;
    }
}

/*!
 * Absorb rbe3 element from a records structure
 * \param[in] records vector of string records
 * \param[out] ID ID of the RBE3 element
 * \param[out] PID PID of the RBE3 element
 * \param[out] connectivity connectivity of the RBE3 element
 * \param[out] components reaction components of the element
 */
void
NastranInterface::absorbRBE3(std::vector<std::string> & records, long & ID, long & PID, std::vector<long> & connectivity, std::string & components)
{
    connectivity.clear();
    ID = std::stoi(records[1]);
    records[2] = bitpit::utils::string::trim(records[2]);
    if(!records[2].empty()){
        PID = std::stoi(records[2]);
    }
    else{
        PID = 0;
    }
    components = records[3];
    std::size_t startpos = 4;
    std::size_t endpos = records.size();
    for (std::size_t pos = startpos; pos < endpos; pos++){
        // If it is not integer it is the string UM or the weight factor
        // It can be even a new line tag (ANSA format)
        if(!isInteger(records[pos]) || records[pos].empty()){
            if(records[pos] == "UM" || records[pos].empty()){
                return;
            }
            else{
                // If not integer it is weight factor
                if(!isInteger(records[pos])){
                    if (pos < endpos - 1 && records[pos][0] != '+'){
                        // Skip even components numbers after weight factor
                        pos++;
                    }
                }
            }
        }
        else{
            long item = std::stoi(records[pos]);
            if (!std::count(connectivity.begin(), connectivity.end(), item)){
                connectivity.push_back(item);
            }
        }
    }
}

/*!
 * Absorb rbe2 element from a records structure
 * \param[in] records vector of string records
 * \param[out] ID ID of the RBE2 element
 * \param[out] PID PID of the RBE2 element
 * \param[out] connectivity connectivity of the RBE2 element
 * \param[out] components reaction components of the element
 */
void
NastranInterface::absorbRBE2(std::vector<std::string> & records, long & ID, long & PID, std::vector<long> & connectivity, std::string & components)
{
    connectivity.clear();
    ID = std::stoi(records[1]);
    PID = std::stoi(records[2]);
    components = records[3];
    std::size_t startpos = 4;
    std::size_t endpos = records.size();
    for (std::size_t pos = startpos; pos < endpos; pos++){
        if(isInteger(records[pos])){
            long item = std::stoi(records[pos]);
            if (!std::count(connectivity.begin(), connectivity.end(), item))
                connectivity.push_back(item);
        }
    }
}

/*!
 * Check if a string contains an integer
 * \param[in] str string to test
 * \return true if the string has only digits
 */
bool
NastranInterface::isInteger(std::string & str){
    for (auto &val : str){
        if (!std::isdigit(val))
            return false;
    }
    return true;
}

/*!
 * Get if an element type is enbaled
 * \param[in] type element type
 * \return true if the element type is enbaled for the current read session
 */
bool
NastranInterface::isEnabled(NastranElementType type){
    // If not in the map return false
    if (!m_enabled.count(type)){
        return false;
    }
    return m_enabled[type];
}

/*!
 * Enable an element type for the current read session
 * \param[in] type element type
 */
void
NastranInterface::enable(NastranElementType type){
    m_enabled[type] = true;
}

/*!
 * Disable an element type for the current read session
 * \param[in] type element type
 */
void
NastranInterface::disable(NastranElementType type){
    m_enabled[type] = false;
}

}
