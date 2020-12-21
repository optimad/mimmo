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
#ifndef __MIMMOGEOMETRY_HPP__
#define __MIMMOGEOMETRY_HPP__

#include "BaseManipulation.hpp"
#include "enum.hpp"
#include <typeinfo>
#include <type_traits>

BETTER_ENUM(FileType, int, STL = 0, SURFVTU = 1, VOLVTU = 2, NAS = 3, PCVTU = 4, CURVEVTU = 5, MIMMO = 99);

BETTER_ENUM(NastranElementType, int, GRID = 0, CTRIA = 1, CQUAD = 2, CBAR = 3, RBE2 = 4, RBE3 = 5);

namespace mimmo{

class NastranInterface;

/*!
 * \ingroup iogeneric
 * Format of data to read/write the Nastran geometries.
 */
enum WFORMAT{
    Short = 0, /**<Single precision data.*/
    Long =  1  /**<Double precision data.*/
};


/*!
 * \class MimmoGeometry
 * \ingroup iogeneric
 * \brief MimmoGeometry is an executable block class wrapping(linking or internally instantiating) a Mimmo Object, handling geometry.
 *
 *  MimmoGeometry is the object to manage the import/export/substantial modifications of geometries. It returns a
 *  standard MimmoObject class as product of its execution;
 *  The mesh to import/export/convert has to be a mesh with constant type elements.
 *  The valid format are: binary .stl, ascii .vtu (triangle/quadrilateral elements) and
 *  ascii .nas (triangle elements) for surface mesh; ascii .vtu (tetra/hexa elements)
 *  for volume mesh.
 *
 * On distributed archs, MimmoGeometry can write and read in parallel in VTU formats and in mimmo dumping format.
 * In case of VTU formats the collective file extension .pvtu is used for the input/output filename.
 * When reading, the number of processors used to write the input file had to be equal
 * to the current number of processors used to read.
 * If MPI support is enabled but the number of processors is 1, MimmoGeometry reads/writes VTU file as serial one,
 * by postponing the usual file extension .vtu.
 * Nastran files are always written in serial mode, i.e. only the master rank 0 writes on file its partition. In order
 * to write a distributed mesh in a Nastran format, the user has to serialize it before to pass the MimmoObject
 * to the MimmoGeometry block.
 *
 * When MPI is enabled a read geometry in mimmo dumping can be distributed among the processes (parallel)
 * or repeated serially among the processes (not parallel).
 * The different behavior can be controlled by a set method (setParallelRestore).
 * Note that the parallel/not parallel property set on block
 * has to be coherent with the restored geometry, i.e. it has to be known a-priori by the user.
 * Default value of the feature is true (parallel) in case of MPI enabled when reading a mimmo dumping format geometry.
 *
 *  \n
 *  It can be used in three modes reader/writer/converter. To set the mode it uses an enum
 *  IOMode list:
 *  - <B>READ    = 0</B>
 *  - <B>WRITE   = 1</B>
 *  - <B>CONVERT = 2</B>
 *
 * \n
 * It uses smart enums FileType list of available geometry formats, which are:
 *
 * - <B>STL     = 0</B>    Ascii/Binary triangulation stl.
 * - <B>SURFVTU = 1</B> Surface mesh vtu, of any 2D bitpit::ElementType.
 * - <B>VOLVTU  = 2</B> Volume mesh VTU, of any 3D bitpit::ElementType.
 * - <B>NAS     = 3</B> Nastran surface triangular/quad surface meshes.
 * - <B>PCVTU   = 4</B> Point Cloud VTU, of only VERTEX elements
 * - <B>CURVEVTU= 5</B> 3D Curve in VTU, of only LINE elements
 * - <B>MIMMO   = 99</B> mimmo dump/restore format *.geomimmo
 *
 * Outside this list of options, the class cannot hold any other type of formats for now.
 * The smart enum can be recalled in every moment in your code, just using <tt>mimmo::FileType</tt>
 * and including MimmoGeometry header.
 *
 * \n
 * Ports available in MimmoGeometry Class :
 *
 *    =========================================================

     |Port Input|||
     ||||
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_GEOM   | m_geometry        | (MC_SCALAR, MD_MIMMO_)      |


     |Port Output |||
     ||||
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_GEOM    | getGeometry       | (MC_SCALAR, MD_MIMMO_)      |

 *    =========================================================
 * \n
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.Geometry</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>IOMode</B>: activate READ, WRITE or CONVERT mode;
 * - <B>Dir</B>: directory path;
 * - <B>Filename</B>: name of file for reading/writing;
 * - <B>FileType</B>: file type identifier;
 * - <B>ReadFileType</B>: file type identifier for reader (to be used in case of converter mode);
 * - <B>WriteFileType</B>: file type identifier for writer (to be used in case of converter mode);
 * - <B>WriteMultiSolidSTL</B>: 0-false/1-true, if WriteFileType is STL, MultiSolid STL writing can be activated or not
 * - <B>ReadDir</B>: directory path (to be used in case of converter mode and different paths);
 * - <B>ReadFilename</B>: name of file for reading/writing (to be used in case of converter mode and different filenames);
 * - <B>WriteDir</B>: directory path (to be used in case of converter mode and different paths);
 * - <B>WriteFilename</B>: name of file for reading/writing (to be used in case of converter mode and different filenames);
 * - <B>Codex</B>: boolean to write ascii/binary;
 * - <B>SkdTree</B>: evaluate SkdTree true 1/false 0;
 * - <B>KdTree</B>: evaluate kdTree true 1/false 0.
 * - <B>AssignRefPID</B>: assign a reference PID on the whole geometry, after reading or just before writing. If the geometry is already pidded,
 *                     translate all existent PIDs w.r.t. the reference PID assigned. Default value is RefPID = 0.
 * - <B>Tolerance</B>:value of the geometric tolerance to be used;
 * - <B>Clean</B>: clean the geometry after read true 1/false 0;
 * - <B>FormatNAS</B>: 0-singlePrecision 1-doubleprecision for writing nas files;
 * - <B>ParallelRestore</B>: set if the read geometry is parallel true 1/false 0;

 *
 * In case of writing mode Geometry has to be mandatorily passed through port.
 *
 */
// TODO ADD CUSTOM NUMBER OF PROCESSORS FOR INPUT/OUTPUT OF VTU FILES
class MimmoGeometry: public BaseManipulation{

protected:
    bool         m_read;         /**<If true it reads the geometry from file during the execution.*/
    FileDataInfo m_rinfo;       /**< Info on the external file to read */

    bool         m_write;         /**<If true it writes the geometry on file during the execution.*/
    FileDataInfo m_winfo;       /**< Info on the external file to write */

    bool        m_codex;                    /**< Set codex format for writing true binary, false ascii */
    WFORMAT        m_wformat;                    /**<Format for .nas import/export. (Short/Long).*/

    bool        m_buildSkdTree;             /**<If true the simplex ordered SkdTree of the geometry is built in execution, whenever geometry support simplicies. */
    bool        m_buildKdTree;                /**<If true the vertex ordered KdTree of the geometry is built in execution*/
    long        m_refPID;                     /**<Reference PID, to be assigned on all cells of geometry in read/convert mode*/
    bool        m_multiSolidSTL;            /**< activate or not MultiSolid STL writing if STL writing Filetype is selected */

    double		m_tolerance;				/**<Geometric tolerance of the related geometry. */
    bool        m_clean;                    /**<Set if the geometry has to cleaned after reading. */

    bool        m_parallelRestore;               /**<Set if the geometry to read is parallel. */

public:
    /*!
     * \ingroup iogeneric
     * enum for MimmoGeometry class mode
     */
    enum class IOMode{
        READ = 0,    /**<reading mode */
        WRITE = 1,   /**< writing mode */
        CONVERT = 2 /**< convert mode, i.e. first read from, then write to */
    };

    MimmoGeometry(IOMode mode = IOMode::READ);
    MimmoGeometry(const bitpit::Config::Section & rootXML);
    virtual ~MimmoGeometry();

    MimmoGeometry(const MimmoGeometry & other);
    MimmoGeometry & operator=(MimmoGeometry other);

    void buildPorts();

    const MimmoGeometry *    getCopy();

public:
    void        setReadDir(std::string dir);
    void        setWriteDir(std::string dir);
    void        setReadFilename(std::string filename);
    void        setWriteFilename(std::string filename);

    void        setDir(std::string dir);
    void        setFilename(std::string filename);
    void        setReadFileType(FileType type);
    void        setReadFileType(int type);
    void        setWriteFileType(FileType type);
    void        setWriteFileType(int type);
    void        setFileType(FileType type);
    void        setFileType(int type);
    void        setCodex(bool binary = true);
    void        setMultiSolidSTL(bool multi = true);

    void 		setTolerance(double tol);
    void        setClean(bool clean = true);
    void        setParallelRestore(bool parallelRestore = MIMMO_ENABLE_MPI);

    using BaseManipulation::setGeometry;
    void        setGeometry(int type=1);

    bitpit::PiercedVector<bitpit::Vertex> *     getVertices();
    bitpit::PiercedVector<bitpit::Cell> *         getCells();


    void        setPID(livector1D pids);
    void        setPID(std::unordered_map<long,long> pidsMap);
    void        setReferencePID(long pid);

    void        setFormatNAS(WFORMAT wform);

    void        setBuildSkdTree(bool build);
    void        setBuildKdTree(bool build);

    bool         isEmpty();
    void         clear();
    bool        write();
    bool        read();

    void         execute();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

protected:
    void swap(MimmoGeometry & x) noexcept;
    void        setIOMode(IOMode mode);
    void        setIOMode(int mode);
    void    setDefaults();
    void    _setRead(bool read = true);
    void    _setWrite(bool write = true);
    bool   fileExist(const std::string & filename);

};

/*!
 * \class NastranInterface
 * \ingroup iogeneric
 * \brief NastranInterface is an interface class for I/O handling of BDF bulk nastran format *.nas.
 */
class NastranInterface{
    static const char nl = '\n'; /**< static const declaration of new line command. */
public:

    WFORMAT    m_wformat; /**< member storing the file type format Short/Long */

    std::map<NastranElementType, bool> m_enabled; /** element type enabled for the current session */

    NastranInterface();

    void setWFormat(WFORMAT);
    void writeKeyword(std::string key, std::ofstream& os);
    void writeCoord(const darray3E & p, const long& pointI, std::ofstream& os);
    void writeFace(std::string faceType, const livector1D& facePts, const long& nFace, std::ofstream& os, long PID);
    bool writeGeometry(dvecarr3E &points,livector1D &pointsID, livector2D &faces, livector1D &facesID, std::ofstream &os, livector1D *PIDS = nullptr);
    void writeFooter(std::ofstream& os, std::unordered_set<long>* PIDSSET = nullptr);
    void write(std::string& outputDir, std::string& surfaceName, dvecarr3E& points, livector1D& pointsID, livector2D& faces, livector1D& facesID, livector1D* PIDS = nullptr, std::unordered_set<long>* PIDSSET = nullptr);
    void read(std::string& inputDir, std::string& surfaceName, dvecarr3E& points, livector1D& pointsID, livector2D& faces, livector1D& facesID, livector1D& PIDS);

    std::string trim(std::string in);
    std::string convertVertex(std::string in);

    void appendLineRecords(std::string & sread, std::size_t recordlength, std::vector<std::string> & records);
    void absorbRBE2(std::vector<std::string> & records, long & ID, long & PID, std::vector<long> & connectivity, std::string & components);
    void absorbRBE3(std::vector<std::string> & records, long & ID, long & PID, std::vector<long> & connectivity, std::string & components);
    bool isInteger(std::string & str);

    bool isEnabled(NastranElementType type);

protected:
    void enable(NastranElementType type);
    void disable(NastranElementType type);


    /*!
     * Write a template value according to WFORMAT chosen by the User
     * \param[in]     value    value of class Type
     * \param[in,out] os    ofstream where the value is written
     */
    template<class Type>
    void writeValue (Type& value, std::ofstream& os){

        std::size_t offset = 5;
        if (m_wformat == WFORMAT::Long) offset = 13;

        std::stringstream manip;
        manip <<std::fixed<<std::setprecision(offset+3);
        manip<<value;
        std::string manips = manip.str();
        std::string mantissa, expon;
        std::size_t pos = manips.find("E");
        if (pos < manips.size()){
            mantissa = manips.substr(0,std::min(offset,pos));
            expon = manips.substr(pos+1,manips.size());
            manips = mantissa + expon;
        }
        pos = manips.find("e");
        if (pos < manips.size()){
            mantissa = manips.substr(0,std::min(offset,pos));
            expon = manips.substr(pos+1,manips.size());
            manips = mantissa + expon;
        }
        manips = manips.substr(0, offset+3);

        switch (m_wformat)
        {
        case Short:
        {
            os << std::setw(8) << manips;
            break;
        }
        case Long:
        {
            os << std::setw(16) << manips;
            break;
        }
        }
    }
};

REGISTER_PORT(M_GEOM, MC_SCALAR, MD_MIMMO_, __MIMMOGEOMETRY_HPP__)
REGISTER(BaseManipulation, MimmoGeometry, "mimmo.Geometry")

};

#endif /* __MIMMOGEOMETRY_HPP__ */
