/*
    Get a License Header here- Temporary POC for Click-Ins
*/

#ifndef __IOWAVEFRONTOBJ__HPP__
#define __IOWAVEFRONTOBJ__HPP__

#include <BaseManipulation.hpp>

//data definition
#define XD_WOBJDATA_ "XD_WOBJDATA_"     /**< pointer to Wavefront OBJ data structure */

//port definition
#define X_WDATA         "X_WDATA"          /**< port to pass a XD_WOBJDATA_ data (pointer to WavefrontObjData)*/

namespace mimmo{

/*!
    \class WavefrontObjData
    \ingroup iogeneric
    \brief struct for storing cell data attached to Wavefront OBJ polygonal mesh
*/
class WavefrontObjData{

    friend class IOWavefrontOBJ;

public:

    MimmoPiercedVector<std::string> materials; /**< label identifying material group which a cell can possibly belong to */
    MimmoPiercedVector<long> smoothids; /**< marker identifying smooth group which a cell can possibly belong to */

    std::unordered_map<std::string, long> materialsList; /**< list of all materials inside the class, argument marks long id */
    std::unordered_set<long> smoothidsList; /**< list of all smoothing group ids inside the class */

    std::string materialfile; /**< path to materials file associated to the obj file */

    MimmoPiercedVector<std::array<double,3>> vertexTexture; /**< replaced vertex normals of the object */
    MimmoPiercedVector<std::array<double,3>> vertexNormal; /**< replaced vertex normals of the object */

    WavefrontObjData();
    /*! destructor */
    virtual ~WavefrontObjData(){};

    void swap(WavefrontObjData &x) noexcept;
    void syncListsOnData();

protected:

    void dump(std::ostream & out);
    void restore(std::istream & out);

};

/*!
\class IOWavefrontOBJ
\ingroup iogeneric
\brief Executable block handling io of 3D surface polygonal mesh in *.obj format.

IOWavefrontOBJ manages reading/writing of surface mesh from/to WaveFront ASCII obj format to/from a
MimmoObject surface mesh container (type=1).
More information at https://en.wikipedia.org/wiki/Wavefront_.obj_file
The class does not cover all the features supported by the format, but is restricted to
description of 3D surface tesselleted meshes. Points and Lines cell elements are not supported
in the current class.

Materials file *.mtl attached to the inner *.obj file is not needed by the current class.

Thus, the obj file must have:
- Polygonal geometry statement only, No free form statement features (NURBS-CAD) are supported.
- Polygonal cells (only facets, no lines or vertex cells)

Mesh Object names/parts are absorbed/flushed through PID/PIDNames mechanism
of MimmoObject surface mesh cells.

Textures information (if any) are managed as a MimmoObject mesh which shares
the same cell ids as the original polygonal Mesh.

Materials assigned on cells and smoothing group ids are managed as special
fields attached to the MimmoObject mesh cells-ids.

Group naming (g key entry) is ignored and not supported in the current class;

\n
Ports available in IOWaveFrontOBJ Class :

========================================================

|Port Input | | |
|-|-|-|
| <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
| M_GEOM       | setGeometry          | (MC_SCALAR, MD_MIMMO_)          |
| X_WDATA      | setData              | (MC_SCALAR, XD_WOBJDATA_) |
| M_STRINGFIELD      | setMaterials              | (MC_SCALAR, MD_MPVECSTRING_) |
| M_LONGFIELD      | setSmoothIds              | (MC_SCALAR, MD_MPVECLONG_) |
| M_VECTORFIELD      | setNormals              | (MC_SCALAR, MD_MPVECARR3FLOAT_) |
| M_VECTORFIELD2    | setTexture           | (MC_SCALAR, MD_MPVECARR3FLOAT_) |

|Port Output | | |
|-|-|-|
| <B>PortType</B> | <B>variable/function</B> |<B>DataType</B>|
| M_GEOM       | getGeometry          | (MC_SCALAR,MD_MIMMO_)      |
| X_WDATA      | getData              | (MC_SCALAR, XD_WOBJDATA_) |
| M_STRINGFIELD      | getMaterials              | (MC_SCALAR, MD_MPVECSTRING_) |
| M_LONGFIELD      | getSmoothIds              | (MC_SCALAR, MD_MPVECLONG_) |
| M_VECTORFIELD      | getNormals              | (MC_SCALAR, MD_MPVECARR3FLOAT_) |
| M_VECTORFIELD2    | getTexture           | (MC_SCALAR, MD_MPVECARR3FLOAT_) |

=========================================================
    \n

     The xml available parameters, sections and subsections are the following :

     Inherited from BaseManipulation:
     - <B>ClassName</B>: name of the class as <tt>mimmo.IOWavefrontOBJ</tt>;
     - <B>Priority</B>: uint marking priority in multi-chain execution;

     Proper of the class:
     - <B>IOMode</B>: 0-read from obj file, 1-restore from class dump file, 2- write to obj file, 3-write class dump file.
     - <B>Dir</B> directory path to read from/write to
     - <B>Filename</B> name of the file to be read/written (without tag .obj)
     - <B>GeomTolerance</B> (float) geometric tolerance for double nodes collapsing.
     - <B>TextureTolerance</B> (float) texture geometric tolerance for double nodes collapsing.
     - <B>SkipTexture</B>: 0/1 if 1, force to skip reading/writing of any texture associated to mesh.
     - <B>PrintResumeFile</B>: 0/1 print a resume file of the mesh contents after execution.

     Geometry and additional fields have to be mandatorily passed through ports.
*/
class IOWavefrontOBJ: public mimmo::BaseManipulation{

public:
    /*!
        \ingroup iogeneric
        \brief Working mode for class IOWavefrontOBJ
    */
    enum class IOMode{
        READ = 0, /**< Read from file *.obj */
        RESTORE = 1, /**< Restore from dump file *.dump */
        WRITE = 2, /**< write to file *.obj*/
        DUMP = 3 /**< write to dump file *.dump */
    };

    IOWavefrontOBJ(IOMode mode = IOMode::READ);
    virtual ~IOWavefrontOBJ();
    IOWavefrontOBJ(const bitpit::Config::Section & rootXML);

    IOMode         whichMode();
    int            whichModeInt();

    MimmoObject*                          	getGeometry();
    std::string								getMaterialFile();
    WavefrontObjData*                     	getData();
    std::unordered_map<long, std::string>	getSubParts();

    MimmoPiercedVector<std::string>*			getMaterials();
    MimmoPiercedVector<long>*					getSmoothIds();
    MimmoPiercedVector<std::array<double,3>>*	getNormals();
    MimmoPiercedVector<std::array<double,3>>*	getTexture();

    void	setGeometry(MimmoObject * geo);
    void    setMaterialFile(std::string materialfile);
    void	setData(WavefrontObjData* data);

    void	setMaterials(MimmoPiercedVector<std::string>* materials);
    void	setSmoothIds(MimmoPiercedVector<long>* smoothids);
    void	setNormals(MimmoPiercedVector<std::array<double,3>>* normals);
    void	setTexture(MimmoPiercedVector<std::array<double,3>>* texture);
    void	setGeometryDisplacements(MimmoPiercedVector<std::array<double,3>>* displacements);

    void	setDir(const std::string & pathdir);
    void    setFilename(const std::string & name);
    void    setGeomTolerance(double tolerance);
    void    setTextureTolerance(double tolerance);
    void    printResumeFile(bool print);

    void    execute();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name = "");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name= "");

protected:
    void swap(IOWavefrontOBJ & x) noexcept;
    void buildPorts();

    void read(const std::string & fullpath);
    void write(const std::string & fullpath);
    void restore(std::istream &);
    void dump(std::ostream &);
    void writeResumeFile();

    std::string searchMaterialFile(std::ifstream & in);

    void searchObjectPosition(std::ifstream & in,
                              std::vector<std::streampos> & mapPos,
                              std::vector<std::string>& mapNames,
                              long &nVertTot, long &nCellTot);

    void readObjectData(std::ifstream & in, const std::streampos &begObjectStream, const long &PID,
                        long &vOffset, long &vnOffset, long &vTxtOffset, long &cOffset);

    void writeObjectData(WavefrontObjData* objData, std::ofstream & out, const livector1D &vertList, const livector1D &txtVertList, const livector1D & vnVertList,
                         const livector1D & cellList, long & vOffset, long & vnOffset, long & vTxtOffset, long& cOffset);

    std::unordered_map<std::string, std::vector<long>> regroupCellsByMaterials(const WavefrontObjData* objData, const livector1D & cellList);

    void computeMovedNormals();

private:

    IOMode m_mode;      /**< working mode */
    std::unique_ptr<MimmoObject> m_intPatch; /**< internal mesh  */
    std::unique_ptr<WavefrontObjData> m_intData; /**< internal data  */
    WavefrontObjData * m_extData; /**< externally linked data*/

    std::string m_dir; /**< io directory path  */
    std::string m_filename; /**< io name of the file  */
    bool m_resume; /**< boolean to print resume file */
    double m_tol; /**< geometric tolerance for duplicate vertex collapsing */

    std::size_t m_totalVertexCount;
    std::size_t m_totalCellCount;

    MimmoPiercedVector<std::array<double,3>> m_displacements; /**< Geometry vertex displacements. If filled is used to interpolate new normals during write.*/

    int convertKeyEntryToInt(const std::string & key);
};

REGISTER_PORT(M_GEOM, MC_SCALAR, MD_MIMMO_, __IOWAVEFRONTOBJ__HPP__)
REGISTER_PORT(X_WDATA, MC_SCALAR, XD_WOBJDATA_, __IOWAVEFRONTOBJ__HPP__)
REGISTER_PORT(M_STRINGFIELD, MC_SCALAR, MD_MPVECSTRING_, __IOWAVEFRONTOBJ__HPP__)
REGISTER_PORT(M_LONGFIELD, MC_SCALAR, MD_MPVECLONG_, __IOWAVEFRONTOBJ__HPP__)
REGISTER_PORT(M_VECTORFIELD, MC_SCALAR, MD_MPVECARR3FLOAT_, __IOWAVEFRONTOBJ__HPP__)
REGISTER_PORT(M_VECTORFIELD2, MC_SCALAR, MD_MPVECARR3FLOAT_, __IOWAVEFRONTOBJ__HPP__)
REGISTER_PORT(M_GDISPLS, MC_SCALAR, MD_MPVECARR3FLOAT_, __IOWAVEFRONTOBJ__HPP__)

REGISTER(BaseManipulation, IOWavefrontOBJ, "mimmo.IOWaveFrontOBJ")

}

#endif /* __IOWAVEFRONTOBJ__HPP__ */
