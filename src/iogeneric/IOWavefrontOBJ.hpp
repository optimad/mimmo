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

#ifndef __IOWAVEFRONTOBJ__HPP__
#define __IOWAVEFRONTOBJ__HPP__

#include <BaseManipulation.hpp>

namespace mimmo{

/*!
    \class WavefrontOBJData
    \ingroup iogeneric
    \brief struct for storing cell data attached to Wavefront OBJ polygonal mesh

    BEWARE MimmoPiercedVector for cellgroups, materials and smoothids must be referred to
    Mesh linked in refGeometry. Textures and Normals are treated as indipendent MimmoObject w.r.t.
    Mesh MimmoObject( each one has its own set of nodes, and local connectivities), but they must always share
    the same number and cell-id labels for element/simplex cross referencing between the three MimmoObjects.

*/
class WavefrontOBJData{

    friend class IOWavefrontOBJ;

public:

    //data
    MimmoPiercedVector<std::string> materials; /**< label identifying material group which a cell can possibly belong to */
    MimmoPiercedVector<long> smoothids; /**< marker identifying smooth group which a cell can possibly belong to */
    MimmoPiercedVector<std::string> cellgroups; /**< label identifying cell group label which a cell can possibly belong to */

    MimmoSharedPointer<MimmoObject> textures; /**< MimmoObject container for texture properties*/
    MimmoSharedPointer<MimmoObject> normals; /**<MimmoObject container for vnormals properties*/

    //list
    std::unordered_map<std::string, long> materialsList; /**< list of all materials inside the class, argument marks long id */
    std::unordered_map<std::string, long> smoothidsList; /**< list of all smoothing group ids inside the class, arguments marks long id */
    std::unordered_map<std::string, long> cellgroupsList; /**< list of all cell group labels inside the class, arguments marks long id */

    std::unordered_map<long, std::string> inv_materialsList; /**< inverse of materialList*/
    std::unordered_map<long, std::string> inv_smoothidsList; /**< inverse of smoothidslist*/
    std::unordered_map<long, std::string> inv_cellgroupsList; /**< inverse of cellgroupslist*/

    // accessory info
    std::string materialfile; /**< path to materials file associated to the obj file */
    MimmoSharedPointer<MimmoObject> refGeometry; /**<reference geometry which the class belongs to */

    WavefrontOBJData();
    /*! destructor */
    virtual ~WavefrontOBJData(){};

    void swap(WavefrontOBJData &x) noexcept;
    void syncListsOnData();
    void autoCompleteCellFields();

protected:

    void dump(std::ostream & out);
    void restore(std::istream & out);

};

/*!
\class ManipulateWFOBJData
\ingroup iogeneric
\brief Executable block manipulating optional data of WavefrontOBJ mesh

The class performs manipulation of Wavefront mesh optional data such as:

- cell annotation (annotating a group of cells throughout the sub-objects of the mesh).
  Different strategies for multiple annotations on a single cell are available.
  See OverlapAnnotationMode enum. OverlapAnnotationMode::HARD is class default.

- recompute vertex normals on nodes of a list of candidate cells of the mesh. Once a list of
  mesh cell-ids is provided data on normals (stored in WavefrontOBJData) are recomputed
  according to NormalsComputeMode enum. NormalsComputeMode::FLAT_BITPIT is the default.
  Please note that recomputing normals will modify MimmoObject internal structures,
  so its properties (kdtrees, adjacencies etc..) will be invalidated.

- recover list of cell-ids referring to target Wavefront data cell properties such as
  materials, cellgroups, smoothids (f.e. extract all mesh cell-ids with material
  property "body").
\n
Ports available in ManipulateWFOBJData Class :

========================================================

|Port Input | | |
|-|-|-|
| <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
| M_WAVEFRONTDATA   | setData                   | (MC_SCALAR, MD_WOBJDATA_) |
| M_LONGFIELD       | addAnnotation             | (MC_SCALAR, MD_MPVECLONG_)
| M_LONGFIELD2      | setRecomputeNormalsCells  | (MC_SCALAR, MD_MPVECLONG_) |

|Port Output | | |
|-|-|-|
| <B>PortType</B> | <B>variable/function</B>  |<B>DataType</B>|
| M_WAVEFRONTDATA     | getData               | (MC_SCALAR, MD_WOBJDATA_) |
| M_VECTORLI          | getPinnedMaterialGroup| (MC_VECTOR, MD_LONG) |
| M_VECTORLI2         | getPinnedCellGroup    | (MC_VECTOR, MD_LONG) |
| M_VECTORLI3         | getPinnedSmoothGroup  | (MC_VECTOR, MD_LONG) |
| M_VECTORLI4         | getPinnedObjectGroup  | (MC_VECTOR, MD_LONG) |


=========================================================
    \n

     The xml available parameters, sections and subsections are the following :

     Inherited from BaseManipulation:
     - <B>ClassName</B>: name of the class as <tt>mimmo.ManipulateWFOBJData</tt>;
     - <B>Priority</B>: uint marking priority in multi-chain execution;

     Proper of the class:
     - <B>CheckNormalsMagnitude</B>: 0/1 boolean, if 1, check and force Wavefront Data to have normals magnitude equal to 1. If no normals are available in wavefront data it does nothing.
     - <B>MultipleAnnotationStrategy</B>: 0/1/2 choose strategy to deal with multiple annotations concurring into a target cell. see enum OverlapAnnotationMode. Default is 0 - FAWIN.
     - <B>NormalsComputeStrategy<B>: 0/1/2 choose strategy to recompute normals on candidate cells, see NormalsComputeMode enum.
     - <B>PinObjects</B> list of objects/subparts names (blank separated) to pin mesh cells belonging to them.
                           Cells-ids will be available with method getPinnedObjectGroup method/M_VECTORLI4 port.

     - <B>PinMaterials</B> list of materials names (blank separated) to pin mesh cells owning these properties.
                           Cells-ids will be available with method getPinnedMaterialGroup method/M_VECTORLI port.
     - <B>PinCellGroups</B> list of cellgroups names (blank separated) to pin mesh cells owning these properties.
                           Cells-ids will be available with method getPinnedCellGroup method/M_VECTORLI2 port.
     - <B>PinSmoothIds</B> list of smooth-ids int (blank separated) to pin mesh cells owning these properties.
                           Cells-ids will be available with method getPinnedSmoothGroup method/M_VECTORLI3 port.

     Data and additional input fields have to be mandatorily passed through ports.
*/
class ManipulateWFOBJData: public mimmo::BaseManipulation{

public:
    /*!
        \ingroup iogeneric
        \brief Strategy Mode to deal with Multiple Annotations on a target cell in ManipulateWFOBJData
    */
    enum class OverlapAnnotationMode{
        HARD = 0,  /**< Cell is continuously marked by annotations which refer to it. The last annotation in the list wins.*/
        SOFT = 1,  /**< Cell is marked with an annotation value, if and only if no other
                        annotations are present. Successive referring annotations in the list are ignored */
        GETALL = 2,    /**< Multiple annotations values on a cell are chained toghether (blank spaced) to form a new unique marker */
        GETALLNOBLANKS = 3    /**< do exactly as GETALL without blank space separators*/
    };

    /*!
        \ingroup iogeneric
        \brief Strategy Mode to recompute Wavefront vertex normals on a set of candidate cells
    */
    enum class NormalsComputeMode{
        FLAT_BITPIT = 0,  /**< vertex normals recomputed on cell nodes as standard average of 1-Ring facet cells */
        AREA_WEIGHTED = 1,  /**< vertex normals recomputed on cell nodes as area-weighted average of 1-Ring facet cells */
        ANGLE_WEIGHTED = 2  /**< vertex normals recomputed on cell nodes as incident edges angle-weighted average of 1-Ring facet cells */
    };

    ManipulateWFOBJData();
    virtual ~ManipulateWFOBJData();
    ManipulateWFOBJData(const bitpit::Config::Section & rootXML);

    WavefrontOBJData*       getData();
    bool                    getCheckNormalsMagnitude();
    OverlapAnnotationMode   getMultipleAnnotationStrategy();
    NormalsComputeMode      getNormalsComputeStrategy();

    std::vector<long>      getPinnedObjectGroup();
    std::vector<long>      getPinnedMaterialGroup();
    std::vector<long>      getPinnedCellGroup();
    std::vector<long>      getPinnedSmoothGroup();

    void    setRecomputeNormalsCells(MimmoPiercedVector<long>* cellList);
    void    setData(WavefrontOBJData* data);
    void    addAnnotation(MimmoPiercedVector<long>* annotation);

    void    setCheckNormalsMagnitude(bool flag);
    void    setMultipleAnnotationStrategy(OverlapAnnotationMode mode);
    void    setNormalsComputeStrategy(NormalsComputeMode mode);

    void    setPinObjects(const std::vector<std::string> & objectsList);
    void    setPinMaterials(const std::vector<std::string> & materialsList);
    void    setPinCellGroups(const std::vector<std::string> & cellgroupsList);
    void    setPinSmoothIds(const std::vector<int> & smoothidsList);

    void clearAnnotations();
    void clearRecomputeNormalsCells();
    void clearPinLists();
    void clear();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name = "");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name= "");

    void    execute();

protected:
    void swap(ManipulateWFOBJData & x) noexcept;
    void buildPorts();

    void computeAnnotations();
    void computeNormals();
    void checkNormalsMagnitude();
    void extractPinnedLists();

    std::array<double,3> evalVNormal(bitpit::SurfaceKernel * mesh, long idCell, int locVertex);

    WavefrontOBJData * m_extData; /**< externally linked data*/
    MimmoPiercedVector<long> m_normalsCells; /**< candidate cell list for normals recomputation.*/
    std::vector<MimmoPiercedVector<long> > m_annotations; /**< list of annotations*/
    bool m_checkNormalsMag; /**< boolean to set check on normals magnitude */
    OverlapAnnotationMode m_annMode; /**< mode to deal with multiple concurrent annotations on a cell*/
    NormalsComputeMode m_normalsMode; /**< mode to deal with normals recompute on a list of candidate cells*/

    std::unordered_set<std::string> m_pinObjects; /**< list of object/parts names for cell list pinning*/
    std::unordered_set<std::string> m_pinMaterials; /**< list of material names for cell list pinning*/
    std::unordered_set<std::string> m_pinCellGroups; /**< list of cellgroup names for cell list pinning*/
    std::unordered_set<int> m_pinSmoothIds; /**< list of smooth-ids int entries for cell list pinning*/

    std::array<std::vector<long>, 4> m_pinnedCellLists; /**< list of cell pinned 0-material, 1-cellgroups, 2-smoothids, 3-objects/subparts */

    bool checkEntry(const std::string& entry, const std::string& root);

private:
    //make useless base class methods private;
    MimmoSharedPointer<MimmoObject> getGeometry(){return nullptr;};
    void setGeometry(MimmoSharedPointer<MimmoObject> geo){BITPIT_UNUSED(geo);};
};

/*!
\class IOWavefrontOBJ
\ingroup iogeneric
\brief Executable block handling io of 3D surface polygonal mesh in *.obj format.

IOWavefrontOBJ manages reading/writing of surface mesh from/to WaveFront ASCII obj format to/from a
MimmoObject surface mesh container (type=1).
More information at https://en.wikipedia.org/wiki/Wavefront_.obj_file
The class does not cover all the features supported by the format, but is restricted to
description of 3D surface tessellated meshes. Points and Lines cell elements are not supported
in the current class.

Materials file *.mtl attached to the inner *.obj file is not needed by the current class.

Thus, the obj file must have:
- Polygonal geometry statement only, No free form statement features (NURBS-CAD) are supported.
- Polygonal cells (only facets, no lines or vertex cells)

Mesh Object sub-parts are absorbed/flushed through PID/PIDNames mechanism
of MimmoObject surface mesh cells. The class identifies as parts only sub-objects defined
with o key entry (grouping with g key entry is ignored.)

Mesh accessory data are stored in a WaveFrontObjData class, containing:

- Textures and Normals properties (if any) are managed in independent MimmoObjects
(their node coordinates list and their own connectivities).The link is that cell-ids
are the same in all 3 objects (reference mesh, textures and normals).

 - Materials assigned on cells, smoothing group ids and cellgroups id are managed as special
fields attached to the MimmoObject mesh cells-ids.


WARNING 1: for reading/writing of v,vt,vn data chunks, the class supposes that the three of them
are organized in the file in blocks of consecutive entries, i.e. after 'o' sub-object declaration
all the v's of the object are written, then the vt's if any, and finally the vn's if any.
This may be not true for the general Wavefront format(in which declarations can be random and not organized at all),
but meet quite well with the practice of well known OBJ file exporter like the Blender's one.

WARNING 2: During OBJ writing if empty cellgroup entry is encountered, it will be assigned the
object name (o) by default.

\n
Ports available in IOWaveFrontOBJ Class :

========================================================

|Port Input | | |
|-|-|-|
| <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
| M_GEOM            | setGeometry               | (MC_SCALAR, MD_MIMMO_)          |
| M_WAVEFRONTDATA   | setData                   | (MC_SCALAR, MD_WOBJDATA_) |

|Port Output | | |
|-|-|-|
| <B>PortType</B> | <B>variable/function</B> |<B>DataType</B>|
| M_GEOM              | getGeometry          | (MC_SCALAR,MD_MIMMO_)     |
| M_WAVEFRONTDATA     | getData              | (MC_SCALAR, MD_WOBJDATA_) |
| M_STRINGFIELD       | getMaterials         | (MC_SCALAR, MD_MPVECSTRING_) |
| M_STRINGFIELD2      | getCellGroups        | (MC_SCALAR, MD_MPVECSTRING_) |
| M_LONGFIELD         | getSmoothIds         | (MC_SCALAR, MD_MPVECLONG_) |
| M_GEOM2             | getTextures          | (MC_SCALAR,MD_MIMMO_)     |
| M_GEOM3             | getNormals           | (MC_SCALAR,MD_MIMMO_)     |
| M_NAME              | getMaterialFile      | (MC_SCALAR, MD_STRING) |

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
     - <B>GeomTolerance</B> (float) geometric tolerance for repeated mesh vertices.
     - <B>CleanDoubleMeshVertices</B>: for READ Mode only, boolean 0/1 if 1, force mesh vertices collapsing after reading with GeomTolerance.
     - <B>IgnoreCellGroups</B>: for READ Mode only, boolean 0/1 if 1, ignore cellgroups labeling of the OBJ mesh.
     - <B>TextureUVMode</B>: for WRITE Mode only, boolean 0/1 if 1, force writing textures with first 2 components only.
     - <B>PrintResumeFile</B>: 0/1 print a resume file of the mesh contents after execution.

     Geometry and additional fields have to be mandatorily passed through ports.
*/
class IOWavefrontOBJ: public mimmo::BaseManipulation{

public:
    typedef std::map<long,std::map<long, std::vector<long>>> TreeGroups;
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

    WavefrontOBJData*                       getData();
    std::unordered_map<long, std::string>   getSubParts();
    std::string                             getMaterialFile();

    MimmoPiercedVector<std::string>*        getMaterials();
    MimmoPiercedVector<std::string>*        getCellGroups();
    MimmoPiercedVector<long>*               getSmoothIds();

    MimmoSharedPointer<MimmoObject>                           getNormals();
    MimmoSharedPointer<MimmoObject>                           getTextures();

    void    setGeometry(MimmoSharedPointer<MimmoObject> geo);
    void    setData(WavefrontOBJData* data);
    void    setMaterialFile(std::string materialfile);
    void    setDir(const std::string & pathdir);
    void    setFilename(const std::string & name);
    void    setGeometryTolerance(double tolerance);
    void    setCleanDoubleMeshVertices(bool clean);
    void    setIgnoreCellGroups(bool ignore);
    void    setTextureUVMode(bool UVmode);
    void    printResumeFile(bool print);

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name = "");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name= "");

    void    execute();

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
                              std::vector<std::array<long,3>>& mapVCountObject,
                              std::array<long,3>& mapVCountTotal,
                              long &nCellTot);

    void readObjectData(std::ifstream & in, const std::streampos &begObjectStream, const long &PID,
                        const std::string & defaultGroup, const std::array<long,3> & vCounter,
                        long &vOffset, long &vnOffset, long &vTxtOffset, long &cOffset);

    void writeObjectData(WavefrontOBJData* objData, std::ofstream & out,
                         const std::array<std::vector<long>,3> & vertexLists,
                         const std::vector<long> & cellList,
                         std::array<long,3> &vOffsets,
                         long &cOffset,
                         std::array<std::unordered_map<long,long>,3> & vinsertion_maps,
                         const std::string & defaultGroup,
                         long & activeGroup, long & activeMaterial, long &activeSmoothId);

    TreeGroups regroupCells(const WavefrontOBJData* objData, const livector1D & cellList);


private:

    IOMode m_mode;      /**< working mode */
    std::unique_ptr<WavefrontOBJData> m_intData; /**< internal data  */
    WavefrontOBJData * m_extData; /**< externally linked data*/

    std::string m_dir; /**< io directory path  */
    std::string m_filename; /**< io name of the file  */
    bool m_resume; /**< boolean to print resume file */
    double m_tol; /**< geometric tolerance for duplicate mesh vertex/ vnormals collapsing */
    bool m_cleanDoubleVertices; /**< valid for read mode only, clean repeated mesh vertices if required */
    bool m_ignoringCellGroups; /**< boolean to skip absorbing cell groups g while reading */
    bool m_textureUVMode;  /**< true to force writing textures in UV mode (first two components only) */

    //INTERNAL METHODS
    int convertKeyEntryToInt(char key);
    long pushCell(MimmoSharedPointer<MimmoObject> geo, std::vector<long> &conn, long PID, long id, int rank = -1 );
    int checkFacetDefinition(const std::string & str);
};

REGISTER_PORT(M_GEOM, MC_SCALAR, MD_MIMMO_, __IOWAVEFRONTOBJ__HPP__)
REGISTER_PORT(M_WAVEFRONTDATA, MC_SCALAR, MD_WOBJDATA_, __IOWAVEFRONTOBJ__HPP__)
REGISTER_PORT(M_STRINGFIELD, MC_SCALAR, MD_MPVECSTRING_, __IOWAVEFRONTOBJ__HPP__)
REGISTER_PORT(M_STRINGFIELD2, MC_SCALAR, MD_MPVECSTRING_, __IOWAVEFRONTOBJ__HPP__)
REGISTER_PORT(M_LONGFIELD, MC_SCALAR, MD_MPVECLONG_, __IOWAVEFRONTOBJ__HPP__)
REGISTER_PORT(M_LONGFIELD2, MC_SCALAR, MD_MPVECLONG_, __IOWAVEFRONTOBJ__HPP__)
REGISTER_PORT(M_GEOM2, MC_SCALAR, MD_MIMMO_, __IOWAVEFRONTOBJ__HPP__)
REGISTER_PORT(M_GEOM3, MC_SCALAR, MD_MIMMO_, __IOWAVEFRONTOBJ__HPP__)
REGISTER_PORT(M_GDISPLS, MC_SCALAR, MD_MPVECARR3FLOAT_, __IOWAVEFRONTOBJ__HPP__)
REGISTER_PORT(M_NAME, MC_SCALAR, MD_STRING, __IOWAVEFRONTOBJ__HPP__)
REGISTER_PORT(M_VECTORLI, MC_VECTOR, MD_LONG, __IOWAVEFRONTOBJ__HPP__)
REGISTER_PORT(M_VECTORLI2, MC_VECTOR, MD_LONG, __IOWAVEFRONTOBJ__HPP__)
REGISTER_PORT(M_VECTORLI3, MC_VECTOR, MD_LONG, __IOWAVEFRONTOBJ__HPP__)
REGISTER_PORT(M_VECTORLI4, MC_VECTOR, MD_LONG, __IOWAVEFRONTOBJ__HPP__)


REGISTER(BaseManipulation, IOWavefrontOBJ, "mimmo.IOWaveFrontOBJ")
REGISTER(BaseManipulation, ManipulateWFOBJData, "mimmo.ManipulateWFOBJData")

}

#endif /* __IOWAVEFRONTOBJ__HPP__ */
