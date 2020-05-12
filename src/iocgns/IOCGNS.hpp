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

#ifndef __IOCGNS_HPP__
#define __IOCGNS_HPP__

#include "BaseManipulation.hpp"
#include <map>
#include <unordered_map>
#include <memory>

namespace mimmo{

/*! \ingroup typedefs{ */
typedef int M_CG_BCType_t;          /**<custom typedef for IOCGNS .*/
/*! } */

/*!
 * \class BCCGNS
 * \ingroup iocgns
 * \brief Class to store the Info related to boundary conditions of CGNS grid.
 *
 *  It works only if the boundary conditions are previously read from
    a CGNS file and passed by the related input port.
 *
 */
class BCCGNS{
    friend class IOCGNS;

protected:
    std::map<int, M_CG_BCType_t >       mcg_pidtobc;      /**<Converter of boundary conditions types (PID->conditionType).*/
    std::map<int, std::vector<int> >    mcg_zonetobndpid;  /**<map that associates volume zone to interested boundary pid*/
    std::map<int, int >               mcg_pidtolisttype;  /**< Track the list type for BC 0-PointList, 1-ElementList */
    std::map<int, std::string >       mcg_bcpidnames;     /**< map boundary conditions pids and their names if any*/
    std::map<int, std::string >       mcg_zonepidnames;   /**< map zone pids and their names if any*/

public:
    BCCGNS(){};
    ~BCCGNS(){
        clear();
    }

    void dump(std::ostream & out);
    void restore(std::istream & in);

private:

    void clear(){
        mcg_pidtobc.clear();
        mcg_zonetobndpid.clear();
        mcg_pidtolisttype.clear();
        mcg_bcpidnames.clear();
        mcg_zonepidnames.clear();
    }
};

/*!
 * \class IOCGNS
 * \ingroup iocgns
 * \brief IOCGNS is the class to import/ export an unstructured volume mesh written in CNGS.
 *
 * The object provides an interface to retrieve volume mesh and its boundary surfaces and
 * importing them in mimmo's native Data structures of type MimmoObject
 *
 * The class is in beta version and has the following limitations:
 * - Only unstructured mesh with single base supported.
 * - Singular of multi-zone meshes are supported in reading native cgns.
 * - Singular zone meshes are supported in writing cgns format.
 * - Multi-zone meshes cgns writing will be supported in future, but not yet available.
 * - Data attached to boundary conditions patches are not read/written, except for
 *   their names and their cgns type.
   - CGNS meshes exported from Pointwise16 and StarCCM++ are still unreadable with
     the current class. Errors are known and will be fixed in later versions.
   - Reading/writing partitioned mesh in cgns format is not yet available.
 *
 * Dependencies : cgns libraries.
 *
 * Ports available in IOCGNS Class :
 *
 *  =========================================================

   |                     Port Input  ||                                      |
   |------------------|------------------------|-------------------|
   | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
   | M_GEOM           | setGeometry            | (MC_SCALAR, MD_MIMMO_)  |
   | M_GEOM2          | setSurfaceBoundary     | (MC_SCALAR, MD_MIMMO_)  |
   | M_BCCGNS         | setBoundaryConditions  | (MC_SCALAR, MD_BCCGNS_)  |


   |               Port Output   ||                                          |
   |------------------|------------------------|-------------------|
   | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
   | M_GEOM           | getGeometry            | (MC_SCALAR, MD_MIMMO_)  |
   | M_GEOM2          | getSurfaceBoundary     | (MC_SCALAR, MD_MIMMO_)  |
   | M_BCCGNS         | getBoundaryConditions  | (MC_SCALAR, MD_BCCGNS_)  |

 *    =========================================================
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.IOCGNS</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;

 * Proper of the class:
 * - <B>IOCGNS_Mode</B>: working mode of the class, see IOCGNS_Mode enum.
 * - <B>Dir</B>: directory path for reading/writing;
 * - <B>Filename</B>: name of file to read/write without .cgns/dump tag;
 * - <B>WriteInfo</B>: boolean (1/0) write on file zoneNames, bcNames, either in reading and writing mode. The save directory path is specified with Dir.
 * - <B>WriteFormat</B>: writing format supported by the class, see IOCGNS_WriteType enum
 * - <B>WriteMultiZone</B>: 0- write single zone, 1- write multizone(if multi zone are available in the mesh).
 *
 * Geometry has to be mandatorily read or passed through port.
 *
 */
class IOCGNS: public BaseManipulation{

public:

    /*!
        \ingroup iocgns
        \brief enumeration of IOCGNS class I/O modes.
    */
    enum IOCGNS_Mode{
        READ = 0        , /**< 0 - 0rank only, Read the cgns from file */
        RESTORE=1       , /**< 1 - Read the mesh from a previous <>.dump file*/
        WRITE=2         , /**< 2 - 0rank only, Write the mesh on a cgns file */
        DUMP=3            /**< 3 - Write the mesh to dump file*/
    };

    /*!
        \ingroup iocgns
        \brief enumeration of IOCGNS class format type for writing only.
    */
    enum IOCGNS_WriteType{
        HDF5 = 1        , /**< 1-write cgns in default HDF5 mode */
        ADF  = 2        , /**< 2-write cgns in ADF mode*/
        ADF2 = 3          /**< 3-write cgns in old ADF2 mode for 2.5 cgns on 32 bit system only*/
    };

    IOCGNS(IOCGNS_Mode mode= IOCGNS_Mode::READ);
    IOCGNS(const bitpit::Config::Section & rootXML);
    ~IOCGNS();

    IOCGNS(const IOCGNS & other);
    IOCGNS & operator=(IOCGNS other);

    void            setDefaults();
    void            buildPorts();

    MimmoSharedPointer<MimmoObject>    getSurfaceBoundary();
    MimmoSharedPointer<MimmoObject>    getGeometry();
    BCCGNS*         getBoundaryConditions();

    IOCGNS_Mode     whichMode();
    int             whichModeInt();

    IOCGNS_WriteType  whatWritingFormat();
    bool              isWritingMultiZone();

    void            setDir(const std::string &dir);
    void            setFilename(const std::string &filename);

    void            setGeometry(MimmoSharedPointer<MimmoObject>);
    void            setSurfaceBoundary(MimmoSharedPointer<MimmoObject>);
    void            setBoundaryConditions(BCCGNS*);

    void            setWritingFormat(IOCGNS_WriteType type);
    void            setWritingMultiZone(bool multizone);
    void            setWriteOnFileMeshInfo(bool write);

    void            execute();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

protected:
    void            swap(IOCGNS &) noexcept;
    bool            write(const std::string & file);
    bool            read(const std::string & file);
    bool            dump(std::ostream & stream);
    bool            restore(std::istream & stream);
    bool            belongToPool(const bitpit::ConstProxyVector<long> & elementconn, const std::set<long> &pool);

    void            unpackMixedConns( const ivector1D & connsArray,
                                     MimmoSharedPointer<MimmoObject> patchVol,
                                     std::unordered_map<long, std::vector<long>> surfElem,
                                     const long & idVertexOffset,
                                     const long & idCellOffset,
                                     const long & PIDZoneVolume,
                                     long & idwork);

    void            writeInfoFile();
    std::map<int, std::vector<std::size_t> >
                    getZoneConn(const livector1D& cellIds,
                                const std::unordered_map<long,int> & mapToLocVert,
                                std::map<int, std::size_t> & ncells);
    std::map<int, std::vector<std::size_t> >
                    getBCElementsConn(const livector1D& cellIds,
                                      const std::unordered_map<long,int> & mapToLocVert,
                                      std::map<int, std::size_t> &ncells,
                                      std::unordered_map<long, int> & surfCellGlobToLoc);

#if MIMMO_ENABLE_MPI
    void communicateAllProcsStoredBC();
#endif

private:
    IOCGNS_Mode      m_mode;       /**<Mode of execution.See setMode configuration.*/
    IOCGNS_WriteType m_wtype;      /**<Writing type HDF5 or ADF*/
    bool             m_multizone;  /**<Writing multizone true, or single zone false*/

    std::string     m_dir;         /**<Name of directory path*/
    std::string     m_filename;    /**<Name of file */

    MimmoSharedPointer<MimmoObject> m_volmesh;          /**<Original volume mesh, instantiated in reading */
    MimmoSharedPointer<MimmoObject> m_surfmesh;         /**<Original boundary mesh, instantiated in reading */
    MimmoSharedPointer<MimmoObject> m_surfmesh_not;     /**<Pointed external boundary mesh*/

    std::unique_ptr<BCCGNS>             m_storedBC;         /**<Information of boundary conditions of a CGNS read mesh.*/

    bool    m_writeOnFile;                                  /**! write mesh info on file */
    std::unordered_map<int, std::string> m_elementsSectionName;  /**< facility to giv name to group of elements in writing mode*/
};


REGISTER_PORT(M_GEOM, MC_SCALAR, MD_MIMMO_, __IOCGNS_HPP__)
REGISTER_PORT(M_GEOM2, MC_SCALAR, MD_MIMMO_,__IOCGNS_HPP__)
REGISTER_PORT(M_BCCGNS, MC_SCALAR, MD_BCCGNS_,__IOCGNS_HPP__)

REGISTER(BaseManipulation, IOCGNS,"mimmo.IOCGNS")

}

#endif /* __IOCGNS_HPP__ */
