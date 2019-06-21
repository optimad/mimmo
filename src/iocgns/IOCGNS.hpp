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


namespace mimmo{

typedef int M_CG_ElementType_t;     /**<Casting of CG_ElementType_t to int.*/
typedef int M_CG_BCType_t;          /**<Casting of CG_BCType_t to int.*/

// /*!
//  * \class InfoCGNS
//  * \ingroup iocgns
//  *  \brief Class to store the Info related to MimmoObject in order to export in CGNS grid format. Only one zone supported.
//  *
//  */
// class InfoCGNS{
//     friend class IOCGNS;
// private:
//
//     std::map<bitpit::ElementType, M_CG_ElementType_t>   mcg_typeconverter; /**<Element type conversion bitpit->CGNS.*/
//     std::map<M_CG_ElementType_t, std::string>                 mcg_typetostring;  /**<Element type conversion CGNS->sting.*/
//
//     std::map<M_CG_ElementType_t, int>                         mcg_number;     /**<Number of elements in each volume section per type.*/
//     std::map<M_CG_ElementType_t, std::vector<long int> >      mcg_typetoid;   /**<Map of elements (bitpit ID vector) in each volume section per type.*/
//     std::map<M_CG_ElementType_t, std::vector<int> >           mcg_typetoconn; /**<Map of connectivity of elements (bitpit vertex ID vector) in each volume section per type.*/
//     std::map<long int, int>                                 mcg_idtoindex;  /**<Map of elements bitpit ID to local CGNS index.*/
//     std::vector<long int>                                   mcg_indextoid;  /**<Vector of elements local CGNS index to bitpit ID.*/
//
//     std::map<M_CG_ElementType_t, int>                         mcg_bndnumber;      /**<Number of elements in each boundary section per type.*/
//     std::map<M_CG_ElementType_t, std::vector<long int> >      mcg_bndtypetoid;    /**<Map of boundary elements (bitpit ID vector) in each boundary section per type.*/
//     std::map<M_CG_ElementType_t, std::vector<int> >           mcg_bndtypetoconn;  /**<Map of connectivity of boundary elements (bitpit vertex ID vector) in each boundary section per type.*/
//     std::map<long int, int>                                 mcg_bndidtoindex;   /**<Map of boundary elements bitpit ID to local CGNS index.*/
//     std::vector<long int>                                   mcg_bndindextoid;   /**<Vector of boundary elements local CGNS index to bitpit ID.*/
//
// public:
//     InfoCGNS(){};
//     ~InfoCGNS()
//     {
//         clear();
//     }
//
//     void dump(std::ostream & out);
//     void restore(std::istream & in);
//
// private:
//     /*!
//      * Clear members of the class.
//      */
//     void
//     clear(){
//         mcg_number.clear();
//         mcg_typetoid.clear();
//         mcg_typetoconn.clear();
//         mcg_idtoindex.clear();
//         mcg_indextoid.clear();
//         mcg_bndnumber.clear();
//         mcg_bndtypetoid.clear();
//         mcg_bndtypetoconn.clear();
//         mcg_bndidtoindex.clear();
//         mcg_bndindextoid.clear();
//     }
//
// };

/*!
 * \class BCCGNS
 * \ingroup iocgns
 * \brief Class to store the Info related to boundary conditions of CGNS grid.
 *
 *  It works only if the boundary conditions are previously read from a CGNS file and passed by the related input port.
 *
 */
class BCCGNS{
    friend class IOCGNS;
private:
    std::map<int, M_CG_BCType_t >       mcg_pidtobc;      /**<Converter of boundary conditions types (PID->conditionType).*/
    std::map<int, std::string >         mcg_pidtoname;      /**<Converter of boundary conditions types (PID->nameBC).*/
    std::map<int, std::vector<int> >    mcg_zonetobndpid;  /**<map that associates volume zone to interested boundary pid*/

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
        mcg_pidtoname.clear();
        mcg_zonetobndpid.clear();
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
 * - Singular of multi-zone meshes are supported.
 * - Data attached to boundary conditions patches are not read/written, except for
 *   their names and their cgns type.
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
 * - <B>IOCGNS_Mode</B>: 0-read from file cgns, 1-read from file cgns and save dump file,
 *                       2-read from dump, 3-write grid on cgns file, 4-write on dump, 5-convert dump to cgns file.
 * - <B>Dir</B>: directory path for reading/writing; In case of IOCGNS_Mode 1, cgns file is read here and the dump file will be saved here
 * - <B>Filename</B>: name of file to read/write without .cgns tag; In case of IOCGNS_Mode 1, this is the name of the cgns file. The dump file will be saved with the same name + "_MDUMP.dump"
 * - <B>WriteInfo</B>: boolean (1/0) write on file zoneNames, bcNames, either in reading and writing mode. The save directory path is specified with Dir.
 *
 * Geometry has to be mandatorily read or passed through port.
 *
 */
class IOCGNS: public BaseManipulation{

public:

    enum IOCGNS_Mode{
        READ = 0        , /**< Read the cgns from file */
        READANDDUMP = 1 , /**< Read the cgns from file and dump it in file <>.dump*/
        RESTORE = 2     , /**< Read the mesh from a previous <>.dump file*/
        WRITE=3         , /**< Write the mesh on a cgns file */
        WRITEONDUMP =4  , /**< Write the mesh on a <>.dump file */
        CONVERTFROMDUMP=5 /**< Convert a mesh read from a dump file in the cgns native format*/
    };

    IOCGNS(IOCGNS_Mode mode= IOCGNS_Mode::READ);
    IOCGNS(const bitpit::Config::Section & rootXML);
    ~IOCGNS();

    IOCGNS(const IOCGNS & other);
    IOCGNS & operator=(IOCGNS other);

    void            setDefaults();
    void            buildPorts();

    MimmoObject*    getSurfaceBoundary();
    MimmoObject*    getGeometry();
    BCCGNS*         getBoundaryConditions();

    IOCGNS_Mode     whichMode();
    int             whichModeInt();

    void            setDir(const std::string &dir);
    void            setMode(IOCGNS_Mode mode);
    void            setFilename(const std::string &filename);

    void            setGeometry(MimmoObject*);
    void            setSurfaceBoundary(MimmoObject*);
    void            setBoundaryConditions(BCCGNS*);

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

//    void    unpack3DElementsMixedConns(MimmoObject*,MimmoObject*, ivector1D &, long &startId);
//    void    recoverCGNSInfo();
    void writeInfoFile(); 
    std::unordered_map<int, std::vector<std::size_t> >
                getZoneConn(const livector1D& cellIds, const std::unordered_map<long,int> & mapToLocVert);
private:
    IOCGNS_Mode    m_mode;         /**<Mode of execution.See setMode configuration.*/

    std::string     m_dir;         /**<Name of directory path*/
    std::string     m_filename;    /**<Name of file */

    std::unique_ptr<MimmoObject>        m_volmesh;          /**<Original volume mesh, instantiated in reading */
    std::unique_ptr<MimmoObject>        m_surfmesh;         /**<Original boundary mesh, instantiated in reading */
    MimmoObject *                       m_surfmesh_not;     /**<Pointed external boundary mesh*/

//    std::unique_ptr<InfoCGNS>           m_storedInfo;       /**<Information of a CGNS read mesh.*/
    std::unique_ptr<BCCGNS>             m_storedBC;         /**<Information of boundary conditions of a CGNS read mesh.*/

    bool    m_writeOnFile; /**! write mesh info on file */

};


REGISTER_PORT(M_GEOM, MC_SCALAR, MD_MIMMO_, __IOCGNS_HPP__)
REGISTER_PORT(M_GEOM2, MC_SCALAR, MD_MIMMO_,__IOCGNS_HPP__)
REGISTER_PORT(M_BCCGNS, MC_SCALAR, MD_BCCGNS_,__IOCGNS_HPP__)

REGISTER(BaseManipulation, IOCGNS,"mimmo.IOCGNS")

}

#endif /* __IOCGNS_HPP__ */
