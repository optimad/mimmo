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
//#include <cgnslib.h>


namespace mimmo{

typedef int M_CG_ElementType_t;     /**<Casting of CG_ElementType_t to int.*/
typedef int M_CG_BCType_t;          /**<Casting of CG_BCType_t to int.*/

/*!
 * \class InfoCGNS
 * \ingroup iocgns
 *  \brief Class to store the Info related to MimmoObject in order to export in CGNS grid format. Only one zone supported.
 *
 */
class InfoCGNS{
    friend class IOCGNS;
private:

    std::map<bitpit::ElementInfo::Type, M_CG_ElementType_t>   mcg_typeconverter; /**<Element type conversion bitpit->CGNS.*/
    std::map<M_CG_ElementType_t, std::string>                 mcg_typetostring;  /**<Element type conversion CGNS->sting.*/

    std::map<M_CG_ElementType_t, int>                         mcg_number;     /**<Number of elements in each volume section per type.*/
    std::map<M_CG_ElementType_t, std::vector<long int> >      mcg_typetoid;   /**<Map of elements (bitpit ID vector) in each volume section per type.*/
    std::map<M_CG_ElementType_t, std::vector<int> >           mcg_typetoconn; /**<Map of connectivity of elements (bitpit vertex ID vector) in each volume section per type.*/
    std::map<long int, int>                                 mcg_idtoindex;  /**<Map of elements bitpit ID to local CGNS index.*/
    std::vector<long int>                                   mcg_indextoid;  /**<Vector of elements local CGNS index to bitpit ID.*/

    std::map<M_CG_ElementType_t, int>                         mcg_bndnumber;      /**<Number of elements in each boundary section per type.*/
    std::map<M_CG_ElementType_t, std::vector<long int> >      mcg_bndtypetoid;    /**<Map of boundary elements (bitpit ID vector) in each boundary section per type.*/
    std::map<M_CG_ElementType_t, std::vector<int> >           mcg_bndtypetoconn;  /**<Map of connectivity of boundary elements (bitpit vertex ID vector) in each boundary section per type.*/
    std::map<long int, int>                                 mcg_bndidtoindex;   /**<Map of boundary elements bitpit ID to local CGNS index.*/
    std::vector<long int>                                   mcg_bndindextoid;   /**<Vector of boundary elements local CGNS index to bitpit ID.*/

public:
    InfoCGNS(){};
    ~InfoCGNS()
    {
        clear();
    }
    
private:
    /*!
     * Clear members of the class.
     */
    void
    clear(){
        mcg_number.clear();
        mcg_typetoid.clear();
        mcg_typetoconn.clear();
        mcg_idtoindex.clear();
        mcg_indextoid.clear();
        mcg_bndnumber.clear();
        mcg_bndtypetoid.clear();
        mcg_bndtypetoconn.clear();
        mcg_bndidtoindex.clear();
        mcg_bndindextoid.clear();
    }

};

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
    std::map<int, M_CG_BCType_t >         mcg_pidtobc;        /**<Converter of boundary conditions types (PID->conditionType).*/
    std::map<int, std::string >         mcg_pidtoname;       /**<Converter of boundary conditions types (PID->nameBC).*/
    std::map<int, std::vector<int> >    mcg_bndpidtoindex;  /**<Map of boundary elements bitpit PID to vector of local CGNS index.*/

public:
    BCCGNS(){};
    ~BCCGNS()
    {
        mcg_pidtobc.clear();
        mcg_pidtoname.clear();
        mcg_bndpidtoindex.clear();
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
 * - When imported, a Volume mesh will be absorbed as Point Cloud MimmoObject;
 * - When imported boundary surface meshes will be absorbed as surface triangulations, even if the original
 *  CGNS allows mixed element types.
 * - Only unstructured mesh with single base and single zone supported.
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
 *
 * Proper of the class:
 * - <B>ReadFlag</B>: activate reading mode boolean 1-reading mode, 0-writing mode;
 * - <B>ReadDir</B>: reading directory path;
 * - <B>ReadFilename</B>: name of file for reading;
 * - <B>WriteDir</B>: writing directory path;
 * - <B>WriteFilename</B>: name of file for writing;
 *
 * Geometry has to be mandatorily read or passed through port.
 *
 */
class IOCGNS: public BaseManipulation{
private:
    bool            m_read;         /**<If true it reads the geometry from file during the execution.*/
    bool            m_write;        /**<If true it writes the geometry on file during the execution.*/

    std::string     m_rdir;         /**<Name of directory to read the geometry (without final "/").*/
    std::string     m_rfilename;    /**<Name of file to read the geometry (without extension).*/
    std::string     m_wdir;         /**<Name of directory to write the geometry (without final "/").*/
    std::string     m_wfilename;    /**<Name of file to write the geometry.*/

    std::unique_ptr<MimmoObject>        m_volmesh;          /**<Original volume mesh, instantiated in reading */
    std::unique_ptr<MimmoObject>        m_surfmesh;         /**<Original boundary mesh, instantiated in reading */
    MimmoObject *                       m_surfmesh_not;     /**<Pointed external boundary mesh*/

    std::unique_ptr<InfoCGNS>           m_storedInfo;       /**<Information of a CGNS read mesh.*/
    std::unique_ptr<BCCGNS>             m_storedBC;         /**<Information of boundary conditions of a CGNS read mesh.*/

public:
    IOCGNS(bool read = false);
    IOCGNS(const bitpit::Config::Section & rootXML);
    ~IOCGNS();

    IOCGNS(const IOCGNS & other);
    IOCGNS & operator=(IOCGNS other);

    void            setDefaults();
    void            buildPorts();

    MimmoObject*    getSurfaceBoundary();
    MimmoObject*    getGeometry();
    BCCGNS*         getBoundaryConditions();

    bool            isReadingMode();
    bool            isWritingMode();

    void            setReadDir(std::string dir);
    void            setRead(bool read);
    void            setWriteDir(std::string dir);
    void            setReadFilename(std::string filename);
    void            setWrite(bool write);
    void            setWriteFilename(std::string filename);

    void            setGeometry(MimmoObject*);
    void            setSurfaceBoundary(MimmoObject*);
    void            setBoundaryConditions(BCCGNS*);
    void            execute();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

protected:
    void            swap(IOCGNS &) noexcept;
    bool            write();
    bool            read();

private:
    void    unpack3DElementsMixedConns(MimmoObject*,MimmoObject*, ivector1D &, long &startId);
    void    recoverCGNSInfo();

};

REGISTER(BaseManipulation, IOCGNS,"mimmo.IOCGNS")

}

#endif /* __IOCGNS_HPP__ */
