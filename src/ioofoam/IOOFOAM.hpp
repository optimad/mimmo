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
#ifndef __IOOFOAM_HPP__
#define __IOOFOAM_HPP__

#include "BaseManipulation.hpp"

namespace mimmo{

/*!
 * \class IOOFOAM
 * \ingroup ioofoam
 * \brief IOOFOAM is the class to import/export and modify an OpenFOAM volume mesh as mimmo data structure MimmoObject.
 *
 * Data Field associated to the mesh and boundary information can be read or written/modified, directy using the class interface.
 *
 * Dependencies : OpenFOAM libraries (tested with OpenFOAM 2.4.x version).
 *
 * \n
 *
 * Ports available in IOOFOAM Class :
 *
 *    =========================================================

   |                     Port Input   ||                                     |
   |------------------|---------------------|----------------------|
   | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
   | M_GEOM           | setGeometry         | (MC_SCALAR, MD_MIMMO_)     |
   
 
   | M_SCALARFIELD    | setField            | (MC_MPVECTOR, MD_FLOAT)      |
   | M_GEOM2          | setSurfaceBoundary  | (MC_SCALAR, MD_MIMMO_)     |


   |               Port Output    ||                                         |
   |------------------|--------------------|-----------------------|
   | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
   | M_GEOM           | getGeometry        | (MC_SCALAR, MD_MIMMO_)      |
   
   | M_SCALARFIELD    | getField           | (MC_MPVECTOR, MD_FLOAT)       |
   | M_GEOM2          | getSurfaceBoundary | (MC_SCALAR, MD_MIMMO_)      |

 * =========================================================
 * \n
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.IOOFOAM</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 *
 * Proper of the class:
 * - <B>ReadFlag</B>: activate reading mode boolean;
 * - <B>ReadDir</B>: path to the current OpenFOAM mesh for reading purposes;
 * - <B>WriteFlag</B>: activate writing mode boolean;
 * - <B>WriteDir</B>: path to the current OpenFOAM mesh for writing purposes;
 *
 * In case of writing mode Geometries have to be mandatorily passed by port.
 *
 */
class IOOFOAM: public BaseManipulation{
private:
    bool                            m_read;         /**<If true it reads the geometry from files during the execution.*/
    std::string                     m_readPath;     /**< path to current mesh for reading */
    bool                            m_write;         /**<If true it writes the geometries on file during the execution*/
    std::string                     m_writePath;     /**< path to current mesh for writing */

    std::unique_ptr<MimmoObject>    m_volmesh;      /**<Original volume mesh, instantiated in reading */
//     std::unique_ptr<MimmoObject>    m_surfmesh;     /**<Original boundary mesh, instantiated in reading */
//     MimmoObject*                    m_surfmesh_ext; /**<Pointed external boundary mesh */

//     dmpvector1D                     m_field;        /**<Scalar field related to polydata mesh (pint located values).*/
//     bool                            m_normalize;    /**<If true the field is normalized with the maximum absolute value after reading.*/
//     double                          m_maxf;         /**<Max value of the read scalar field.*/
//     double                          m_scaling;      /**<Value used to scale the scalar field (default m_scaling = 1). */
    std::unordered_map<std::string, bitpit::ElementType>    m_OFE_supp;     /**<list of openfoam shapes actually supported as it is, and not as generic polyhedron*/
 
public:
    IOOFOAM();
    IOOFOAM(const bitpit::Config::Section & rootXML);
    ~IOOFOAM();

    IOOFOAM(const IOOFOAM & other);
    IOOFOAM & operator=(IOOFOAM other);
    
    void            buildPorts();
    void            setDefaults();

    void            setRead(bool read);
    void            setReadDir(std::string dir);
    void            setWrite(bool write);
    void            setWriteDir(std::string dir);
    void            setGeometry(MimmoObject* geom);

//     MimmoObject*    getSurfaceBoundary();
    MimmoObject*    getGeometry();

    std::unique_ptr<MimmoObject> cloneInternalMesh();

    bool            write();
    bool            read();

    void             execute();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

protected:
    void swap(IOOFOAM & x) noexcept;
    livector1D mapEleVConnectivity(const livector1D &, const bitpit::ElementType &);
    
};

REGISTER_PORT(M_GEOM, MC_SCALAR, MD_MIMMO_, __IOOFOAM_HPP__)
// REGISTER_PORT(M_GEOM2, MC_SCALAR, MD_MIMMO_, __IOOFOAM_HPP__)
// REGISTER_PORT(M_SCALARFIELD, MC_MPVECTOR, MD_FLOAT, __IOOFOAM_HPP__)


REGISTER(BaseManipulation, IOOFOAM,"mimmo.IOOFOAM")

}

#endif /* __IOOFOAM_HPP__ */
