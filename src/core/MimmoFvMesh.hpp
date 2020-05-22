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
#ifndef __MIMMOFVMESH_HPP__
#define __MIMMOFVMESH_HPP__

#include "MimmoObject.hpp"
#include "BaseManipulation.hpp"

namespace mimmo{

/*!
 * \struct InfoBoundaryPatch
 * \brief struct containing basic information on boundary patch defined in MimmoFvMesh class.
 * \ingroup core
 */
struct InfoBoundaryPatch{

    std::string name;  /**< name of boundary patch */
    std::string type;  /**< type of boundary condition string identifier */

    InfoBoundaryPatch();
    ~InfoBoundaryPatch();
    /*! Copy constructor*/
    InfoBoundaryPatch(const InfoBoundaryPatch &) = default;
    /*! Assignment Operator */
    InfoBoundaryPatch & operator=(const InfoBoundaryPatch &) = default;
};

/*!
* \class MimmoFvMesh
* \brief MimmoFvMesh is an abstract executable class for handling
* 3D and 2D mesh with their boundary information
* \ingroup core

* Basically it is a wrapping container of two MimmoObject,
  one holding the bulk volume/surface mesh information,
* and the other holding esplicitly the mesh boundaries.
* The link between bulk and boundaries is guaranteed with a one-to-one
  correspondance between ids of bulk Interfaces at mesh border and boundary
  mesh Cells. Multi-patch subdivision of the boundary mesh is achieved
  assigning PID on boundary mesh cells.
* External Additional info can be appended to each boundary patch.
*
* Ports available in MimmoFvMesh class:
*
* =========================================================
*  |Port Input  | | |
*  |-|-|-|
*  |<B>PortType</B>|<B>variable/function</B>|<B>DataType</B>|
*  | M_GEOM      | setGeometry                           | (MC_SCALAR, MD_MIMMO_)      |
*  | M_GEOM2     | setBoundaryGeometry                   | (MC_SCALAR, MD_MIMMO_)      |
*
*
*  | Port Output | | |
*  |-|-|-|
*  |<B>PortType</B>|<B>variable/function</B>|<B>DataType</B>|
*  | M_GEOM      | getGeometry                           | (MC_SCALAR, MD_MIMMO_)      |
*  | M_GEOM2     | getBoundaryGeometry                   | (MC_SCALAR, MD_MIMMO_)      |
*
* =========================================================
*
* The xml available parameters, sections and subsections  are the following:
*
* Inherited from BaseManipulation:
* - <B>ClassName</B>: name of the class as <tt>mimmo.MimmoFvMesh</tt>;
* - <B>Priority</B>: uint marking priority in multi-chain execution;
*
* Geometries has to be mandatorily passed through port.
*/
class MimmoFvMesh: public BaseManipulation{

protected:
    MimmoSharedPointer<MimmoObject>            m_bulk;              /**<Reference to bulk MimmoObject. */
    MimmoSharedPointer<MimmoObject>            m_boundary;          /**<Reference to boundary MimmoObject. */

    //members
    std::unordered_map<long, InfoBoundaryPatch>   m_infoBoundary; /**< list of additional info on boundary patches */

public:
    MimmoFvMesh();
    MimmoFvMesh(MimmoSharedPointer<MimmoObject> bulk, MimmoSharedPointer<MimmoObject> boundary);
    virtual ~MimmoFvMesh();

    MimmoFvMesh(const MimmoFvMesh & other);

    void    setGeometry(MimmoSharedPointer<MimmoObject> bulk);
    void    setBoundaryGeometry(MimmoSharedPointer<MimmoObject> boundary);
    MimmoSharedPointer<MimmoObject> getGeometry();
    MimmoSharedPointer<MimmoObject> getBoundaryGeometry();

    void                addInfoBoundaryPatch(const long & PID, const InfoBoundaryPatch & info);
    InfoBoundaryPatch   getInfoBoundaryPatch(const long & PID);

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name = "");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name= "");

protected:
    void    swap(MimmoFvMesh & ) noexcept;
    void    createBoundaryMesh();
    bool    checkMeshCoherence();
    void buildPorts();

};

REGISTER_PORT(M_GEOM, MC_SCALAR, MD_MIMMO_, __MIMMOFVMESH_HPP__)
REGISTER_PORT(M_GEOM2, MC_SCALAR, MD_MIMMO_, __MIMMOFVMESH_HPP__)

}

#endif /* __MIMMOFVMESH_HPP__*/
