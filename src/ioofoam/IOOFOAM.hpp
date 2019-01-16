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

#include "MimmoGeometry.hpp"
#include "MimmoFvMesh.hpp"
#include "enum.hpp"

BETTER_ENUM(IOOFMode, int, READ = 0, WRITE = 1, WRITEPOINTSONLY = 2);

namespace mimmo{

/*!
 * \class IOOFOAM
 * \ingroup ioofoam
 * \brief IOOFOAM is the class to import/export an OpenFOAM volume mesh alongside its boundary information.
 * 
 * The class is derived from MimmoFvMesh interface.
 *
 * It works as reader and writer. Its mode types are summarized in the following better enum
 * IOOFMode list:
 *  - <B>READ           </B> : import bulk and boundary mesh from an OpenFOAM case
 *  - <B>WRITE          </B> : export bulk and boundary mesh to an OpenFOAM case from scratch
 *  - <B>WRITEPOINTSONLY</B> : export only bulk mesh points coordinates to a pre-existent and compatible OpenFOAM case.
 *
 * Warning: WRITE mode is not available yet. It triggers the WRITEPOINTSONLY mode, 
 * which is a temporary mode to write only the coordinates of OpenFOAM mesh points.
 * 
 * Dependencies : OpenFOAM libraries (tested with OpenFOAM 2.4.x version).
 *
 * \n
 *
 * Ports available in IOOFOAM Class :
 *
 * Inherited from MimmoFvMesh:
 *
 *    =========================================================

   |                     Port Input   ||                                     |
   |------------------|---------------------|----------------------|
   | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
   | M_GEOM           | setGeometry         | (MC_SCALAR, MD_MIMMO_)     |
   | M_GEOM2          | setBoundaryGeometry | (MC_SCALAR, MD_MIMMO_)     |


   |               Port Output    ||                                         |
   |------------------|--------------------|-----------------------|
   | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
   | M_GEOM           | getGeometry        | (MC_SCALAR, MD_MIMMO_)      |
   | M_GEOM2          | getBoundaryGeometry| (MC_SCALAR, MD_MIMMO_)      |

 * Proper of the class:
 *
 |                     Port Input   ||                                     |
 |------------------|---------------------|----------------------|
 | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |


 |               Port Output    ||                                         |
 |------------------|--------------------|-----------------------|
 | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
 | M_GEOMOFOAM       | getGeometry               | (MC_SCALAR, MD_MIMMO_)      |
 | M_GEOMOFOAM2      | getBoundaryGeometry       | (MC_SCALAR, MD_MIMMO_)      |

 * =========================================================
 * \n
 *
 *
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation/MimmoFvMesh:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.IOOFOAM</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 *
 * Proper of the class:
 * - <B>IOMode</B>: activate mode of the class: READ, WRITE, WRITEPOINTSONLY;
 * - <B>Dir</B>: path to the current OpenFOAM mesh for reading/writing purposes;
 * - <B>Overwrite</B>: option valid only in WRITEPOINTSONLY mode: if 1-true overwrite 
                       points in the current OpenFoam case time of the mesh at WriteDir. 
                       If 0-false (DEFAULT) save them in a newly created case time at current time + 1;
                       
 * In case of writing mode Geometries have to be mandatorily passed by port.
 *
 */
class IOOFOAM: public MimmoFvMesh{
protected:
    int                             m_type;         /**<mode type of the class 0-read, 1-write, 2-writepointsonly*/
    std::string                     m_path;         /**< path to current mesh for reading/writing */
    bool                            m_overwrite;    /**< Overwrite in time case when in mode WRITEPOINTSONLY */

    std::unordered_map<std::string, bitpit::ElementType>    m_OFE_supp;     /**<list of openfoam shapes actually supported as it is, and not as generic polyhedron*/

public:
    IOOFOAM(int type = IOOFMode::READ);
    IOOFOAM(std::unique_ptr<MimmoObject> & bulk, std::unique_ptr<MimmoObject> &boundary, int type = IOOFMode::WRITE);
    IOOFOAM(const bitpit::Config::Section & rootXML);
    virtual ~IOOFOAM();

    IOOFOAM(const IOOFOAM & other);
    IOOFOAM & operator=(IOOFOAM other);
    
    void            buildPorts();

    void            setDir(std::string dir);
    void            setOverwrite(bool flag);
    bool            getOverwrite();
    void            execute();

    void            setGeometry(MimmoObject * bulk);
    void            setBoundaryGeometry(MimmoObject * boundary);

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

protected:
    void 			swap(IOOFOAM & x) noexcept;
    void            setDefaults();
    bool            write();
    bool            read();
    bool            writePointsOnly();

};


/*!
 * \class IOOFOAMField
 * \ingroup ioofoam
 * \brief IOOFOAMField is the base class to import/export an OpenFOAM field defined on a volume mesh and/or the related boundary fields.
 *
 * It works as reader and writer. Its mode types are summarized in the following better enum
 * IOOFMode list:
 *  - <B>READ           </B> : import bulk and/or boundary field from an OpenFOAM case
 *  - <B>WRITE          </B> : export bulk and/or boundary field to an OpenFOAM case from scratch.
 *
 * Warning: WRITE mode is not available yet.
 *
 * The read field can be an internal volume field or/and a surface field, in the latter case all the patches are read together.
 * The field is read in function of which geometry is set.
 *
 * Dependencies : OpenFOAM libraries (tested with OpenFOAM 2.4.x version).
 *
 * \n
 *
 * Ports available in IOOFOAMField Class :
 *
 * Inherited from MimmoFvMesh:
 *
 *    =========================================================

   |                     Port Input   ||                                     |
   |------------------|---------------------|----------------------|
   | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
   | M_GEOM           | setGeometry         | (MC_SCALAR, MD_MIMMO_)     |
   | M_GEOM2          | setBoundaryGeometry | (MC_SCALAR, MD_MIMMO_)     |


   |               Port Output    ||                                         |
   |------------------|--------------------|-----------------------|
   | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
   | M_GEOM           | getGeometry        | (MC_SCALAR, MD_MIMMO_)      |
   | M_GEOM2          | getBoundaryGeometry| (MC_SCALAR, MD_MIMMO_)      |


 * Proper of the class:
 *
 |                     Port Input   ||                                     |
 |------------------|---------------------|----------------------|
 | <B>PortType</B>  | <B>variable/function</B>  |<B>DataType</B> |
 | M_GEOMOFOAM       | setGeometry               | (MC_SCALAR, MD_MIMMO_)      |
 | M_GEOMOFOAM2      | setBoundaryGeometry       | (MC_SCALAR, MD_MIMMO_)      |


 |               Port Output    ||                                         |
 |------------------|--------------------|-----------------------|
 | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |

 * =========================================================
 * \n
 *
 *
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation/MimmoFvMesh:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.IOOFOAM</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 *
 * Proper of the class:
 * - <B>IOMode</B>: activate mode of the class: READ, WRITE;
 * - <B>Dir</B>: path to the current OpenFOAM field for reading/writing purposes;
 * - <B>Name</B>: name of the OpenFOAM field for reading/writing purposes;
 * - <B>Overwrite</B>: if 1-true overwrite the field in the current OpenFoam case time at WriteDir.
                       If 0-false (DEFAULT) save them in a newly created case time at current time + 1;

 * One between M_GEOMOF or M_GEOMOF2 (Geometries) is a mandatory port.
 *
 */
class IOOFOAMField: public MimmoFvMesh{
protected:
    int                             m_type;         /**<mode type of the class 0-read, 1-write*/
    std::string                     m_path;         /**< path to current mesh for reading field */
    std::string                     m_fieldname;    /**< name of current field for reading/writing */
    bool                            m_overwrite;    /**< Overwrite in time */

    std::unordered_map<std::string, bitpit::ElementType>    m_OFE_supp;     /**<list of openfoam shapes actually supported as it is, and not as generic polyhedron*/

public:
    IOOFOAMField(int type = IOOFMode::READ);
    virtual ~IOOFOAMField();

//    IOOFOAMField(const IOOFOAMField & other);
//    IOOFOAMField & operator=(IOOFOAMField other);

    void            buildPorts();

    void            setDir(std::string dir);
    void            setField(std::string fieldname);
    void            setOverwrite(bool flag);
    bool            getOverwrite();
    void			setType(int type);
    int				getType();
    void            execute();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

protected:
//    void swap(IOOFOAMField & x) noexcept;
    void            setDefaults();
    virtual bool    write() = 0;
    virtual bool    read() = 0;

};



/*!
 * \class IOOFOAMScalarField
 * \ingroup ioofoam
 * \brief IOOFOAMScalarField is the class derived from IOOFOAMField to import/export an OpenFOAM scalar field defined on a volume mesh and/or the related boundary fields.
 *
 * The read field can be an internal volume field or/and a surface field, in the latter case all the patches are read together.
 * The field is read in function of which geometry is set.
 *
 * Dependencies : OpenFOAM libraries (tested with OpenFOAM 2.4.x version).
 *
 * \n
 *
 * Ports available in IOOFOAMScalarField Class :
 *
 * Inherited from MimmoFvMesh:
 *
 *    =========================================================

   |                     Port Input   ||                                     |
   |------------------|---------------------|----------------------|
   | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
   | M_GEOM           | setGeometry         | (MC_SCALAR, MD_MIMMO_)     |
   | M_GEOM2          | setBoundaryGeometry | (MC_SCALAR, MD_MIMMO_)     |


   |               Port Output    ||                                         |
   |------------------|--------------------|-----------------------|
   | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
   | M_GEOM           | getGeometry        | (MC_SCALAR, MD_MIMMO_)      |
   | M_GEOM2          | getBoundaryGeometry| (MC_SCALAR, MD_MIMMO_)      |


 * Inherited from IOOFOAMField:
 *
 |                     Port Input   ||                                     |
 |------------------|---------------------|----------------------|
 | <B>PortType</B>  | <B>variable/function</B>  |<B>DataType</B> |
 | M_GEOMOFOAM       | setGeometry               | (MC_SCALAR, MD_MIMMO_)      |
 | M_GEOMOFOAM2      | setBoundaryGeometry       | (MC_SCALAR, MD_MIMMO_)      |


 |               Port Output    ||                                         |
 |------------------|--------------------|-----------------------|
 | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |

 * Proper of the class:
 *
 |                     Port Input   ||                                     |
 |------------------|---------------------|----------------------|
 | <B>PortType</B>  | <B>variable/function</B>  |<B>DataType</B> |


 |               Port Output    ||                                         |
 |------------------|--------------------|-----------------------|
 | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
 | M_SCALARFIELD  | getField      		| (MC_MPVECTOR, MD_FLOAT)       |
 | M_SCALARFIELD2 | getBoundaryField    | (MC_MPVECTOR, MD_FLOAT)       |

 * =========================================================
 * \n
 *
 *
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation/MimmoFvMesh:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.IOOFOAM</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 *
 * Proper of the class:
 * - <B>IOMode</B>: activate mode of the class: READ, WRITE;
 * - <B>Dir</B>: path to the current OpenFOAM mesh for reading/writing purposes;
 * - <B>Name</B>: name of the OpenFOAM field for reading/writing purposes;
 * - <B>Overwrite</B>: if 1-true overwrite the field in the current OpenFoam case time at WriteDir.
                       If 0-false (DEFAULT) save them in a newly created case time at current time + 1;
 *
 */
class IOOFOAMScalarField: public IOOFOAMField{
protected:
    dmpvector1D m_field;        		/**<Internal volume field. */
    dmpvector1D m_boundaryField;        /**<Surface boundary field. */

public:
    IOOFOAMScalarField(int type = IOOFMode::READ);
    IOOFOAMScalarField(const bitpit::Config::Section & rootXML);
    virtual ~IOOFOAMScalarField();

    IOOFOAMScalarField(const IOOFOAMScalarField & other);
    IOOFOAMScalarField & operator=(IOOFOAMScalarField other);

    void            buildPorts();

    dmpvector1D     getField();
    dmpvector1D     getBoundaryField();

protected:
    void swap(IOOFOAMScalarField & x) noexcept;
    void            setDefaults();
    virtual bool    write();
    virtual bool    read();

};





livector1D mapEleVConnectivity(const livector1D &, const bitpit::ElementType &);


REGISTER_PORT(M_GEOMOFOAM, MC_SCALAR, MD_MIMMO_, __IOOFOAM_HPP__)
REGISTER_PORT(M_GEOMOFOAM2, MC_SCALAR, MD_MIMMO_, __IOOFOAM_HPP__)
REGISTER_PORT(M_SCALARFIELD, MC_MPVECTOR, MD_FLOAT, __IOOFOAM_HPP__)
REGISTER_PORT(M_SCALARFIELD2, MC_MPVECTOR, MD_FLOAT, __IOOFOAM_HPP__)

REGISTER(BaseManipulation, IOOFOAM, "mimmo.IOOFOAM")

}

#endif /* __IOOFOAM_HPP__ */
