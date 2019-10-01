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

#include "MimmoFvMesh.hpp"

#include "enum.hpp"
#include <unordered_map>

BETTER_ENUM(IOOFMode, int, READ = 0, WRITE = 1, WRITEPOINTSONLY = 2);

namespace mimmo{

/*!
 * \class IOOFOAM_Kernel
 * \ingroup ioofoam
 * \brief Abstract class to import/export an OpenFOAM volume mesh alongside its boundary information
 *        or field attached to the mesh.
 *
 * The class is derived from MimmoFvMesh interface.
 *
 * Dependencies : OpenFOAM libraries (tested with OpenFOAM Fundation official releases from 2.4.x up to 6)
 *
 * \n
 *
 * Ports available in IOOFOAM_Kernel Class :
 *
 |                     Port Input   ||                                 |
 |------------------|---------------------|----------------------------|
 | <B>PortType</B>  | <B>variable/function</B>  |<B>DataType</B>       |
 | M_GEOMOFOAM      | setGeometry         | (MC_SCALAR, MD_MIMMO_)     |
 | M_GEOMOFOAM2     | setBoundaryGeometry | (MC_SCALAR, MD_MIMMO_)     |
 | M_UMAPIDS        | setFacesMap         | (MC_UMAP, MD_LONG) |


 |               Port Output    ||                                           |
 |-------------------|---------------------------|---------------------------|
 | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B>            |
 | M_GEOMOFOAM       | getGeometry               | (MC_SCALAR, MD_MIMMO_)    |
 | M_GEOMOFOAM2      | getBoundaryGeometry       | (MC_SCALAR, MD_MIMMO_)    |
 | M_UMAPIDS         | getFacesMap               | (MC_UMAP, MD_LONG)|

 * =========================================================
 *
 * \n
 *
 * The xml available parameters, sections and subsections are specified each time
 * in the proper derived classes.
 *
 */
class IOOFOAM_Kernel: public MimmoFvMesh{

protected:
    int                      m_type;         /**<mode type of the class see IOOFMode enum */
    std::string              m_path;         /**< path to current mesh for reading/writing */
    std::unordered_map<std::string, bitpit::ElementType>    m_OFE_supp;     /**<list of openfoam shapes actually supported as it is, and not as generic polyhedron*/
	std::unordered_map<long,long> 	m_OFbitpitmapfaces; /**< OpenFoam faces -> bitpit Interfaces map. Used to to detect boundaries correspondence. */
    std::string                     m_fieldname;    /**< name of current field for reading/writing */

public:
    IOOFOAM_Kernel(int type = IOOFMode::READ);
    IOOFOAM_Kernel(std::unique_ptr<MimmoObject> & bulk, std::unique_ptr<MimmoObject> &boundary, int type = IOOFMode::WRITE);
    IOOFOAM_Kernel(const IOOFOAM_Kernel & other);
    virtual ~IOOFOAM_Kernel();

    virtual void                    buildPorts();
    std::unordered_map<long,long>   getFacesMap();
    int                             getType();


    void            setDir(const std::string &dir);
    void            setGeometry(MimmoObject * bulk);
    void            setBoundaryGeometry(MimmoObject * boundary);
    void            setFacesMap(std::unordered_map<long,long> mapFaces);
    void            setFieldName(const std::string & fieldname);
    void            setType(int type);
    void            setType(IOOFMode type);

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

protected:
    void 			swap(IOOFOAM_Kernel & x) noexcept;
    virtual void    setDefaults();
    virtual bool    write() = 0;
    virtual bool    read() = 0;

};


/*!
 \class IOOFOAM
 \ingroup ioofoam
 \brief IOOFOAM_Kernel specialization to import/export an OpenFOAM volume mesh only alongside its boundary information.

 The class is derived from IOOFOAM_Kernel interface.

 It works as reader and writer. Its mode types are summarized in the following better enum
 IOOFMode list:
  - <B>READ           </B> : import bulk and boundary mesh from an OpenFOAM case
  - <B>WRITE          </B> : export bulk and boundary mesh to an OpenFOAM case from scratch
  - <B>WRITEPOINTSONLY</B> : export only bulk mesh points coordinates to a pre-existent and compatible OpenFOAM case.

 Warning: WRITE mode is not available yet. It triggers the WRITEPOINTSONLY mode,
 which is a temporary mode to write only the coordinates of OpenFOAM mesh points.

 Ports inherited from IOOFOAM_Kernel Class :
*
|                     Port Input   ||                                 |
|------------------|---------------------|----------------------------|
| <B>PortType</B>  | <B>variable/function</B>  |<B>DataType</B>       |
| M_GEOMOFOAM      | setGeometry         | (MC_SCALAR, MD_MIMMO_)     |
| M_GEOMOFOAM2     | setBoundaryGeometry | (MC_SCALAR, MD_MIMMO_)     |
| M_UMAPIDS        | setFacesMap         | (MC_UMAP, MD_LONG) |


|               Port Output    ||                                           |
|-------------------|---------------------------|---------------------------|
| <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B>            |
| M_GEOMOFOAM       | getGeometry               | (MC_SCALAR, MD_MIMMO_)    |
| M_GEOMOFOAM2      | getBoundaryGeometry       | (MC_SCALAR, MD_MIMMO_)    |
| M_UMAPIDS         | getFacesMap               | (MC_UMAP, MD_LONG)|

* =========================================================
* The xml available parameters, sections and subsections are the following :
*
* - <B>ClassName</B>: name of the class as <tt>mimmo.IOOFOAM</tt>;
* - <B>Priority</B>: uint marking priority in multi-chain execution;
* - <B>IOMode</B>: activate mode of the class: READ, WRITE, WRITEPOINTSONLY;
* - <B>Dir</B>: path to the current OpenFOAM mesh for reading/writing purposes;
* - <B>Overwrite</B>: option valid only in WRITEPOINTSONLY mode: if 1-true overwrite
                      points in the current OpenFoam case time of the mesh at WriteDir.
                      If 0-false (DEFAULT) save them in a newly created case time at current time + 1;

* In case of writing mode Geometries have to be mandatorily passed by port.
*
*/
class IOOFOAM: public IOOFOAM_Kernel{

protected:
    bool        m_overwrite;    /**< Overwrite in time case when in mode WRITEPOINTSONLY */

public:
   IOOFOAM(int type = IOOFMode::READ);
   IOOFOAM(std::unique_ptr<MimmoObject> & bulk, std::unique_ptr<MimmoObject> &boundary, int type = IOOFMode::WRITE);
   IOOFOAM(const bitpit::Config::Section & rootXML);
   virtual ~IOOFOAM();

   IOOFOAM(const IOOFOAM & other);
   IOOFOAM & operator=(IOOFOAM other);

   void            execute();

   bool            getOverwrite();
   void            setOverwrite(bool flag);

   virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
   virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

protected:
   void swap(IOOFOAM & x) noexcept;
   virtual void setDefaults();
   virtual bool read();
   virtual bool write();
   virtual bool writePointsOnly();

private:
    // hide the method from interface.
    void    setFieldName(const std::string &fieldname){ BITPIT_UNUSED(fieldname);};

};

/*!
  \class IOOFOAMScalarField
  \ingroup ioofoam
  \brief IOOFOAM_Kernel specialization to import/export an OpenFOAM scalar field defined on a volume mesh
  and/or the related boundary fields.

  The class presumes mesh is already read by the User with a IOOFOAM specialized class. Thus it requires
  a valid link for both bulk and boundary geometry in input (M_GEOMOFOAM ports), coherent with the mesh defined the in
  OpenFoam directory where the field is defined.

  WARNING: WRITE mode for fields non available yet.

 * \n
 *
 Ports inherited from IOOFOAM_Kernel Class :
*
|                     Port Input   ||                                 |
|------------------|---------------------|----------------------------|
| <B>PortType</B>  | <B>variable/function</B>  |<B>DataType</B>       |
| M_GEOMOFOAM      | setGeometry         | (MC_SCALAR, MD_MIMMO_)     |
| M_GEOMOFOAM2     | setBoundaryGeometry | (MC_SCALAR, MD_MIMMO_)     |
| M_UMAPIDS        | setFacesMap         | (MC_UMAP, MD_LONG) |


|               Port Output    ||                                           |
|-------------------|---------------------------|---------------------------|
| <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B>            |
| M_GEOMOFOAM       | getGeometry               | (MC_SCALAR, MD_MIMMO_)    |
| M_GEOMOFOAM2      | getBoundaryGeometry       | (MC_SCALAR, MD_MIMMO_)    |
| M_UMAPIDS         | getFacesMap               | (MC_UMAP, MD_LONG)|

 Proper of the class:

|                     Port Input   ||                                     |
|------------------|---------------------|----------------------|
| <B>PortType</B>  | <B>variable/function</B>  |<B>DataType</B> |


|               Port Output    ||                                         |
|-------------------|--------------------|-----------------------|
| <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
| M_SCALARFIELD     | getField      		| (MC_MPVECTOR, MD_FLOAT)       |
| M_SCALARFIELD2    | getBoundaryField    | (MC_MPVECTOR, MD_FLOAT)       |

 =========================================================
 \n



 The xml available parameters, sections and subsections are the following :

  - <B>ClassName</B>: name of the class as <tt>mimmo.IOOFOAMScalarField</tt>;
  - <B>Priority</B>: uint marking priority in multi-chain execution;
  - <B>IOMode</B>: activate mode of the class: READ, WRITE;
  - <B>Dir</B>: path to the current OpenFOAM mesh for reading/writing purposes;
  - <B>FieldName</B>: name of the OpenFOAM field for reading/writing purposes;

  Geometries must be passed by ports.
 */
class IOOFOAMScalarField: public IOOFOAM_Kernel{

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

    void            execute();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

protected:
    void swap(IOOFOAMScalarField & x) noexcept;
    virtual bool    write();
    virtual bool    read();

};


/*!
  \class IOOFOAMVectorField
  \ingroup ioofoam
  \brief IOOFOAM_Kernel specialization to import/export an OpenFOAM vector field defined on a volume mesh
  and/or the related boundary fields.

  The class presumes mesh is already read by the User with a IOOFOAM specialized class. Thus it requires
  a valid link for bulk or boundary geometry in input (M_GEOMOFOAM ports), coherent with the mesh defined the in
  OpenFoam directory where the field is defined.

  WARNING: WRITE mode for fields non available yet.

 * \n
 *
 Ports inherited from IOOFOAM_Kernel Class :
*
|                     Port Input   ||                                 |
|------------------|---------------------|----------------------------|
| <B>PortType</B>  | <B>variable/function</B>  |<B>DataType</B>       |
| M_GEOMOFOAM      | setGeometry         | (MC_SCALAR, MD_MIMMO_)     |
| M_GEOMOFOAM2     | setBoundaryGeometry | (MC_SCALAR, MD_MIMMO_)     |
| M_UMAPIDS        | setFacesMap         | (MC_UMAP, MD_LONG) |


|               Port Output    ||                                           |
|-------------------|---------------------------|---------------------------|
| <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B>            |
| M_GEOMOFOAM       | getGeometry               | (MC_SCALAR, MD_MIMMO_)    |
| M_GEOMOFOAM2      | getBoundaryGeometry       | (MC_SCALAR, MD_MIMMO_)    |
| M_UMAPIDS         | getFacesMap               | (MC_UMAP, MD_LONG)|

 Proper of the class:

|                     Port Input   ||                                     |
|------------------|---------------------|----------------------|
| <B>PortType</B>  | <B>variable/function</B>  |<B>DataType</B> |


|               Port Output    ||                                         |
|-------------------|--------------------|-----------------------|
| <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
| M_VECTORFIELD     | getField      		| (MC_MPVECTOR, MD_FLOAT)       |
| M_VECTORFIELD2    | getBoundaryField    | (MC_MPVECTOR, MD_FLOAT)       |

 =========================================================
 \n



 The xml available parameters, sections and subsections are the following :

  - <B>ClassName</B>: name of the class as <tt>mimmo.IOOFOAMVectorField</tt>;
  - <B>Priority</B>: uint marking priority in multi-chain execution;
  - <B>IOMode</B>: activate mode of the class: READ, WRITE;
  - <B>Dir</B>: path to the current OpenFOAM mesh for reading/writing purposes;
  - <B>FieldName</B>: name of the OpenFOAM field for reading/writing purposes;

  Geometries must be passed by ports.
 */
class IOOFOAMVectorField: public IOOFOAM_Kernel{
protected:
    dmpvecarr3E m_field;        		/**<Internal volume field. */
    dmpvecarr3E m_boundaryField;        /**<Surface boundary field. */

public:
    IOOFOAMVectorField(int type = IOOFMode::READ);
    IOOFOAMVectorField(const bitpit::Config::Section & rootXML);
    virtual ~IOOFOAMVectorField();

    IOOFOAMVectorField(const IOOFOAMVectorField & other);
    IOOFOAMVectorField & operator=(IOOFOAMVectorField other);

    void            buildPorts();

    dmpvecarr3E     getField();
    dmpvecarr3E     getBoundaryField();

    void            execute();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

protected:
    void swap(IOOFOAMVectorField & x) noexcept;
    virtual bool    write();
    virtual bool    read();

};

REGISTER_PORT(M_GEOMOFOAM, MC_SCALAR, MD_MIMMO_, __IOOFOAM_HPP__)
REGISTER_PORT(M_GEOMOFOAM2, MC_SCALAR, MD_MIMMO_, __IOOFOAM_HPP__)
REGISTER_PORT(M_UMAPIDS, MC_UMAP, MD_LONG, __IOOFOAM_HPP__)
REGISTER_PORT(M_SCALARFIELD, MC_MPVECTOR, MD_FLOAT, __IOOFOAM_HPP__)
REGISTER_PORT(M_SCALARFIELD2, MC_MPVECTOR, MD_FLOAT, __IOOFOAM_HPP__)
REGISTER_PORT(M_VECTORFIELD, MC_MPVECARR3, MD_FLOAT, __IOOFOAM_HPP__)
REGISTER_PORT(M_VECTORFIELD2, MC_MPVECARR3, MD_FLOAT, __IOOFOAM_HPP__)

REGISTER(BaseManipulation, IOOFOAM, "mimmo.IOOFOAM")
REGISTER(BaseManipulation, IOOFOAMScalarField, "mimmo.IOOFOAMScalarField")
REGISTER(BaseManipulation, IOOFOAMVectorField, "mimmo.IOOFOAMVectorField")

}

#endif /* __IOOFOAM_HPP__ */
