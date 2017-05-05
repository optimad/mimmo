/*---------------------------------------------------------------------------*\
 *
 *  mimmo
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
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
#ifndef __IOVTKScalar_HPP__
#define __IOVTKScalar_HPP__

#include "BaseManipulation.hpp"
#include <vtkPolyData.h>

namespace mimmo{

/*!
 *	\date			28/apr/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief IOVTKScalar is the class to import/export and modify a surface vtk Polydata geometry.
 *
 * The object provides an interface to retrieve/modify a scalar field evaluated on
 * the points of the surface polymesh.
 *
 * When imported a triangulation filter is applied to the
 * polymesh and insert in the linked MimmoObject (if not linked a MimmoObject is locally
 * instantiated).
 *
 * When exported a Polydata linked geometry is used if linked, otherwise the
 * linked MimmoObject is used to write a Polydata vtk surface.
 *
 * Dependencies : vtk libraries (tested with vtk DataFile Version 4.0).
 *

 *  =========================================================
 * ~~~
 *  |-----------------------------------------------------------------------|
 *  |                     Port Input                                        |
 *  |-------|------------------|---------------------|----------------------|
 *  |PortID | PortType         | variable/function   | DataTypes            |
 *  |-------|------------------|---------------------|----------------------|
 *  | 19    | M_SCALARFIELD    | setField            | (VECTOR, FLOAT)      |
 *  | 30    | M_VALUED         | setScaling          | (SCALAR, FLOAT)      |
 *  | 99    | M_GEOM           | setGeometry         | (SCALAR, MIMMO_)     |
 *  | 1100  | M_POLYDATA_      | setPolyData         | (SCALAR, POLYDATA_)  |
 *  |-------|------------------|---------------------|----------------------|
 *
 *
 *  |-----------------------------------------------------------------------|
 *  |               Port Output                                             |
 *  |-------|------------------|-------------------|------------------------|
 *  |PortID | PortType         | variable/function | DataTypes              |
 *  |-------|------------------|-------------------|------------------------|
 *  | 19    | M_SCALARFIELD    | getField          | (VECTOR, FLOAT)        |
 *  | 99    | M_GEOM           | getGeometry       | (SCALAR, MIMMO_)       |
 *  | 1100  | M_POLYDATA_      | getPolyData       | (SCALAR, POLYDATA_)    |
 *	|-------|------------------|-------------------|------------------------|
 *
 */
class IOVTKScalar: public BaseManipulation{
private:
    bool            m_read;         /**<If true it reads the geometry from file during the execution.*/
    std::string     m_rdir;         /**<Name of directory to read the geometry (without final "/").*/
    std::string     m_rfilename;    /**<Name of file to read the geometry (without extension).*/

    bool            m_write;        /**<If true it writes the geometry on file during the execution.*/
    std::string     m_wdir;         /**<Name of directory to write the geometry (without final "/").*/
    std::string     m_wfilename;    /**<Name of file to write the geometry.*/

    bool            m_local;        /**<Is the geometry locally instantiated?.*/

    vtkPolyData*    m_polydata;     /**<VTK Polydata geometry member.*/

    dvector1D       m_field;        /**<Scalar field related to polydata mesh (pint located values).*/
    bool            m_normalize;    /**<If true the field is normalized with the maximum absolute value after reading.*/
    double          m_scaling;      /**<Value used to scale the scalar field (default m_scaling = 1). */

public:
    IOVTKScalar();
    IOVTKScalar(const bitpit::Config::Section & rootXML);
    ~IOVTKScalar();

    IOVTKScalar(const IOVTKScalar & other);
    IOVTKScalar & operator=(const IOVTKScalar & other);

    void            buildPorts();

    vtkPolyData*    getPolyData();
    double          getScaling();
    dvector1D       getField();

    void            setReadDir(std::string dir);
    void            setRead(bool read);
    void            setWriteDir(std::string dir);
    void            setReadFilename(std::string filename);
    void            setWrite(bool write);
    void            setWriteFilename(std::string filename);
    void            setPolyData(vtkPolyData* polydata);
    void            setNormalize(bool normalize);
    void            setScaling(double scaling);
    void            setField(dvector1D field);

    bool            write();
    bool            read();

    void            execute();
    
    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

};

REGISTER(BaseManipulation, IOVTKScalar, "mimmo.IOVTKScalar")
}

#endif /* __IOVTKScalar_HPP__ */
