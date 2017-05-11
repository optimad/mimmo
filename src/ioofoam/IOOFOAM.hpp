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
#include <vtkPolyData.h>

namespace mimmo{

/*!
 * \class IOOFOAM
 * \ingroup ioofoam
 *    \brief IOOFOAM is the class to import/export and modify an OpenFOAM mesh in format of surface vtk Polydata for surface meshes and points cloud for volume mesh.
 *
 * The surfaces patch (boundary of volume mesh) have to be imported as .vtk files, while the volume mesh has to be
 * read as points file with OpenFOAM format.
 *
 * When imported a triangulation filter is applied to the each polymesh of surface
 * meshes (input of multiple patches allowed) and insert together
 * in a unique pointer to MimmoObject objects (locally instantiated).
 *
 * A scalar field related to each surface mesh is read if it is present.
 *
 * The volume points cloud (import of one "points" file allowed)
 * is insert in pointer to a MimmoObject (locally instantiated).
 *
 * Dependencies : vtk libraries (tested with vtk DataFile Version 4.0).
 *
 * \n
 *
 * Ports available in IOOFOAM Class :
 *
 *    =========================================================

   |                     Port Input   |||                                     |
   |-------|------------------|---------------------|----------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
   | 19    | M_SCALARFIELD    | setField            | (VECTOR, FLOAT)      |
   | 99    | M_GEOM           | setGeometry         | (SCALAR, MIMMO_)     |
   | 98    | M_GEOM2          | setSurfaceBoundary  | (SCALAR, MIMMO_)     |


   |               Port Output    |||                                         |
   |-------|------------------|--------------------|-----------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
   | 19    | M_SCALARFIELD    | getField           | (VECTOR, FLOAT)       |
   | 99    | M_GEOM           | getGeometry        | (SCALAR, MIMMO_)      |
   | 98    | M_GEOM2          | getSurfaceBoundary | (SCALAR, MIMMO_)      |

 * =========================================================
 * \n
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * - <B>ClassName</B>: name of the class as <tt>mimmo.IOOFOAM</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>ReadFlag</B>: activate reading mode boolean;
 * - <B>VTKReadDirs</B>: VTK reading directories path; there can be more than one with the following sub-structure : \n
 *              <tt> \<VTKReadDirs\> \n
 *                  \<dir0\> \n
 *                      \<dir\> path to the directory of VTK file \</dir\> \n
 *                  \</dir0\> \n
 *                  \<dir1\> \n
 *                      \<dir\> path to the directory of VTK file \</dir\> \n
 *                  \</dir1\> \n
 *                  ... \n
 *                  ... \n
 *              \</VTKReadDirs\> </tt> \n
 * - <B>VTKReadFilenames</B>: names of files VTK for reading there can be more than one with the following sub-structure : \n
 *           <tt>  \<VTKReadFilenames\> \n
 *                  \<file0\> \n
 *                      \<filename\> path to the directory of VTK file \</dir\> \n
 *                  \</file0\> \n
 *                  \<file1\> \n
 *                      \<filename\> path to the directory of VTK file \</dir\> \n
 *                  \</file1\> \n
 *                  ... \n
 *                  ... \n
 *              \</VTKReadFilenames\> </tt> \n
 * - <B>PointsReadDir</B>: points file reading directory path (there can be only one);
 * - <B>PointsReadFilename</B>: name of points file for reading (there can be only one);
 * - <B>WriteFlag</B>: activate writing mode boolean;
 * - <B>VTKWriteDir</B>: VTK file (one in output) writing directory path;
 * - <B>VTKWriteFilename</B>: VTK name of file for writing;
 * - <B>Scalar</B>: value to scale the scalar field eventually present on volume mesh (defined on points) [default 1.0];
 * - <B>Normalize</B>: bool to define if the scalar field has to be normalize after read [default false].
 *
 * In case of writing mode Geometries have to be mandatorily passed by port.
 *
 */
class IOOFOAM: public BaseManipulation{
private:
    bool                            m_read;         /**<If true it reads the geometry from files during the execution.*/
    std::vector<std::string>        m_rdirS;        /**<Name of directories to read the surface geometries (without final "/").*/
    std::vector<std::string>        m_rfilenameS;   /**<Name of files to read the surface geometries (without ".vtk" extension).*/
    std::string                     m_rdirV;        /**<Name of directories to read the volume point cloud (without final "/").*/
    std::string                     m_rfilenameV;   /**<Name of files to read the volume point cloud (with extension) in openFOAM points format.*/

    bool                            m_write;         /**<If true it writes the geometries on file during the execution. [WARNING : only point cloud file written for this version.]*/
    std::string                     m_wdirS;        /**<Name of directory to write the surface geometry (without final "/").*/
    std::string                     m_wfilenameS;    /**<Name of file to write the geometry.*/
    std::string                     m_wdirV;        /**<Name of directory to write the point cloud (without final "/").*/
    std::string                     m_wfilenameV;   /**<Name of file to write the point cloud (with extension) in openFOAM points format.*/

    short                           m_stopat;       /**<Index of first patch not read (OpenFOAM points = MAX_SHORT) or, during write, flag of patch not linked (0=surface, 1=volume, 2=both). */

    std::unique_ptr<MimmoObject>    m_volmesh;      /**<Original volume mesh, instantiated in reading */
    std::unique_ptr<MimmoObject>    m_surfmesh;     /**<Original boundary mesh, instantiated in reading */
    MimmoObject*                    m_surfmesh_ext; /**<Pointed external boundary mesh */

    dvector1D                       m_field;        /**<Scalar field related to polydata mesh (pint located values).*/
    bool                            m_normalize;    /**<If true the field is normalized with the maximum absolute value after reading.*/
    double                          m_maxf;         /**<Max value of the read scalar field.*/
    double                          m_scaling;      /**<Value used to scale the scalar field (default m_scaling = 1). */


public:
    IOOFOAM();
    IOOFOAM(const bitpit::Config::Section & rootXML);
    ~IOOFOAM();

    IOOFOAM(const IOOFOAM & other);
    IOOFOAM & operator=(const IOOFOAM & other);

    void            buildPorts();
    void            setDefaults();

    void            setRead(bool read);
    void            setVTKReadDir(std::vector<std::string> dir);
    void            setVTKReadFilename(std::vector<std::string> filename);
    void            addVTKReadDir(std::string dir);
    void            addVTKReadFilename(std::string filename);
    void            setPointsReadDir(std::string dir);
    void            setPointsReadFilename(std::string filename);
    void            setWrite(bool write);
    void            setVTKWriteDir(std::string dir);
    void            setVTKWriteFilename(std::string filename);
    void            setPointsWriteDir(std::string dir);
    void            setPointsWriteFilename(std::string filename);
    void            setGeometry(MimmoObject* geom);
    void            setSurfaceBoundary(MimmoObject* surf);
    void            setNormalize(bool normalize);
    void            setScaling(double scaling);
    void            setField(dvector1D field);

    MimmoObject*    getSurfaceBoundary();
    MimmoObject*    getGeometry();
    dvector1D       getField();

    void            readOFP(std::string& inputDir, std::string& pointsName, dvecarr3E& points);
    void            writeOFP(std::string& outputDir, std::string& pointsName, bitpit::PiercedVector<bitpit::Vertex> &vertices);
    bool            readVTK(std::string& inputDir, std::string& surfaceName, short PID, MimmoObject* patchBnd);

    bool            write();
    bool            read();

    void             execute();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

};

REGISTER(BaseManipulation, IOOFOAM,"mimmo.IOOFOAM")

}

#endif /* __IOOFOAM_HPP__ */
