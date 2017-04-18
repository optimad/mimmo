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
#ifndef __IOOFOAM_HPP__
#define __IOOFOAM_HPP__

#include "BaseManipulation.hpp"
#include <vtkPolyData.h>

namespace mimmo{

/*!
 *	\date			20/oct/2016
 *	\authors		Edoardo Lombardi
 *
 *	\brief IOOFOAM is the class to import/export and modify an OpenFOAM mesh in format of surface vtk Polydata for surface meshes and points cloud for volume mesh.
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
 * Even the volume points cloud (import of one "points" file allowed)
 * is insert in pointer to a MimmoObject (locally instantiated).
 *
 * vtk libraries are needed.
 *

 *	=========================================================
 * ~~~
 *	|-----------------------------------------------------------------------|
 *	|                     Port Input                                    	|
 *	|-------|------------------|---------------------|----------------------|
 *	|PortID | PortType         | variable/function   | DataTypes			|
 *	|-------|------------------|---------------------|----------------------|
 *  | 19    | M_SCALARFIELD    | setField            | (VECTOR, FLOAT)      |
 *  | 99	| M_GEOM 		   | setGeometry  		 | (SCALAR, MIMMO_)		|
 *  | 98    | M_GEOM2          | setSurfaceBoundary  | (SCALAR, MIMMO_)     |
 *	|-------|------------------|---------------------|----------------------|
 *
 *
 *	|-----------------------------------------------------------------------|
 *	|               Port Output                                      		|
 *	|-------|------------------|--------------------|-----------------------|
 *	|PortID | PortType         | variable/function  | DataTypes				|
 *	|-------|------------------|--------------------|-----------------------|
 *  | 19    | M_SCALARFIELD    | getField           | (VECTOR, FLOAT)       |
 *  | 99	| M_GEOM 		   | getGeometry  	    | (SCALAR, MIMMO_) 		|
 *  | 98    | M_GEOM2          | getSurfaceBoundary | (SCALAR, MIMMO_)      |
 *	|-------|------------------|--------------------|-----------------------|
 *
 */
class IOOFOAM: public BaseManipulation{
private:
	bool			                m_read;         /**<If true it reads the geometry from files during the execution.*/
    std::vector<std::string>        m_rdirS;        /**<Name of directories to read the surface geometries (without final "/").*/
    std::vector<std::string>        m_rfilenameS;   /**<Name of files to read the surface geometries (without ".vtk" extension).*/
    std::string                     m_rdirV;        /**<Name of directories to read the volume point cloud (without final "/").*/
    std::string                     m_rfilenameV;   /**<Name of files to read the volume point cloud (with extension) in openFOAM points format.*/

	bool			                m_write; 		/**<If true it writes the geometries on file during the execution. [WARNING : only point cloud file written for this version.]*/
	std::string                     m_wdirS;		/**<Name of directory to write the surface geometry (without final "/").*/
	std::string                     m_wfilenameS;	/**<Name of file to write the geometry.*/
	std::string                     m_wdirV;        /**<Name of directory to write the point cloud (without final "/").*/
	std::string                     m_wfilenameV;   /**<Name of file to write the point cloud (with extension) in openFOAM points format.*/

    short                           m_stopat;       /**<Index of first patch not read (OpenFOAM points = MAX_DHORT) or, during write, flag of patch not linked (0=surface, 1=volume, 2=both). */

    std::unique_ptr<MimmoObject>    m_volmesh;      /**<Original volume mesh, instantiated in reading */
    std::unique_ptr<MimmoObject>    m_surfmesh;     /**<Original boundary mesh, instantiated in reading */
    MimmoObject*                    m_surfmesh_ext; /**<Pointed external boundary mesh */

    dvector1D                       m_field;        /**<Scalar field related to polydata mesh (pint located values).*/
    bool                            m_normalize;    /**<If true the field is normalized with the maximum absolute value after reading.*/
    double                          m_maxf;         /**<Max value of the read scalar field.*/
    double                          m_scaling;      /**<Value used to scale the scalar field (default m_scaling = 1). */


public:
	IOOFOAM();
	~IOOFOAM();

	IOOFOAM(const IOOFOAM & other);
	IOOFOAM & operator=(const IOOFOAM & other);

	void			buildPorts();

	void			setRead(bool read);
    void            setVTKReadDir(std::vector<std::string> dir);
    void            setVTKReadFilename(std::vector<std::string> filename);
    void            addVTKReadDir(std::string dir);
    void            addVTKReadFilename(std::string filename);
    void            setPointsReadDir(std::string dir);
    void            setPointsReadFilename(std::string filename);
	void			setWrite(bool write);
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

    void            readOFP(std::string& inputDir, std::string& surfaceName, dvecarr3E& points);
    void            writeOFP(std::string& outputDir, std::string& surfaceName, bitpit::PiercedVector<bitpit::Vertex> &vertices);
    bool            readVTK(std::string& inputDir, std::string& surfaceName, short PID, MimmoObject* patchBnd);

	bool			write();
	bool			read();

	void 			execute();

};

}

#endif /* __IOOFOAM_HPP__ */
