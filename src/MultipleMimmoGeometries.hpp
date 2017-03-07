/*---------------------------------------------------------------------------*\
 *
 *  MiMMO
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of MiMMO.
 *
 *  MiMMO is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  MiMMO is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with MiMMO. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/
#ifndef __MULTIPLEMIMMOGEOMETRIES_HPP__
#define __MULTIPLEMIMMOGEOMETRIES_HPP__

#include "MimmoGeometry.hpp"

namespace mimmo{

/*!
 *	\date			30/nov/2016
 *	\authors		Rocco Arpa
 *	\authors		
 *
 *	\brief MultipleMimmoGeometries is an executable block class capable of  
 *         I/O hanfling of multiple geometries files of the same topology 
 *
 *	MultipleMultipleMimmoGeometries is the object to manage the reading/writing
 *  of one or more geometry files of the same topology to/from a list of standard MimmoObjects. 
 *  It works as a MimmoGeometry class, and retains almost the same features.
 *	The meshes to import/export have to be meshes with constant type elements.
 *	The valid format are: binary .stl, ascii .vtu (triangle/quadrilateral elements) and
 *	ascii .nas (triangle elements) for surface mesh; ascii .vtu (tetra/hexa elements)
 *	for volume mesh. ofp or pcvtu point cloud formats for generic point cloud structures.
 *  WARNING The current class can be used in Reading or Writing modes once, not both at the same time 
 * 
 * 
 * It uses smart enums FileType list of available geometry formats, which are:
 * 
 *	1)	Ascii/Binary triangulation stl.	STL 	= 0,
 *  2)	Surface triangulation vtu.		STVTU 	= 1,
 *	3)	Surface quadrilateral vtu.		SQVTU 	= 2,
 *	4)	Volume tetrahedral VTU.			VTVTU 	= 3,
 *	5)	Volume hexahedral VTU.			VHVTU 	= 4,
 *	6)	Nastran triangulation nas.		NAS 	= 5,
 *	7)	Ascii OpenFoam point cloud.		OFP 	= 6,
 *  8)	VTU point cloud.				PCVTU 	= 7,
 *
 * Outside this list of options, the class cannot hold any other type of formats for now.
 * The smart enum can be recalled in every moment in your code, just using mimmo::FileType.
 * and including MultipleMimmoGeometries header. 
 *	=========================================================
 * ~~~
 *	|-------------------------------------------------------------------------------|
 *	|                 Port Input                           					        |
 *	|-------|----------|------------------------------------|-----------------------|
 *	|PortID | PortType | variable/function 					| DataType		        |
 *	|-------|----------|------------------------------------|-----------------------|
 *	| 100   | M_VECGEOM| setGeometry      					| (VECTOR, MIMMO_)	    |
 *  | 101   | M_MAPGEOM| setObjMap  	      				| (UN_MAP, STRINGPAIRINTMIMMO_)|
 *  | 102   | M_FINFO  | setReadListFDI 			        | (VECTOR, FILEINFODATA)|
 *  | 103   | M_FINFO2 | setWriteListFDI     			    | (VECTOR, FILEINFODATA)| 
 *	|-------|----------|------------------------------------|-----------------------|
 *
 *
 *	|--------------------------------------------------------|-----------------------|
 *	|            Port Output                				 |                       |
 *	|-------|-----------|------------------------------------|-----------------------|
 *	|PortID | PortType  | variable/function 				 | DataType		      	 | 
 *	|-------|-----------|------------------------------------|-----------------------|
 *	| 100   | M_VECGEOM | getGeometry      					 | (VECTOR, MIMMO_)	     |
 *  | 101   | M_MAPGEOM | getObjMap		       				 | (UN_MAP, STRINGPAIRINTMIMMO_)|
 *  | 102   | M_FINFO   | getReadListFDI			         | (VECTOR, FILEINFODATA)|
 *  | 103   | M_FINFO2  | getWriteListFDI			         | (VECTOR, FILEINFODATA)| 
 *	|-------|-----------|------------------------------------|-----------------------|
 * ~~~
 *	=========================================================
 *
 */
class MultipleMimmoGeometries: public BaseManipulation{
	
private:
	int 						m_topo;			/**<Mark topology of your multi-file I/O geometry 1-surface, 2-volume, 3-pointcloud*/
	std::set<int>				m_ftype_allow;  /**< list of file type currnetly allowed for the class, according to its topology */
	bool						m_read;			/**<If true it reads the geometry from file during the execution.*/
	std::vector<FileDataInfo>	m_rinfo;		/**<List of filenames data given to the reader */
	
	bool						m_write;		/**<If true it writes the geometry from file during the execution.*/
	std::vector<FileDataInfo>	m_winfo;		/**<List of filenames data where the geometry must be written */
	
	bool		m_isInternal;								/**< flag for internal instantiated MimmoObjects */
	std::vector< std::unique_ptr<MimmoObject>> m_intgeo;	/**< pointers to internal allocated geometries, if any */
	std::vector<MimmoObject*>	m_extgeo;					/**< pointers to external allocated geometries, if any */

	bool		m_codex;					/**< Set codex format for writing true binary, false ascii */
	WFORMAT		m_wformat;					/**<Format for .nas import/export. (Short/Long).*/

	bool		m_buildBvTree;				/**<If true the simplex ordered BvTree of every geometries is built in execution, whenever geometry support simplicies. */
	bool		m_buildKdTree;				/**<If true the vertex ordered KdTree of every geometries is built in execution*/

	
public:
	MultipleMimmoGeometries(int topo);
	MultipleMimmoGeometries(int topo, bool IOMode);
	MultipleMimmoGeometries(const bitpit::Config::Section & rootXML);
	
	virtual ~MultipleMimmoGeometries();

	MultipleMimmoGeometries(const MultipleMimmoGeometries & other);
	MultipleMimmoGeometries & operator=(const MultipleMimmoGeometries & other);

	void buildPorts();

	std::vector<int>		getFileTypeAllowed();
	
	const MultipleMimmoGeometries *				getCopy();
	std::vector< MimmoObject *> 				getGeometry();
	const std::vector< MimmoObject *> 			getGeometry()const ;
	
	std::vector<FileDataInfo> 		getReadListFDI();
	std::vector<FileDataInfo> 		getWriteListFDI();

	std::unordered_map<std::string, std::pair<int, MimmoObject*> >	getObjMAP();
	
	void		setAddReadFile(std::string dir, std::string name, FileType ftype);
	void		setAddWriteFile(std::string dir, std::string name, FileType ftype);
	
	void 		setReadListFDI(std::vector<FileDataInfo> data);
	void 		setWriteListFDI(std::vector<FileDataInfo> data);
		
	void 		setObjMAP(std::unordered_map<std::string,std::pair<int, MimmoObject*> > map);
	
	void		setRead(bool read);
	void		setWrite(bool write);
	void		setCodex(bool binary = true);
	
	void		setHARDCopy( const MultipleMimmoGeometries * other);		
	void		setSOFTCopy( const MultipleMimmoGeometries * other);
	
	void		setGeometry( std::vector<MimmoObject *> external);
	
	
	void		setFormatNAS(WFORMAT wform);
	
	void		setBuildBvTree(bool build);
	void		setBuildKdTree(bool build);
	
	
	bool 		isEmpty();
	bool		isInternal();
	bool 		whichMode();
	
	void 		clear();
	void 		clearReadPaths();
	void 		clearWritePaths();
	bool		write();
	bool		read();

	void 		execute();
	
	virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
	virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");
	
private:
	void 	setDefaults();
	void 	initializeClass(int topo, bool IOMode);
	void	setGeometry();
};

REGISTER(BaseManipulation, MultipleMimmoGeometries, "MiMMO.MultipleGeometries")

};

#endif /* __MULTIPLEMIMMOGEOMETRIES_HPP__ */
