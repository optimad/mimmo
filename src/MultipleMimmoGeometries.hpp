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
 *         handling multiple geometries files of the same topology into a 
 *         unique MimmoObject container
 *
 *	MultipleMultipleMimmoGeometries is the object to manage the import/export/substantial modifications 
 *  of one or more geometry files of the same topology in a unique, standard MimmoObject structure as 
 *  product of its execution. It works as a  MultipleMimmoGeometries class, and retains almost the same features.
 *	The meshes to import/export have to be meshes with constant type elements.
 *	The valid format are: binary .stl, ascii .vtu (triangle/quadrilateral elements) and
 *	ascii .nas (triangle elements) for surface mesh; ascii .vtu (tetra/hexa elements)
 *	for volume mesh.
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
 *	| 99    | M_GEOM   | setGeometry      					| (SCALAR, MIMMO_)	    |
 *  | 101   | M_MAPGEOM| setDivisionMap        				| (UN_MAP, LONGSHORT)   |
 *  | 102   | M_FINFO  | setReadListAbsolutePathFiles       | (VECTOR, FILEINFODATA)|
 *  | 103   | M_FINFO2 | setWriteListAbsolutePathFiles      | (VECTOR, FILEINFODATA)| 
 *	|-------|----------|------------------------------------|-----------------------|
 *
 *
 *	|--------------------------------------------------------|-----------------------|
 *	|            Port Output                				 |                       |
 *	|-------|-----------|------------------------------------|-----------------------|
 *	|PortID | PortType  | variable/function 				 | DataType		      	 | 
 *	|-------|-----------|------------------------------------|-----------------------|
 *	| 99    | M_GEOM    | getGeometry      					 | (SCALAR, MIMMO_)	     |
 *  | 101   | M_MAPGEOM | getDivisionMap       				 | (UN_MAP, LONGSHORT)   |
 *  | 102   | M_FINFO   | getReadListAbsolutePathFiles       | (VECTOR, FILEINFODATA)|
 *  | 103   | M_FINFO2  | getWriteListAbsolutePathFiles      | (VECTOR, FILEINFODATA)| 
 *	|-------|-----------|------------------------------------|-----------------------|
 * ~~~
 *	=========================================================
 *
 */
class MultipleMimmoGeometries: public BaseManipulation{
	
private:
	int 						m_topo;			/**<Mark topology of your multi-file I/O geometry 1-surface, 2-volume, 3-pointcloud*/
	int 						m_tagtype; 		/**<Mark type of tag format admissible. See FileType enum */
	bool						m_read;			/**<If true it reads the geometry from file during the execution.*/
	std::vector<FileDataInfo>	m_rinfo;		/**<List of filenames data given to the reader */
	
	bool						m_write;		/**<If true it writes the geometry from file during the execution.*/
	std::vector<FileDataInfo>	m_winfo;		/**<List of filenames data where the geometry must be written */
	
	std::unordered_map<long, short>	m_mapMGeo;	/**<mapping cells ids of the MimmoObject with short ids linking to each sub-geometry filename il m_rinfo, m_winfo lists. 
													The map is compiled automatically in reading mode.It is required in writing mode to divide geometry
													in the various sub-files*/	
			
	bool		m_isInternal;				/**< flag for internal instantiated MimmoObject */
	std::unique_ptr<MimmoObject> m_intgeo;	/**< pointer to internal allocated geometry, if any */
	bool		m_codex;					/**< Set codex format for writing true binary, false ascii */
	WFORMAT		m_wformat;					/**<Format for .nas import/export. (Short/Long).*/

	bool		m_buildBvTree;				/**<If true the simplex ordered BvTree of the geometry is built in execution, whenever geometry support simplicies. */
	bool		m_buildKdTree;				/**<If true the vertex ordered KdTree of the geometry is built in execution*/

	
public:
	MultipleMimmoGeometries(int topo, int formattype);
	MultipleMimmoGeometries(int topo, FileType formattype);
	virtual ~MultipleMimmoGeometries();

	MultipleMimmoGeometries(const MultipleMimmoGeometries & other);
	MultipleMimmoGeometries & operator=(const MultipleMimmoGeometries & other);

	void buildPorts();

	int 	getFormatTypeAllowed();
	const MultipleMimmoGeometries *	getCopy();
	MimmoObject * 	getGeometry();
	const MimmoObject * 	getGeometry()const ;
	
	std::unordered_map<long, short> getDivisionMap();
	std::unordered_set<short>		getHowManySubDivisions();
	std::vector<FileDataInfo> 		getReadListAbsolutePathFiles();
	std::vector<FileDataInfo> 		getWriteListAbsolutePathFiles();
	
	std::unordered_set<short> &	 getPIDList();
	
	void		setAddReadAbsolutePathFile(std::string dir, std::string name);
	void		setAddWriteAbsolutePathFile(std::string dir, std::string name);
	
	void 		setReadListAbsolutePathFiles(std::vector<FileDataInfo> data);
	void 		setWriteListAbsolutePathFiles(std::vector<FileDataInfo> data);
	
	void		setRead(bool read);
	void		setWrite(bool write);
	void		setCodex(bool binary = true);
	
	void		setHARDCopy( const MultipleMimmoGeometries * other);		
	void		setSOFTCopy( const MultipleMimmoGeometries * other);
	
	void		setGeometry( MimmoObject * external);
	void		setGeometry();
	
	void 		setDivisionMap(std::unordered_map<long,short> mapMG);

	bitpit::PiercedVector<bitpit::Vertex> * 	getVertices();
	bitpit::PiercedVector<bitpit::Cell> * 		getCells();
	
	void		setVertices(bitpit::PiercedVector<bitpit::Vertex> * vertices);
	void		setCells(bitpit::PiercedVector<bitpit::Cell> * cells);
	
	void		setPID(shivector1D pids);
	void		setPID(std::unordered_map<long,short> pidsMap);
	
	void		setFormatNAS(WFORMAT wform);
	
	void		setBuildBvTree(bool build);
	void		setBuildKdTree(bool build);
	
	bool 		isEmpty();
	bool		isInternal();
	void 		clear();
	void 		clearReadPaths();
	void 		clearWritePaths();
	bool		write();
	bool		read();

	void 		execute();
	
	virtual void absorbSectionXML(bitpit::Config::Section & slotXML, std::string name="");
	virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");
	
private:
	void 	setDefaults();
	void 	initializeClass(int topo, int formattype);
	livector1D    cellExtractor(short id);
	void		  fillSubStructure(livector1D &cellList, MimmoObject* subG);
	bool checkCoherentFormat(MimmoObject *);
	bool checkReadingFiles(FileDataInfo & filedata);
	
	
};


}

#endif /* __MULTIPLEMIMMOGEOMETRIES_HPP__ */
