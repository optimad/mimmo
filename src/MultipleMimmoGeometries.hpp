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
#ifndef __MULTIPLEMIMMOGEOMETRY_HPP__
#define __MULTIPLEMIMMOGEOMETRY_HPP__

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
 *	MultipleMimmoGeometry is the object to manage the import/export/substantial modifications 
 *  of one or more geometry files of the same topology in a unique, standard MimmoObject structure as 
 *  product of its execution. It works as a  MimmoGeometry class, and retains almost the same features.
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
 *	7)	Ascii OpenFoam point cloud.		OFP 	= 6
 *
 * Outside this list of options, the class cannot hold any other type of formats for now.
 * The smart enum can be recalled in every moment in your code, just using mimmo::FileType.
 * and including MultipleMimmoGeometries header. 
 *	=========================================================
 * ~~~
 *	|--------------------------------------------------------------|
 *	|                 Port Input                                   |
 *	|-------|----------|-------------------|-----------------------|
 *	|PortID | PortType | variable/function | DataType		       |
 *	|-------|----------|-------------------|-----------------------|
 *	| 99    | M_GEOM   | m_geometry        | (SCALAR, MIMMO_)	   |
 *	|-------|----------|-------------------|-----------------------|
 *
 *
 *	|---------------------------------------|-----------------------|
 *	|            Port Output                |                       |
 *	|-------|-----------|-------------------|-----------------------|
 *	|PortID | PortType  | variable/function | DataType		      	|
 *	|-------|-----------|-------------------|-----------------------|
 *	| 99    | M_GEOM    | getGeometry       | (SCALAR, MIMMO_)		|
 *	|-------|-----------|-------------------|-----------------------|
 * ~~~
 *	=========================================================
 *
 */
class MultipleMimmoGeometries: public BaseManipulation{
	
private:
	int  		m_rtype;		/**<Extension of file to read the geometry.*/
	bool		m_read; 		/**<If true it reads the geometry from file during the execution.*/
	std::string	m_rdir;			/**<Name of directory to read the geometry (without final "/").*/
	std::string	m_rfilename;	/**<Name of file to read the geometry (without extension).*/

	int 		m_wtype;		/**<Extension of file to write the geometry.*/
	bool		m_write; 		/**<If true it writes the geometry on file during the execution.*/
	std::string	m_wdir;			/**<Name of directory to write the geometry (without final "/").*/
	std::string	m_wfilename;	/**<Name of file to write the geometry.*/

	bool		m_isInternal;				/**< flag for internal instantiated MimmoObject */
	std::unique_ptr<MimmoObject> m_intgeo;	/**< pointer to internal allocated geometry, if any */
	bool		m_codex;					/**< Set codex format for writing true binary, false ascii */
	WFORMAT		m_wformat;					/**<Format for .nas import/export. (Short/Long).*/

	bool		m_buildBvTree;				/**<If true the simplex ordered BvTree of the geometry is built in execution, whenever geometry support simplicies. */
	bool		m_buildKdTree;				/**<If true the vertex ordered KdTree of the geometry is built in execution*/

public:
	MimmoGeometry();
	virtual ~MimmoGeometry();

	MimmoGeometry(const MimmoGeometry & other);
	MimmoGeometry & operator=(const MimmoGeometry & other);

	void buildPorts();

	const MimmoGeometry *	getCopy();
	MimmoObject * 	getGeometry();
	const MimmoObject * 	getGeometry()const ;
	
	void		setReadDir(std::string dir);
	void		setReadFileType(FileType type);
	void		setReadFileType(int type);
	void		setWriteFileType(FileType type);
	void		setWriteFileType(int type);
	void		setRead(bool read);
	void		setWriteDir(std::string dir);
	void		setReadFilename(std::string filename);
	void		setWrite(bool write);
	void		setWriteFilename(std::string filename);
	void		setCodex(bool binary = true);
	
	void		setHARDCopy( const MimmoGeometry * other);		
	void		setSOFTCopy( const MimmoGeometry * other);
	
	void		setGeometry( MimmoObject * external);
	void		setGeometry(int type=1);

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
	bool		write();
	bool		read();

	void 		execute();

	void 		readOFP(std::string& inputDir, std::string& surfaceName, dvecarr3E& points); 
	void 		writeOFP(std::string& outputDir, std::string& surfaceName, dvecarr3E& points); 
	
	virtual void absorbSectionXML(bitpit::Config::Section & slotXML, std::string name="");
	virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");
	
private:
	void 	setDefaults();
	
};


class NastranInterface{
	static const char nl = '\n';
public:

	WFORMAT	wformat;

	void setWFormat(WFORMAT);
	void writeKeyword(std::string key, std::ofstream& os);
	void writeCoord(darray3E & p, int& pointI, std::ofstream& os);
	void writeFace(std::string faceType, ivector1D& facePts, int& nFace, std::ofstream& os, int PID);
	void writeGeometry(dvecarr3E& points, ivector2D& faces, std::ofstream& os, shivector1D* PIDS = NULL);
	void writeFooter(std::ofstream& os, std::unordered_set<short>* PIDSSET = NULL);
	void write(std::string& outputDir, std::string& surfaceName, dvecarr3E& points, ivector2D& faces, shivector1D* PIDS = NULL, std::unordered_set<short>* PIDSSET = NULL);
	void read(std::string& inputDir, std::string& surfaceName, dvecarr3E& points, ivector2D& faces, shivector1D& PIDS);

	std::string trim(std::string in);
	std::string convertVertex(std::string in);

	template<class Type>
	void writeValue (Type& value, std::ofstream& os){

		int offset = 5;
		if (wformat == Long) offset = 13;

		std::stringstream manip;
		manip << value;
		std::string manips = manip.str();
		std::string mantissa, expon;
		int pos = manips.find("E");
		if (pos < manips.size()){
			mantissa = manips.substr(0,std::min(offset,pos));
			expon = manips.substr(pos+1,manips.size());
			manips = mantissa + expon;
		}
		pos = manips.find("e");
		if (pos < manips.size()){
			mantissa = manips.substr(0,std::min(offset,pos));
			expon = manips.substr(pos+1,manips.size());
			manips = mantissa + expon;
		}
		manips = manips.substr(0, offset+3);

		switch (wformat)
		{
		case Short:
		{
			os << std::setw(8) << manips;
			break;
		}
		case Long:
		{
			os << std::setw(16) << manips;
			break;
		}
		}
	}
};


}

#endif /* __MIMMOGEOMETRY_HPP__ */
