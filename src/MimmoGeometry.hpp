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
#ifndef __MIMMOGEOMETRY_HPP__
#define __MIMMOGEOMETRY_HPP__

#include "BaseManipulation.hpp"

namespace mimmo{

class NastranInterface;

enum FileType{/*!Ascii/Binary stl.*/ STL = 0, /*!Surface vtu.*/ SVTU = 1,
			/*!Volume VTU.*/ VVTU = 2, /*!Nastran .*/ NAS = 3,
			/*!Ascii OpenFoam point cloud.*/ OFP = 4};					/**Type of file to read
																	the geometry.*/
enum WFORMAT{Short, Long};

/*!
 *	\date			30/mar/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief MimmoGeometry is the class that wrap a geometry Mimmo Object .
 *
 *	The parameter of linked geometry are given by wrapper functions of this
 *	BaseManipulation derived class.
 *	MimmoGeometry is the object to manage the import/export of geometry file.
 *	The valid format are: binary .stl, ascii .vtu (triangle elements) and
 *	ascii .nas (triangle elements) for surface mesh; ascii .vtu (tetra elements) for volume mesh.
 *
 */
class MimmoGeometry: public BaseManipulation{
private:
	FileType	m_rtype;		/**<Extension of file to read the geometry.*/
	bool		m_read; 		/**<If true it reads the geometry from file during the execution.*/
	std::string	m_rdir;			/**<Name of directory to read the geometry (without final "/").*/
	std::string	m_rfilename;	/**<Name of file to read the geometry (without extension).*/

	FileType	m_wtype;		/**<Extension of file to write the geometry.*/
	bool		m_write; 		/**<If true it writes the geometry on file during the execution.*/
	std::string	m_wdir;			/**<Name of directory to write the geometry (without final "/").*/
	std::string	m_wfilename;	/**<Name of file to write the geometry.*/

	bool		m_local;		/**<Is the geometry locally instantiated?.*/

	WFORMAT		m_wformat;		/**<Format for .nas import/export. (Short/Long).*/
	ivector1D	m_pids;			/**<Cells' PID.*/

public:
	MimmoGeometry();
	~MimmoGeometry();

	MimmoGeometry(const MimmoGeometry & other);
	MimmoGeometry & operator=(const MimmoGeometry & other);

	int			getType();
	long		getNVertex();
	long		getNCells();
	dvecarr3E	getVertex();
	darray3E	getVertex(long i);
	ivector1D	getConnectivity(long i);
	ivector2D	getConnectivity();
	ivector1D	getPID();
	int			getPID(long i);

	bool		setVertex(dvecarr3E & vertex);
	bool		setVertex(darray3E & vertex);
	bool		modifyVertex(darray3E & vertex, long id);
	bool		setConnectivity(ivector2D * connectivity);
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

	void		activatePID();
	void		setPID(ivector1D pids);
	void		setPID(long i, int pid);
	void		setPIDforce(long i, int pid);

	void		setFormatNAS(WFORMAT wform);

	bool		write();
	bool		read();

	void 		execute();


	void 		readOFP(std::string& inputDir, std::string& surfaceName, dvecarr3E& points);
	void 		writeOFP(std::string& outputDir, std::string& surfaceName, dvecarr3E& points);

};


class NastranInterface{
	static const char nl = '\n';
public:

	WFORMAT	wformat;

	void setWFormat(WFORMAT);
	void writeKeyword(std::string key, std::ofstream& os);
	void writeCoord(darray3E & p, int& pointI, std::ofstream& os);
	void writeFace(std::string faceType, ivector1D& facePts, int& nFace, std::ofstream& os, int PID);
	void writeGeometry(dvecarr3E& points, ivector2D& faces, std::ofstream& os, ivector1D* PIDS = NULL);
	void writeFooter(std::ofstream& os);
	void write(std::string& outputDir, std::string& surfaceName, dvecarr3E& points, ivector2D& faces, ivector1D* PIDS = NULL);
	void read(std::string& inputDir, std::string& surfaceName, dvecarr3E& points, ivector2D& faces, ivector1D& PIDS);

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
