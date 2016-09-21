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
 \ *---------------------------------------------------------------------------*/
#ifndef __IOVIOLATIONFILE_HPP__
#define __IOVIOLATIONFILE_HPP__

#include "BaseManipulation.hpp"
#include <unordered_map>

namespace mimmo{

/*! 
 * Class storing violation data info of one or more mimmo::ControlDeformMaxDistance and mimmo::ControlDeformExtSurface classes
 * (see their documentation).
 * Write data received to a file with the following scheme:
 * 
 * $MAXIMUM_VIOLATION	<double value>  
 * 
 * $SLOT	<geometry source name> 	<local value>
 * $SLOT ...
 * $SLOT ...
 * 
 * File has always an extension .dat.
 * Each $SLOT raw reports the name of the object who send the reference geometry to ControlDeformXXX class and the value
 * of deformation violation detected by it. 
 * $MAXIMUM_VIOLATION reports the maximum violation value, amongst those scheduled in SLOT list.
 * 
 * PORTS AVAILABLE IN the current class 
 *	=========================================================
 * ~~~
 *	|-----------------------------------------------------------------------------------------------------|
 *	|                   Port Input                                                                        |
 *	|-------|----------------|---------------------------------------|------------------------------------|
 *	|PortID | PortType       | variable/function                     | DataTypes	                      |
 *	|-------|----------------|---------------------------------------|------------------------------------|
 *	| 82    | M_VIOLATION    | addViolationData	                     | (PAIR, PAIRMIMMO_OBJFLOAT_)		  |
 *  |-------|----------------|---------------------------------------|------------------------------------|
 *
 *
 *	|-----------------------------------------------------------------------|
 *	|             Port Output                     							|
 *	|-------|----------------|--------------------|-------------------------|
 *	|PortID | PortType       | variable/function  | DataTypes	         	|
 *	|-------|----------------|--------------------|-------------------------|
 *	|-------|----------------|--------------------|-------------------------|
 * ~~~
 * 
 * *** this option works only with an active geometry set, and the geometry must be a superficial tessellation;
 *	=========================================================
 */
	class IOViolationFile: public mimmo::BaseManipulation{
		
	private:
		std::vector<std::pair<std::string, double> > m_list;		/**< list of violation pair values */
		std::string m_dir; 								/**<folder output */
		std::string m_namefile; 						/**<filename */
	public:
		
		IOViolationFile();
		virtual ~IOViolationFile();
		IOViolationFile(const IOViolationFile & other);
		IOViolationFile & operator=(const IOViolationFile & other);
		
		void	buildPorts();

		//get methods
		double											getMaximumViolation();
		std::vector<std::pair<std::string, double> >	getViolationMap();
		std::string										getDir();
		std::string										getFileName();
		
		//set methods
		void  	setDir(std::string dir);
		void 	setFileName(std::string name);
		void	setViolationMap(std::vector<std::pair<std::string, double> > list);
		void	addViolationData(std::pair<BaseManipulation *, double> data);
		void	addViolationDataString(std::pair<std::string, double> data);
		//cleaners
		void clearViolationMap();
		void clear();
		//execute
		void		execute();
		
		//XML utilities from reading writing settings to file
		virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");
 };
	
}
#endif /* __IOVIOLATIONFILE_HPP__ */
