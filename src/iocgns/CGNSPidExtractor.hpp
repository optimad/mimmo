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
#ifndef __CGNSPIDEXTRACTOR_HPP__
#define __CGNSPIDEXTRACTOR_HPP__

#include "BaseManipulation.hpp"
#include <set>

namespace mimmo{

/*!
 *	\date			28/apr/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief CGNSPidExtractor is the class to extract a pidded patch from a surface boundary mesh readed form cgns file.
 *
 * The object return the pidded patch as an independent surface mesh. It triangulates the patch if forced to

 *	=========================================================
 * ~~~
 *	|-----------------------------------------------------------------------|
 *	|                     Port Input                                    	|
 *	|-------|------------------|---------------------|----------------------|
 *	|PortID | PortType         | variable/function   | DataTypes			|
 *	|-------|------------------|---------------------|----------------------|
 *  | 99	| M_GEOM 		   | setGeometry  		 | (SCALAR, MIMMO_)		|
 *	|-------|------------------|---------------------|----------------------|
 *
 *
 *	|-------------------------------------------------------------------------|
 *	|               Port Output                           					  |
 *	|-------|------------------|---------------------|------------------------|
 *	|PortID | PortType         | variable/function   | DataTypes			  |
 *	|-------|------------------|---------------------|------------------------|
 *  | 99	| M_GEOM 		   | getPatch()  		 | (SCALAR, MIMMO_)		  |
 *	|-------|------------------|---------------------|------------------------|
 *
 */
class CGNSPidExtractor: public BaseManipulation{
private:
	bool					m_force; 		/**<force triangulation of the extracted patch if true.*/
	std::set<short>			m_targetpid;	/**<If true it writes the geometry on file during the execution.*/
	
	std::unique_ptr<MimmoObject> m_patch; /*!extracted patch */
	
public:
	CGNSPidExtractor();
	~CGNSPidExtractor();

	CGNSPidExtractor(const CGNSPidExtractor & other);
	CGNSPidExtractor & operator=(const CGNSPidExtractor & other);

	void			buildPorts();
	
	MimmoObject*			getPatch();
	std::set<short>			whatPIDActive();
	bool					isForcedToTriangulate();
	
	void			addPID(short val);
	void			setPID(std::vector<short int> vval);
	void 			setForcedToTriangulate(bool flag);
	void			setGeometry(MimmoObject*);
	
	void 			execute();

protected:
	
	void plotOptionalResults();
	
};

}

#endif /* __CGNSPIDEXTRACTOR_HPP__ */
