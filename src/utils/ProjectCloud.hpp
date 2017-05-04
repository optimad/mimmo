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
#ifndef __PROJECTCLOUD_HPP__
#define __PROJECTCLOUD_HPP__

#include "BaseManipulation.hpp"

namespace mimmo{

/*!
 *  \class ProjectCloud
 *  \brief ProjectCloud is the class that project a 3D cloud of points on a target geometry
 *
 *
 *  =========================================================
 * ~~~
 *  |------------------------------------------------------------|
 *  |                 Port Input                                 |
 *  |-------|----------|-------------------|---------------------|
 *  |PortID | PortType | variable/function | DataType            |
 *  |-------|----------|-------------------|---------------------|
 *  | 0     | M_COORDS | setCoords         | (VECARR3, FLOAT)    |
 *  | 99    | M_GEOM   | setGeometry       | (SCALAR, MIMMO_)    |
 *  |-------|----------|-------------------|---------------------|
 *
 *
 *  |------------------------------------------------------------|
 *  |            Port Output                                     |
 *  |-------|----------|-------------------|---------------------|
 *  |PortID | PortType | variable/function | DataType            |
 *  |-------|----------|-------------------|---------------------|
 *  | 0     | M_COORDS | getCloudResult    | (VECARR3, FLOAT)    |
 *  |-------|----------|-------------------|---------------------|
 * ~~~
 *	=========================================================
 *
 */

class ProjectCloud: public BaseManipulation{
protected:
	dvecarr3E			m_points;	/**<Coordinates of 3D points in the cloud.*/
	dvecarr3E			m_proj;	 /**<Projected points coordinates.*/

public:
	ProjectCloud();
	ProjectCloud(const bitpit::Config::Section & rootXML);
	~ProjectCloud();

	ProjectCloud(const ProjectCloud & other);
	ProjectCloud & operator=(const ProjectCloud & other);

	void	buildPorts();

	dvecarr3E	getCoords();
	dvecarr3E	getCloudResult();

	void	setCoords(dvecarr3E coords);
	void 	execute();
	
	void clear();
	
	virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name= "");
	virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name= "");
	
protected:
	virtual void plotOptionalResults();

};

REGISTER(BaseManipulation, ProjectCloud,"mimmo.ProjectCloud")

};

#endif /* __PROJECTCLOUD_HPP__ */
