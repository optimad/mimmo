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
#include "Mask.hpp"

namespace mimmo{

/*!Default constructor of Mask
*/
Mask::Mask(){
	m_name = "mimmo.Mask";
	m_thres.fill({{0,0}});
	m_inside = {{true, true, true}};
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
Mask::Mask(const bitpit::Config::Section & rootXML){
	
	m_name = "mimmo.Mask";
	m_thres.fill({{0,0}});
	m_inside = {{true, true, true}};
	
	std::string fallback_name = "ClassNONE";	
	std::string input = rootXML.get("ClassName", fallback_name);
	input = bitpit::utils::trim(input);
	if(input == "mimmo.Mask"){
		absorbSectionXML(rootXML);
	}else{	
		std::cout<<"Warning in custom xml mimmo::Mask constructor. No valid xml data found"<<std::endl;
	};
}

/*!Default destructor of Mask
 */
Mask::~Mask(){};

/*!Copy constructor of Mask.
 */
Mask::Mask(const Mask & other):BaseManipulation(other){};

/*!Assignement operator of Mask.
 */
Mask & Mask::operator=(const Mask & other){
	*(static_cast<BaseManipulation*> (this)) = *(static_cast<const BaseManipulation*> (&other));
	return(*this);
};

/*! It builds the input/output ports of the object
 */
void
Mask::buildPorts(){
	bool built = true;
	built = (built && createPortIn<dvecarr3E, Mask>(&m_coords, PortType::M_COORDS, mimmo::pin::containerTAG::VECARR3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<dvecarr3E, Mask>(&m_displ, PortType::M_DISPLS, mimmo::pin::containerTAG::VECARR3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<dmatrix32E, Mask>(&m_thres, PortType::M_RANGE, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<std::array<bool,3>, Mask>(&m_inside, PortType::M_BOOLS3, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::BOOL));
	built = (built && createPortOut<dvecarr3E, Mask>(this, &mimmo::Mask::getCoords, PortType::M_COORDS, mimmo::pin::containerTAG::VECARR3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortOut<dvecarr3E, Mask>(this, &mimmo::Mask::getDisplacements, PortType::M_DISPLS, mimmo::pin::containerTAG::VECARR3, mimmo::pin::dataTAG::FLOAT));
	m_arePortsBuilt = built;
};

/*!It gets the coordinates of points stored in the object.
 * \return Coordinates of points stored in the object.
 */
dvecarr3E
Mask::getCoords(){
	return(m_coords);
};

/*!It gets the displacements of points stored in the object.
 * \return Displacements of points stored in the object.
 */
dvecarr3E
Mask::getDisplacements(){
	return(m_displ);
};

/*!It sets the coordinates of points used by the masking.
 * \param[in] coords Coordinates of points used by the masking.
 */
void
Mask::setCoords(dvecarr3E coords){
	m_coords = coords;
};

/*!It sets the limits of coordinates to apply the masking.
 * \param[in] thres Limits of coordinates to apply the masking.
 */
void
Mask::setThresholds(dmatrix32E thres){
	m_thres = thres;
};

/*!It sets the lower limits of the three coordinates to apply the masking.
 * \param[in] thres Minimum limits of coordinates to apply the masking.
 */
void
Mask::setMinThresholds(darray3E thres){
	for (int dir = 0; dir < 3; ++dir) m_thres[dir][0] = thres[dir];
};

/*!It sets the greater limits of the three coordinates to apply the masking.
 * \param[in] thres Maximum limits of coordinates to apply the masking.
 */
void
Mask::setMaxThresholds(darray3E thres){
	for (int dir = 0; dir < 3; ++dir) m_thres[dir][1] = thres[dir];
};

/*!It sets the limits of one coordinate to apply the masking.
 * \param[in] thres Limits of coordinate to apply the masking.
 * \param[in] dir Index of component.
 */
void
Mask::setThresholds(darray2E thres, int dir){
	m_thres[dir] = thres;
};

/*!It sets the limits of x-coordinate to apply the masking.
 * \param[in] thres Limits of x-coordinate to apply the masking.
 */
void
Mask::setThresholdx(darray2E thres){
	m_thres[0] = thres;
};

/*!It sets the limits of y-coordinate to apply the masking.
 * \param[in] thres Limits of y-coordinate to apply the masking.
 */
void
Mask::setThresholdy(darray2E thres){
	m_thres[1] = thres;
};

/*!It sets the limits of z-coordinate to apply the masking.
 * \param[in] thres Limits of z-coordinate to apply the masking.
 */
void
Mask::setThresholdz(darray2E thres){
	m_thres[2] = thres;
};

/*!It sets the condition to apply the masking
 * (true/false to set to zero the displacements inside/outside the thresholds).
 * \param[in] inside Condition to apply the mask for all the components.
 */
void
Mask::setInside(bool inside){
	for (int i=0; i<3; i++){
		m_inside[i] = inside;
	}
};

/*!It sets the condition to apply the masking
 * (true/false to set to zero the displacements inside/outside the thresholds).
 * \param[in] inside Condition to apply the mask for the three components.
 */
void
Mask::setInside(std::array<bool,3> inside){
	for (int i=0; i<3; i++){
		m_inside[i] = inside[i];
	}
};

/*!It sets the condition to apply the mask
 * (true/false to set to zero the displacements inside/outside the thresholds).
 * \param[in] i Index of component.
 * \param[in] inside Condition to apply the mask for i-th component.
 */
void
Mask::setInside(int i, bool inside){
	if (i >= 0 && i < 3) m_inside[i] = inside;
};

/*!Execution command. It modifies the initial values of m_displ with the masking conditions.
 * After exec() the modified values are stored in the m_displ member.
 */
void
Mask::execute(){
	int	ndispl = m_displ.size();
	ndispl = std::min(ndispl, int(m_coords.size()));
	for (int i=0; i<ndispl; i++){
		if (m_coords[i][0]>m_thres[0][0] && m_coords[i][0]<m_thres[0][1] &&
				m_coords[i][1]>m_thres[1][0] && m_coords[i][1]<m_thres[1][1] &&
				m_coords[i][2]>m_thres[2][0] && m_coords[i][2]<m_thres[2][1]){
			for (int j=0; j<3; j++){
				m_displ[i][j] = (1-m_inside[j])*m_displ[i][j];
			}
		}
		else{
			for (int j=0; j<3; j++){
				m_displ[i][j] = (m_inside[j])*m_displ[i][j];
			}
		}
	}
	return;
};

/*!
 * Get settings of the class from bitpit::Config::Section slot. Reimplemented from
 * BaseManipulation::absorbSectionXML.The class read essential parameters to mask a point cloud 
 * with a vector field attached. Point coordinates and fields are meant to be passed through ports.
 * 
 * --> Absorbing data:
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>LowerThreshold</B>: array of 3 float elements containing the coordinates of the lower threshold point. 
 * - <B>UpperThreshold</B>: array of 3 float elements containing the coordinates of the upper threshold point. 
 * - <B>Inside</B>: array of 3 booleans (0/1), each for space coordinate, to perform mask inside the threshold(0) or not (0) 
 * 
 * 
 * \param[in] slotXML bitpit::Config::Section which reads from
 * \param[in] name   name associated to the slot
 */
void Mask::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){
	
	BITPIT_UNUSED(name);
	
	std::string input; 

	if(slotXML.hasOption("Priority")){
		input = slotXML.get("Priority");
		int value =0;
		if(!input.empty()){
			std::stringstream ss(bitpit::utils::trim(input));
			ss>>value;
		}
		setPriority(value);
	};
	
	if(slotXML.hasOption("LowerThreshold")){
		std::string input = slotXML.get("LowerThreshold");
		input = bitpit::utils::trim(input);
		darray3E temp = {{-1.E18,-1.E18,-1.E18}};
		if(!input.empty()){
			std::stringstream ss(input);
			for(auto &val : temp) ss>>val;
		}
		setMinThresholds(temp);
	}; 
	
	if(slotXML.hasOption("UpperThreshold")){
		std::string input = slotXML.get("UpperThreshold");
		input = bitpit::utils::trim(input);
		darray3E temp = {{1.E18,1.E18,1.E18}};
		if(!input.empty()){
			std::stringstream ss(input);
			for(auto &val : temp) ss>>val;
		}
		setMaxThresholds(temp);
	}; 
	
	
	if(slotXML.hasOption("Inside")){
		std::string input = slotXML.get("Dimension");
		input = bitpit::utils::trim(input);
		std::array<bool,3> temp = {{false,false,false}};
		if(!input.empty()){
			std::stringstream ss(input);
			for(auto &val : temp) ss>>val;
		}
		setInside(temp);
	};
}

/*!
 * Write settings of the class from bitpit::Config::Section slot. Reimplemented from
 * BaseManipulation::absorbSectionXML.The class read essential parameters to build lattice.
 * 
 * --> Flushing data// how to write it on XML:
 * - <B>ClassName</B>: name of the class as "mimmo.Mask"
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>LowerThreshold</B>: array of 3 float elements containing the coordinates of the lower threshold point. 
 * - <B>UpperThreshold</B>: array of 3 float elements containing the coordinates of the upper threshold point. 
 * - <B>Inside</B>: array of 3 booleans (0/1), each for space coordinate, to perform mask inside the threshold(0) or not (0) 
 * 
 * \param[in] slotXML bitpit::Config::Section which writes to
 * \param[in] name   name associated to the slot
 */
void Mask::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
	BITPIT_UNUSED(name);
	
	slotXML.set("ClassName", m_name);
	slotXML.set("Priority", std::to_string(getPriority()));
	
	
	{
		std::stringstream ss;
		ss<<std::scientific<<m_thres[0][0]<<'\t'<<m_thres[1][0]<<'\t'<<m_thres[2][0];
		slotXML.set("LowerThreshold", ss.str());
	}
	
	{
		std::stringstream ss;
		ss<<std::scientific<<m_thres[0][1]<<'\t'<<m_thres[1][1]<<'\t'<<m_thres[2][1];
		slotXML.set("UpperThreshold", ss.str());
	}
	
	{
		std::stringstream ss;
		ss<<std::scientific<<m_inside[0]<<'\t'<<m_inside[1]<<'\t'<<m_inside[2];
		slotXML.set("Inside", ss.str());
	}
	
};	

}

