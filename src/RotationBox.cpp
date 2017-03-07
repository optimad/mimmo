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
#include "RotationBox.hpp"

using namespace mimmo;


/*!
 * Default constructor of RotationBox
 */
RotationBox::RotationBox(darray3E origin, darray3E direction){
	m_origin = origin;
	m_direction = direction;
	m_name = "MiMMO.RotationBox";
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
RotationBox::RotationBox(const bitpit::Config::Section & rootXML){
	
	m_origin.fill(0.0);
	m_direction.fill(0.0);
	m_name = "MiMMO.RotationBox";
	
	std::string fallback_name = "ClassNONE";	
	std::string input = rootXML.get("ClassName", fallback_name);
	input = bitpit::utils::trim(input);
	if(input == "MiMMO.RotationBox"){
		absorbSectionXML(rootXML);
	}else{	
		std::cout<<"Warning in custom xml MiMMO::RotationBox constructor. No valid xml data found"<<std::endl;
	};
}

/*!Default destructor of RotationBox
 */
RotationBox::~RotationBox(){};

/*!Copy constructor of RotationBox.
 */
RotationBox::RotationBox(const RotationBox & other):BaseManipulation(other){
	m_origin = other.m_origin;
	m_direction = other.m_direction;
};

/*!Assignement operator of RotationBox.
 */
RotationBox & RotationBox::operator=(const RotationBox & other){
	*(static_cast<BaseManipulation*> (this)) = *(static_cast<const BaseManipulation*> (&other));
	m_origin = other.m_origin;
	m_direction = other.m_direction;
	return(*this);
};


/*! It builds the input/output ports of the object
 */
void RotationBox::buildPorts(){
	bool built = true;
	built = (built && createPortIn<darray3E, RotationBox>(&m_origin, PortType::M_POINT, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<darray3E, RotationBox>(&m_direction, PortType::M_AXIS, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<double, RotationBox>(&m_alpha, PortType::M_VALUED, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<darray3E, RotationBox>(&m_axes_origin, PortType::M_POINT2, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<dmatrix33E, RotationBox>(&m_axes, PortType::M_AXES, mimmo::pin::containerTAG::ARR3ARR3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortOut<darray3E, RotationBox>(this, &mimmo::RotationBox::getRotatedOrigin, PortType::M_POINT, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortOut<dmatrix33E, RotationBox>(this, &mimmo::RotationBox::getRotatedAxes, PortType::M_AXES, mimmo::pin::containerTAG::ARR3ARR3, mimmo::pin::dataTAG::FLOAT));
	m_arePortsBuilt = built;
};

/*!It sets the origin and direction of the rotation axis.
 * \param[in] origin Origin of rotation axis.
 * \param[in] direction Direction of rotation axis.
 */
void
RotationBox::setAxis(darray3E origin, darray3E direction){
	m_origin = origin;
	m_direction = direction;
}

/*!It sets the origin of the rotation axis.
 * \param[in] origin Origin of rotation axis.
 */
void
RotationBox::setOrigin(darray3E origin){
	m_origin = origin;
}

/*!It sets the direction of the rotation axis.
 * \param[in] direction Direction of rotation axis.
 */
void
RotationBox::setDirection(darray3E direction){
	m_direction = direction;
	double L = norm2(m_direction);
	for (int i=0; i<3; i++)
		m_direction[i] /= L;
}

/*!It sets the value of the rotation.
 * \param[in] alpha Value of rotation axis.
 */
void
RotationBox::setRotation(double alpha){
	m_alpha = alpha;
}

/*!It sets the reference system to be rotated.
 * \param[in] axes Original reference system.
 */
void
RotationBox::setAxes(dmatrix33E axes){
	m_axes = axes;
}

/*!It sets the origin of the reference system to be rotated.
 * \param[in] axes_origin Origin of reference system.
 */
void
RotationBox::setAxesOrigin(darray3E axes_origin){
	m_axes_origin = axes_origin;
}

/*!It gets the rotated reference system.
 * \return Rotated reference system.
 */
dmatrix33E
RotationBox::getRotatedAxes(){
	return(m_rotax);
}

/*!It gets the rotated origin of the reference system.
 * \return Rotated origin of reference system.
 */
darray3E
RotationBox::getRotatedOrigin(){
	return(m_rotax_origin);
}

/*!Execution command. It saves in "rot"-terms the modified axes and origin, by the
 * rotation conditions. This terms can be recovered and passed by a pin to a child object
 * by the related get-methods.
 */
void
RotationBox::execute(){

	//Rotation of origin
	m_rotax_origin = {{0,0,0}};
	m_axes_origin -= m_origin;
	//rodrigues formula
	m_rotax_origin = m_axes_origin * cos(m_alpha) +
			dotProduct(m_direction, m_axes_origin) * (1 - cos(m_alpha)) * m_direction +
			crossProduct(m_direction, m_axes_origin) * sin(m_alpha);

	m_rotax_origin += m_origin;
	m_axes_origin += m_origin;

	//rotation of axes
	m_rotax.fill(darray3E{{0,0,0}});
	//rodrigues formula
	for (int i=0; i<3; i++){
		m_rotax[i] = m_axes[i] * cos(m_alpha) +
				dotProduct(m_direction, m_axes[i]) * (1 - cos(m_alpha)) * m_direction +
				crossProduct(m_direction, m_axes[i]) * sin(m_alpha);
	}

	return;
};


/*!
 * Get settings of the class from bitpit::Config::Section slot. Reimplemented from
 * BaseManipulation::absorbSectionXML.The class read essential parameters to perform rotation of an axes ref system
 * 
 * --> Absorbing data:
 * 		Priority  : uint marking priority in multi-chain execution; 
 * 		Origin: rotation axis origin
 * 		Direction: axis direction coordinates
 * 		Rotation : rotation angle in radians. Positive on counterclockwise rotations around reference axis
 * 		RefSystem: current reference system to be rotated
 * 		OriginRS: origin of the target reference system to be rotated
 * 
 * \param[in] slotXML bitpit::Config::Section which reads from
 * \param[in] name   name associated to the slot
 */
void RotationBox::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){
	
	
	if(slotXML.hasOption("Priority")){
		std::string input = slotXML.get("Priority");
		int value =0;
		if(!input.empty()){
			std::stringstream ss(bitpit::utils::trim(input));
			ss>>value;
		}
		setPriority(value);
	};
	
	if(slotXML.hasOption("Origin")){
		std::string input = slotXML.get("Origin");
		input = bitpit::utils::trim(input);
		darray3E temp = {{0.0,0.0,0.0}};
		if(!input.empty()){
			std::stringstream ss(input);
			for(auto &val : temp) ss>>val;
		}
		setOrigin(temp);
	} 
	
	if(slotXML.hasOption("Direction")){
		std::string input = slotXML.get("Direction");
		input = bitpit::utils::trim(input);
		darray3E temp = {{0.0,0.0,0.0}};
		if(!input.empty()){
			std::stringstream ss(input);
			for(auto &val : temp) ss>>val;
		}
		setDirection(temp);
	} 
	
	if(slotXML.hasOption("Rotation")){
		std::string input = slotXML.get("Rotation");
		input = bitpit::utils::trim(input);
		double temp = 0.0;
		if(!input.empty()){
			std::stringstream ss(input);
			ss>>temp;
		}
		setRotation(temp);
	} 
	
	if(slotXML.hasSection("RefSystem")){
		const bitpit::Config::Section & rfXML = slotXML.getSection("RefSystem");
		std::string rootAxis = "axis";
		std::string axis;
		dmatrix33E temp;
		temp[0].fill(0.0); temp[0][0] = 1.0;
		temp[1].fill(0.0); temp[1][1] = 1.0;
		temp[2].fill(0.0); temp[2][2] = 1.0;
		for(int i=0; i<3; ++i){			
			axis = rootAxis + std::to_string(i);
			std::string input = rfXML.get(axis);
			input = bitpit::utils::trim(input);
			if(!input.empty()){
				std::stringstream ss(input);
				for(auto &val : temp[i]) ss>>val;
			}
		}
		setAxes(temp);
	} 
	
	if(slotXML.hasOption("OriginRS")){
		std::string input = slotXML.get("OriginRS");
		input = bitpit::utils::trim(input);
		darray3E temp = {{0.0,0.0,0.0}};
		if(!input.empty()){
			std::stringstream ss(input);
			for(auto &val : temp) ss>>val;
		}
		setAxesOrigin(temp);
	} 
};	
/*!
 * Write settings of the class from bitpit::Config::Section slot. Reimplemented from
 * BaseManipulation::absorbSectionXML.The class read essential parameters to perform rotation of an axes ref system.
 * 
 * --> Flushing data// how to write it on XML:
 * 		ClassName : name of the class as "MiMMO.RotationBox"
 * 		Priority  : uint marking priority in multi-chain execution; 
 * 		Origin: rotation axis origin
 * 		Direction: axis direction coordinates
 * 		Rotation : rotation angle in radians. Positive on counterclockwise rotations around reference axis
 * 		RefSystem: axes of current shape reference system. written in XML as:
 * 					<RefSystem>
 * 						<axis0>	1.0 0.0 0.0 </axis0>
 * 						<axis1>	0.0 1.0 0.0 </axis1>
 * 						<axis2>	0.0 0.0 1.0 </axis2>
 * 					</RefSystem>
 * 		OriginRS: origin of the target reference system to be rotated
 * 
 * \param[in] slotXML bitpit::Config::Section which writes to
 * \param[in] name   name associated to the slot
 */
void RotationBox::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
	slotXML.set("ClassName", m_name);
	slotXML.set("Priority", std::to_string(getPriority()));
	
	
	{
		std::stringstream ss;
		ss<<std::scientific<<m_origin[0]<<'\t'<<m_origin[1]<<'\t'<<m_origin[2];
		slotXML.set("Origin", ss.str());
	}
	
	{
		std::stringstream ss;
		ss<<std::scientific<<m_direction[0]<<'\t'<<m_direction[1]<<'\t'<<m_direction[2];
		slotXML.set("Direction", ss.str());
	}
	
	slotXML.set("Rotation", std::to_string(m_alpha));	
	
	{
		bitpit::Config::Section & rsXML = slotXML.addSection("RefSystem");
		std::string rootAxis = "axis";
		std::string localAxis;
		int counter=0;
		for(auto &axis : m_axes){
			localAxis = rootAxis+std::to_string(counter);
			std::stringstream ss;
			ss<<std::scientific<<axis[0]<<'\t'<<axis[1]<<'\t'<<axis[2];
			rsXML.set(localAxis, ss.str());
			++counter;
		}
	}
	
	{
		std::stringstream ss;
		ss<<std::scientific<<m_axes_origin[0]<<'\t'<<m_axes_origin[1]<<'\t'<<m_axes_origin[2];
		slotXML.set("OriginRS", ss.str());
	}

};	



