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
#include "Bend.hpp"

namespace mimmo{


/*!Default constructor of Bend
 */
Bend::Bend(){
	m_name = "MiMMO.Bend";
	m_degree.fill({{0,0,0}});
};


/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
Bend::Bend(const bitpit::Config::Section & rootXML){
	
	m_name = "MiMMO.Bend";
	m_degree.fill({{0,0,0}});
	
	std::string fallback_name = "ClassNONE";	
	std::string input = rootXML.get("ClassName", fallback_name);
	input = bitpit::utils::trim(input);
	if(input == "MiMMO.Bend"){
		absorbSectionXML(rootXML);
	}else{	
		std::cout<<"Warning in custom xml MiMMO::Bend constructor. No valid xml data found"<<std::endl;
	};
}

/*!Default destructor of Bend
 */
Bend::~Bend(){};

/*!Copy constructor of Bend.
 */
Bend::Bend(const Bend & other):BaseManipulation(other){
	m_coords = other.m_coords;
	m_degree = other.m_degree;
	m_coeffs = other.m_coeffs;
};

/*!Assignement operator of Bend.
 */
Bend & Bend::operator=(const Bend & other){
	*(static_cast<BaseManipulation*> (this)) = *(static_cast<const BaseManipulation*> (&other));
	m_coords = other.m_coords;
	m_degree = other.m_degree;
	m_coeffs = other.m_coeffs;
	return	*this;
};


/*! It builds the input/output ports of the object
 */
void
Bend::buildPorts(){
	bool built = true;
	built = (built && createPortIn<dvecarr3E, Bend>(&m_displ, PortType::M_DISPLS, mimmo::pin::containerTAG::VECARR3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<dvecarr3E, Bend>(&m_coords, PortType::M_COORDS, mimmo::pin::containerTAG::VECARR3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<umatrix33E, Bend>(&m_degree, PortType::M_BMATRIX, mimmo::pin::containerTAG::ARR3ARR3, mimmo::pin::dataTAG::INT));
	built = (built && createPortIn<dmat33Evec, Bend>(&m_coeffs, PortType::M_BCOEFFS, mimmo::pin::containerTAG::ARR3ARR3VEC, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortOut<dvecarr3E, Bend>(this, &mimmo::Bend::getDisplacements, PortType::M_DISPLS, mimmo::pin::containerTAG::VECARR3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortOut<dvecarr3E, Bend>(this, &mimmo::Bend::getCoords, PortType::M_COORDS, mimmo::pin::containerTAG::VECARR3, mimmo::pin::dataTAG::FLOAT));
	m_arePortsBuilt = built;
};

/*!It gets the coordinates of the degrees of freedom.
 * \return Coordinates of the degrees of freedom.
 */
dvecarr3E
Bend::getCoords(){
	return(m_coords);
};

/*!It gets the degrees of polynomial law for each component of displacements of degrees of freedom.
 * \return Degrees of polynomial laws (degree[i][j] = degree of displacement function si = f(xj)).
 */
umatrix33E
Bend::getDegree(){
	return(m_degree);
};

/*!It gets the coefficients of the polynomial laws.
 * \return Coefficients of the polynomial laws. (coeffs[i][j][k] = coefficients aijk of term si = aij * xj^k).
 */
dmat33Evec
Bend::getCoeffs(){
	return(m_coeffs);
};

/*!It gets the displacements of the degrees of freedom (modified after the execution).
 * \return Displacements of the degrees of freedom.
 */
dvecarr3E
Bend::getDisplacements(){
	return(m_displ);
};

/*!It sets the initial displacements of the degrees of freedom.
 * \param[in] displ Displacements of the degrees of freedom.
 */
void
Bend::setDisplacements(dvecarr3E displ){
	m_displ = displ;
};


/*!It sets the coordinates of the degrees of freedom.
 * \param[in] coords Coordinates of the degrees of freedom.
 */
void
Bend::setCoords(dvecarr3E coords){
	m_coords = coords;
};

/*!It sets the degrees of polynomial law for each component of displacements of degrees of freedom.
 * \param[in] degree Degrees of polynomial laws (degree[i][j] = degree of displacement function si = f(xj)).
 */
void
Bend::setDegree(umatrix33E degree){
	for(int i=0; i<3; ++i){
		for(int j=0; j<3; ++j){
			setDegree(i,j, degree[i][j]);
		}
	}
};

/*!It sets the degrees of a term of a polynomial law for a component of displacements of degrees of freedom.
 * \param[in] i Components of displacement.
 * \param[in] j Coordinate of the function related to input degree.
 * \param[in] degree Degrees of polynomial laws (degree[i][j] = degree of displacement function si = f(xj)).
 */
void
Bend::setDegree(int i, int j, uint32_t degree){
	m_degree[i][j] = degree;
	m_coeffs[i][j].resize(degree+1, 0.0);
};

/*!It sets the coefficients of the polynomial laws.
 * \param[in] coeffs Coefficients of the polynomial laws. (coeffs[i][j][k] = coefficients aijk of term si = aij * xj^k).
 */
void
Bend::setCoeffs(dmat33Evec coeffs){
	m_coeffs = coeffs;
};

/*!It sets the coefficients of the polynomial laws.
 * \param[in] i Components of displacement.
 * \param[in] j Coordinate of the function related to input degree.
 * \param[in] coeffs Coefficients of the polynomial laws. (coeffs[i][j][k] = coefficients aijk of term si = aij * xj^k).
 */
void
Bend::setCoeffs(int i, int j, dvector1D coeffs){
	m_coeffs[i][j] = coeffs;
};

/*!Execution command. It modifies the displacements given by the input manipulation object
 * with the polynomial law.
 * The input has to be set with a dvecarr3E variable (bend it casts the template method
 * getInput to this type) and the result will be of the same type.
 * After exec() the modified displacements are stored in result of base class.
 */
void
Bend::execute(){
	
	//check coherence of degrees and coeffs;
	for(int i=0; i<3; ++i){
		for(int j=0; j<3; ++j){
			m_coeffs[i][j].resize(m_degree[i][j]+1, 0.0);
		}
	}
	
	int	ndispl = m_displ.size();
	ndispl = std::min(ndispl, int(m_coords.size()));
	for (int j=0; j<3; j++){
		for (int i=0; i<ndispl; i++){
			for (int z=0; z<3; z++){
				if (m_degree[j][z] > 0){
					for (int k=0; k<(int)m_degree[j][z]+1; k++){
						m_displ[i][j] += pow(m_coords[i][z],(double)k)*m_coeffs[j][z][k];
					}
				}
			}
		}
	}
	return;
};

/*!
 * Get settings of the class from bitpit::Config::Section slot. Reimplemented from
 * BaseManipulation::absorbSectionXML.The class read only RefreshGeometryTrees parameter, 
 * while Input and Geometry parameters are meant to be passed only through Port linking.
 * 
 * Assuming that x,y,z are the reference directions in space, as well as the preference directions
 * for nodal displacements, we can assume that:
 * 
 * --> Absorbing data:
 * - <B>Priority</B>  : uint marking priority of class execution in multichain frame	
 * - <B>DegreesMatrix(3x3)</B> : degrees of each polynomial function referred to a displacement 
 *      in direction i (x,y,z) and modulating displacement in direction j (x,y,z). Degree 0
 *      marks a constant function 
 * - <B>PolyCoefficients</B>: coefficients of each 9 bending polynomial functions.
 * 
 * \param[in] slotXML bitpit::Config::Section which reads from
 * \param[in] name   name associated to the slot
 */
void Bend::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){
	
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

	if(slotXML.hasSection("DegreesMatrix")){
		auto & subslot = slotXML.getSection("DegreesMatrix");
		umatrix33E temp;
		temp.fill({{0,0,0}});
		
		if(subslot.hasOption("xDispl")){
			input = subslot.get("xDispl");
			std::stringstream ss(bitpit::utils::trim(input));
			for(int i=0; i<3; ++i)	ss>>temp[0][i];
		}
		
		if(subslot.hasOption("yDispl")){
			input = subslot.get("yDispl");
			std::stringstream ss(bitpit::utils::trim(input));
			for(int i=0; i<3; ++i)	ss>>temp[1][i];
		}
		
		if(subslot.hasOption("yDispl")){
			input = subslot.get("yDispl");
			std::stringstream ss(bitpit::utils::trim(input));
			for(int i=0; i<3; ++i)	ss>>temp[2][i];
		}
		
		setDegree(temp);
	}; 
	
	if(slotXML.hasSection("PolyCoefficients")){
		auto & subslot = slotXML.getSection("PolyCoefficients");
		dmat33Evec temp = getCoeffs();
		std::string rootPoly = "Poly";
		std::string locPoly;
		int ik,jk;
		for (int k=0; k<9; ++k){
			locPoly = rootPoly + std::to_string(k);
			if(subslot.hasOption(locPoly)){
				input = subslot.get(locPoly);
				std::stringstream ss(bitpit::utils::trim(input));
				ik = (int)(k/3);
				jk = k%3;
				for(auto & val: temp[ik][jk])	ss>>val;
			}
		}
		setCoeffs(temp);
	}; 
};

/*!
 * Write settings of the class to bitpit::Config::Section slot. Reimplemented from
 * BaseManipulation::flushSectionXML;
 * The class write only RefreshGeometryTrees parameter, if it is different from its default value, 
 * while Input and Geometry parameters are meant to be passed only through Port linking.
 * 
 *    --> Flushing data// how to write it on XML:
 * - <B>ClassName</B> : name of the class as "MiMMO.Bend"
 * - <B>Priority</B>  : uint marking priority of class execution in multichain frame	
 *
 * - <B>DegreesMatrix(3x3)</B>: degrees of each polynomial function referred to a displacement 
 * 					        in direction i (x,y,z) and modulating displacement in direction j (x,y,z). Degree 0
 * 						    marks a constant function.
 *							Written in XML as:
 * 							<DegreesMatrix>
 *								<xDispl> 1 0 0 </xDispl> (linear x-displacement distribution in x bending direction)
 *								<yDispl> 2 0 0 </yDispl> (quadratic y-displacement distribution in x bending direction)
 *								<zDispl> 0 3 0  </zDispl> (cubic z-displacement in y bending direction)
 *							</DegreesMatrix>
 *
 * - <B>PolyCoefficients</B>: coefficients of each 9 bending polynomial functions. Writing following 
 * 						  the enumeration n = i*3 + j, where i is the displacement direction and j the bending direction.
 * 						  For example n=7 corresponds to i=2, j=1, reflecting in a z displacement distribution in y-bending direction. 
 * 						  Please note the number of coefficients for a ij-bending function is equal to the degree
 * 						  of freedom DegreesMatrix(i,j), ordered as c0 + c1*x +c2*x^2+...
 * 						  Written in xml as :
 * 						  <PolyCoefficients>
 * 							<Poly0> 1.0 1.5 </Poly0>
 * 							<Poly3> -0.1 0.2 -0.01 </Poly3>
 * 							<Poly7> 1.5 0.0 0.1 0.2 </Poly7>
 * 						  </PolyCoefficients>
 * 
 * \param[in]	slotXML bitpit::Config::Section which writes to
 * \param[in] name   name associated to the slot
 */
void Bend::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
	BITPIT_UNUSED(name);
	
	slotXML.set("ClassName", m_name);
	slotXML.set("Priority", std::to_string(getPriority()));
	
	bitpit::Config::Section & degXML = slotXML.addSection("DegreesMatrix");
	
	svector1D d_str(3);
	int counter = 0;
	for(auto & val: m_degree){
		std::stringstream ss;
		for(auto & loc: val)	ss<<loc<<'\t';
		d_str[counter] = ss.str();
		++counter;
	}
	
	degXML.set("xDispl", d_str[0]);
	degXML.set("yDispl", d_str[1]);
	degXML.set("zDispl", d_str[2]);
	
	bitpit::Config::Section & polyXML = slotXML.addSection("PolyCoefficients");
	std::string rootPoly = "Poly";
	std::string locPoly;
	for(int i=0; i<3; i++){
		for(int j =0; j<3; j++){
			if(m_coeffs[i][j].empty()) continue;
			locPoly = rootPoly + std::to_string(int(i*3+j));
			std::stringstream ss;
			for(auto & loc: m_coeffs[i][j])	ss<<loc<<'\t';
			polyXML.set(locPoly, ss.str());
		}
	}
	
};	

}

