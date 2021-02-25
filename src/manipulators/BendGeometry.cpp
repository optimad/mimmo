/*---------------------------------------------------------------------------*\
 *
 *  mimmo
 *
 *  Copyright (C) 2015-2021 OPTIMAD engineering Srl
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
#include "BendGeometry.hpp"

namespace mimmo{


/*!
    Default constructor of BendGeometry
 */
BendGeometry::BendGeometry(){
    m_name = "mimmo.BendGeometry";
    m_degree.fill({{0,0,0}});
    m_origin.fill(0.0);
    m_system = { 1.0 , 0.0 , 0.0 , 0.0, 1.0, 0.0 , 0.0 , 0.0 , 1.0 };
    m_local = false;
};


/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
BendGeometry::BendGeometry(const bitpit::Config::Section & rootXML){

    m_name = "mimmo.BendGeometry";
    m_degree.fill({{0,0,0}});
    m_origin.fill(0.0);
    m_system = { 1.0 , 0.0 , 0.0 , 0.0, 1.0, 0.0 , 0.0 , 0.0 , 1.0 };
    m_local = false;

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::string::trim(input);
    if(input == "mimmo.BendGeometry"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*!
    Default destructor of BendGeometry
 */
BendGeometry::~BendGeometry(){};

/*!
    Copy constructor of BendGeometry. Resulting displacements are not copied.
 */
BendGeometry::BendGeometry(const BendGeometry & other):BaseManipulation(other){
    m_origin = other.m_origin;
    m_system = other.m_system;
    m_local  = other.m_local;
    m_filter = other.m_filter;
    m_degree = other.m_degree;
    m_coeffs = other.m_coeffs;
};

/*!
 * Assignment operator.
 */
BendGeometry & BendGeometry::operator=(BendGeometry other){
    swap(other);
    return *this;
}

/*!
 * Swap function
 * \param[in] x object to be swapped
 */
void BendGeometry::swap(BendGeometry & x) noexcept
{
   std::swap(m_origin, x.m_origin);
   std::swap(m_system, x.m_system);
   std::swap(m_local , x.m_local);
   m_filter.swap(x.m_filter);
   std::swap(m_degree, x.m_degree);
   std::swap(m_coeffs, x.m_coeffs);
   m_displ.swap(x.m_displ);
   BaseManipulation::swap(x);
}
/*!
It builds the input/output ports of the object
 */
void
BendGeometry::buildPorts(){
    bool built = true;
    built = (built && createPortIn<MimmoSharedPointer<MimmoObject>, BendGeometry>(&m_geometry, M_GEOM, true));
    built = (built && createPortIn<dmpvector1D*, BendGeometry>(this, &mimmo::BendGeometry::setFilter, M_FILTER));
    built = (built && createPortIn<umatrix33E, BendGeometry>(&m_degree, M_BMATRIX));
    built = (built && createPortIn<dmat33Evec*, BendGeometry>(this, &BendGeometry::setCoeffs, M_BCOEFFS));
    built = (built && createPortIn<dmatrix33E, BendGeometry>(&m_system, M_AXES));
    built = (built && createPortIn<darray3E, BendGeometry>(&m_origin, M_POINT));
    built = (built && createPortOut<dmpvecarr3E*, BendGeometry>(this, &mimmo::BendGeometry::getDisplacements, M_GDISPLS));
    built = (built && createPortOut<MimmoSharedPointer<MimmoObject>, BendGeometry>(this, &BaseManipulation::getGeometry, M_GEOM));
    m_arePortsBuilt = built;
};

/*!
    \return current origin of reference system
 */
darray3E
BendGeometry::getOrigin(){
    return(m_origin);
}

/*!
    \return current reference system axes
 */
dmatrix33E
BendGeometry::getRefSystem(){
    return(m_system);
}

/*!It gets the degrees of polynomial law for each component of displacements of degrees of freedom.
 * \return Degrees of polynomial laws (degree[i][j] = degree of displacement function si = f(xj)).
 */
umatrix33E
BendGeometry::getDegree(){
    return(m_degree);
};

/*!
    It gets the matrix of coefficients of the polynomial functions.
    coeffs[i][j] = \f${a_{ij0},a_{ij1},a_{ij2},...,a_{ijN}}\f$
 * \return polynomials coefficents matrix.
 */
dmat33Evec *
BendGeometry::getCoeffs(){
    return  &m_coeffs;
};

/*!
    \return bending displacements computed on geometry nodes
 */
dmpvecarr3E*
BendGeometry::getDisplacements(){
    return &m_displ;
};

/*!
    It sets the degrees matrix (3x3) of each polynomial law.
 * \param[in] degree matrix of degrees \f$d_{ij}\f$
 */
void
BendGeometry::setDegree(umatrix33E degree){
    for(int i=0; i<3; ++i){
        for(int j=0; j<3; ++j){
            setDegree(i,j, degree[i][j]);
        }
    }
};

/*!
    It sets the degree \f$d_{ij}\f$ of polynomial \f$P_{ij}\f$.
 * \param[in] i component of displacement vector.
 * \param[in] j coordinate of the axis
 * \param[in] degree \f$d_{ij}\f$
 */
void
BendGeometry::setDegree(int i, int j, uint32_t degree){
    m_degree[i][j] = degree;
    m_coeffs[i][j].resize(degree+1, 0.0);
};

/*!
    It sets the matrix of coefficients of the polynomial functions.
    Singular entry ij of the matrix is a vector
    coeffs[i][j] = \f${a_{ij0},a_{ij1},a_{ij2},...,a_{ijN}}\f$
 * \param[in] coeffs pointer to matrix of coefficients
 */
void
BendGeometry::setCoeffs(dmat33Evec * coeffs){
    if(!coeffs) return;
    m_coeffs = *coeffs;
};

/*!
   It sets the coefficients of polynomial \f$P_{ij}\f$
 * \param[in] i component of displacement vector
 * \param[in] j coordinate of the axis
 * \param[in] coeffs vector of polynomial coefficients coeffs[i][j] = \f${a_{ij0},a_{ij1},a_{ij2},...,a_{ijN}}\f$
 */
void
BendGeometry::setCoeffs(int i, int j, dvector1D coeffs){
    m_coeffs[i][j] = coeffs;
};

/*!
It sets the origin of the local reference system.
 * \param[in] origin  of local reference system.
 */
void
BendGeometry::setOrigin(darray3E origin){
    m_origin = origin;
    m_local = true;
}

/*!
    It sets a new axes orientation frame for bending s.d.r.
 * \param[in] axes matrix.
 */
void
BendGeometry::setRefSystem(dmatrix33E axes){
    m_system[0] = axes[0];
    m_system[1] = axes[1];
    m_system[2] = axes[2];
    m_local        = true;
}

/*!
   It sets a filter field to modulate the displacementsof the target geometry nodes.
 * \param[in] filter field defined on geometry vertices.
 */
void
BendGeometry::setFilter(dmpvector1D * filter){
    m_filter = *filter;
}

/*!
    Execution command. It computes the displacements on geometry nodes
    with the bending polynomial set.
    The displacements are stored in local m_displ variable.
 */
void
BendGeometry::execute(){

    if(getGeometry() == nullptr){
        (*m_log)<<m_name + " : nullptr pointer to linked geometry found"<<std::endl;
        throw std::runtime_error(m_name + "nullptr pointer to linked geometry found");
    }

    m_displ.clear();
    m_displ.setDataLocation(mimmo::MPVLocation::POINT);
    m_displ.reserve(getGeometry()->getNVertices());
    m_displ.setGeometry(getGeometry());


    //check coherence of degrees and coeffs;
    for(int i=0; i<3; ++i){
        for(int j=0; j<3; ++j){
            m_coeffs[i][j].resize(m_degree[i][j]+1, 0.0);
        }
    }

    checkFilter();

    long ID;
    darray3E value;
    darray3E point, point0;
    for (const auto & vertex : m_geometry->getVertices()){
        point = vertex.getCoords();
        if (m_local){
            point0 = point;
            point = toLocalCoord(point);
        }
        ID = vertex.getId();
        value.fill(0.0);
        for (int j=0; j<3; j++){
            for (int z=0; z<3; z++){
                if (m_degree[j][z] > 0){
                    for (int k=0; k<(int)m_degree[j][z]+1; k++){
                        value[j] += pow(point[z],(double)k)*m_coeffs[j][z][k]*m_filter[ID];
                    }
                }
            }
        }
        if (m_local){
            point += value;
            point = toGlobalCoord(point);
            value = point - point0;
        }
        m_displ.insert(ID, value);
    }
};

/*!
 * Directly apply deformation field to target geometry.
 */
void
BendGeometry::apply(){
    _apply(m_displ);
}

/*!
 * Check if the filter is related to the target geometry.
 * If not create a unitary filter field.
 */
void
BendGeometry::checkFilter(){
    bool check = m_filter.getDataLocation() == mimmo::MPVLocation::POINT;
    check = check && m_filter.completeMissingData(0.0);
    check = check && m_filter.getGeometry() == getGeometry();

    if (!check){
        m_log->setPriority(bitpit::log::Verbosity::DEBUG);
        (*m_log)<<"Not valid filter found in "<<m_name<<". Proceeding with default unitary field"<<std::endl;
        m_log->setPriority(bitpit::log::Verbosity::NORMAL);

        m_filter.clear();
        m_filter.setGeometry(m_geometry);
        m_filter.setDataLocation(mimmo::MPVLocation::POINT);
        m_filter.reserve(getGeometry()->getNVertices());
        for (const auto & vertex : getGeometry()->getVertices()){
            m_filter.insert(vertex.getId(), 1.0);
        }
    }
}

/*! It gets the local coordinates of a point w.r.t. the bending local reference system
 * \param[in] point Input point global coordinates
 * \return Local point coordinates
 */
darray3E
BendGeometry::toLocalCoord(darray3E  point){
    darray3E work;
    //unapply origin translation
    work = point - m_origin;
    //apply change to local sdr transformation
    return  matmul(m_system, work);
};

/*!
 * Transform point from local reference system of the bending,
 * to world reference system.
 * \param[in] point target
 * \return transformed point
 */
darray3E    BendGeometry::toGlobalCoord(darray3E  point){

    darray3E work, work2;
    //unscale your local point
    for(int i =0; i<3; ++i){
        work[i] = point[i];
    }

    //unapply change to local sdr transformation
    work2 = matmul(work, m_system);

    //unapply origin translation
    return (work2 + m_origin);
};

/*!
  Matrix M-vector V multiplication V*M (V row dot M columns).
  \param[in] vec vector V with 3 elements
  \param[in] mat matrix M 3x3 square
  \return 3 element vector result of multiplication
*/
darray3E   BendGeometry::matmul(const darray3E & vec, const dmatrix33E & mat){

    darray3E out;
    for (std::size_t i = 0; i < 3; i++) {
        out[i] = 0.0;
        for (std::size_t j = 0; j <3; j++) {
            out[i] += vec[j]*mat[j][i];
        } //next j
    } //next i

    return out;
}

/*!
  Matrix M-vector V multiplication M*V (M rows dot V column).
  \param[in] vec vector V with 3 elements
  \param[in] mat matrix M 3x3 square
  \return 3 element vector result of multiplication
*/
darray3E   BendGeometry::matmul(const dmatrix33E & mat, const darray3E & vec){

    darray3E out;
    for (std::size_t i = 0; i < 3; i++) {
        out[i] = 0.0;
        for (std::size_t j = 0; j <3; j++) {
            out[i] += vec[j]*mat[i][j];
        } //next j
    } //next i

    return out;
}

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
BendGeometry::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    std::string input;

    BaseManipulation::absorbSectionXML(slotXML, name);

    if(slotXML.hasSection("DegreesMatrix")){
        auto & subslot = slotXML.getSection("DegreesMatrix");
        umatrix33E temp;
        temp.fill({{0,0,0}});

        if(subslot.hasOption("xDispl")){
            input = subslot.get("xDispl");
            std::stringstream ss(bitpit::utils::string::trim(input));
            for(int i=0; i<3; ++i)    ss>>temp[0][i];
        }

        if(subslot.hasOption("yDispl")){
            input = subslot.get("yDispl");
            std::stringstream ss(bitpit::utils::string::trim(input));
            for(int i=0; i<3; ++i)    ss>>temp[1][i];
        }

        if(subslot.hasOption("zDispl")){
            input = subslot.get("zDispl");
            std::stringstream ss(bitpit::utils::string::trim(input));
            for(int i=0; i<3; ++i)    ss>>temp[2][i];
        }

        setDegree(temp);
    };

    // set Degree already set correctly the matrix of coeffs..
    if(slotXML.hasSection("PolyCoefficients")){
        auto & subslot = slotXML.getSection("PolyCoefficients");
        dmat33Evec &temp = m_coeffs;
        std::string rootPoly = "Poly";
        std::string locPoly;
        int ik,jk;
        for (int k=0; k<9; ++k){
            locPoly = rootPoly + std::to_string(k);
            if(subslot.hasOption(locPoly)){
                input = subslot.get(locPoly);
                std::stringstream ss(bitpit::utils::string::trim(input));
                ik = (int)(k/3);
                jk = k%3;
                for(double & val: temp[ik][jk])    ss>>val;
            }
        }
    };

    if(slotXML.hasOption("Origin")){
        std::string input = slotXML.get("Origin");
        input = bitpit::utils::string::trim(input);
        darray3E temp = {{0.0,0.0,0.0}};
        if(!input.empty()){
            std::stringstream ss(input);
            for(auto &val : temp) ss>>val;
        }
        setOrigin(temp);
    };

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
            input = bitpit::utils::string::trim(input);
            if(!input.empty()){
                std::stringstream ss(input);
                for(auto &val : temp[i]) ss>>val;
            }
        }
        setRefSystem(temp);
    };
};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
BendGeometry::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    BaseManipulation::flushSectionXML(slotXML, name);

    bitpit::Config::Section & degXML = slotXML.addSection("DegreesMatrix");

    svector1D d_str(3);
    int counter = 0;
    for(auto & val: m_degree){
        std::stringstream ss;
        for(auto & loc: val)    ss<<loc<<'\t';
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
            for(auto & loc: m_coeffs[i][j])    ss<<loc<<'\t';
            polyXML.set(locPoly, ss.str());
        }
    }

    {
        std::stringstream ss;
        ss<<std::scientific<<getOrigin()[0]<<'\t'<<getOrigin()[1]<<'\t'<<getOrigin()[2];
        slotXML.set("Origin", ss.str());
    }

    {
        auto rs = getRefSystem();
        bitpit::Config::Section & rsXML = slotXML.addSection("RefSystem");
        std::string rootAxis = "axis";
        std::string localAxis;
        int counter=0;
        for(auto &axis : rs){
            localAxis = rootAxis+std::to_string(counter);
            std::stringstream ss;
            ss<<std::scientific<<axis[0]<<'\t'<<axis[1]<<'\t'<<axis[2];
            rsXML.set(localAxis, ss.str());
            ++counter;
        }
    }

};

}
