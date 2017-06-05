/*---------------------------------------------------------------------------*\
 *
 *  mimmo
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
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
\ *---------------------------------------------------------------------------*/

#include "Lattice.hpp"
#include "Operators.hpp"
#include "customOperators.hpp"


using namespace std;
namespace mimmo{

/*! 
 * Basic Constructor
 */
Lattice::Lattice(){
    m_np = 0;
    m_intMapDOF.clear();
    m_name = "mimmo.Lattice";
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
Lattice::Lattice(const bitpit::Config::Section & rootXML){
    m_np = 0;
    m_intMapDOF.clear();
    m_name = "mimmo.Lattice";
    std::string fallback_name = "ClassNONE";	
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::trim(input);
    if(input == "mimmo.Lattice"){
        absorbSectionXML(rootXML);
    }else{	
        warningXML(m_log, m_name);
    };
}

/*! Destructor */
Lattice::~Lattice(){};

/*! Copy Constructor
 *\param[in] other Lattice object
 */
Lattice::Lattice(const Lattice & other):BaseManipulation(), UStructMesh(){
    *this = other;
};

/*! Copy Operator
 * \param[in] other Lattice object
 */
Lattice & Lattice::operator=(const Lattice & other){

    *(static_cast<UStructMesh *>(this))  = *(static_cast<const UStructMesh *>(&other));
    m_intMapDOF = other.m_intMapDOF;
    m_np = other.m_np;
    return(*this);
};

/*!
 * It builds the input/output ports of the object
 */
void Lattice::buildPorts(){

    bool built = true;
    built = (built && createPortIn<MimmoObject*, Lattice>(&m_geometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_, true));
    built = (built && createPortIn<iarray3E, Lattice>(this, &mimmo::Lattice::setDimension, PortType::M_DIMENSION, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::INT));
    built = (built && createPortIn<darray3E, Lattice>(this, &mimmo::Lattice::setInfLimits, PortType::M_INFLIMITS, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::FLOAT));
    built = (built && createPortIn<dmatrix33E, Lattice>(this, &mimmo::Lattice::setRefSystem, PortType::M_AXES, mimmo::pin::containerTAG::ARR3ARR3, mimmo::pin::dataTAG::FLOAT));
    built = (built && createPortIn<darray3E, Lattice>(this, &mimmo::Lattice::setSpan, PortType::M_SPAN, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::FLOAT));
    built = (built && createPortIn<darray3E, Lattice>(this, &mimmo::Lattice::setOrigin, PortType::M_POINT, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::FLOAT));
    built = (built && createPortIn<mimmo::ShapeType, Lattice>(this, &mimmo::Lattice::setShape, PortType::M_SHAPE, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::SHAPET));
    built = (built && createPortIn<const BasicShape *, Lattice>(this, &mimmo::Lattice::setShape, PortType::M_COPYSHAPE, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::SHAPE_));
    built = (built && createPortIn<int, Lattice>(this, &mimmo::Lattice::setShape, PortType::M_SHAPEI, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::INT));

    // creating output ports
    built = (built && createPortOut<dvecarr3E, Lattice>(this, &mimmo::Lattice::getGlobalCoords, PortType::M_GLOBAL, mimmo::pin::containerTAG::VECARR3, mimmo::pin::dataTAG::FLOAT));
    built = (built && createPortOut<dvecarr3E, Lattice>(this, &mimmo::Lattice::getLocalCoords, PortType::M_LOCAL, mimmo::pin::containerTAG::VECARR3, mimmo::pin::dataTAG::FLOAT));
    built = (built && createPortOut<darray3E, Lattice>(this, &mimmo::Lattice::getOrigin, PortType::M_POINT, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::FLOAT));
    built = (built && createPortOut<dmatrix33E, Lattice>(this, &mimmo::Lattice::getRefSystem, PortType::M_AXES, mimmo::pin::containerTAG::ARR3ARR3, mimmo::pin::dataTAG::FLOAT));
    built = (built && createPortOut<darray3E, Lattice>(this, &mimmo::Lattice::getSpan, PortType::M_SPAN, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::FLOAT));
    built = (built && createPortOut<iarray3E, Lattice>(this, &mimmo::Lattice::getDimension, PortType::M_DIMENSION, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::INT));
    built = (built && createPortOut<BasicShape *, Lattice>(this, &mimmo::Lattice::getShape, PortType::M_COPYSHAPE, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::SHAPE_));
    built = (built && createPortOut<MimmoObject*, Lattice>(this, &mimmo::Lattice::getGeometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));

    m_arePortsBuilt = built;
};

/*!
 * Clean every setting and data hold by the class
 */
void Lattice::clearLattice(){
    clear(); //base manipulation stuff clear
    clearMesh(); // structured mesh cleaned
    m_intMapDOF.clear();
};

/*! 
 * Get the total number of control nodes.
 * \return number of control nodes
 */
int Lattice::getNNodes(){
    return(m_np);
}

/*! 
 * Return position of effective mesh nodes in the lattice, in absolute reference system.
 *  Reimplemented from UstructMesh::getGlobalCoords.
 * \return effective mesh nodes position
 */
dvecarr3E
Lattice::getGlobalCoords(){
    int np = (getNNodes());
    dvecarr3E coords(np);
    int index, i0, i1, i2;
    for (int i=0; i<np; i++){
        index = accessGridFromDOF(i);
        accessPointIndex(index,i0,i1,i2);
        coords[i] = getGlobalPoint(i0,i1,i2);
    }
    return(coords);
};

/*! 
 * Return position of effective mesh nodes in the lattice, in local shape reference system.
 *  Reimplemented from UstructMesh::getLocalCoords.
 * \return effective mesh nodes position
 */
dvecarr3E
Lattice::getLocalCoords(){
    int np = (getNNodes());
    dvecarr3E coords(np);
    int index, i0, i1, i2;
    for (int i=0; i<np; i++){
        index = accessGridFromDOF(i);
        accessPointIndex(index,i0,i1,i2);
        coords[i] = getLocalPoint(i0,i1,i2);
    }
    return(coords);
};

/*! 
 * Find a corrispondent degree of freedom index of a lattice grid node.
 * Not found indices are marked as -1.
 * \param[in] index lattice grid global index
 * \return corrispondent DOF global index
 */
int
Lattice::accessDOFFromGrid(int index){
    return(m_intMapDOF[index]);
}

/*! 
 * Find a corrispondent lattice grid index of a degree of freedom node.
 * Not found indices are marked as -1.
 * \param[in] index DOF global index
 * \return corrispondent lattice grid global index
 */
int
Lattice::accessGridFromDOF(int index){
    return(posVectorFind(m_intMapDOF, index));
}

/*! 
 * Find a corrispondent degree of freedom index list of a lattice grid node list.
 * Not found indices are marked as -1.
 * \param[in] gNindex lattice grid global index list
 * \return corrispondent DOF global index list
 */
ivector1D
Lattice::accessDOFFromGrid(ivector1D gNindex) {
    ivector1D result(gNindex.size());
    int counter = 0;
    for(auto &val : gNindex){
        result[counter] = accessDOFFromGrid(val);
        ++counter;
    }
    return(result);
};

/*! 
 * Find a corrispondent lattice grid index list of a degree of freedom node list.
 * Not found indices are marked as -1.
 * \param[in] dofIndex DOF global index list
 * \return corrispondent lattice grid global index
 */
ivector1D
Lattice::accessGridFromDOF(ivector1D dofIndex){
    ivector1D result(dofIndex.size());
    int counter = 0;
    for(auto &val : dofIndex){
        result[counter] = accessGridFromDOF(val);
        ++counter;
    }
    return(result);
};

/*! 
 * Plot your current lattice as a structured grid to a VTK *.vtu file. 
 * Wrapped method of plotGrid of its base class UStructMesh.
 * \param[in] directory output directory
 * \param[in] filename  output filename w/out tag
 * \param[in] counter   integer identifier of the file
 * \param[in] binary     boolean flag for 0-"ascii" or 1-"appended" writing
 */
void
Lattice::plotGrid(std::string directory, std::string filename,int counter, bool binary){
    dvecarr3E* pnull = NULL;
    UStructMesh::plotGrid(directory, filename, counter, binary,  pnull);
};

/*! 
 * Plot your current lattice as a point cloud to VTK *.vtu file. 
 * Wrapped method of plotCloud of its base class UStructMesh.
 * \param[in] directory output directory
 * \param[in] filename  output filename w/out tag
 * \param[in] counter   integer identifier of the file
 * \param[in] binary     boolean flag for 0-"ascii" or 1-"appended" writing
 */
void
Lattice::plotCloud(std::string directory, std::string filename, int counter, bool binary){
    dvecarr3E* pnull = NULL;
    UStructMesh::plotCloud(directory, filename, counter, binary, pnull);
};

/*! 
 * Build the structured mesh and create a wrapped map of effective degree of freedom 
 * of the current lattice mesh. Map is stored in internal member m_intMapDOF and accessible 
 * through methods accessDOFFromGrid and accessGridFromDOF. Reimplemented from UstructMesh::build();
 */
void 
Lattice::build(){
    UStructMesh::build();
    resizeMapDof();
};

/*!
 * Execute your object, that is, recall the method Lattice::build(). 
 */
void
Lattice::execute(){
    build();
};

/*! 
 * Resize map of effective nodes of the lattice grid to fit 
 * a total number od degree of freedom nx*ny*nz. Old structure 
 * is deleted and reset to zero.
 */
void
Lattice::resizeMapDof(){
    //reallocate your displacement node
    m_intMapDOF.clear();
    m_intMapDOF.resize((m_nx+1)*(m_ny+1)*(m_nz+1), -1);
    ivector1D::iterator itMapBegin = m_intMapDOF.begin();
    ivector1D::iterator itMap = itMapBegin;
    ivector1D::iterator itMapEnd = m_intMapDOF.end();
    bvector1D info;
    m_np = reduceDimToDOF(m_nx+1,m_ny+1,m_nz+1, info);

    //set m_intMapDOF

    int target;
    int index;
    ivector1D dummy;

    int i0,i1,i2;
    switch(getShapeType()){

    case ShapeType::CYLINDER :
        target=0;
        while(itMap != itMapEnd){

            *itMap = target;
            index = std::distance(itMapBegin, itMap);
            accessPointIndex(index,i0,i1,i2);

            if(info[0] && i0 == 0){
                for(int k=0; k<=m_ny;++k){
                    m_intMapDOF[accessPointIndex(i0,k,i2)] = target;
                }
            }
            if(info[1] && i1 == 0){
                m_intMapDOF[accessPointIndex(i0,m_ny,i2)] = target;
            }

            itMap = find(m_intMapDOF.begin(), itMapEnd,-1);
            target++;
        }
        break;

    case ShapeType::SPHERE :

        target = 0;
        while(itMap != itMapEnd){

            *itMap = target;
            index = std::distance(itMapBegin, itMap);
            accessPointIndex(index,i0,i1,i2);

            if(info[0] && i0 == 0){
                for(int k1=0; k1<=m_ny;++k1){
                    for(int k2=0; k2<=m_nz; ++k2){
                        m_intMapDOF[accessPointIndex(i0,k1,k2)] = target;
                    }
                }
            }

            if(info[1] && i1 == 0){
                m_intMapDOF[accessPointIndex(i0,m_ny-1,i2)] = target;
            }

            if(info[2] && i2 == 0){
                for(int k1=0; k1<=m_ny; ++k1){
                    m_intMapDOF[accessPointIndex(i0,k1,i2)] = target;
                }
            }

            if(info[3] && i2 == (m_nz)){
                for(int k1=0; k1<=m_ny;++k1){
                    m_intMapDOF[accessPointIndex(i0,k1,i2)] = target;
                }
            }

            itMap = find(m_intMapDOF.begin(), itMapEnd,-1);
            target++;
        }
        break;


    case ShapeType::CUBE :
        target = 0;
        while(itMap != itMapEnd){

            *itMap = target;
            itMap = find(m_intMapDOF.begin(), itMapEnd,-1);
            target++;
        }
        break;

    default: //doing nothing
    break;
    }//end switch

}

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void Lattice::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

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

    if(slotXML.hasOption("Shape")){
        std::string input = slotXML.get("Shape");
        input = bitpit::utils::trim(input);

        if(input == "CYLINDER"){
            setShape(ShapeType::CYLINDER);
        }else if(input =="SPHERE"){
            setShape(ShapeType::SPHERE);
        }else{
            setShape(ShapeType::CUBE);
        }
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
    }; 

    if(slotXML.hasOption("Span")){
        std::string input = slotXML.get("Span");
        input = bitpit::utils::trim(input);
        darray3E temp = {{0.0,0.0,0.0}};
        if(!input.empty()){
            std::stringstream ss(input);
            for(auto &val : temp) ss>>val;
        }
        setSpan(temp);
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
            input = bitpit::utils::trim(input);
            if(!input.empty()){
                std::stringstream ss(input);
                for(auto &val : temp[i]) ss>>val;
            }
        }
        setRefSystem(temp);
    }; 

    if(slotXML.hasOption("InfLimits")){
        std::string input = slotXML.get("InfLimits");
        input = bitpit::utils::trim(input);
        darray3E temp = {{0.0,0.0,0.0}};
        if(!input.empty()){
            std::stringstream ss(input);
            for(auto &val : temp) ss>>val;
        }
        setSpan(temp);
    }; 

    if(slotXML.hasOption("Dimension")){
        std::string input = slotXML.get("Dimension");
        input = bitpit::utils::trim(input);
        iarray3E temp = {{2,2,2}};
        if(!input.empty()){
            std::stringstream ss(input);
            for(auto &val : temp) ss>>val;
        }
        setDimension(temp);
    };

    if(slotXML.hasOption("PlotInExecution")){
        std::string input = slotXML.get("PlotInExecution");
        input = bitpit::utils::trim(input);
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
        }
        setPlotInExecution(value);
    }

    if(slotXML.hasOption("OutputPlot")){
        std::string input = slotXML.get("OutputPlot");
        input = bitpit::utils::trim(input);
        std::string temp = ".";
        if(!input.empty())	setOutputPlot(input);
        else			  	setOutputPlot(temp);
    }

}

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void Lattice::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    slotXML.set("ClassName", m_name);
    slotXML.set("Priority", std::to_string(getPriority()));


    std::string towrite = "CUBE";

    if(getShapeType() == ShapeType::CYLINDER){
        towrite = "CYLINDER";
    } else if(getShapeType() == ShapeType::SPHERE){
        towrite = "SPHERE";
    }	
    slotXML.set("Shape", towrite);

    {
        std::stringstream ss;
        ss<<std::scientific<<getOrigin()[0]<<'\t'<<getOrigin()[1]<<'\t'<<getOrigin()[2];
        slotXML.set("Origin", ss.str());
    }

    {
        std::stringstream ss;
        ss<<std::scientific<<getSpan()[0]<<'\t'<<getSpan()[1]<<'\t'<<getSpan()[2];
        slotXML.set("Span", ss.str());
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

    {
        std::stringstream ss;
        ss<<std::scientific<<getInfLimits()[0]<<'\t'<<getInfLimits()[1]<<'\t'<<getInfLimits()[2];
        slotXML.set("InfLimits", ss.str());
    }

    {
        std::stringstream ss;
        ss<<getDimension()[0]<<'\t'<<getDimension()[1]<<'\t'<<getDimension()[2];
        slotXML.set("Dimension", ss.str());
    }

    if(isPlotInExecution()){
        slotXML.set("PlotInExecution", std::to_string(1));
    }

    if(m_outputPlot != "."){
        slotXML.set("OutputPlot", m_outputPlot);
    }

};



/*!
 * Plot Optional results of the class, that is, the lattice grid as an 
 * unstructured mesh of hexahedrons in a VTK *.vtu format.
 */
void 	Lattice::plotOptionalResults(){
    std::string dir = m_outputPlot;
    std::string nameGrid  = m_name+"GRID";
    plotGrid(dir, nameGrid, getClassCounter(), true );
}

/*!
 * \return the effective dof size of the lattice according to its shape.
 * \param[in] nx number of points in the first coordinate
 * \param[in] ny number of points in the second coordinate
 * \param[in] nz number of points in the third coordinate
 * \param[out] info boolean list, if all true, dof map of the lattice can be built safely. 
 */
int
Lattice::reduceDimToDOF(int nx, int ny, int nz, bvector1D & info){

    int delta = 0;
    double dval;
    switch(getShapeType()){

    case ShapeType::CYLINDER :
        delta += nz;
        nx--;
        if(getCoordType(1) == CoordType::PERIODIC)	ny--;

        info.push_back(true);
        info.push_back(getCoordType(1) == CoordType::PERIODIC);
        break;

    case ShapeType::SPHERE :
        delta ++;
        nx--;
        if(getCoordType(1) == CoordType::PERIODIC)	ny--;
        dval = getInfLimits()[2];
        if(dval == 0.0)	{
            nz--;
            delta += nx;
        }
        if((dval + getLocalSpan()[2]) == M_PI){
            nz--;
            delta += nx;
        }

        info.push_back(true);
        info.push_back(getCoordType(1) == CoordType::PERIODIC);
        info.push_back(dval==0.0);
        info.push_back((dval + getLocalSpan()[2]) == M_PI);
        break;

    default:
        //doing nothing
        break;
    }

    int result = nx*ny*nz + delta;
    return(result);
};

}
