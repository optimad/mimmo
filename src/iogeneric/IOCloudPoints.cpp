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
\*---------------------------------------------------------------------------*/

#include "IOCloudPoints.hpp"

namespace mimmo {

/*!
 * Default constructor of IOCloudPoints.
 * \param[in] readMode True if the object is in read mode, false if in Write mode.
 */
IOCloudPoints::IOCloudPoints(bool readMode){
    m_name         = "mimmo.IOCloudPoints";
    m_read         = readMode;
    m_template     = false;
    m_dir       = ".";
    m_filename     = m_name+"_source.dat";
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
IOCloudPoints::IOCloudPoints(const bitpit::Config::Section & rootXML){

    m_name         = "mimmo.IOCloudPoints";
    m_template     = false;
    m_dir       = ".";
    m_filename     = m_name+"_source.dat";
    m_read         = true;

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::string::trim(input);

    std::string fallback_name2 = "1";
    std::string input2 = rootXML.get("IOmode", fallback_name2);
    input2 = bitpit::utils::string::trim(input2);

    m_read = bool(std::atoi(input2.c_str()));

    if(input == "mimmo.IOCloudPoints"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}


IOCloudPoints::~IOCloudPoints(){};

/*!
 * Copy constructor of IOCloudPoints. Labels, points and data attached are no copied.
 */
IOCloudPoints::IOCloudPoints(const IOCloudPoints & other):BaseManipulation(other){
    m_read         = other.m_read;
    m_dir       = other.m_dir;
    m_filename     = other.m_filename;
    m_template = other.m_template;

};

/*!Assignement operator of IOCloudPoints. Labels, points and data attached are no copied.
 */
IOCloudPoints & IOCloudPoints::operator=(IOCloudPoints other){
    swap(other);
    return *this;
};

/*!
 * Swap method
 * \param[in] x object to be swapped
 */
void IOCloudPoints::swap(IOCloudPoints & x) noexcept
{
    std::swap(m_read        , x.m_read);
    std::swap(m_dir         , x.m_dir);
    std::swap(m_filename    , x.m_filename);
    std::swap(m_template    , x.m_template);
    std::swap(m_labels      , x.m_labels);
    std::swap(m_points      , x.m_points);
    std::swap(m_scalarfield , x.m_scalarfield);
    std::swap(m_vectorfield , x.m_vectorfield);
    BaseManipulation::swap(x);
}

/*!
 * It builds the input/output ports of the object
 */
void
IOCloudPoints::buildPorts(){
    bool built = true;

    built = (built && createPortIn<dvecarr3E, IOCloudPoints>(this, &mimmo::IOCloudPoints::setPoints, M_COORDS));
    built = (built && createPortIn<dvecarr3E, IOCloudPoints>(this, &mimmo::IOCloudPoints::setVectorField, M_DISPLS));
    built = (built && createPortIn<dvector1D, IOCloudPoints>(this, &mimmo::IOCloudPoints::setScalarField, M_DATAFIELD));
    built = (built && createPortIn<livector1D, IOCloudPoints>(this, &mimmo::IOCloudPoints::setLabels, M_VECTORLI));

    built = (built && createPortOut<dvecarr3E, IOCloudPoints>(this, &mimmo::IOCloudPoints::getPoints, M_COORDS));
    built = (built && createPortOut<dvecarr3E, IOCloudPoints>(this, &mimmo::IOCloudPoints::getVectorField, M_DISPLS));
    built = (built && createPortOut<dvector1D, IOCloudPoints>(this, &mimmo::IOCloudPoints::getScalarField, M_DATAFIELD));
    built = (built && createPortOut<livector1D, IOCloudPoints>(this, &mimmo::IOCloudPoints::getLabels, M_VECTORLI));

    m_arePortsBuilt = built;
}

/*!
 * Return the points coordinates stored in the class
 * \return the points coordinates
 */
dvecarr3E
IOCloudPoints::getPoints(){
    return m_points;
};

/*!
 * Return the vector field stored in the class. Field size is always checked and automatically
 * fit to point list during the class execution
 * \return vector field stored in the class
 *
 */
dvecarr3E
IOCloudPoints::getVectorField(){
    return m_vectorfield;
};

/*!
 * Return the scalar field stored in the class. Field size is always checked and automatically
 * fit to point list during the class execution
 * \return scalar field stored in the class
 *
 */
dvector1D
IOCloudPoints::getScalarField(){
    return m_scalarfield;
};


/*!
 * Return the labels attached to displacements actually stored into the class.
 * Labels are always checked and automatically fit to point list during the class execution
 * \return labels of displacements
 */
livector1D
IOCloudPoints::getLabels(){
    return m_labels;
};

/*!
 * Return if template option is active. This method is only meant in class working in Write mode
 * \return template flag
 */
bool
IOCloudPoints::isTemplate(){
    return m_template;
};

/*!
 * It sets the name of the input directory. Active only in read mode.
 * \param[in] dir directory path
 */
void
IOCloudPoints::setReadDir(std::string dir){
    if(!m_read)    return;
    m_dir = dir;
};

/*!
 * It sets the name of the input file. Active only in read mode.
 * \param[in] filename filename with tag extension included.
 */
void
IOCloudPoints::setReadFilename(std::string filename){
    if(!m_read)    return;
    m_filename = filename;
};

/*!
 * It sets the name of the output directory. Active only in write mode.
 * \param[in] dir directory path
 */
void
IOCloudPoints::setWriteDir(std::string dir){
    if(m_read)    return;
    m_dir = dir;
};

/*!
 * It sets the name of the output file. Active only in write mode.
 * \param[in] filename filename with tag extension included.
 */
void
IOCloudPoints::setWriteFilename(std::string filename){
    if(m_read)    return;
    m_filename = filename;
};

/*!
 * Set the point list into the class
 * The method is not active in Read mode.
 * \param[in] points list of 3D point
 */
void
IOCloudPoints::setPoints(dvecarr3E points){
    if(m_read) return;
    m_points = points;
};

/*!
 * It sets the labels attached to each point.
 * The method is not active in Read mode.
 * \param[in] labels list of label ids
 */
void
IOCloudPoints::setLabels(livector1D labels){
    if(m_read) return;
    m_labels = labels;
};

/*!
 * It sets the vector field associated to point.
 * The method is not active in Read mode.
 * \param[in] vectorfield vector field
 */
void
IOCloudPoints::setVectorField(dvecarr3E vectorfield){
    if(m_read) return;
    m_vectorfield = vectorfield;
};

/*!
 * It sets the scalar field associated to point.
 * The method is not active in Read mode.
 * \param[in] scalarfield scalar field
 */
void
IOCloudPoints::setScalarField(dvector1D scalarfield){
    if(m_read) return;
    m_scalarfield = scalarfield;
};

/*!
 * Enables the template writing mode. The method is not active in Read mode.
 * \param[in] flag true to enable the template writing
 */
void
IOCloudPoints::setTemplate(bool flag){
    if(m_read) return;
    m_template = flag;
};

/*!
 * Clear all data stored in the class
 */
void
IOCloudPoints::clear(){
    m_vectorfield.clear();
    m_scalarfield.clear();
    m_points.clear();
    m_labels.clear();
    m_template     = false;
    m_filename     = m_name+"_source.dat";
}

/*!
 * Execution command.
 * Read data from or Write data on linked filename
 */
void
IOCloudPoints::execute(){

    if(m_read){
        read();
    }else{
        write();
    }
};

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
IOCloudPoints::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    std::string input;

    //checking IOmode

    if(slotXML.hasOption("IOmode")){
        input = slotXML.get("IOmode");
        bool value = true;
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::string::trim(input));
            ss>>value;
        }
        if (m_read != value){
            (*m_log)<< "warning in class "<<m_name<<": cannot absorb xml parameters for class IOmode mismatching"<<std::endl;
            throw std::runtime_error (m_name + " : xml absorbing failed ");
        }
    };

    BaseManipulation::absorbSectionXML(slotXML, name);

    if(m_read){
        if(slotXML.hasOption("ReadDir")){
            std::string input = slotXML.get("ReadDir");
            input = bitpit::utils::string::trim(input);
            setReadDir(input);
        };

        if(slotXML.hasOption("ReadFilename")){
            std::string input = slotXML.get("ReadFilename");
            input = bitpit::utils::string::trim(input);
            setReadFilename(input);
        };
    }else{
        if(slotXML.hasOption("WriteDir")){
            std::string input = slotXML.get("WriteDir");
            input = bitpit::utils::string::trim(input);
            setWriteDir(input);
        };

        if(slotXML.hasOption("WriteFilename")){
            std::string input = slotXML.get("WriteFilename");
            input = bitpit::utils::string::trim(input);
            setWriteFilename(input);
        };
    }


    if(slotXML.hasOption("Template")){
        std::string input = slotXML.get("Template");
        input = bitpit::utils::string::trim(input);
        bool temp = false;
        if(!input.empty()){
            std::stringstream ss(input);
            ss>>temp;
        }
        setTemplate(temp);
    };

}

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
IOCloudPoints::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    BaseManipulation::flushSectionXML(slotXML, name);

    slotXML.set("IOmode", std::to_string(int(m_read)));

    if(m_read){
        slotXML.set("ReadDir", m_dir);
        slotXML.set("ReadFilename", m_filename);
    }else{
        slotXML.set("WriteDir", m_dir);
        slotXML.set("WriteFilename", m_filename);
    }

    slotXML.set("Template", std::to_string(int(m_template)));
};


/*!
 * Plot cloud point and store it in *.vtu file
 */
void
IOCloudPoints::plotOptionalResults(){

    std::string path = m_outputPlot;
    std::string name = m_name +".Cloud_" + std::to_string(getId());
    bitpit::VTKUnstructuredGrid output( path, name, bitpit::VTKElementType::VERTEX);

    int size = m_points.size();
    ivector1D conn(size);
    for(int i=0; i<size; ++i)    conn[i] = i;
    output.setGeomData(bitpit::VTKUnstructuredField::POINTS, m_points);
    output.setGeomData(bitpit::VTKUnstructuredField::CONNECTIVITY, conn);
    output.setDimensions(size, size);

    std::string sfield = "scalarfield";
    dvector1D scafield = m_scalarfield;
    scafield.resize(size, 0.0);

    std::string vfield = "vectorfield";
    dvecarr3E vecfield = m_vectorfield;
    vecfield.resize(size, {{0.0,0.0,0.0}});

    output.addData( sfield, bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, scafield ) ;
    output.addData( vfield, bitpit::VTKFieldType::VECTOR, bitpit::VTKLocation::POINT, vecfield ) ;
    output.write() ;
}

/*!
 * Read displacement data from file
 */
void
IOCloudPoints::read(){

    std::unordered_map<long, int> mapP;
    std::ifstream reading;
    std::string source = m_dir+"/"+m_filename;
    reading.open(source.c_str());
    if(reading.is_open()){

        m_points.clear();
        m_labels.clear();
        m_scalarfield.clear();
        m_vectorfield.clear();

        std::string keyword, line;
        long label;
        darray3E dtrial;
        double temp;

        do{
            line.clear();
            keyword.clear();
            std::getline(reading,line);
            line = bitpit::utils::string::trim(line);
            std::stringstream ss(line);
            ss>>keyword;
            keyword = bitpit::utils::string::trim(keyword);
            if(keyword == "$POINT")    {

                ss>>label;
                m_labels.push_back(label);

                dtrial.fill(0.0);
                ss>>dtrial[0]>>dtrial[1]>>dtrial[2];
                m_points.push_back(dtrial);
            }

        }while(!reading.eof());

        int counter = 0;
        for(auto &lab :m_labels){
            mapP[lab] = counter;
            ++counter;
        }

        m_scalarfield.resize(m_points.size(),0.0);
        m_vectorfield.resize(m_points.size(),{{0.0,0.0,0.0}});

        reading.clear();
        reading.seekg(0, std::ios::beg);

        do{
            line.clear();
            keyword.clear();
            std::getline(reading,line);
            line = bitpit::utils::string::trim(line);
            std::stringstream ss(line);
            ss>>keyword;
            keyword = bitpit::utils::string::trim(keyword);
            if(keyword == "$SCALARF")    {

                ss>>label>>temp;
                m_scalarfield[mapP[label]] = temp;
            }
            if(keyword == "$VECTORF")    {

                dtrial.fill(0.0);
                ss>>label>>dtrial[0]>>dtrial[1]>>dtrial[2];
                m_vectorfield[mapP[label]] = dtrial;
            }

        }while(!reading.eof());

    }else{
        (*m_log)<<"error of "<<m_name<<" : cannot open "<<m_filename<< " requested. Exiting... "<<std::endl;
        throw std::runtime_error (m_name + " : cannot open " + m_filename + " requested. Exiting... ");
    }

    reading.close();
};

/*!
 * Write displacement data to file
 */
void
IOCloudPoints::write(){

    //assessing data;
    int sizelabels = m_labels.size();
    long maxlabel = 0;
    for(auto & val: m_labels){
        maxlabel = std::max(maxlabel,val);
    }

    int sizedispl = m_points.size();

    m_labels.resize(m_points.size(), -1);
    if(sizelabels < sizedispl){
        for(int i = sizelabels; i<sizedispl; ++i){
            m_labels[i] = maxlabel+1;
            ++maxlabel;
        }
    }

    m_scalarfield.resize(m_points.size(), 0.0);
    m_vectorfield.resize(m_points.size(), {{0.0,0.0,0.0}});

    std::string filename;
    if(m_template){
        filename = "TEMPLATE_"+m_filename;
    }else{
        filename = m_filename;
    }

    std::ofstream writing;
    std::string source = m_dir+"/"+m_filename;
    writing.open(source.c_str());
    std::string keyT1 = "{", keyT2 = "}";
    if(writing.is_open()){

        int counter = 0;
        for(auto & dd : m_points){
            writing<<"$POINT"<<'\t'<<m_labels[counter]<<'\t'<<dd[0]<<'\t'<<dd[1]<<'\t'<<dd[2]<<std::endl;
            ++counter;
        }
        writing<<""<<std::endl;

        counter = 0;
        for(auto & dd : m_scalarfield){
            if(m_template){
                std::string str1 = keyT1+"s"+std::to_string(m_labels[counter])+keyT2;
                writing<<"$SCALARF"<<'\t'<<m_labels[counter]<<'\t'<<str1<<std::endl;
            }else{
                writing<<"$SCALARF"<<'\t'<<m_labels[counter]<<'\t'<<dd<<std::endl;
            }
            ++counter;
        }
        writing<<""<<std::endl;

        counter = 0;
        for(auto & dd : m_vectorfield){
            if(m_template){
                std::string str1 = keyT1+"x"+std::to_string(m_labels[counter])+keyT2;
                std::string str2 = keyT1+"y"+std::to_string(m_labels[counter])+keyT2;
                std::string str3 = keyT1+"z"+std::to_string(m_labels[counter])+keyT2;

                writing<<"$VECTORF"<<'\t'<<m_labels[counter]<<'\t'<<str1<<'\t'<<str2<<'\t'<<str3<<std::endl;
            }else{
                writing<<"$VECTORF"<<'\t'<<m_labels[counter]<<'\t'<<dd[0]<<'\t'<<dd[1]<<'\t'<<dd[2]<<std::endl;
            }
            ++counter;
        }
        writing<<""<<std::endl;
    }else{
        (*m_log)<<"error of "<<m_name<<" : cannot open "<<m_filename<< " requested. Exiting... "<<std::endl;
        throw std::runtime_error (m_name + " : cannot open " + m_filename + " requested. Exiting... ");
    }

    writing.close();
};


}
