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
#include "Apply.hpp"


namespace mimmo{

/*!
 * Default constructor of Apply
 */
Apply::Apply():BaseManipulation(){
    m_name = "mimmo.Apply";
    m_input.clear();
    m_scalarinput.clear();
    m_factor = 1.;
    m_annotation = false;
    m_annotationThres = 1.0E-18;
    m_annCellLabel = "DeformedCells";
    m_annVertexLabel = "DeformedVertices";
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
Apply::Apply(const bitpit::Config::Section & rootXML){

	m_name = "mimmo.Apply";
	m_factor = 1.;
    m_annotation = false;
    m_annotationThres = 1.0E-18;
    m_annCellLabel = "DeformedCells";
    m_annVertexLabel = "DeformedVertices";

	std::string fallback_name = "ClassNONE";
	std::string input = rootXML.get("ClassName", fallback_name);
	input = bitpit::utils::string::trim(input);
	if(input == "mimmo.Apply"){
		absorbSectionXML(rootXML);
	}else{
		warningXML(m_log, m_name);
	};
}

/*!Default destructor of Apply
 */
Apply::~Apply(){};

/*!
    Copy constructor of Apply.
 */
Apply::Apply(const Apply & other):BaseManipulation(other){
	m_input = other.m_input;
	m_scalarinput = other.m_scalarinput;
	m_factor = other.m_factor;
    m_annotation = other.m_annotation;
    m_annotationThres = other.m_annotationThres;
    m_annCellLabel = other.m_annCellLabel;
    m_annVertexLabel = other.m_annVertexLabel;
    //annotation structures are not copied, they are filled each time in execution if m_annotation is true.
};

/*!
 * Swap function
 * \param[in] x object to be swapped
 */
void Apply::swap(Apply & x) noexcept
{
	//std::swap(m_input, x.m_input);
	m_input.swap(x.m_input);
	m_scalarinput.swap(x.m_scalarinput);
	std::swap(m_factor, x.m_factor);
    std::swap(m_annotation, x.m_annotation);
    std::swap(m_annotationThres, x.m_annotationThres);
    std::swap(m_annCellLabel, x.m_annCellLabel);
    std::swap(m_annVertexLabel, x.m_annVertexLabel);

    m_cellAnnotation.swap(x.m_cellAnnotation);
    m_vertexAnnotation.swap(x.m_vertexAnnotation);

	BaseManipulation::swap(x);
}

/*! It builds the input/output ports of the object
 */
void
Apply::buildPorts(){
	bool built = true;
	built = (built && createPortIn<dmpvecarr3E*, Apply>(this, &Apply::setInput, M_GDISPLS, true, 1));
    built = (built && createPortIn<dmpvector1D*, Apply>(this, &Apply::setScalarInput, M_SCALARFIELD, true, 1));
    built = (built && createPortIn<dmpvector1D*, Apply>(this, &Apply::setFilter, M_FILTER));
	built = (built && createPortIn<MimmoSharedPointer<MimmoObject>, Apply>(this, &BaseManipulation::setGeometry, M_GEOM, true));

	built = (built && createPortOut<MimmoSharedPointer<MimmoObject>, Apply>(this, &BaseManipulation::getGeometry, M_GEOM));
	built = (built && createPortOut<dmpvecarr3E*, Apply>(this, &Apply::getOutput, M_GDISPLS));
    built = (built && createPortOut<MimmoPiercedVector<long>*, Apply>(this, &Apply::getAnnotatedVertices, M_LONGFIELD));
    built = (built && createPortOut<MimmoPiercedVector<long>*, Apply>(this, &Apply::getAnnotatedCells, M_LONGFIELD2));

	m_arePortsBuilt = built;
};


/*!It sets the displacements input.
 * \param[in] input Input displacements of the geometry vertices.
 */
void
Apply::setInput(dmpvecarr3E *input){
    if(!input)  return;
	m_input = *input;
};

/*!It sets the displacements given as scalar input. It makes sense only for surface geometries.
 * The displacements will be directed as the surface normals.
 * \param[in] input Input displacements of the geometry vertices given as module of normal directed vectors.
 */
void
Apply::setScalarInput(dmpvector1D *input){
    if(!input) return;
    m_scalarinput = *input;
};

/*!It sets the filter to be applied during the deformation of the geometry.
 * The deformation field will be multiplied by the scalar filter field passed as input.
 * \param[in] input Input filter field used during the deformation.
 */
void
Apply::setFilter(dmpvector1D *input){
    if(!input) return;
    m_filter = *input;
};

/*!It sets the displacements scalar factor.
 * The displacements will be scaled by using this factor.
 * \param[in] alpha Scale factor of displacements field.
 */
void
Apply::setScaling(double alpha){
	m_factor = alpha;
};

/*!
   If true, annotate vertices and cells involved into mesh
   deformation and store it into the class. Annotations are acessible through
   methods getAnnotatedVertices/getAnnotatedCells.
 * \param[in] activate  true/false to enable annotation creation.
 */
void
Apply::setAnnotation(bool activate){
    m_annotation = activate;
}

/*!
   Set threshold for deformation field norm, over which a cell/vertex interested by it is
   eligible for annotation.
 * \param[in] threshold value (double, >0.0)
 */
void
Apply::setAnnotationThreshold(double threshold){
    m_annotationThres = std::max(1.0E-18, threshold);
}

/*!
   Provide the string label to uniquely name the deformed cells annotation.
   Class default label is DeformedCells.
 * \param[in] label cells annotation name.
 */
void
Apply::setCellsAnnotationName(const std::string & label){
    if(label.empty()){
        m_annCellLabel = "DeformedCells";
    }else{
        m_annCellLabel = label;
    }
}

/*!
   Provide the string label to uniquely name the deformed vertices annotation.
   Class default label is DeformedVertices.
 * \param[in] label vertices annotation name.
 */
void
Apply::setVerticesAnnotationName(const std::string & label){
    if(label.empty()){
        m_annVertexLabel = "DeformedVertices";
    }else{
        m_annVertexLabel = label;
    }

}

/*!It get the displacements output.
 * \return Output displacements of the geometry vertices.
 */
dmpvecarr3E*
Apply::getOutput(){
    return	&m_output;
};

/*!
    Return annotation structure for deformed vertices.
    The structure is filled during class execution and if and only Annotation
    is active (see setAnnotation method).
    \return vertex annotation structure
*/
MimmoPiercedVector<long> * Apply::getAnnotatedVertices(){
    return &m_vertexAnnotation;
}
/*!
    Return annotation structure for deformed cells.
    The structure is filled during class execution and if and only Annotation
    is active (see setAnnotation method).
    \return cell annotation structure
*/
MimmoPiercedVector<long> * Apply::getAnnotatedCells(){
    return &m_cellAnnotation;
}

/*!Execution command.
 * It applies the deformation stored in the input of base class (casting the input
 * for apply object to dvecarr3E) to the linked geometry.
 * After exec() the original geometry will be permanently modified.
 */
void
Apply::execute(){
    if(!getGeometry()){
        (*m_log)<<m_name + " : nullptr pointer to linked geometry found"<<std::endl;
        return;
    }

    if(getGeometry()->isEmpty()){
#if MIMMO_ENABLE_MPI
        (*m_log)<<"WARNING "<<m_name <<"on rank "<<m_rank<<" : empty linked geometry found"<<std::endl;
#else
        (*m_log)<<"WARNING "<<m_name <<" : empty linked geometry found"<<std::endl;
#endif
    }

	checkInput();

	m_output = m_input;

	darray3E vertexcoords;
	long int ID;
	for (const auto & vertex : m_geometry->getVertices()){
		vertexcoords = vertex.getCoords();
		ID = vertex.getId();
		std::array<double,3> val = m_factor*m_input[ID];
		m_output[ID] = val;
		vertexcoords += val;
		getGeometry()->modifyVertex(vertexcoords, ID);
	}

    //step 2: produce annotations.
    if(m_annotation){
        //-->get the list of vertices involved, whose deformation norm is greater
        //   then prescribed tolerance;
        double tol = m_annotationThres;
        m_vertexAnnotation.clear();
        m_vertexAnnotation.setGeometry(getGeometry());
        m_vertexAnnotation.setDataLocation(MPVLocation::POINT);
        m_vertexAnnotation.setName(m_annVertexLabel);
        m_vertexAnnotation.reserve(m_output.size());

        for(auto it = m_output.begin(); it!= m_output.end(); ++it){
            if(norm2(*it) < tol) continue;
            m_vertexAnnotation.insert(it.getId(),it.getId());
        }

        m_vertexAnnotation.shrinkToFit();

        //cell round
        livector1D celllist = getGeometry()->getCellFromVertexList(m_vertexAnnotation.getIds(), false); //false! any cell sharing a moved vertex
        m_cellAnnotation.clear();
        m_cellAnnotation.setGeometry(getGeometry());
        m_cellAnnotation.setDataLocation(MPVLocation::CELL);
        m_cellAnnotation.setName(m_annCellLabel);
        m_cellAnnotation.reserve(celllist.size());

        for (long id: celllist){
            m_cellAnnotation.insert(id,id);
        }
        //no need to shrink cellAnnotation.
    }

    // Update geometry
    getGeometry()->update();

};


/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
Apply::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

	BITPIT_UNUSED(name);
	BaseManipulation::absorbSectionXML(slotXML, name);

    if(slotXML.hasOption("Scaling")){
        std::string input = slotXML.get("Scaling");
        input = bitpit::utils::string::trim(input);
        double val = 1.;
        if(!input.empty()){
            std::stringstream ss(input);
            ss>>val;
        }
        setScaling(val);
    }

    if(slotXML.hasOption("Annotation")){
        std::string input = slotXML.get("Annotation");
        input = bitpit::utils::string::trim(input);
        bool val = false;
        if(!input.empty()){
            std::stringstream ss(input);
            ss>>val;
        }
        setAnnotation(val);
    }

    if(slotXML.hasOption("AnnotationThreshold")){
        std::string input = slotXML.get("AnnotationThreshold");
        input = bitpit::utils::string::trim(input);
        double val = 1.0E-18;
        if(!input.empty()){
            std::stringstream ss(input);
            ss>>val;
        }
        setAnnotationThreshold(val);
    }

    if(slotXML.hasOption("CellsAnnotationName")){
        std::string input = slotXML.get("CellsAnnotationName");
        input = bitpit::utils::string::trim(input);
        setCellsAnnotationName(input);
    }

    if(slotXML.hasOption("VerticesAnnotationName")){
        std::string input = slotXML.get("VerticesAnnotationName");
        input = bitpit::utils::string::trim(input);
        setVerticesAnnotationName(input);
    }


};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
Apply::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

	BITPIT_UNUSED(name);
	BaseManipulation::flushSectionXML(slotXML, name);

    slotXML.set("Scaling", std::to_string(m_factor));
    slotXML.set("Annotation", std::to_string(m_annotation));
    slotXML.set("AnnotationThreshold", std::to_string(m_annotationThres));
    slotXML.set("CellsAnnotationName", m_annCellLabel);
    slotXML.set("VerticesAnnotationName", m_annVertexLabel);
};

/*!
 * Check if the input is related to the target geometry.
 * If not create a zero input field.
 */
void
Apply::checkInput(){

	bool check = false;
	if (m_input.size()){
		check = m_input.getDataLocation() == mimmo::MPVLocation::POINT;
		check = check && m_input.getGeometry() == m_geometry;
		check = check && m_input.completeMissingData({{0.0,0.0,0.0}});
	}
	if (m_scalarinput.size()){
		check = m_scalarinput.getDataLocation() == mimmo::MPVLocation::POINT;
		check = check && m_geometry->getType() == 1;
		check = check && m_scalarinput.getGeometry() == m_geometry;
		check = check && m_scalarinput.completeMissingData(0.0);
		if (check){
			bitpit::SurfaceKernel* skernel = static_cast<bitpit::SurfaceKernel*>(m_geometry->getPatch());
			bitpit::PiercedVector<darray3E> vNormals;
			bitpit::ConstProxyVector<long> verts;
			std::size_t size;
			long idN;
			for(const auto & cell: m_geometry->getCells()){
				verts = cell.getVertexIds();
				size = verts.size();
				for(std::size_t i=0; i<size; ++i){
					idN = verts[i];
					if(!vNormals.exists(idN)){
						vNormals.insert(idN, skernel->evalVertexNormal(cell.getId(), i));
					}
				}
			}
			m_input.clear();
			for(const auto & vertex: m_geometry->getVertices()){
				m_input.setDataLocation(mimmo::MPVLocation::POINT);
				m_input.setGeometry(m_geometry);
				long id = vertex.getId();
				m_input.insert(id, darray3E(m_scalarinput[id]*vNormals[id]));
			}
		}
	}
	if (m_filter.size()){
	    check = m_filter.getDataLocation() == mimmo::MPVLocation::POINT;
	    check = check && m_geometry->getType() == 1;
        // Force filter to link to the morphing geometry
	    bool geometry_check = m_filter.getGeometry() == m_geometry;
	    if (!geometry_check){
	        m_filter.setGeometry(m_geometry);
	    }
	    // Complete missing data force filter to 0.
	    check = check && m_filter.completeMissingData(0.0);
	    if (check){
	        // Apply filter to deformation field
	        for (const auto & vertex : m_geometry->getVertices()){
	            long ID = vertex.getId();
	            m_input[ID] *= m_filter[ID];
	        }
	    }
	}
	if (!check){
		m_log->setPriority(bitpit::log::Verbosity::DEBUG);
		(*m_log)<<"Not valid input found in "<<m_name<<". Proceeding with default zero field"<<std::endl;
		m_log->setPriority(bitpit::log::Verbosity::NORMAL);
		m_input.clear();
		m_input.setGeometry(m_geometry);
		m_input.setDataLocation(mimmo::MPVLocation::POINT);
		m_input.reserve(m_geometry->getNVertices());
		for (const auto & vertex : m_geometry->getVertices()){
			m_input.insert(vertex.getId(), {{0.0,0.0,0.0}});
		}
	}
};

}
