#include <fstream>
#include "Operators.hpp"

/*!Overloaded function of base class setInput.
 * It sets the input/result and write on file at the same time.
 * \param[in] data Pointer to data to be written and to be used to set the input/result.
 */
template<typename T>
void
mimmo::GenericOutput::setInput(T* data){
	BaseManipulation::setInput(data);
	BaseManipulation::setResult(data);
	std::ofstream file;
	file.open(m_filename);
	if (file.is_open()){
		file << *data;
		file.close();
	}
}

/*!Overloaded function of base class setInput.
 * It sets the input/result and write on file at the same time.
 * \param[in] data Data to be written and to be used to set the input/result.
 */
template<typename T>
void
mimmo::GenericOutput::setInput(T& data){
	BaseManipulation::setInput(data);
	BaseManipulation::setResult(data);
	std::ofstream file;
	file.open(m_filename);
	if (file.is_open()){
		file << data;
		file.close();
	}
}


