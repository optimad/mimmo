#include <fstream>
#include "Operators.hpp"

/*!It sets the input/result and write on file at the same time.
 */
template<typename T>
void
GenericOutput::setInput(T* data){
	BaseManipulation::setInput(data);
	BaseManipulation::setResult(data);
	std::ofstream file;
	file.open(m_filename);
	if (file.is_open()){
		file << *data;
		file.close();
	}
}

/*!It sets the input/result and write on file at the same time.
 */
template<typename T>
void
GenericOutput::setInput(T& data){
	BaseManipulation::setInput(data);
	BaseManipulation::setResult(data);
	std::ofstream file;
	file.open(m_filename);
	if (file.is_open()){
		file << data;
		file.close();
	}
}


