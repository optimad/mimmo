template<typename T, typename C, typename O>
int
mimmo::Chain::addDofInput(O* obj, T (C::*fget) ()){
	DofOutT<T,C,O>* dofin = new DofOutT<T,C,O>(obj, fget);
	m_dof.push_back(dofin);
	return (m_dof.size()-1);
}

template<typename T, typename C, typename O>
int
mimmo::Chain::addDofInput(O* obj, T* (C::*fget) ()){
	DofOutT<T,C,O>* out = new DofOutT<T,C,O>(obj, fget);
	m_dof.push_back(out);
	return (m_dof.size()-1);
}

template<typename T, typename C, typename O>
int
mimmo::Chain::addDofInput(O* obj, T (C::*fget) (), int nglob, int nuse, std::vector<bool> map){
	DofOutT<T,C,O>* dofin = new DofOutT<T,C,O>(obj, fget);
	m_dof.push_back(dofin);
	setNdofs(m_dof.size()-1, nglob, nuse, map);
	return (m_dof.size()-1);
}

template<typename T, typename C, typename O>
int
mimmo::Chain::addDofInput(O* obj, T* (C::*fget) (), int nglob, int nuse, std::vector<bool> map){
	DofOutT<T,C,O>* dofin = new DofOutT<T,C,O>(obj, fget);
	m_dof.push_back(dofin);
	setNdofs(m_dof.size()-1, nglob, nuse, map);
	return (m_dof.size()-1);
}


template<typename T, typename C, typename O>
int
mimmo::Chain::addDofInput(O* obj, T (C::*fget) (), std::vector<bool> map){
	DofOutT<T,C,O>* dofin = new DofOutT<T,C,O>(obj, fget);
	m_dof.push_back(dofin);
	int nglob = map.size();
	int nuse = std::accumulate(map.begin(), map.end(), 0);
	setNdofs(m_dof.size()-1, nglob, nuse, map);
	return (m_dof.size()-1);
}

template<typename T, typename C, typename O>
int
mimmo::Chain::addDofInput(O* obj, T* (C::*fget) (), std::vector<bool> map){
	DofOutT<T,C,O>* dofin = new DofOutT<T,C,O>(obj, fget);
	m_dof.push_back(dofin);
	int nglob = map.size();
	int nuse = std::accumulate(map.begin(), map.end(), 0);
	setNdofs(m_dof.size()-1, nglob, nuse, map);
	return (m_dof.size()-1);
}


template<typename T, typename C, typename O>
void
mimmo::Chain::activateDof(O* obj, T* (C::*fget) (), int j){
	if (m_dof.size() == 0) return;
	int i = findDofInput(obj,fget);
	if (i >= m_dof.size()) return;
	if (j>m_dof[i]->getNgdof()) return;
	m_dof[i]->m_actives[j] = true;
};

template<typename T, typename C, typename O>
void
mimmo::Chain::activateDofs(O* obj, T* (C::*fget) ()){
	if (m_dof.size() == 0) return;
	int i = findDofInput(obj,fget);
	if (i >= m_dof.size()) return;
	m_dof[i]->m_actives = std::vector<bool>(m_dof[i]->getNgdof(), true);
};


template<typename T, typename C, typename O>
void
mimmo::Chain::activateDof(O* obj, T (C::*fget) (), int j){
	if (m_dof.size() == 0) return;
	int i = findDofInput(obj,fget);
	if (i >= m_dof.size()) return;
	if (j>m_dof[i]->getNgdof()) return;
	m_dof[i]->m_actives[j] = true;
};

template<typename T, typename C, typename O>
void
mimmo::Chain::activateDofs(O* obj, T (C::*fget) ()){
	if (m_dof.size() == 0) return;
	int i = findDofInput(obj,fget);
	if (i >= m_dof.size()) return;
	m_dof[i]->m_actives = std::vector<bool>(m_dof[i]->getNgdof(), true);
};


template<typename T, typename C, typename O>
void
mimmo::Chain::disableDof(O* obj, T* (C::*fget) (), int j){
	if (m_dof.size() == 0) return;
	int i = findDofInput(obj,fget);
	if (i >= m_dof.size()) return;
	if (j>m_dof[i]->getNgdof()) return;
	m_dof[i]->m_actives[j] = false;
};

template<typename T, typename C, typename O>
void
mimmo::Chain::disableDofs(O* obj, T* (C::*fget) ()){
	if (m_dof.size() == 0) return;
	int i = findDofInput(obj,fget);
	if (i >= m_dof.size()) return;
	m_dof[i]->m_actives = std::vector<bool>(m_dof[i]->getNgdof(), false);
};


template<typename T, typename C, typename O>
void
mimmo::Chain::disableDof(O* obj, T (C::*fget) (), int j){
	if (m_dof.size() == 0) return;
	int i = findDofInput(obj,fget);
	if (i >= m_dof.size()) return;
	if (j>m_dof[i]->getNgdof()) return;
	m_dof[i]->m_actives[j] = false;
};

template<typename T, typename C, typename O>
void
mimmo::Chain::disableDofs(O* obj, T (C::*fget) ()){
	if (m_dof.size() == 0) return;
	int i = findDofInput(obj,fget);
	if (i >= m_dof.size()) return;
	m_dof[i]->m_actives = std::vector<bool>(m_dof[i]->getNgdof(), false);
};





template<typename T, typename C, typename O>
int
mimmo::Chain::getNdofs(O* obj, T* (C::*fget) ()){
	int i = findDofInput(obj, fget ());
	return getNdofs(i);
}

template<typename T, typename C, typename O>
dvector1D
mimmo::Chain::getDofs(O* obj, T* (C::*fget) ()){
	int i = findDofInput(obj, fget ());
	return getDofs(i);
}

template<typename T, typename C, typename O>
double
mimmo::Chain::getDof(O* obj, T* (C::*fget) (), int j){
	int i = findDofInput(obj, fget ());
	return getDof(i,j);
}


template<typename T, typename C, typename O>
int
mimmo::Chain::getNdofs(O* obj, T (C::*fget) ()){
	int i = findDofInput(obj, fget ());
	return getNdofs(i);
}

template<typename T, typename C, typename O>
dvector1D
mimmo::Chain::getDofs(O* obj, T (C::*fget) ()){
	int i = findDofInput(obj, fget ());
	return getDofs(i);
}

template<typename T, typename C, typename O>
double
mimmo::Chain::getDof(O* obj, T (C::*fget) (), int j){
	int i = findDofInput(obj, fget ());
	return getDof(i,j);
}




template<typename T, typename C, typename O>
int
mimmo::Chain::findDofInput(O* obj, T* (C::*fget) ()){
	DofOutT<T,C,O>* dofin = new DofOutT<T,C,O>(obj, fget);
	std::vector<DofOut*>::iterator it = findDofInput(m_dof.begin(), m_dof.end(), dofin);
	return(distance(m_dof.begin(), it));
}


template<typename T, typename C, typename O>
int
mimmo::Chain::findDofInput(O* obj, T (C::*fget) ()){
	DofOutT<T,C,O>* dofin = new DofOutT<T,C,O>(obj, fget);
	std::vector<DofOut*>::iterator it = findDofInput(m_dof.begin(), m_dof.end(), dofin);
	return(distance(m_dof.begin(), it));
}



