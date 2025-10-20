:  Vector stream of events

NEURON {
	ARTIFICIAL_CELL VecStim
	POINTER ptr
}

ASSIGNED {
	index
	etime (ms)
	ptr
}


INITIAL {
	index = 0
	element()
	if (index > 0) {
		net_send(etime - t, 1)
	}
}

NET_RECEIVE (w) {
	if (flag == 1) {
		net_event(t)
		element()
		if (index > 0) {
			net_send(etime - t, 1)
		}
	}
}

DESTRUCTOR {
VERBATIM
	IvocVect* vv = (IvocVect*)(_p_ptr);
	if (vv) {
		hoc_obj_unref(*vector_pobj(vv));
	}
ENDVERBATIM
}

PROCEDURE element() {
VERBATIM	
  { void* vv; int i, size; double* px;
	i = (int)index;
	if (i >= 0) {
		vv = (void*)(_p_ptr);
		if (vv) {
			size = vector_capacity((IvocVect*)vv);
			px = vector_vec((IvocVect*)vv);
			if (i < size) {
				etime = px[i];
				index += 1.;
			}else{
				index = -1.;
			}
		}else{
			index = -1.;
		}
	}
  }
ENDVERBATIM
}

PROCEDURE play() {
VERBATIM
	void** pv;
	IvocVect* ptmp = NULL;
	if (ifarg(1)) {
		ptmp = vector_arg(1);
		hoc_obj_ref(*vector_pobj(ptmp));
	}
	pv = (void**)(&_p_ptr);
	if (*pv) {
		hoc_obj_unref(*vector_pobj((IvocVect*)*pv));
	}
	*pv = ptmp;
ENDVERBATIM
}
