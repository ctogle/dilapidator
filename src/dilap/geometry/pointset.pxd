









# dilapidators implementation of a pointset (for quats,vec3s,vec2s,etc.)
cdef class pointset:

    cdef public list ps
    cdef public int pcnt

    cdef list gpscp_c(self)
    cdef list gps_c(self)
    cdef int ap_c(self,np)
    cdef list aps_c(self,list nps)
    cdef int np_c(self,np)
    cdef list nps_c(self,list nps)
    cdef int fp_c(self,p)
    cdef list fps_c(self,list ps)
    cdef bint disjoint_c(self,pointset o)

    cpdef list gpscp(self)
    cpdef list gps(self)
    cpdef int ap(self,np)
    cpdef list aps(self,list nps)
    cpdef int np(self,np)
    cpdef list nps(self,list nps)
    cpdef int fp(self,p)
    cpdef list fps(self,list ps)
    cpdef bint disjoint(self,pointset o)
        

 







