//Hyrodium's POV-Ray format
//Hy_transformations.inc

//#include "Hy_transformations.inc"
#macro Matrix(e1,e2,e3,p0)
	matrix<
	vdot(e1,x), vdot(e1,y), vdot(e1,z),
	vdot(e2,x), vdot(e2,y), vdot(e2,z),
	vdot(e3,x), vdot(e3,y), vdot(e3,z),
	vdot(p0,x), vdot(p0,y), vdot(p0,z)>
#end

#macro RotateMatrix(n,w,p0)
	translate -p0
	Matrix(RotateVec(x,n,w),RotateVec(y,n,w),RotateVec(z,n,w),p0)
#end

#macro RotateMatrix2(n,w,p0)
	translate -p0
	Matrix(vaxis_rotate(x,n,degrees(w)),vaxis_rotate(y,n,degrees(w)),vaxis_rotate(z,n,degrees(w)),p0)
#end

#macro Front(level)
	#if(Pers)
		translate -Loc
		scale 1/level
		translate Loc
	#else
		translate Z*level
	#end
#end
