//Hyrodium's POV-Ray format
//Hy_functions.inc

//#include "Hy_functions.inc"
//lim

#macro Sign(a)
	#if(a>0)
		#local b=1;
	#elseif(a<0)
		#local b=-1;
	#else
		#local b=0;
	#end
	b
#end


#macro SGP(a,w)
	(2*a*<1,1,cos(w)>+z*(vdot(a,a)-1)*sin(w))/((vdot(a,a)-1)*(1-cos(w))+2*(vdot(a,z)*sin(w)+1))
#end

#macro SCS(lng,lat)
	<cos(radians(lat))*cos(radians(lng)),cos(radians(lat))*sin(radians(lng)),sin(radians(lat))>
#end

#macro Area3(a,b,c)
	vlength((vcross(a,b)+vcross(b,c)+vcross(c,a))/2)
#end

#macro Deg3(a,b,c)
	vdot(vcross(a,b)+vcross(b,c)+vcross(c,a),vcross(a,b)+vcross(b,c)+vcross(c,a))*(vdot(a-b,a-b)+vdot(b-c,b-c)+vdot(c-a,c-a))/(vdot(a-b,a-b)*vdot(b-c,b-c)*vdot(c-a,c-a))
#end

#macro DEG3(a,b,c)
	#if(Deg3(a,b,c)>0.01)
		#local b=0;
	#else
		#if(vdot(c-b,b-a)>0)
			#local b=1;
		#elseif(vdot(c-b,a-c)<0)
			#local b=2;
		#else
			#local b=3;
		#end
	#end
	b
#end

#macro CircumCenter(a,b,c)
	#local A=vdot(b-c,b-c)*(vdot(c-a,c-a)+vdot(a-b,a-b)-vdot(b-c,b-c));
	#local B=vdot(c-a,c-a)*(vdot(a-b,a-b)+vdot(b-c,b-c)-vdot(c-a,c-a));
	#local C=vdot(a-b,a-b)*(vdot(b-c,b-c)+vdot(c-a,c-a)-vdot(a-b,a-b));
	(A*a+B*b+C*c)/(A+B+C)
#end

#macro CircumRadius(a,b,c)
	vlength(a-CircumCenter(a,b,c))
#end

#macro CircumRadius2(a,b,c)
	sqrt(vdot(a-b,a-b)*vdot(b-c,b-c)*vdot(c-a,c-a))/(4*Area3(a,b,c))
#end

#macro CircumCenter4(a,b,c,d)
	#local a1=vdot(a,x);
	#local b1=vdot(b,x);
	#local c1=vdot(c,x);
	#local d1=vdot(d,x);
	#local a2=vdot(a,y);
	#local b2=vdot(b,y);
	#local c2=vdot(c,y);
	#local d2=vdot(d,y);
	#local a3=vdot(a,z);
	#local b3=vdot(b,z);
	#local c3=vdot(c,z);
	#local d3=vdot(d,z);
	#local A=vdot(a,a);
	#local B=vdot(b,b);
	#local C=vdot(c,c);
	#local D=vdot(d,d);
	#local K=<B-A,C-A,D-A>;
	#local k1=<b1-a1,c1-a1,d1-a1>;
	#local k2=<b2-a2,c2-a2,d2-a2>;
	#local k3=<b3-a3,c3-a3,d3-a3>;
	<Det(K,k2,k3),Det(k1,K,k3),Det(k1,k2,K)>/(Det(b,c,d)-Det(a,c,d)+Det(a,b,d)-Det(a,b,c))/2
#end

#macro CircumRadius4(a,b,c,d)
	vlength(a-CircumCenter4(a,b,c,d))
#end

#macro Normal(a,b,c)
	vnormalize(vcross(a,b)+vcross(b,c)+vcross(c,a))
#end

#macro Det(a,b,c)
	vdot(a,vcross(b,c))
#end

#macro Orth_vec(a)
	#local e1=vnormalize(a);
	#if(abs(vdot(e1,y))>0.9)
		#local e2=vnormalize(z-e1*vdot(e1,z));
	#else
		#local e2=vnormalize(y-e1*vdot(e1,y));
	#end
	e2
#end

#macro smootha(a)
	#if(a>0)
		#local b=exp(-1/a);
	#else
		#local b=0;
	#end
	b
#end

#macro usmooth(a,b,s)
	smootha((2/sqrt(3))*((s-a)/(b-a)))/(smootha((2/sqrt(3))*((s-a)/(b-a)))+smootha((2/sqrt(3))*(1-(s-a)/(b-a))))
#end

#macro RotateVec(w,n,s)
	#local nn=vnormalize(n);
	vdot(w,nn)*nn+(w-vdot(w,nn)*nn)*cos(s)-vcross(w,nn)*sin(s)
#end

#macro Exp(w)
	<exp(vdot(w,u))*cos(vdot(w,v)),exp(vdot(w,u))*sin(vdot(w,v))>
#end
//Inverse
