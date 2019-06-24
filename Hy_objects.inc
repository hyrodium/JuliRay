//Hyrodium's POV-Ray format
//Hy_objects.inc

//#include "Hy_objects.inc"
#declare Empty=difference{sphere{0,1} sphere{0,2}};

#macro Sphere(a,r)
	#if(r=0)
		Empty
	#else
		object{sphere{a,abs(r)}}
	#end
#end

#macro Cylinder(a,b,r)
	#if(r=0)
		Empty
	#elseif(vlength(a-b)=0)
		Empty
	#else
		object{cylinder{a,b,abs(r)}}
	#end
#end

#macro Cone(a,b,r)
	cone{a,abs(r),b,0}
#end

#macro Cone2(p1,r1,p2,r2)
	#local d=vlength(p2-p1);
	#local P1=p1-r1*(r2-r1)*vnormalize(p2-p1)/d;
	#local P2=p2-r2*(r2-r1)*vnormalize(p2-p1)/d;
	#local R1=sqrt(r1*r1-r1*r1*(r2-r1)*(r2-r1)/d/d);
	#local R2=sqrt(r2*r2-r2*r2*(r2-r1)*(r2-r1)/d/d);
	#if(r1*r2<0)
		union{
			cone{P1,R1,(R2*P1+R1*P2)/(R1+R2),0}
			cone{(R2*P1+R1*P2)/(R1+R2),0,P2,R2}
		}
	#else
		cone{P1,R1,P2,R2}
	#end
#end

#macro Arrow(a,b,r)
	union{
		difference{
			Cylinder(a,b,r)
			Cylinder(b-6*r*vnormalize(b-a),2*b-a,2*r)
		}
		Cone(b-6*r*vnormalize(b-a),b,2*r)
	}
#end

#macro Disc(a,n,r)
	disc{a,n,r}
#end

#macro Torus(a,n,r1,r2)
	#local e1=vnormalize(n);
	#local e2=Orth_vec(e1);
	#local e3=vcross(e1,e2);
	torus{r1,r2 Matrix(e3,e1,e2,a)}
#end

#macro Torus3(a,b,c,r)
	Torus(CircumCenter(a,b,c),Normal(a,b,c),CircumRadius(a,b,c),r)
#end

#macro Disc3(a,b,c)
	Disc(CircumCenter(a,b,c),Normal(a,b,c),CircumRadius(a,b,c))
#end

#macro Cylinder3(a,b,c,d)
	Cylinder(CircumCenter(a,b,c)+d*Normal(a,b,c),CircumCenter(a,b,c)-d*Normal(a,b,c),CircumRadius(a,b,c))
#end

#macro Blocks3(a,b,c,w)
	#local n=Normal(a,b,c);
	#if(Det(n,a-CircumCenter(a,b,c),c-CircumCenter(a,b,c))>0)
		union
	#else
		intersection
	#end
		{
			#local e1=a-CircumCenter(a,b,c);
			#local e3=vcross(e1,n);
			box{2*<1,1,1>,2*<-1,-1,0> Matrix(e1,w*n,e3,CircumCenter(a,b,c))}
			#local e1=c-CircumCenter(a,b,c);
			#local e3=vcross(e1,n);
			box{2*<1,1,-1>,2*<-1,-1,0> Matrix(e1,w*n,e3,CircumCenter(a,b,c))}
		}
#end

#macro Arc3(a,b,c,r)
	difference{
		Torus3(a,b,c,r)
		Blocks3(a,b,c,r)
	}
#end

#macro Pie3(a,b,c,d)
	difference{
		Cylinder3(a,b,c,d)
		Blocks3(a,b,c,d)
	}
#end

#macro Pac3(a,b,c)
	intersection{
		Disc3(a,b,c)
		Cylinder3(a,b,c,0.00001)
		Pie3(a,b,c,0.00001)
	}
#end

#macro Horn0(k,a,N)
	#local d=1/N;
	#local r=d/sqrt(1+k*k);
	#local i=1;
	union{
		#while(i<N+1)
			Sphere(<0,1-cos(i*d*a),-sin(i*d*a)>/a,r*i)
			Cone2(<0,1-cos((i-1)*d*a),-sin((i-1)*d*a)>/a,r*(i-1),<0,1-cos(i*d*a),-sin(i*d*a)>/a,r*i)
			#local i=i+1;
		#end
	}
#end

#macro Horn(k,a,N)
	#if(a)
		intersection{
			object{Horn0(k,a*1.5,N) scale 1.5}
			object{Pie3(<0,1/a-(1+2/k)/a,0>,<0,1/a-(1+2/k)*cos(a/2)/a,-(1+2/k)*sin(a/2)/a>,<0,1/a-(1+2/k)*cos(a)/a,-(1+2/k)*sin(a)/a>,2/k)}
		}
	#else
		object{Cone(-z,0,1/k)}
	#end
#end

#macro Arrow3(a,b,c,r)
	#local e1=Normal(a,b,c);
	#local e2=vnormalize(c-CircumCenter(a,b,c));
	#local e3=vcross(e1,e2);
	union{
		difference{
			object{Arc3(a,b,c,r)}
			intersection{
				object{Sphere(c,4.5*r)}
				cylinder{c+r*e3/10,c-5*r*e3,5*r}
			}
		}
		object{Horn(3,6*r/CircumRadius(a,b,c),10) scale 6*r Matrix(-e1,-e2,e3,c)}
	}
#end

#macro Tile(I,k)
	#local i=-floor(I);
	union{
		#while(i<floor(I))
			#local j=-floor(I);
			#while(j<floor(I))
				#if(Mod(i+j,2)=k)
					polygon{5,<i,j,0>,<i+1,j,0>,<i+1,j+1,0>,<i,j+1,0>,<i,j,0>}
				#end
				#local j=j+1;
			#end
			#local i=i+1;
		#end
	}
#end

#macro UnitSphericalTriangularPrism(a,b,c)
	#local A=vnormalize(a);
	#local B=vnormalize(b);
	#local C=vnormalize(c);
	#if(Det(A,B,C)>0)
		intersection{
			sphere{0,1}
			box{0,1 scale 1+1/vlength(A+B+C) Matrix(a,b,c,0)}
		}
	#else
		difference{
			sphere{0,1}
			box{0,1 scale 1+1/vlength(A+B+C) Matrix(a,b,c,0)}
		}
	#end
#end

#macro Sphere4(a,b,c,d)
	Sphere(CircumCenter4(a,b,c,d),vlength(a-CircumCenter4(a,b,c,d)))
#end

//Textit
//Textrm
//Display
//TTFObj

//DEG3
