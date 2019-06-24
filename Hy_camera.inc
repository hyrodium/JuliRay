//Hyrodium's POV-Ray format
//Hy_transformations.inc

//#include "Hy_camera.inc"
#ifndef(AspectRatio) #declare AspectRatio=image_width/image_height; #end
#ifndef(Lng) #declare Lng=-90; #end
#ifndef(Lat) #declare Lat=90; #end
#ifndef(Tilt) #declare Tilt=0; #end
#ifndef(Pers) #declare Pers=0; #end
#ifndef(Zoom) #declare Zoom=1; #end
#ifndef(LookAt) #declare LookAt=<0,0,0>; #end

#declare Z=SCS(Lng,Lat);
#declare X=vaxis_rotate(<-sin(radians(Lng)),cos(radians(Lng)),0>,Z,Tilt);
#declare Y=vcross(Z,X);

#if(Pers)
	#declare Loc=LookAt+SCS(Lng,Lat)/(Zoom*Pers);
	camera{
		perspective
		location Loc
		right -2*X*sqrt(AspectRatio)/Zoom
		up 2*Y/(sqrt(AspectRatio)*Zoom)
		direction Z/(Zoom*Pers)
	 	sky Y
		look_at LookAt
	}
	light_source{
		Loc
		color White
	}
#else
	#declare Loc=SCS(Lng,Lat);
	camera{
		orthographic
		location Loc*100
		right -2*X*sqrt(AspectRatio)/Zoom
		up 2*Y/(sqrt(AspectRatio)*Zoom)
	 	sky Y
		look_at LookAt
	}
	light_source{
		SCS(Lng,Lat)
		color White
		parallel
		point_at 0
	}
#end

background{White}




/*
#declare AspectRatio=4/3;
#declare Lng=0;
#declare Lat=30;
#declare Pers=0.2;
#declare Zoom=0.8;
#declare LookAt=<0,0,0>;
#local loc=LookAt+<cos(radians(Lat))*cos(radians(Lng)),cos(radians(Lat))*sin(radians(Lng)),sin(radians(Lat))>/(Zoom*Pers);
camera{
	location loc
	right <-AspectRatio,0,0>
	sky <0,0,1>
	look_at LookAt
	angle 180*Pers
}
*/
