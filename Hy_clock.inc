//Hyrodium's POV-Ray format
//Hy_color.inc

//#include "Hy_color.inc"
#ifndef(Dir_Time) #declare Dir_Time=1; #end
#ifndef(Time) #declare Time=clock; #end
#if(Dir_Time=-1)
	#declare Time=1-Time;
#end
