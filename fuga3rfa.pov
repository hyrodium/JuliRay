#version 3.7;
global_settings{assumed_gamma 1.0}

#include "Hy_constants.inc"
#include "Hy_functions.inc"
#include "Hy_colors.inc"

//#declare Time=1;
/* #declare N_Time=5; */
#declare Dir_Time=1;
#include "Hy_clock.inc"

#declare Lng=30;
#declare Lat=30;
#declare Pers=0.1;
#declare Zoom=0.9;
#declare LookAt=<0,0,0>;
#include "Hy_camera.inc"


object{sphere{<1.0, 0.0, 0.0>,0.1}}
