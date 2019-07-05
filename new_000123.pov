#version 3.7;
global_settings{assumed_gamma 1.0}

#include "Hy_constants.inc"
#include "Hy_functions.inc"
#include "Hy_colors.inc"

#declare Lng=30;
#declare Lat=30;
#declare Pers=0.1;
#declare Zoom=0.9;
#declare LookAt=<0,0,0>;
#include "Hy_camera.inc"

merge{object{object{sphere{<1.0, 0.0, 0.0>,0.5} pigment{rgbft<1.0,0.0,0.0,0.1,0.2>}} matrix<1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0>}object{sphere{<0.0, 1.0, 0.0>,0.5} pigment{rgbft<0.0,1.0,0.0,0.1,0.2>}}}
