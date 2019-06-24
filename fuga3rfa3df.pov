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

union{sphere{<0.0, 0.0, 0.0>,0.05}object{union{cylinder{<0.0, 0.0, 0.0>,<0.94, 0.0, 0.0>,0.01}cone{<0.94, 0.0, 0.0>,0.02,<1.0, 0.0, 0.0>,0}} pigment{rgb<1.0, 0.0, 0.0>}}object{union{cylinder{<0.0, 0.0, 0.0>,<0.0, 0.88, 0.0>,0.02}cone{<0.0, 0.88, 0.0>,0.04,<0.0, 1.0, 0.0>,0}} pigment{rgb<0.0, 1.0, 0.0>}}object{union{cylinder{<0.0, 0.0, 0.0>,<0.0, 0.0, 0.82>,0.03}cone{<0.0, 0.0, 0.82>,0.06,<0.0, 0.0, 1.0>,0}} pigment{rgb<0.0, 0.0, 1.0>}}}
