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

union{sphere{<0.0, 0.0, 0.0>,0.05}object{union{cylinder{<0.0, 0.0, 0.0>,<0.82, 0.0, 0.0>,0.03}cone{<0.82, 0.0, 0.0>,0.06,<1.0, 0.0, 0.0>,0}} pigment{rgb<1.0, 0.0, 0.0>}}object{union{cylinder{<0.0, 0.0, 0.0>,<0.0, 0.91, 0.0>,0.015000000000000003}cone{<0.0, 0.91, 0.0>,0.030000000000000006,<0.0, 1.0, 0.0>,0}} pigment{rgb<0.0, 1.0, 0.0>}}object{union{cylinder{<0.0, 0.0, 0.0>,<0.0, 0.0, 0.91>,0.014999999999999996}cone{<0.0, 0.0, 0.91>,0.029999999999999992,<0.0, 0.0, 1.0>,0}} pigment{rgb<0.0, 0.0, 1.0>}}object{cone{<-1.0, 0.0, 0.0>,0.1,<-0.8, 0.5, 1.0>,0} pigment{rgb<1.0, 1.0, 0.5>}}object{object{cylinder{<0.0, 0.0, 0.0>,<0.0, 0.0, 1.0>,0.1} pigment{rgb<1.0, 0.2, 0.3>}} matrix<1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, -0.7, 0.0>}object{box{<0.9, 0.9, 0.9>,<0.0, 0.7, 0.4>} pigment{rgb<0.2, 0.3, 0.5>}}object{object{torus{0.6,0.15} pigment{rgb<0.9, 0.1, 0.15>}} matrix<1.0, 0.0, 0.0, 0.0, 6.12323e-17, 1.0, 0.0, -1.0, 6.12323e-17, 0.0, 0.0, 0.0>}object{object{torus{0.6,0.1} pigment{rgb<0.2, 0.8, 0.1>}} matrix<6.12323e-17, 1.0, 0.0, -1.0, 6.12323e-17, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0>}object{disc{<1.0, -1.0, 1.0>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.9, 0.9>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.8, 0.8>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.7, 0.7>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.6, 0.6>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.5, 0.5>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.4, 0.4>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.3, 0.3>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.2, 0.2>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.1, 0.1>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.0, -0.0>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.1, -0.1>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.2, -0.2>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.3, -0.3>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.4, -0.4>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.5, -0.5>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.6, -0.6>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.7, -0.7>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.8, -0.8>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.9, -0.9>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 1.0, -1.0>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -1.0, 1.0>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.9, 0.9>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.8, 0.8>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.7, 0.7>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.6, 0.6>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.5, 0.5>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.4, 0.4>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.3, 0.3>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.2, 0.2>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.1, 0.1>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.0, -0.0>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.1, -0.1>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.2, -0.2>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.3, -0.3>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.4, -0.4>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.5, -0.5>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.6, -0.6>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.7, -0.7>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.8, -0.8>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.9, -0.9>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 1.0, -1.0>} pigment{rgb<0.0, 0.8, 0.8>}}object{polygon{4,<0.0, 0.0, 0.5>,<1.0, 0.0, 0.0>,<1.0, 1.0, -0.5>,<0.0, 1.0, 0.0>} pigment{rgb<0.0, 1.0, 1.0>}}}
