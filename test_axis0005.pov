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

union{sphere{<0.0, 0.0, 0.0>,0.05}object{union{cylinder{<0.0, 0.0, 0.0>,<0.94, 0.0, 0.0>,0.01}cone{<0.94, 0.0, 0.0>,0.02,<1.0, 0.0, 0.0>,0}} pigment{rgb<1.0, 0.0, 0.0>}}object{union{cylinder{<0.0, 0.0, 0.0>,<0.0, 0.85, 0.0>,0.024999999999999994}cone{<0.0, 0.85, 0.0>,0.04999999999999999,<0.0, 1.0, 0.0>,0}} pigment{rgb<0.0, 1.0, 0.0>}}object{union{cylinder{<0.0, 0.0, 0.0>,<0.0, 0.0, 0.85>,0.025000000000000015}cone{<0.0, 0.0, 0.85>,0.05000000000000003,<0.0, 0.0, 1.0>,0}} pigment{rgb<0.0, 0.0, 1.0>}}object{cone{<-1.0, 0.0, 0.0>,0.1,<-0.8, 0.5, 1.0>,0} pigment{rgb<1.0, 1.0, 0.5>}}object{object{cylinder{<0.0, 0.0, 0.0>,<0.0, 0.0, 1.0>,0.1} pigment{rgb<1.0, 0.2, 0.3>}} matrix<1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, -0.7, -0.6>}object{box{<0.9, 0.9, 0.9>,<0.0, 0.7, 0.4>} pigment{rgb<0.2, 0.3, 0.5>}}object{object{torus{0.6,0.15} pigment{rgb<0.9, 0.1, 0.15>}} matrix<1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, -0.0, 1.0, 0.0, 0.0, 0.0>}object{object{torus{0.6,0.1} pigment{rgb<0.2, 0.8, 0.1>}} matrix<1.0, 0.0, 0.0, -0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0>}object{disc{<1.0, -0.95, 0.95>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.85, 0.85>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.75, 0.75>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.65, 0.65>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.55, 0.55>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.45, 0.45>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.35, 0.35>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.25, 0.25>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.15, 0.15>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.05, 0.05>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.05, -0.05>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.15, -0.15>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.25, -0.25>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.35, -0.35>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.45, -0.45>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.55, -0.55>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.65, -0.65>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.75, -0.75>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.85, -0.85>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.95, -0.95>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 1.05, -1.05>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.95, 0.95>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.85, 0.85>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.75, 0.75>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.65, 0.65>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.55, 0.55>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.45, 0.45>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.35, 0.35>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.25, 0.25>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.15, 0.15>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.05, 0.05>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.05, -0.05>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.15, -0.15>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.25, -0.25>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.35, -0.35>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.45, -0.45>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.55, -0.55>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.65, -0.65>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.75, -0.75>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.85, -0.85>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.95, -0.95>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 1.05, -1.05>} pigment{rgb<0.0, 0.8, 0.8>}}object{polygon{4,<0.0, 0.0, 0.5>,<1.0, 0.0, 0.0>,<1.0, 1.0, -0.5>,<0.0, 1.0, 0.0>} pigment{rgb<0.0, 1.0, 1.0>}}}
