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

union{sphere{<0.0, 0.0, 0.0>,0.05}object{union{cylinder{<0.0, 0.0, 0.0>,<0.831459, 0.0, 0.0>,0.028090169943749473}cone{<0.831459, 0.0, 0.0>,0.056180339887498945,<1.0, 0.0, 0.0>,0}} pigment{rgb<1.0, 0.0, 0.0>}}object{union{cylinder{<0.0, 0.0, 0.0>,<0.0, 0.873728, 0.0>,0.021045284632676534}cone{<0.0, 0.873728, 0.0>,0.04209056926535307,<0.0, 1.0, 0.0>,0}} pigment{rgb<0.0, 1.0, 0.0>}}object{union{cylinder{<0.0, 0.0, 0.0>,<0.0, 0.0, 0.934813>,0.010864545423573992}cone{<0.0, 0.0, 0.934813>,0.021729090847147983,<0.0, 0.0, 1.0>,0}} pigment{rgb<0.0, 0.0, 1.0>}}object{cone{<-1.0, 0.0, 0.0>,0.1,<-0.8, 0.5, 1.0>,0} pigment{rgb<1.0, 1.0, 0.793893>}}object{object{cylinder{<0.0, 0.0, 0.0>,<0.0, 0.0, 1.0>,0.1} pigment{rgb<1.0, 0.2, 0.3>}} matrix<1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, -0.7, -0.0572949>}object{box{<0.9, 0.9, 0.9>,<0.0, 0.611832, 0.4>} pigment{rgb<0.2, 0.3, 0.5>}}object{object{torus{0.6,0.15} pigment{rgb<0.9, 0.1, 0.15>}} matrix<1.0, 0.0, 0.0, 0.0, 0.149436, 0.988771, 0.0, -0.988771, 0.149436, 0.0, 0.0, 0.0>}object{object{torus{0.6,0.1} pigment{rgb<0.2, 0.8, 0.1>}} matrix<0.149436, 0.988771, 0.0, -0.988771, 0.149436, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0>}object{disc{<1.0, -0.91, 0.91>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.81, 0.81>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.71, 0.71>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.61, 0.61>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.51, 0.51>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.41, 0.41>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.31, 0.31>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.21, 0.21>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.11, 0.11>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.01, 0.01>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.09, -0.09>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.19, -0.19>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.29, -0.29>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.39, -0.39>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.49, -0.49>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.59, -0.59>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.69, -0.69>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.79, -0.79>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.89, -0.89>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.99, -0.99>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 1.09, -1.09>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.91, 0.91>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.81, 0.81>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.71, 0.71>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.61, 0.61>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.51, 0.51>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.41, 0.41>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.31, 0.31>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.21, 0.21>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.11, 0.11>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.01, 0.01>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.09, -0.09>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.19, -0.19>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.29, -0.29>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.39, -0.39>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.49, -0.49>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.59, -0.59>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.69, -0.69>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.79, -0.79>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.89, -0.89>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.99, -0.99>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 1.09, -1.09>} pigment{rgb<0.0, 0.8, 0.8>}}intersection{object{box{<0.737107, 0.737107, 0.03>,<0.737107, 0.0, -0.03>} matrix<0.0, 0.707107, -0.707107, -0.0, 0.707107, 0.707107, 1.0, 0.0, 0.0, 0.0, 0.5, 0.5>}object{box{<0.737107, 0.737107, 0.03>,<0.737107, 0.0, -0.03>} matrix<-0.0, 0.707107, -0.707107, -0.0, 0.707107, 0.707107, 1.0, 0.0, 0.0, 0.0, 0.5, 0.5>}}}
