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

union{sphere{<0.0, 0.0, 0.0>,0.05}object{union{cylinder{<0.0, 0.0, 0.0>,<0.928541, 0.0, 0.0>,0.011909830056250527}cone{<0.928541, 0.0, 0.0>,0.023819660112501053,<1.0, 0.0, 0.0>,0}} pigment{rgb<1.0, 0.0, 0.0>}}object{union{cylinder{<0.0, 0.0, 0.0>,<0.0, 0.825187, 0.0>,0.02913545457642601}cone{<0.0, 0.825187, 0.0>,0.05827090915285202,<0.0, 1.0, 0.0>,0}} pigment{rgb<0.0, 1.0, 0.0>}}object{union{cylinder{<0.0, 0.0, 0.0>,<0.0, 0.0, 0.886272>,0.01895471536732347}cone{<0.0, 0.0, 0.886272>,0.03790943073464694,<0.0, 0.0, 1.0>,0}} pigment{rgb<0.0, 0.0, 1.0>}}object{cone{<-1.0, 0.0, 0.0>,0.1,<-0.8, 0.5, 1.0>,0} pigment{rgb<1.0, 1.0, 0.793893>}}object{object{cylinder{<0.0, 0.0, 0.0>,<0.0, 0.0, 1.0>,0.1} pigment{rgb<1.0, 0.2, 0.3>}} matrix<1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, -0.7, -0.542705>}object{box{<0.9, 0.9, 0.9>,<0.0, 0.611832, 0.4>} pigment{rgb<0.2, 0.3, 0.5>}}object{object{torus{0.6,0.15} pigment{rgb<0.9, 0.1, 0.15>}} matrix<1.0, 0.0, 0.0, 0.0, 0.988771, 0.149436, 0.0, -0.149436, 0.988771, 0.0, 0.0, 0.0>}object{object{torus{0.6,0.1} pigment{rgb<0.2, 0.8, 0.1>}} matrix<0.988771, 0.149436, 0.0, -0.149436, 0.988771, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0>}object{disc{<1.0, -0.94, 0.94>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.84, 0.84>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.74, 0.74>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.64, 0.64>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.54, 0.54>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.44, 0.44>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.34, 0.34>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.24, 0.24>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.14, 0.14>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.04, 0.04>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.06, -0.06>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.16, -0.16>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.26, -0.26>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.36, -0.36>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.46, -0.46>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.56, -0.56>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.66, -0.66>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.76, -0.76>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.86, -0.86>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.96, -0.96>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 1.06, -1.06>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.94, 0.94>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.84, 0.84>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.74, 0.74>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.64, 0.64>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.54, 0.54>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.44, 0.44>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.34, 0.34>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.24, 0.24>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.14, 0.14>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.04, 0.04>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.06, -0.06>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.16, -0.16>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.26, -0.26>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.36, -0.36>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.46, -0.46>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.56, -0.56>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.66, -0.66>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.76, -0.76>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.86, -0.86>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.96, -0.96>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 1.06, -1.06>} pigment{rgb<0.0, 0.8, 0.8>}}intersection{object{box{<0.737107, 0.737107, 0.03>,<0.737107, 0.0, -0.03>} matrix<0.0, 0.707107, -0.707107, -0.0, 0.707107, 0.707107, 1.0, 0.0, 0.0, 0.0, 0.5, 0.5>}object{box{<0.737107, 0.737107, 0.03>,<0.737107, 0.0, -0.03>} matrix<-0.0, 0.707107, -0.707107, -0.0, 0.707107, 0.707107, 1.0, 0.0, 0.0, 0.0, 0.5, 0.5>}}}
