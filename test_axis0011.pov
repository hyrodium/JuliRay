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

union{sphere{<0.0, 0.0, 0.0>,0.05}object{union{cylinder{<0.0, 0.0, 0.0>,<0.920148, 0.0, 0.0>,0.01330869393641142}cone{<0.920148, 0.0, 0.0>,0.02661738787282284,<1.0, 0.0, 0.0>,0}} pigment{rgb<1.0, 0.0, 0.0>}}object{union{cylinder{<0.0, 0.0, 0.0>,<0.0, 0.898541, 0.0>,0.016909830056250526}cone{<0.0, 0.898541, 0.0>,0.03381966011250105,<0.0, 1.0, 0.0>,0}} pigment{rgb<0.0, 1.0, 0.0>}}object{union{cylinder{<0.0, 0.0, 0.0>,<0.0, 0.0, 0.821311>,0.02978147600733806}cone{<0.0, 0.0, 0.821311>,0.05956295201467612,<0.0, 0.0, 1.0>,0}} pigment{rgb<0.0, 0.0, 1.0>}}object{cone{<-1.0, 0.0, 0.0>,0.1,<-0.8, 0.5, 1.0>,0} pigment{rgb<1.0, 1.0, 0.128428>}}object{object{cylinder{<0.0, 0.0, 0.0>,<0.0, 0.0, 1.0>,0.1} pigment{rgb<1.0, 0.2, 0.3>}} matrix<1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, -0.7, -0.500739>}object{box{<0.9, 0.9, 0.9>,<0.0, 0.811472, 0.4>} pigment{rgb<0.2, 0.3, 0.5>}}object{object{torus{0.6,0.15} pigment{rgb<0.9, 0.1, 0.15>}} matrix<1.0, 0.0, 0.0, 0.0, 0.966425, 0.256949, 0.0, -0.256949, 0.966425, 0.0, 0.0, 0.0>}object{object{torus{0.6,0.1} pigment{rgb<0.2, 0.8, 0.1>}} matrix<0.966425, 0.256949, 0.0, -0.256949, 0.966425, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0>}object{disc{<1.0, -0.963333, 0.963333>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.863333, 0.863333>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.763333, 0.763333>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.663333, 0.663333>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.563333, 0.563333>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.463333, 0.463333>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.363333, 0.363333>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.263333, 0.263333>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.163333, 0.163333>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, -0.0633333, 0.0633333>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.0366667, -0.0366667>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.136667, -0.136667>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.236667, -0.236667>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.336667, -0.336667>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.436667, -0.436667>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.536667, -0.536667>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.636667, -0.636667>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.736667, -0.736667>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.836667, -0.836667>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 0.936667, -0.936667>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{disc{<1.0, 1.03667, -1.03667>,<0.2, -0.4, 1.0>,0.1} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.963333, 0.963333>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.863333, 0.863333>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.763333, 0.763333>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.663333, 0.663333>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.563333, 0.563333>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.463333, 0.463333>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.363333, 0.363333>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.263333, 0.263333>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.163333, 0.163333>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, -0.0633333, 0.0633333>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.0366667, -0.0366667>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.136667, -0.136667>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.236667, -0.236667>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.336667, -0.336667>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.436667, -0.436667>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.536667, -0.536667>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.636667, -0.636667>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.736667, -0.736667>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.836667, -0.836667>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 0.936667, -0.936667>} pigment{rgb<0.0, 0.8, 0.8>}}object{object{torus{0.1,0.02} matrix<0.0, 0.928477, 0.371391, 0.182574, -0.365148, 0.912871, 0.983192, 0.0678064, -0.169516, 1.0, 1.03667, -1.03667>} pigment{rgb<0.0, 0.8, 0.8>}}}
