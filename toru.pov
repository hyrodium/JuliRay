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

object{torus{0.816496580927726,0.01} matrix<0.0, 0.707107, -0.707107, 0.57735, 0.57735, 0.57735, 0.816497, -0.408248, -0.408248, 0.333333, 0.333333, 0.333333>}
