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

object{intersection{object{torus{0.7071067811865476,0.03} matrix<0.0, 0.0, -1.0, -1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.5, 0.5>}merge{object{box{<0.767107, 0.767107, 0.06>,<-0.767107, -0.0, -0.06>} matrix<0.0, 0.707107, -0.707107, -0.0, -0.707107, -0.707107, -1.0, 0.0, 0.0, 0.0, 0.5, 0.5>}object{box{<0.767107, 0.767107, 0.06>,<-0.767107, -0.0, -0.06>} matrix<-0.0, -0.707107, -0.707107, 0.0, -0.707107, 0.707107, -1.0, 0.0, 0.0, 0.0, 0.5, 0.5>}}} pigment{rgbft<1.0,0.0,0.0,0.5,0.5>}}
