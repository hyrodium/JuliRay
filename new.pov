#version 3.7;
global_settings{assumed_gamma 1.0}
camera{perspective location <0.0, 0.0, 0.0> right <1.0, 0.0, 0.0> up <0.0, 1.0, 0.0> direction <0.0, 0.0, 1.0> sky <0.0, 1.0, 0.0> look_at <0.0, 0.0, -1.0>}
light_source{<0.0, 0.0, 0.0> rgb<1.0,1.0,1.0>}
background{rgb<1,1,1>}
merge{object{sphere{<0.666667, -0.333333, -0.333333>,0.6} pigment{rgbft<1.0,0.0,0.0,0.1,0.2>}}object{sphere{<-0.333333, 0.666667, -0.333333>,0.6} pigment{rgbft<0.5,0.0,0.5,0.1,0.2>}}object{sphere{<-0.333333, -0.333333, 0.666667>,0.6} pigment{rgbft<0.0,1.0,0.0,0.1,0.2>}}}
