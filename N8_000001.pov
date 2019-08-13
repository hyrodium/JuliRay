#version 3.7;
global_settings{assumed_gamma 1.0}
camera{perspective location <18.750000000000000000000000, 10.825317547305481014063844, 12.499999999999996447286321> right <0.133333333333333331482962, -0.230940107675850353485814, -0.000000000000000000000000> up <-0.064951905283832905846353, -0.037499999999999998612221, 0.129903810567665811692706> direction <0.750000000000000111022302, 0.433012701892219298294151, 0.499999999999999944488849> sky <-0.433012701892219298294151, -0.249999999999999944488849, 0.866025403784438707610605> look_at <18.000000000000000000000000, 10.392304845413262270881205, 11.999999999999996447286321>}
light_source{<18.750000000000000000000000, 10.825317547305481014063844, 12.499999999999996447286321> rgb<1.0,1.0,1.0>}
background{rgb<1,1,1>}
union{union{object{union{sphere{<-0.408248290463862795185435, 0.707106781186547461715008, -0.577350269189625842081171>,0.029999999999999998889777}sphere{<-0.408248290463863350296947, -0.707106781186547350692706, -0.577350269189625842081171>,0.029999999999999998889777}sphere{<0.816496580927725923437777, -0.000000000000000199983985, -0.577350269189625842081171>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<-0.408248290463862795185435, 0.707106781186547461715008, -0.577350269189625842081171>,<-0.408248290463863350296947, -0.707106781186547350692706, -0.577350269189625842081171>,0.010000000000000000208167}cylinder{<-0.408248290463863350296947, -0.707106781186547350692706, -0.577350269189625842081171>,<0.816496580927725923437777, -0.000000000000000199983985, -0.577350269189625842081171>,0.010000000000000000208167}cylinder{<0.816496580927725923437777, -0.000000000000000199983985, -0.577350269189625842081171>,<-0.408248290463862795185435, 0.707106781186547461715008, -0.577350269189625842081171>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{3.000000000000000000000000,<-0.408248290463862795185435, 0.707106781186547461715008, -0.577350269189625842081171>,<-0.408248290463863350296947, -0.707106781186547350692706, -0.577350269189625842081171>,<0.816496580927725923437777, -0.000000000000000199983985, -0.577350269189625842081171>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<-0.408248290463862850696586, 0.707106781186547461715008, -0.577350269189625731058868>,0.029999999999999998889777}sphere{<-0.408248290463863350296947, -0.707106781186547350692706, -0.577350269189625953103473>,0.029999999999999998889777}sphere{<-0.816496580927726034460079, 0.000000000000000410638679, 0.577350269189625731058868>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<-0.408248290463862850696586, 0.707106781186547461715008, -0.577350269189625731058868>,<-0.408248290463863350296947, -0.707106781186547350692706, -0.577350269189625953103473>,0.010000000000000000208167}cylinder{<-0.408248290463863350296947, -0.707106781186547350692706, -0.577350269189625953103473>,<-0.816496580927726034460079, 0.000000000000000410638679, 0.577350269189625731058868>,0.010000000000000000208167}cylinder{<-0.816496580927726034460079, 0.000000000000000410638679, 0.577350269189625731058868>,<-0.408248290463862850696586, 0.707106781186547461715008, -0.577350269189625731058868>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{3.000000000000000000000000,<-0.408248290463862850696586, 0.707106781186547461715008, -0.577350269189625731058868>,<-0.408248290463863350296947, -0.707106781186547350692706, -0.577350269189625953103473>,<-0.816496580927726034460079, 0.000000000000000410638679, 0.577350269189625731058868>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<-0.408248290463862850696586, 0.707106781186547350692706, -0.577350269189625731058868>,0.029999999999999998889777}sphere{<0.408248290463863350296947, 0.707106781186548016826521, 0.577350269189625509014263>,0.029999999999999998889777}sphere{<-0.816496580927726034460079, 0.000000000000000527355937, 0.577350269189625842081171>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<-0.408248290463862850696586, 0.707106781186547350692706, -0.577350269189625731058868>,<0.408248290463863350296947, 0.707106781186548016826521, 0.577350269189625509014263>,0.010000000000000000208167}cylinder{<0.408248290463863350296947, 0.707106781186548016826521, 0.577350269189625509014263>,<-0.816496580927726034460079, 0.000000000000000527355937, 0.577350269189625842081171>,0.010000000000000000208167}cylinder{<-0.816496580927726034460079, 0.000000000000000527355937, 0.577350269189625842081171>,<-0.408248290463862850696586, 0.707106781186547350692706, -0.577350269189625731058868>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{3.000000000000000000000000,<-0.408248290463862850696586, 0.707106781186547350692706, -0.577350269189625731058868>,<0.408248290463863350296947, 0.707106781186548016826521, 0.577350269189625509014263>,<-0.816496580927726034460079, 0.000000000000000527355937, 0.577350269189625842081171>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<0.408248290463862906207737, -0.707106781186547350692706, 0.577350269189626619237288>,0.029999999999999998889777}sphere{<0.408248290463863294785796, 0.707106781186548127848823, 0.577350269189625731058868>,0.029999999999999998889777}sphere{<-0.816496580927726034460079, 0.000000000000000610622664, 0.577350269189625953103473>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<0.408248290463862906207737, -0.707106781186547350692706, 0.577350269189626619237288>,<0.408248290463863294785796, 0.707106781186548127848823, 0.577350269189625731058868>,0.010000000000000000208167}cylinder{<0.408248290463863294785796, 0.707106781186548127848823, 0.577350269189625731058868>,<-0.816496580927726034460079, 0.000000000000000610622664, 0.577350269189625953103473>,0.010000000000000000208167}cylinder{<-0.816496580927726034460079, 0.000000000000000610622664, 0.577350269189625953103473>,<0.408248290463862906207737, -0.707106781186547350692706, 0.577350269189626619237288>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{3.000000000000000000000000,<0.408248290463862906207737, -0.707106781186547350692706, 0.577350269189626619237288>,<0.408248290463863294785796, 0.707106781186548127848823, 0.577350269189625731058868>,<-0.816496580927726034460079, 0.000000000000000610622664, 0.577350269189625953103473>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<0.408248290463862906207737, -0.707106781186547461715008, 0.577350269189626619237288>,0.029999999999999998889777}sphere{<-0.408248290463863572341552, -0.707106781186548349893428, -0.577350269189625286969658>,0.029999999999999998889777}sphere{<-0.816496580927726034460079, 0.000000000000000610622664, 0.577350269189625953103473>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<0.408248290463862906207737, -0.707106781186547461715008, 0.577350269189626619237288>,<-0.408248290463863572341552, -0.707106781186548349893428, -0.577350269189625286969658>,0.010000000000000000208167}cylinder{<-0.408248290463863572341552, -0.707106781186548349893428, -0.577350269189625286969658>,<-0.816496580927726034460079, 0.000000000000000610622664, 0.577350269189625953103473>,0.010000000000000000208167}cylinder{<-0.816496580927726034460079, 0.000000000000000610622664, 0.577350269189625953103473>,<0.408248290463862906207737, -0.707106781186547461715008, 0.577350269189626619237288>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{3.000000000000000000000000,<0.408248290463862906207737, -0.707106781186547461715008, 0.577350269189626619237288>,<-0.408248290463863572341552, -0.707106781186548349893428, -0.577350269189625286969658>,<-0.816496580927726034460079, 0.000000000000000610622664, 0.577350269189625953103473>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<0.408248290463863239274644, -0.707106781186547905804218, 0.577350269189625509014263>,0.029999999999999998889777}sphere{<-0.408248290463863350296947, -0.707106781186547572737311, -0.577350269189625953103473>,0.029999999999999998889777}sphere{<0.816496580927725812415474, -0.000000000000000027755576, -0.577350269189625953103473>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<0.408248290463863239274644, -0.707106781186547905804218, 0.577350269189625509014263>,<-0.408248290463863350296947, -0.707106781186547572737311, -0.577350269189625953103473>,0.010000000000000000208167}cylinder{<-0.408248290463863350296947, -0.707106781186547572737311, -0.577350269189625953103473>,<0.816496580927725812415474, -0.000000000000000027755576, -0.577350269189625953103473>,0.010000000000000000208167}cylinder{<0.816496580927725812415474, -0.000000000000000027755576, -0.577350269189625953103473>,<0.408248290463863239274644, -0.707106781186547905804218, 0.577350269189625509014263>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{3.000000000000000000000000,<0.408248290463863239274644, -0.707106781186547905804218, 0.577350269189625509014263>,<-0.408248290463863350296947, -0.707106781186547572737311, -0.577350269189625953103473>,<0.816496580927725812415474, -0.000000000000000027755576, -0.577350269189625953103473>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<0.408248290463863128252342, -0.707106781186548016826521, 0.577350269189625731058868>,0.029999999999999998889777}sphere{<0.408248290463863627852703, 0.707106781186547461715008, 0.577350269189626175148078>,0.029999999999999998889777}sphere{<0.816496580927725812415474, 0.000000000000000055511151, -0.577350269189626175148078>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<0.408248290463863128252342, -0.707106781186548016826521, 0.577350269189625731058868>,<0.408248290463863627852703, 0.707106781186547461715008, 0.577350269189626175148078>,0.010000000000000000208167}cylinder{<0.408248290463863627852703, 0.707106781186547461715008, 0.577350269189626175148078>,<0.816496580927725812415474, 0.000000000000000055511151, -0.577350269189626175148078>,0.010000000000000000208167}cylinder{<0.816496580927725812415474, 0.000000000000000055511151, -0.577350269189626175148078>,<0.408248290463863128252342, -0.707106781186548016826521, 0.577350269189625731058868>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{3.000000000000000000000000,<0.408248290463863128252342, -0.707106781186548016826521, 0.577350269189625731058868>,<0.408248290463863627852703, 0.707106781186547461715008, 0.577350269189626175148078>,<0.816496580927725812415474, 0.000000000000000055511151, -0.577350269189626175148078>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<-0.408248290463863683363854, 0.707106781186548238871126, -0.577350269189625064925053>,0.029999999999999998889777}sphere{<0.408248290463863683363854, 0.707106781186547683759613, 0.577350269189626397192683>,0.029999999999999998889777}sphere{<0.816496580927725923437777, -0.000000000000000027755576, -0.577350269189626397192683>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<-0.408248290463863683363854, 0.707106781186548238871126, -0.577350269189625064925053>,<0.408248290463863683363854, 0.707106781186547683759613, 0.577350269189626397192683>,0.010000000000000000208167}cylinder{<0.408248290463863683363854, 0.707106781186547683759613, 0.577350269189626397192683>,<0.816496580927725923437777, -0.000000000000000027755576, -0.577350269189626397192683>,0.010000000000000000208167}cylinder{<0.816496580927725923437777, -0.000000000000000027755576, -0.577350269189626397192683>,<-0.408248290463863683363854, 0.707106781186548238871126, -0.577350269189625064925053>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{3.000000000000000000000000,<-0.408248290463863683363854, 0.707106781186548238871126, -0.577350269189625064925053>,<0.408248290463863683363854, 0.707106781186547683759613, 0.577350269189626397192683>,<0.816496580927725923437777, -0.000000000000000027755576, -0.577350269189626397192683>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}}
