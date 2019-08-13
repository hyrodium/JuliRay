#version 3.7;
global_settings{assumed_gamma 1.0}
camera{perspective location <18.750000000000000000000000, 10.825317547305481014063844, 12.499999999999996447286321> right <0.133333333333333331482962, -0.230940107675850353485814, -0.000000000000000000000000> up <-0.064951905283832905846353, -0.037499999999999998612221, 0.129903810567665811692706> direction <0.750000000000000111022302, 0.433012701892219298294151, 0.499999999999999944488849> sky <-0.433012701892219298294151, -0.249999999999999944488849, 0.866025403784438707610605> look_at <18.000000000000000000000000, 10.392304845413262270881205, 11.999999999999996447286321>}
light_source{<18.750000000000000000000000, 10.825317547305481014063844, 12.499999999999996447286321> rgb<1.0,1.0,1.0>}
background{rgb<1,1,1>}
union{union{object{union{sphere{<-0.471404520791031456106168, 0.816496580927726034460079, 0.000000000000000000000000>,0.029999999999999998889777}sphere{<-0.471404520791032122239983, -0.816496580927725923437777, 0.000000000000000000000000>,0.029999999999999998889777}sphere{<0.942809041582063356301546, -0.000000000000000230921615, 0.000000000000000000000000>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<-0.471404520791031456106168, 0.816496580927726034460079, 0.000000000000000000000000>,<-0.471404520791032122239983, -0.816496580927725923437777, 0.000000000000000000000000>,0.010000000000000000208167}cylinder{<-0.471404520791032122239983, -0.816496580927725923437777, 0.000000000000000000000000>,<0.942809041582063356301546, -0.000000000000000230921615, 0.000000000000000000000000>,0.010000000000000000208167}cylinder{<0.942809041582063356301546, -0.000000000000000230921615, 0.000000000000000000000000>,<-0.471404520791031456106168, 0.816496580927726034460079, 0.000000000000000000000000>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{3.000000000000000000000000,<-0.471404520791031456106168, 0.816496580927726034460079, 0.000000000000000000000000>,<-0.471404520791032122239983, -0.816496580927725923437777, 0.000000000000000000000000>,<0.942809041582063356301546, -0.000000000000000230921615, 0.000000000000000000000000>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<-0.471404520791031511617319, 0.816496580927726034460079, -0.000000000008906680407126>,0.029999999999999998889777}sphere{<-0.471404520791032066728832, -0.816496580927725923437777, -0.000000000008906680407126>,0.029999999999999998889777}sphere{<-1.885618083164126934647697, 0.000000000000000828354256, 0.000000000017813360814252>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<-0.471404520791031511617319, 0.816496580927726034460079, -0.000000000008906680407126>,<-0.471404520791032066728832, -0.816496580927725923437777, -0.000000000008906680407126>,0.010000000000000000208167}cylinder{<-0.471404520791032066728832, -0.816496580927725923437777, -0.000000000008906680407126>,<-1.885618083164126934647697, 0.000000000000000828354256, 0.000000000017813360814252>,0.010000000000000000208167}cylinder{<-1.885618083164126934647697, 0.000000000000000828354256, 0.000000000017813360814252>,<-0.471404520791031511617319, 0.816496580927726034460079, -0.000000000008906680407126>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{3.000000000000000000000000,<-0.471404520791031511617319, 0.816496580927726034460079, -0.000000000008906680407126>,<-0.471404520791032066728832, -0.816496580927725923437777, -0.000000000008906680407126>,<-1.885618083164126934647697, 0.000000000000000828354256, 0.000000000017813360814252>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<-1.885618083164127156692302, -1.632993161855451402786343, 0.000000000008906680407126>,0.029999999999999998889777}sphere{<-0.471404520791031955706529, -0.816496580927726145482382, -0.000000000044533402035631>,0.029999999999999998889777}sphere{<-1.885618083164126934647697, 0.000000000000000444089210, -0.000000000017813360814252>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<-1.885618083164127156692302, -1.632993161855451402786343, 0.000000000008906680407126>,<-0.471404520791031955706529, -0.816496580927726145482382, -0.000000000044533402035631>,0.010000000000000000208167}cylinder{<-0.471404520791031955706529, -0.816496580927726145482382, -0.000000000044533402035631>,<-1.885618083164126934647697, 0.000000000000000444089210, -0.000000000017813360814252>,0.010000000000000000208167}cylinder{<-1.885618083164126934647697, 0.000000000000000444089210, -0.000000000017813360814252>,<-1.885618083164127156692302, -1.632993161855451402786343, 0.000000000008906680407126>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{3.000000000000000000000000,<-1.885618083164127156692302, -1.632993161855451402786343, 0.000000000008906680407126>,<-0.471404520791031955706529, -0.816496580927726145482382, -0.000000000044533402035631>,<-1.885618083164126934647697, 0.000000000000000444089210, -0.000000000017813360814252>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<-0.471404520791031400595017, 0.816496580927726145482382, -0.000000000044533402035631>,0.029999999999999998889777}sphere{<-1.885618083164126046469278, 1.632993161855453401187788, 0.000000000008906680407126>,0.029999999999999998889777}sphere{<-1.885618083164127156692302, 0.000000000000000666133815, -0.000000000017813360814252>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<-0.471404520791031400595017, 0.816496580927726145482382, -0.000000000044533402035631>,<-1.885618083164126046469278, 1.632993161855453401187788, 0.000000000008906680407126>,0.010000000000000000208167}cylinder{<-1.885618083164126046469278, 1.632993161855453401187788, 0.000000000008906680407126>,<-1.885618083164127156692302, 0.000000000000000666133815, -0.000000000017813360814252>,0.010000000000000000208167}cylinder{<-1.885618083164127156692302, 0.000000000000000666133815, -0.000000000017813360814252>,<-0.471404520791031400595017, 0.816496580927726145482382, -0.000000000044533402035631>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{3.000000000000000000000000,<-0.471404520791031400595017, 0.816496580927726145482382, -0.000000000044533402035631>,<-1.885618083164126046469278, 1.632993161855453401187788, 0.000000000008906680407126>,<-1.885618083164127156692302, 0.000000000000000666133815, -0.000000000017813360814252>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}}
