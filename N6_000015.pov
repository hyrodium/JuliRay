#version 3.7;
global_settings{assumed_gamma 1.0}
camera{perspective location <18.750000000000000000000000, 10.825317547305481014063844, 12.499999999999996447286321> right <0.133333333333333331482962, -0.230940107675850353485814, -0.000000000000000000000000> up <-0.064951905283832905846353, -0.037499999999999998612221, 0.129903810567665811692706> direction <0.750000000000000111022302, 0.433012701892219298294151, 0.499999999999999944488849> sky <-0.433012701892219298294151, -0.249999999999999944488849, 0.866025403784438707610605> look_at <18.000000000000000000000000, 10.392304845413262270881205, 11.999999999999996447286321>}
light_source{<18.750000000000000000000000, 10.825317547305481014063844, 12.499999999999996447286321> rgb<1.0,1.0,1.0>}
background{rgb<1,1,1>}
union{union{object{union{sphere{<0.000000000000000049995996, 0.816496580927725923437777, -0.180469061561687293959722>,0.029999999999999998889777}sphere{<-0.816496580927725923437777, 0.000000000000000099991992, -0.180469061561687293959722>,0.029999999999999998889777}sphere{<-0.000000000000000149987989, -0.816496580927725923437777, -0.180469061561687293959722>,0.029999999999999998889777}sphere{<0.816496580927725923437777, -0.000000000000000199983985, -0.180469061561687293959722>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<0.000000000000000049995996, 0.816496580927725923437777, -0.180469061561687293959722>,<-0.816496580927725923437777, 0.000000000000000099991992, -0.180469061561687293959722>,0.010000000000000000208167}cylinder{<-0.816496580927725923437777, 0.000000000000000099991992, -0.180469061561687293959722>,<-0.000000000000000149987989, -0.816496580927725923437777, -0.180469061561687293959722>,0.010000000000000000208167}cylinder{<-0.000000000000000149987989, -0.816496580927725923437777, -0.180469061561687293959722>,<0.816496580927725923437777, -0.000000000000000199983985, -0.180469061561687293959722>,0.010000000000000000208167}cylinder{<0.816496580927725923437777, -0.000000000000000199983985, -0.180469061561687293959722>,<0.000000000000000049995996, 0.816496580927725923437777, -0.180469061561687293959722>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{4.000000000000000000000000,<0.000000000000000049995996, 0.816496580927725923437777, -0.180469061561687293959722>,<-0.816496580927725923437777, 0.000000000000000099991992, -0.180469061561687293959722>,<-0.000000000000000149987989, -0.816496580927725923437777, -0.180469061561687293959722>,<0.816496580927725923437777, -0.000000000000000199983985, -0.180469061561687293959722>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<0.000000000000000000000000, 0.816496580927725923437777, -0.180469061561687182937419>,0.029999999999999998889777}sphere{<-0.816496580927726034460079, 0.000000000000000111022302, -0.180469061561687293959722>,0.029999999999999998889777}sphere{<-1.233052701161835962295754, 0.416556120234110149880280, 0.812654964290242309843393>,0.029999999999999998889777}sphere{<-0.416556120234110038857978, 1.233052701161835962295754, 0.812654964290242531887998>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<0.000000000000000000000000, 0.816496580927725923437777, -0.180469061561687182937419>,<-0.816496580927726034460079, 0.000000000000000111022302, -0.180469061561687293959722>,0.010000000000000000208167}cylinder{<-0.816496580927726034460079, 0.000000000000000111022302, -0.180469061561687293959722>,<-1.233052701161835962295754, 0.416556120234110149880280, 0.812654964290242309843393>,0.010000000000000000208167}cylinder{<-1.233052701161835962295754, 0.416556120234110149880280, 0.812654964290242309843393>,<-0.416556120234110038857978, 1.233052701161835962295754, 0.812654964290242531887998>,0.010000000000000000208167}cylinder{<-0.416556120234110038857978, 1.233052701161835962295754, 0.812654964290242531887998>,<0.000000000000000000000000, 0.816496580927725923437777, -0.180469061561687182937419>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{4.000000000000000000000000,<0.000000000000000000000000, 0.816496580927725923437777, -0.180469061561687182937419>,<-0.816496580927726034460079, 0.000000000000000111022302, -0.180469061561687293959722>,<-1.233052701161835962295754, 0.416556120234110149880280, 0.812654964290242309843393>,<-0.416556120234110038857978, 1.233052701161835962295754, 0.812654964290242531887998>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<-0.025092553179723961487291, 0.841589134107449621247099, 1.825989008261649448883190>,0.029999999999999998889777}sphere{<-0.841589134107450287380914, 0.025092553179723416784119, 1.825989008261649226838585>,0.029999999999999998889777}sphere{<-1.233052701161836406384964, 0.416556120234110038857978, 0.812654964290242309843393>,0.029999999999999998889777}sphere{<-0.416556120234109983346826, 1.233052701161836184340359, 0.812654964290242531887998>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<-0.025092553179723961487291, 0.841589134107449621247099, 1.825989008261649448883190>,<-0.841589134107450287380914, 0.025092553179723416784119, 1.825989008261649226838585>,0.010000000000000000208167}cylinder{<-0.841589134107450287380914, 0.025092553179723416784119, 1.825989008261649226838585>,<-1.233052701161836406384964, 0.416556120234110038857978, 0.812654964290242309843393>,0.010000000000000000208167}cylinder{<-1.233052701161836406384964, 0.416556120234110038857978, 0.812654964290242309843393>,<-0.416556120234109983346826, 1.233052701161836184340359, 0.812654964290242531887998>,0.010000000000000000208167}cylinder{<-0.416556120234109983346826, 1.233052701161836184340359, 0.812654964290242531887998>,<-0.025092553179723961487291, 0.841589134107449621247099, 1.825989008261649448883190>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{4.000000000000000000000000,<-0.025092553179723961487291, 0.841589134107449621247099, 1.825989008261649448883190>,<-0.841589134107450287380914, 0.025092553179723416784119, 1.825989008261649226838585>,<-1.233052701161836406384964, 0.416556120234110038857978, 0.812654964290242309843393>,<-0.416556120234109983346826, 1.233052701161836184340359, 0.812654964290242531887998>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<-1.233052701161835962295754, -0.416556120234109927835675, 0.812654964290242531887998>,0.029999999999999998889777}sphere{<-0.816496580927726034460079, 0.000000000000000000000000, -0.180469061561687182937419>,0.029999999999999998889777}sphere{<-0.000000000000000111022302, -0.816496580927725923437777, -0.180469061561687404982024>,0.029999999999999998889777}sphere{<-0.416556120234110149880280, -1.233052701161835962295754, 0.812654964290242309843393>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<-1.233052701161835962295754, -0.416556120234109927835675, 0.812654964290242531887998>,<-0.816496580927726034460079, 0.000000000000000000000000, -0.180469061561687182937419>,0.010000000000000000208167}cylinder{<-0.816496580927726034460079, 0.000000000000000000000000, -0.180469061561687182937419>,<-0.000000000000000111022302, -0.816496580927725923437777, -0.180469061561687404982024>,0.010000000000000000208167}cylinder{<-0.000000000000000111022302, -0.816496580927725923437777, -0.180469061561687404982024>,<-0.416556120234110149880280, -1.233052701161835962295754, 0.812654964290242309843393>,0.010000000000000000208167}cylinder{<-0.416556120234110149880280, -1.233052701161835962295754, 0.812654964290242309843393>,<-1.233052701161835962295754, -0.416556120234109927835675, 0.812654964290242531887998>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{4.000000000000000000000000,<-1.233052701161835962295754, -0.416556120234109927835675, 0.812654964290242531887998>,<-0.816496580927726034460079, 0.000000000000000000000000, -0.180469061561687182937419>,<-0.000000000000000111022302, -0.816496580927725923437777, -0.180469061561687404982024>,<-0.416556120234110149880280, -1.233052701161835962295754, 0.812654964290242309843393>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<0.603980073991055688154006, -1.045628747404890201977423, 1.319321986275945990385594>,0.029999999999999998889777}sphere{<1.020536194225166060078891, -0.629072627170780274141748, 0.326197960424016275560177>,0.029999999999999998889777}sphere{<-0.000000000000000111022302, -0.816496580927726034460079, -0.180469061561687404982024>,0.029999999999999998889777}sphere{<-0.416556120234110260902582, -1.233052701161835962295754, 0.812654964290242420865695>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<0.603980073991055688154006, -1.045628747404890201977423, 1.319321986275945990385594>,<1.020536194225166060078891, -0.629072627170780274141748, 0.326197960424016275560177>,0.010000000000000000208167}cylinder{<1.020536194225166060078891, -0.629072627170780274141748, 0.326197960424016275560177>,<-0.000000000000000111022302, -0.816496580927726034460079, -0.180469061561687404982024>,0.010000000000000000208167}cylinder{<-0.000000000000000111022302, -0.816496580927726034460079, -0.180469061561687404982024>,<-0.416556120234110260902582, -1.233052701161835962295754, 0.812654964290242420865695>,0.010000000000000000208167}cylinder{<-0.416556120234110260902582, -1.233052701161835962295754, 0.812654964290242420865695>,<0.603980073991055688154006, -1.045628747404890201977423, 1.319321986275945990385594>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{4.000000000000000000000000,<0.603980073991055688154006, -1.045628747404890201977423, 1.319321986275945990385594>,<1.020536194225166060078891, -0.629072627170780274141748, 0.326197960424016275560177>,<-0.000000000000000111022302, -0.816496580927726034460079, -0.180469061561687404982024>,<-0.416556120234110260902582, -1.233052701161835962295754, 0.812654964290242420865695>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<0.603980073991055688154006, -1.045628747404890201977423, 1.319321986275945768340989>,0.029999999999999998889777}sphere{<1.020536194225166060078891, -0.629072627170780163119446, 0.326197960424016608627085>,0.029999999999999998889777}sphere{<1.245343661728107509389929, 0.378661974440932802554727, 0.843175627835993979353191>,0.029999999999999998889777}sphere{<0.828787541493997026442742, -0.037894145793177291814402, 1.836299653687923250089398>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<0.603980073991055688154006, -1.045628747404890201977423, 1.319321986275945768340989>,<1.020536194225166060078891, -0.629072627170780163119446, 0.326197960424016608627085>,0.010000000000000000208167}cylinder{<1.020536194225166060078891, -0.629072627170780163119446, 0.326197960424016608627085>,<1.245343661728107509389929, 0.378661974440932802554727, 0.843175627835993979353191>,0.010000000000000000208167}cylinder{<1.245343661728107509389929, 0.378661974440932802554727, 0.843175627835993979353191>,<0.828787541493997026442742, -0.037894145793177291814402, 1.836299653687923250089398>,0.010000000000000000208167}cylinder{<0.828787541493997026442742, -0.037894145793177291814402, 1.836299653687923250089398>,<0.603980073991055688154006, -1.045628747404890201977423, 1.319321986275945768340989>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{4.000000000000000000000000,<0.603980073991055688154006, -1.045628747404890201977423, 1.319321986275945768340989>,<1.020536194225166060078891, -0.629072627170780163119446, 0.326197960424016608627085>,<1.245343661728107509389929, 0.378661974440932802554727, 0.843175627835993979353191>,<0.828787541493997026442742, -0.037894145793177291814402, 1.836299653687923250089398>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}}
