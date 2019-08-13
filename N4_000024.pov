#version 3.7;
global_settings{assumed_gamma 1.0}
camera{perspective location <18.750000000000000000000000, 10.825317547305481014063844, 12.499999999999996447286321> right <0.133333333333333331482962, -0.230940107675850353485814, -0.000000000000000000000000> up <-0.064951905283832905846353, -0.037499999999999998612221, 0.129903810567665811692706> direction <0.750000000000000111022302, 0.433012701892219298294151, 0.499999999999999944488849> sky <-0.433012701892219298294151, -0.249999999999999944488849, 0.866025403784438707610605> look_at <18.000000000000000000000000, 10.392304845413262270881205, 11.999999999999996447286321>}
light_source{<18.750000000000000000000000, 10.825317547305481014063844, 12.499999999999996447286321> rgb<1.0,1.0,1.0>}
background{rgb<1,1,1>}
union{union{object{union{sphere{<-0.471404520791031456106168, 0.816496580927726034460079, -0.006158660039668362173870>,0.029999999999999998889777}sphere{<-0.471404520791032122239983, -0.816496580927725923437777, -0.006158660039668362173870>,0.029999999999999998889777}sphere{<0.942809041582063356301546, -0.000000000000000230921615, -0.006158660039668362173870>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<-0.471404520791031456106168, 0.816496580927726034460079, -0.006158660039668362173870>,<-0.471404520791032122239983, -0.816496580927725923437777, -0.006158660039668362173870>,0.010000000000000000208167}cylinder{<-0.471404520791032122239983, -0.816496580927725923437777, -0.006158660039668362173870>,<0.942809041582063356301546, -0.000000000000000230921615, -0.006158660039668362173870>,0.010000000000000000208167}cylinder{<0.942809041582063356301546, -0.000000000000000230921615, -0.006158660039668362173870>,<-0.471404520791031456106168, 0.816496580927726034460079, -0.006158660039668362173870>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{3.000000000000000000000000,<-0.471404520791031456106168, 0.816496580927726034460079, -0.006158660039668362173870>,<-0.471404520791032122239983, -0.816496580927725923437777, -0.006158660039668362173870>,<0.942809041582063356301546, -0.000000000000000230921615, -0.006158660039668362173870>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<-0.471404520791031511617319, 0.816496580927726034460079, -0.006158660039668348296082>,0.029999999999999998889777}sphere{<-0.471404520791032177751134, -0.816496580927725923437777, -0.006158660039668348296082>,0.029999999999999998889777}sphere{<-1.877920224880325728022967, 0.000000000000000959184826, 0.141196605051872592362372>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<-0.471404520791031511617319, 0.816496580927726034460079, -0.006158660039668348296082>,<-0.471404520791032177751134, -0.816496580927725923437777, -0.006158660039668348296082>,0.010000000000000000208167}cylinder{<-0.471404520791032177751134, -0.816496580927725923437777, -0.006158660039668348296082>,<-1.877920224880325728022967, 0.000000000000000959184826, 0.141196605051872592362372>,0.010000000000000000208167}cylinder{<-1.877920224880325728022967, 0.000000000000000959184826, 0.141196605051872592362372>,<-0.471404520791031511617319, 0.816496580927726034460079, -0.006158660039668348296082>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{3.000000000000000000000000,<-0.471404520791031511617319, 0.816496580927726034460079, -0.006158660039668348296082>,<-0.471404520791032177751134, -0.816496580927725923437777, -0.006158660039668348296082>,<-1.877920224880325728022967, 0.000000000000000959184826, 0.141196605051872592362372>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<-1.858738430734123703302885, -1.626326621026945806747221, 0.287348742152464864396677>,0.029999999999999998889777}sphere{<-0.471404520791032344284588, -0.816496580927725590370869, -0.006158660039668445440597>,0.029999999999999998889777}sphere{<-1.877920224880325283933757, 0.000000000000001110223025, 0.141196605051872453584494>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<-1.858738430734123703302885, -1.626326621026945806747221, 0.287348742152464864396677>,<-0.471404520791032344284588, -0.816496580927725590370869, -0.006158660039668445440597>,0.010000000000000000208167}cylinder{<-0.471404520791032344284588, -0.816496580927725590370869, -0.006158660039668445440597>,<-1.877920224880325283933757, 0.000000000000001110223025, 0.141196605051872453584494>,0.010000000000000000208167}cylinder{<-1.877920224880325283933757, 0.000000000000001110223025, 0.141196605051872453584494>,<-1.858738430734123703302885, -1.626326621026945806747221, 0.287348742152464864396677>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{3.000000000000000000000000,<-1.858738430734123703302885, -1.626326621026945806747221, 0.287348742152464864396677>,<-0.471404520791032344284588, -0.816496580927725590370869, -0.006158660039668445440597>,<-1.877920224880325283933757, 0.000000000000001110223025, 0.141196605051872453584494>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<-0.471404520791031678150773, 0.816496580927726478549289, -0.006158660039668834018656>,0.029999999999999998889777}sphere{<-1.858738430734122148990650, 1.626326621026948915371690, 0.287348742152464531329770>,0.029999999999999998889777}sphere{<-1.877920224880325505978362, 0.000000000000001110223025, 0.141196605051872065006435>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<-0.471404520791031678150773, 0.816496580927726478549289, -0.006158660039668834018656>,<-1.858738430734122148990650, 1.626326621026948915371690, 0.287348742152464531329770>,0.010000000000000000208167}cylinder{<-1.858738430734122148990650, 1.626326621026948915371690, 0.287348742152464531329770>,<-1.877920224880325505978362, 0.000000000000001110223025, 0.141196605051872065006435>,0.010000000000000000208167}cylinder{<-1.877920224880325505978362, 0.000000000000001110223025, 0.141196605051872065006435>,<-0.471404520791031678150773, 0.816496580927726478549289, -0.006158660039668834018656>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{3.000000000000000000000000,<-0.471404520791031678150773, 0.816496580927726478549289, -0.006158660039668834018656>,<-1.858738430734122148990650, 1.626326621026948915371690, 0.287348742152464531329770>,<-1.877920224880325505978362, 0.000000000000001110223025, 0.141196605051872065006435>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}}
