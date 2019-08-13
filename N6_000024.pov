#version 3.7;
global_settings{assumed_gamma 1.0}
camera{perspective location <18.750000000000000000000000, 10.825317547305481014063844, 12.499999999999996447286321> right <0.133333333333333331482962, -0.230940107675850353485814, -0.000000000000000000000000> up <-0.064951905283832905846353, -0.037499999999999998612221, 0.129903810567665811692706> direction <0.750000000000000111022302, 0.433012701892219298294151, 0.499999999999999944488849> sky <-0.433012701892219298294151, -0.249999999999999944488849, 0.866025403784438707610605> look_at <18.000000000000000000000000, 10.392304845413262270881205, 11.999999999999996447286321>}
light_source{<18.750000000000000000000000, 10.825317547305481014063844, 12.499999999999996447286321> rgb<1.0,1.0,1.0>}
background{rgb<1,1,1>}
union{union{object{union{sphere{<0.000000000000000049995996, 0.816496580927725923437777, -0.018463387614030324357373>,0.029999999999999998889777}sphere{<-0.816496580927725923437777, 0.000000000000000099991992, -0.018463387614030324357373>,0.029999999999999998889777}sphere{<-0.000000000000000149987989, -0.816496580927725923437777, -0.018463387614030324357373>,0.029999999999999998889777}sphere{<0.816496580927725923437777, -0.000000000000000199983985, -0.018463387614030324357373>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<0.000000000000000049995996, 0.816496580927725923437777, -0.018463387614030324357373>,<-0.816496580927725923437777, 0.000000000000000099991992, -0.018463387614030324357373>,0.010000000000000000208167}cylinder{<-0.816496580927725923437777, 0.000000000000000099991992, -0.018463387614030324357373>,<-0.000000000000000149987989, -0.816496580927725923437777, -0.018463387614030324357373>,0.010000000000000000208167}cylinder{<-0.000000000000000149987989, -0.816496580927725923437777, -0.018463387614030324357373>,<0.816496580927725923437777, -0.000000000000000199983985, -0.018463387614030324357373>,0.010000000000000000208167}cylinder{<0.816496580927725923437777, -0.000000000000000199983985, -0.018463387614030324357373>,<0.000000000000000049995996, 0.816496580927725923437777, -0.018463387614030324357373>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{4.000000000000000000000000,<0.000000000000000049995996, 0.816496580927725923437777, -0.018463387614030324357373>,<-0.816496580927725923437777, 0.000000000000000099991992, -0.018463387614030324357373>,<-0.000000000000000149987989, -0.816496580927725923437777, -0.018463387614030324357373>,<0.816496580927725923437777, -0.000000000000000199983985, -0.018463387614030324357373>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<0.000000000000000000000000, 0.816496580927725812415474, -0.018463387614029630467982>,0.029999999999999998889777}sphere{<-0.816496580927725701393172, 0.000000000000000000000000, -0.018463387614029644345770>,0.029999999999999998889777}sphere{<-1.626353725899101876350983, 0.809857144971375952913206, 0.128492873943470253816912>,0.029999999999999998889777}sphere{<-0.809857144971376397002416, 1.626353725899101876350983, 0.128492873943470309328063>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<0.000000000000000000000000, 0.816496580927725812415474, -0.018463387614029630467982>,<-0.816496580927725701393172, 0.000000000000000000000000, -0.018463387614029644345770>,0.010000000000000000208167}cylinder{<-0.816496580927725701393172, 0.000000000000000000000000, -0.018463387614029644345770>,<-1.626353725899101876350983, 0.809857144971375952913206, 0.128492873943470253816912>,0.010000000000000000208167}cylinder{<-1.626353725899101876350983, 0.809857144971375952913206, 0.128492873943470253816912>,<-0.809857144971376397002416, 1.626353725899101876350983, 0.128492873943470309328063>,0.010000000000000000208167}cylinder{<-0.809857144971376397002416, 1.626353725899101876350983, 0.128492873943470309328063>,<0.000000000000000000000000, 0.816496580927725812415474, -0.018463387614029630467982>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{4.000000000000000000000000,<0.000000000000000000000000, 0.816496580927725812415474, -0.018463387614029630467982>,<-0.816496580927725701393172, 0.000000000000000000000000, -0.018463387614029644345770>,<-1.626353725899101876350983, 0.809857144971375952913206, 0.128492873943470253816912>,<-0.809857144971376397002416, 1.626353725899101876350983, 0.128492873943470309328063>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<-1.599903960749543774255699, 2.416400541677268698492753, 0.420015413536709347752662>,0.029999999999999998889777}sphere{<-2.416400541677269142581963, 1.599903960749542441988069, 0.420015413536709236730360>,0.029999999999999998889777}sphere{<-1.626353725899102098395588, 0.809857144971376063935509, 0.128492873943470087283458>,0.029999999999999998889777}sphere{<-0.809857144971376508024719, 1.626353725899101876350983, 0.128492873943469976261156>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<-1.599903960749543774255699, 2.416400541677268698492753, 0.420015413536709347752662>,<-2.416400541677269142581963, 1.599903960749542441988069, 0.420015413536709236730360>,0.010000000000000000208167}cylinder{<-2.416400541677269142581963, 1.599903960749542441988069, 0.420015413536709236730360>,<-1.626353725899102098395588, 0.809857144971376063935509, 0.128492873943470087283458>,0.010000000000000000208167}cylinder{<-1.626353725899102098395588, 0.809857144971376063935509, 0.128492873943470087283458>,<-0.809857144971376508024719, 1.626353725899101876350983, 0.128492873943469976261156>,0.010000000000000000208167}cylinder{<-0.809857144971376508024719, 1.626353725899101876350983, 0.128492873943469976261156>,<-1.599903960749543774255699, 2.416400541677268698492753, 0.420015413536709347752662>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{4.000000000000000000000000,<-1.599903960749543774255699, 2.416400541677268698492753, 0.420015413536709347752662>,<-2.416400541677269142581963, 1.599903960749542441988069, 0.420015413536709236730360>,<-1.626353725899102098395588, 0.809857144971376063935509, 0.128492873943470087283458>,<-0.809857144971376508024719, 1.626353725899101876350983, 0.128492873943469976261156>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<-1.626353725899102098395588, -0.809857144971375841890904, 0.128492873943470253816912>,0.029999999999999998889777}sphere{<-0.816496580927725812415474, 0.000000000000000000000000, -0.018463387614029630467982>,0.029999999999999998889777}sphere{<-0.000000000000000222044605, -0.816496580927725923437777, -0.018463387614029630467982>,0.029999999999999998889777}sphere{<-0.809857144971376397002416, -1.626353725899101876350983, 0.128492873943470281572488>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<-1.626353725899102098395588, -0.809857144971375841890904, 0.128492873943470253816912>,<-0.816496580927725812415474, 0.000000000000000000000000, -0.018463387614029630467982>,0.010000000000000000208167}cylinder{<-0.816496580927725812415474, 0.000000000000000000000000, -0.018463387614029630467982>,<-0.000000000000000222044605, -0.816496580927725923437777, -0.018463387614029630467982>,0.010000000000000000208167}cylinder{<-0.000000000000000222044605, -0.816496580927725923437777, -0.018463387614029630467982>,<-0.809857144971376397002416, -1.626353725899101876350983, 0.128492873943470281572488>,0.010000000000000000208167}cylinder{<-0.809857144971376397002416, -1.626353725899101876350983, 0.128492873943470281572488>,<-1.626353725899102098395588, -0.809857144971375841890904, 0.128492873943470253816912>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{4.000000000000000000000000,<-1.626353725899102098395588, -0.809857144971375841890904, 0.128492873943470253816912>,<-0.816496580927725812415474, 0.000000000000000000000000, -0.018463387614029630467982>,<-0.000000000000000222044605, -0.816496580927725923437777, -0.018463387614029630467982>,<-0.809857144971376397002416, -1.626353725899101876350983, 0.128492873943470281572488>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<0.013224882574779273092247, -2.422985988295698334127337, 0.274254143740090605696480>,0.029999999999999998889777}sphere{<0.823082027546155226005453, -1.613128843324322492236433, 0.127297882182590693656010>,0.029999999999999998889777}sphere{<0.000000000000000000000000, -0.816496580927726034460079, -0.018463387614028811678502>,0.029999999999999998889777}sphere{<-0.809857144971375952913206, -1.626353725899101654306378, 0.128492873943471058728605>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<0.013224882574779273092247, -2.422985988295698334127337, 0.274254143740090605696480>,<0.823082027546155226005453, -1.613128843324322492236433, 0.127297882182590693656010>,0.010000000000000000208167}cylinder{<0.823082027546155226005453, -1.613128843324322492236433, 0.127297882182590693656010>,<0.000000000000000000000000, -0.816496580927726034460079, -0.018463387614028811678502>,0.010000000000000000208167}cylinder{<0.000000000000000000000000, -0.816496580927726034460079, -0.018463387614028811678502>,<-0.809857144971375952913206, -1.626353725899101654306378, 0.128492873943471058728605>,0.010000000000000000208167}cylinder{<-0.809857144971375952913206, -1.626353725899101654306378, 0.128492873943471058728605>,<0.013224882574779273092247, -2.422985988295698334127337, 0.274254143740090605696480>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{4.000000000000000000000000,<0.013224882574779273092247, -2.422985988295698334127337, 0.274254143740090605696480>,<0.823082027546155226005453, -1.613128843324322492236433, 0.127297882182590693656010>,<0.000000000000000000000000, -0.816496580927726034460079, -0.018463387614028811678502>,<-0.809857144971375952913206, -1.626353725899101654306378, 0.128492873943471058728605>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<0.013224882574779606159154, -2.422985988295698334127337, 0.274254143740090605696480>,0.029999999999999998889777}sphere{<0.823082027546155448050058, -1.613128843324322714281038, 0.127297882182590749167161>,0.029999999999999998889777}sphere{<1.639363529163234689534079, -2.376940973263576672991348, 0.416449872680250576806316>,0.029999999999999998889777}sphere{<0.829506384191859069687780, -3.186798118234952070793042, 0.563406134237750544357937>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<0.013224882574779606159154, -2.422985988295698334127337, 0.274254143740090605696480>,<0.823082027546155448050058, -1.613128843324322714281038, 0.127297882182590749167161>,0.010000000000000000208167}cylinder{<0.823082027546155448050058, -1.613128843324322714281038, 0.127297882182590749167161>,<1.639363529163234689534079, -2.376940973263576672991348, 0.416449872680250576806316>,0.010000000000000000208167}cylinder{<1.639363529163234689534079, -2.376940973263576672991348, 0.416449872680250576806316>,<0.829506384191859069687780, -3.186798118234952070793042, 0.563406134237750544357937>,0.010000000000000000208167}cylinder{<0.829506384191859069687780, -3.186798118234952070793042, 0.563406134237750544357937>,<0.013224882574779606159154, -2.422985988295698334127337, 0.274254143740090605696480>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{4.000000000000000000000000,<0.013224882574779606159154, -2.422985988295698334127337, 0.274254143740090605696480>,<0.823082027546155448050058, -1.613128843324322714281038, 0.127297882182590749167161>,<1.639363529163234689534079, -2.376940973263576672991348, 0.416449872680250576806316>,<0.829506384191859069687780, -3.186798118234952070793042, 0.563406134237750544357937>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}}
