#version 3.7;
global_settings{assumed_gamma 1.0}
camera{perspective location <18.750000000000000000000000, 10.825317547305481014063844, 12.499999999999996447286321> right <0.133333333333333331482962, -0.230940107675850353485814, -0.000000000000000000000000> up <-0.064951905283832905846353, -0.037499999999999998612221, 0.129903810567665811692706> direction <0.750000000000000111022302, 0.433012701892219298294151, 0.499999999999999944488849> sky <-0.433012701892219298294151, -0.249999999999999944488849, 0.866025403784438707610605> look_at <18.000000000000000000000000, 10.392304845413262270881205, 11.999999999999996447286321>}
light_source{<18.750000000000000000000000, 10.825317547305481014063844, 12.499999999999996447286321> rgb<1.0,1.0,1.0>}
background{rgb<1,1,1>}
union{union{object{union{sphere{<0.000000000000000049995996, 0.816496580927725923437777, -0.576638311083890964248155>,0.029999999999999998889777}sphere{<-0.816496580927725923437777, 0.000000000000000099991992, -0.576638311083890964248155>,0.029999999999999998889777}sphere{<-0.000000000000000149987989, -0.816496580927725923437777, -0.576638311083890964248155>,0.029999999999999998889777}sphere{<0.816496580927725923437777, -0.000000000000000199983985, -0.576638311083890964248155>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<0.000000000000000049995996, 0.816496580927725923437777, -0.576638311083890964248155>,<-0.816496580927725923437777, 0.000000000000000099991992, -0.576638311083890964248155>,0.010000000000000000208167}cylinder{<-0.816496580927725923437777, 0.000000000000000099991992, -0.576638311083890964248155>,<-0.000000000000000149987989, -0.816496580927725923437777, -0.576638311083890964248155>,0.010000000000000000208167}cylinder{<-0.000000000000000149987989, -0.816496580927725923437777, -0.576638311083890964248155>,<0.816496580927725923437777, -0.000000000000000199983985, -0.576638311083890964248155>,0.010000000000000000208167}cylinder{<0.816496580927725923437777, -0.000000000000000199983985, -0.576638311083890964248155>,<0.000000000000000049995996, 0.816496580927725923437777, -0.576638311083890964248155>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{4.000000000000000000000000,<0.000000000000000049995996, 0.816496580927725923437777, -0.576638311083890964248155>,<-0.816496580927725923437777, 0.000000000000000099991992, -0.576638311083890964248155>,<-0.000000000000000149987989, -0.816496580927725923437777, -0.576638311083890964248155>,<0.816496580927725923437777, -0.000000000000000199983985, -0.576638311083890964248155>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<0.000000000000000000000000, 0.816496580927725812415474, -0.576638311083891075270458>,0.029999999999999998889777}sphere{<-0.816496580927725812415474, 0.000000000000000222044605, -0.576638311083890853225853>,0.029999999999999998889777}sphere{<-0.816497202497978635804543, 0.000000621570253100944825, 0.578062227295025987672261>,0.029999999999999998889777}sphere{<-0.000000621570252823389069, 0.816497202497978635804543, 0.578062227295025987672261>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<0.000000000000000000000000, 0.816496580927725812415474, -0.576638311083891075270458>,<-0.816496580927725812415474, 0.000000000000000222044605, -0.576638311083890853225853>,0.010000000000000000208167}cylinder{<-0.816496580927725812415474, 0.000000000000000222044605, -0.576638311083890853225853>,<-0.816497202497978635804543, 0.000000621570253100944825, 0.578062227295025987672261>,0.010000000000000000208167}cylinder{<-0.816497202497978635804543, 0.000000621570253100944825, 0.578062227295025987672261>,<-0.000000621570252823389069, 0.816497202497978635804543, 0.578062227295025987672261>,0.010000000000000000208167}cylinder{<-0.000000621570252823389069, 0.816497202497978635804543, 0.578062227295025987672261>,<0.000000000000000000000000, 0.816496580927725812415474, -0.576638311083891075270458>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{4.000000000000000000000000,<0.000000000000000000000000, 0.816496580927725812415474, -0.576638311083891075270458>,<-0.816496580927725812415474, 0.000000000000000222044605, -0.576638311083890853225853>,<-0.816497202497978635804543, 0.000000621570253100944825, 0.578062227295025987672261>,<-0.000000621570252823389069, 0.816497202497978635804543, 0.578062227295025987672261>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<0.816495959356526856964820, 0.000000621571198733406050, 0.578063985361188859002368>,0.029999999999999998889777}sphere{<-0.000000621571199177495259, -0.816495959356526856964820, 0.578063985361188859002368>,0.029999999999999998889777}sphere{<-0.816497202497978635804543, 0.000000621570253045433674, 0.578062227295025987672261>,0.029999999999999998889777}sphere{<-0.000000621570252823389069, 0.816497202497978635804543, 0.578062227295025987672261>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<0.816495959356526856964820, 0.000000621571198733406050, 0.578063985361188859002368>,<-0.000000621571199177495259, -0.816495959356526856964820, 0.578063985361188859002368>,0.010000000000000000208167}cylinder{<-0.000000621571199177495259, -0.816495959356526856964820, 0.578063985361188859002368>,<-0.816497202497978635804543, 0.000000621570253045433674, 0.578062227295025987672261>,0.010000000000000000208167}cylinder{<-0.816497202497978635804543, 0.000000621570253045433674, 0.578062227295025987672261>,<-0.000000621570252823389069, 0.816497202497978635804543, 0.578062227295025987672261>,0.010000000000000000208167}cylinder{<-0.000000621570252823389069, 0.816497202497978635804543, 0.578062227295025987672261>,<0.816495959356526856964820, 0.000000621571198733406050, 0.578063985361188859002368>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{4.000000000000000000000000,<0.816495959356526856964820, 0.000000621571198733406050, 0.578063985361188859002368>,<-0.000000621571199177495259, -0.816495959356526856964820, 0.578063985361188859002368>,<-0.816497202497978635804543, 0.000000621570253045433674, 0.578062227295025987672261>,<-0.000000621570252823389069, 0.816497202497978635804543, 0.578062227295025987672261>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<-0.816497202497978635804543, -0.000000621570252656855615, 0.578062227295025765627656>,0.029999999999999998889777}sphere{<-0.816496580927725812415474, 0.000000000000000000000000, -0.576638311083890853225853>,0.029999999999999998889777}sphere{<-0.000000000000000222044605, -0.816496580927726034460079, -0.576638311083890853225853>,0.029999999999999998889777}sphere{<-0.000000621570252934411371, -0.816497202497978635804543, 0.578062227295025765627656>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<-0.816497202497978635804543, -0.000000621570252656855615, 0.578062227295025765627656>,<-0.816496580927725812415474, 0.000000000000000000000000, -0.576638311083890853225853>,0.010000000000000000208167}cylinder{<-0.816496580927725812415474, 0.000000000000000000000000, -0.576638311083890853225853>,<-0.000000000000000222044605, -0.816496580927726034460079, -0.576638311083890853225853>,0.010000000000000000208167}cylinder{<-0.000000000000000222044605, -0.816496580927726034460079, -0.576638311083890853225853>,<-0.000000621570252934411371, -0.816497202497978635804543, 0.578062227295025765627656>,0.010000000000000000208167}cylinder{<-0.000000621570252934411371, -0.816497202497978635804543, 0.578062227295025765627656>,<-0.816497202497978635804543, -0.000000621570252656855615, 0.578062227295025765627656>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{4.000000000000000000000000,<-0.816497202497978635804543, -0.000000621570252656855615, 0.578062227295025765627656>,<-0.816496580927725812415474, 0.000000000000000000000000, -0.576638311083890853225853>,<-0.000000000000000222044605, -0.816496580927726034460079, -0.576638311083890853225853>,<-0.000000621570252934411371, -0.816497202497978635804543, 0.578062227295025765627656>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<0.816496580927252635362379, -0.000001243140978812069545, 0.578063106328107312315012>,0.029999999999999998889777}sphere{<0.816497202497505236706843, -0.000000621570726155213718, -0.576637432050809528583102>,0.029999999999999998889777}sphere{<-0.000000000000000222044737, -0.816496580927726034460079, -0.576638311083890964248155>,0.029999999999999998889777}sphere{<-0.000000621570252934411477, -0.816497202497978635804543, 0.578062227295025987672261>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<0.816496580927252635362379, -0.000001243140978812069545, 0.578063106328107312315012>,<0.816497202497505236706843, -0.000000621570726155213718, -0.576637432050809528583102>,0.010000000000000000208167}cylinder{<0.816497202497505236706843, -0.000000621570726155213718, -0.576637432050809528583102>,<-0.000000000000000222044737, -0.816496580927726034460079, -0.576638311083890964248155>,0.010000000000000000208167}cylinder{<-0.000000000000000222044737, -0.816496580927726034460079, -0.576638311083890964248155>,<-0.000000621570252934411477, -0.816497202497978635804543, 0.578062227295025987672261>,0.010000000000000000208167}cylinder{<-0.000000621570252934411477, -0.816497202497978635804543, 0.578062227295025987672261>,<0.816496580927252635362379, -0.000001243140978812069545, 0.578063106328107312315012>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{4.000000000000000000000000,<0.816496580927252635362379, -0.000001243140978812069545, 0.578063106328107312315012>,<0.816497202497505236706843, -0.000000621570726155213718, -0.576637432050809528583102>,<-0.000000000000000222044737, -0.816496580927726034460079, -0.576638311083890964248155>,<-0.000000621570252934411477, -0.816497202497978635804543, 0.578062227295025987672261>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<0.816496580927252746384681, -0.000001243140978812069545, 0.578063106328107534359617>,0.029999999999999998889777}sphere{<0.816497202497505347729145, -0.000000621570726155213612, -0.576637432050809639605404>,0.029999999999999998889777}sphere{<0.000001864711231466217408, 0.816497202496559215667560, -0.576637432049471376771521>,0.029999999999999998889777}sphere{<0.000001243140978753850642, 0.816496580926306503300793, 0.578063106329445797193500>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<0.816496580927252746384681, -0.000001243140978812069545, 0.578063106328107534359617>,<0.816497202497505347729145, -0.000000621570726155213612, -0.576637432050809639605404>,0.010000000000000000208167}cylinder{<0.816497202497505347729145, -0.000000621570726155213612, -0.576637432050809639605404>,<0.000001864711231466217408, 0.816497202496559215667560, -0.576637432049471376771521>,0.010000000000000000208167}cylinder{<0.000001864711231466217408, 0.816497202496559215667560, -0.576637432049471376771521>,<0.000001243140978753850642, 0.816496580926306503300793, 0.578063106329445797193500>,0.010000000000000000208167}cylinder{<0.000001243140978753850642, 0.816496580926306503300793, 0.578063106329445797193500>,<0.816496580927252746384681, -0.000001243140978812069545, 0.578063106328107534359617>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{4.000000000000000000000000,<0.816496580927252746384681, -0.000001243140978812069545, 0.578063106328107534359617>,<0.816497202497505347729145, -0.000000621570726155213612, -0.576637432050809639605404>,<0.000001864711231466217408, 0.816497202496559215667560, -0.576637432049471376771521>,<0.000001243140978753850642, 0.816496580926306503300793, 0.578063106329445797193500>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}}
