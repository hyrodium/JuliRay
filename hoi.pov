#version 3.7;
global_settings{assumed_gamma 1.0}
camera{perspective location <6.545084971874737256314347, 4.755282581475768211021204, 5.877852522924730926945358> right <0.293892626146236568551728, -0.404508497187473725631435, -0.000000000000000000000000> up <-0.237764129073788410551060, -0.172745751406263137184283, 0.404508497187473725631435> direction <0.654508497187473725631435, 0.475528258147576821102120, 0.587785252292473137103457> sky <-0.475528258147576821102120, -0.345491502812526274368565, 0.809016994374947451262869> look_at <5.890576474687263086593703, 4.279754323328191389919084, 5.290067270632257567797296>}
light_source{<6.545084971874737256314347, 4.755282581475768211021204, 5.877852522924730926945358> rgb<1.0,1.0,1.0>}
background{rgb<1,1,1>}
object{object{union{union{sphere{<1.000000000000000000000000, 1.000000000000000000000000, 1.000000000000000000000000>,0.050000000000000002775558}sphere{<-1.000000000000000000000000, 1.000000000000000000000000, 1.000000000000000000000000>,0.050000000000000002775558}sphere{<1.000000000000000000000000, -1.000000000000000000000000, 1.000000000000000000000000>,0.050000000000000002775558}sphere{<-1.000000000000000000000000, -1.000000000000000000000000, 1.000000000000000000000000>,0.050000000000000002775558}sphere{<1.000000000000000000000000, 1.000000000000000000000000, -1.000000000000000000000000>,0.050000000000000002775558}sphere{<-1.000000000000000000000000, 1.000000000000000000000000, -1.000000000000000000000000>,0.050000000000000002775558}sphere{<1.000000000000000000000000, -1.000000000000000000000000, -1.000000000000000000000000>,0.050000000000000002775558}sphere{<-1.000000000000000000000000, -1.000000000000000000000000, -1.000000000000000000000000>,0.050000000000000002775558}}union{cylinder{<1.000000000000000000000000, 1.000000000000000000000000, 1.000000000000000000000000>,<1.000000000000000000000000, 1.000000000000000000000000, -1.000000000000000000000000>,0.010000000000000000208167}cylinder{<-1.000000000000000000000000, 1.000000000000000000000000, 1.000000000000000000000000>,<-1.000000000000000000000000, 1.000000000000000000000000, -1.000000000000000000000000>,0.010000000000000000208167}cylinder{<1.000000000000000000000000, -1.000000000000000000000000, 1.000000000000000000000000>,<1.000000000000000000000000, -1.000000000000000000000000, -1.000000000000000000000000>,0.010000000000000000208167}cylinder{<-1.000000000000000000000000, -1.000000000000000000000000, 1.000000000000000000000000>,<-1.000000000000000000000000, -1.000000000000000000000000, -1.000000000000000000000000>,0.010000000000000000208167}cylinder{<1.000000000000000000000000, 1.000000000000000000000000, 1.000000000000000000000000>,<-1.000000000000000000000000, 1.000000000000000000000000, 1.000000000000000000000000>,0.010000000000000000208167}cylinder{<1.000000000000000000000000, -1.000000000000000000000000, 1.000000000000000000000000>,<-1.000000000000000000000000, -1.000000000000000000000000, 1.000000000000000000000000>,0.010000000000000000208167}cylinder{<1.000000000000000000000000, 1.000000000000000000000000, -1.000000000000000000000000>,<-1.000000000000000000000000, 1.000000000000000000000000, -1.000000000000000000000000>,0.010000000000000000208167}cylinder{<1.000000000000000000000000, -1.000000000000000000000000, -1.000000000000000000000000>,<-1.000000000000000000000000, -1.000000000000000000000000, -1.000000000000000000000000>,0.010000000000000000208167}cylinder{<1.000000000000000000000000, 1.000000000000000000000000, 1.000000000000000000000000>,<1.000000000000000000000000, -1.000000000000000000000000, 1.000000000000000000000000>,0.010000000000000000208167}cylinder{<1.000000000000000000000000, 1.000000000000000000000000, -1.000000000000000000000000>,<1.000000000000000000000000, -1.000000000000000000000000, -1.000000000000000000000000>,0.010000000000000000208167}cylinder{<-1.000000000000000000000000, 1.000000000000000000000000, 1.000000000000000000000000>,<-1.000000000000000000000000, -1.000000000000000000000000, 1.000000000000000000000000>,0.010000000000000000208167}cylinder{<-1.000000000000000000000000, 1.000000000000000000000000, -1.000000000000000000000000>,<-1.000000000000000000000000, -1.000000000000000000000000, -1.000000000000000000000000>,0.010000000000000000208167}}union{polygon{4.000000000000000000000000,<1.000000000000000000000000, 1.000000000000000000000000, 1.000000000000000000000000>,<1.000000000000000000000000, -1.000000000000000000000000, 1.000000000000000000000000>,<-1.000000000000000000000000, -1.000000000000000000000000, 1.000000000000000000000000>,<-1.000000000000000000000000, 1.000000000000000000000000, 1.000000000000000000000000>}polygon{4.000000000000000000000000,<1.000000000000000000000000, 1.000000000000000000000000, -1.000000000000000000000000>,<1.000000000000000000000000, -1.000000000000000000000000, -1.000000000000000000000000>,<-1.000000000000000000000000, -1.000000000000000000000000, -1.000000000000000000000000>,<-1.000000000000000000000000, 1.000000000000000000000000, -1.000000000000000000000000>}polygon{4.000000000000000000000000,<1.000000000000000000000000, 1.000000000000000000000000, 1.000000000000000000000000>,<1.000000000000000000000000, 1.000000000000000000000000, -1.000000000000000000000000>,<1.000000000000000000000000, -1.000000000000000000000000, -1.000000000000000000000000>,<1.000000000000000000000000, -1.000000000000000000000000, 1.000000000000000000000000>}polygon{4.000000000000000000000000,<-1.000000000000000000000000, 1.000000000000000000000000, 1.000000000000000000000000>,<-1.000000000000000000000000, 1.000000000000000000000000, -1.000000000000000000000000>,<-1.000000000000000000000000, -1.000000000000000000000000, -1.000000000000000000000000>,<-1.000000000000000000000000, -1.000000000000000000000000, 1.000000000000000000000000>}polygon{4.000000000000000000000000,<1.000000000000000000000000, 1.000000000000000000000000, 1.000000000000000000000000>,<1.000000000000000000000000, 1.000000000000000000000000, -1.000000000000000000000000>,<-1.000000000000000000000000, 1.000000000000000000000000, -1.000000000000000000000000>,<-1.000000000000000000000000, 1.000000000000000000000000, 1.000000000000000000000000>}polygon{4.000000000000000000000000,<1.000000000000000000000000, -1.000000000000000000000000, 1.000000000000000000000000>,<1.000000000000000000000000, -1.000000000000000000000000, -1.000000000000000000000000>,<-1.000000000000000000000000, -1.000000000000000000000000, -1.000000000000000000000000>,<-1.000000000000000000000000, -1.000000000000000000000000, 1.000000000000000000000000>}}} pigment{rgbft<1.0,0.0,0.0,0.2,0.2>}} scale 0.500000000000000000000000}