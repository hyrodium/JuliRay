#version 3.7;
global_settings{assumed_gamma 1.0}
camera{perspective location <18.750000000000000000000000, 10.825317547305481014063844, 12.499999999999996447286321> right <0.133333333333333331482962, -0.230940107675850353485814, -0.000000000000000000000000> up <-0.064951905283832905846353, -0.037499999999999998612221, 0.129903810567665811692706> direction <0.750000000000000111022302, 0.433012701892219298294151, 0.499999999999999944488849> sky <-0.433012701892219298294151, -0.249999999999999944488849, 0.866025403784438707610605> look_at <18.000000000000000000000000, 10.392304845413262270881205, 11.999999999999996447286321>}
light_source{<18.750000000000000000000000, 10.825317547305481014063844, 12.499999999999996447286321> rgb<1.0,1.0,1.0>}
background{rgb<1,1,1>}
union{union{object{union{sphere{<-0.303530999103342913336689, 0.525731112119133481286326, 0.000000000000000000000000>,0.029999999999999998889777}sphere{<-0.303530999103343301914748, -0.525731112119133370264024, 0.000000000000000000000000>,0.029999999999999998889777}sphere{<0.607061998206686048717984, -0.000000000000000148687307, 0.000000000000000000000000>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<-0.303530999103342913336689, 0.525731112119133481286326, 0.000000000000000000000000>,<-0.303530999103343301914748, -0.525731112119133370264024, 0.000000000000000000000000>,0.010000000000000000208167}cylinder{<-0.303530999103343301914748, -0.525731112119133370264024, 0.000000000000000000000000>,<0.607061998206686048717984, -0.000000000000000148687307, 0.000000000000000000000000>,0.010000000000000000208167}cylinder{<0.607061998206686048717984, -0.000000000000000148687307, 0.000000000000000000000000>,<-0.303530999103342913336689, 0.525731112119133481286326, 0.000000000000000000000000>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{3.000000000000000000000000,<-0.303530999103342913336689, 0.525731112119133481286326, 0.000000000000000000000000>,<-0.303530999103343301914748, -0.525731112119133370264024, 0.000000000000000000000000>,<0.607061998206686048717984, -0.000000000000000148687307, 0.000000000000000000000000>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<-0.303530999103342913336689, 0.525731112119133481286326, -0.000000000000000011282802>,0.029999999999999998889777}sphere{<-0.303530999103343412937051, -0.525731112119133370264024, -0.000000000000000011282802>,0.029999999999999998889777}sphere{<-1.214123996413372319480573, 0.000000000000000593715361, 0.000000000000000022565603>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<-0.303530999103342913336689, 0.525731112119133481286326, -0.000000000000000011282802>,<-0.303530999103343412937051, -0.525731112119133370264024, -0.000000000000000011282802>,0.010000000000000000208167}cylinder{<-0.303530999103343412937051, -0.525731112119133370264024, -0.000000000000000011282802>,<-1.214123996413372319480573, 0.000000000000000593715361, 0.000000000000000022565603>,0.010000000000000000208167}cylinder{<-1.214123996413372319480573, 0.000000000000000593715361, 0.000000000000000022565603>,<-0.303530999103342913336689, 0.525731112119133481286326, -0.000000000000000011282802>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{3.000000000000000000000000,<-0.303530999103342913336689, 0.525731112119133481286326, -0.000000000000000011282802>,<-0.303530999103343412937051, -0.525731112119133370264024, -0.000000000000000011282802>,<-1.214123996413372319480573, 0.000000000000000593715361, 0.000000000000000022565603>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<-0.303530999103342857825538, 0.525731112119133592308629, -0.000000000000000056414008>,0.029999999999999998889777}sphere{<-1.214123996413371653346758, 1.051462224238267850751072, 0.000000000000000011282802>,0.029999999999999998889777}sphere{<-1.214123996413372319480573, 0.000000000000000666133815, -0.000000000000000022565603>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<-0.303530999103342857825538, 0.525731112119133592308629, -0.000000000000000056414008>,<-1.214123996413371653346758, 1.051462224238267850751072, 0.000000000000000011282802>,0.010000000000000000208167}cylinder{<-1.214123996413371653346758, 1.051462224238267850751072, 0.000000000000000011282802>,<-1.214123996413372319480573, 0.000000000000000666133815, -0.000000000000000022565603>,0.010000000000000000208167}cylinder{<-1.214123996413372319480573, 0.000000000000000666133815, -0.000000000000000022565603>,<-0.303530999103342857825538, 0.525731112119133592308629, -0.000000000000000056414008>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{3.000000000000000000000000,<-0.303530999103342857825538, 0.525731112119133592308629, -0.000000000000000056414008>,<-1.214123996413371653346758, 1.051462224238267850751072, 0.000000000000000011282802>,<-1.214123996413372319480573, 0.000000000000000666133815, -0.000000000000000022565603>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<-2.124716993723401170512943, 0.525731112119134813553956, -0.000000000000000101545215>,0.029999999999999998889777}sphere{<-1.214123996413371431302153, 1.051462224238267850751072, -0.000000000000000169242025>,0.029999999999999998889777}sphere{<-1.214123996413372319480573, 0.000000000000000634547970, -0.000000000000000203090430>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<-2.124716993723401170512943, 0.525731112119134813553956, -0.000000000000000101545215>,<-1.214123996413371431302153, 1.051462224238267850751072, -0.000000000000000169242025>,0.010000000000000000208167}cylinder{<-1.214123996413371431302153, 1.051462224238267850751072, -0.000000000000000169242025>,<-1.214123996413372319480573, 0.000000000000000634547970, -0.000000000000000203090430>,0.010000000000000000208167}cylinder{<-1.214123996413372319480573, 0.000000000000000634547970, -0.000000000000000203090430>,<-2.124716993723401170512943, 0.525731112119134813553956, -0.000000000000000101545215>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{3.000000000000000000000000,<-2.124716993723401170512943, 0.525731112119134813553956, -0.000000000000000101545215>,<-1.214123996413371431302153, 1.051462224238267850751072, -0.000000000000000169242025>,<-1.214123996413372319480573, 0.000000000000000634547970, -0.000000000000000203090430>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<-2.124716993723401614602153, 0.525731112119134813553956, -0.000000000000000383615257>,0.029999999999999998889777}sphere{<-1.214123996413371653346758, 1.051462224238267850751072, -0.000000000000000451312067>,0.029999999999999998889777}sphere{<-2.124716993723400726423733, 1.577193336357401998171213, -0.000000000000000315918447>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<-2.124716993723401614602153, 0.525731112119134813553956, -0.000000000000000383615257>,<-1.214123996413371653346758, 1.051462224238267850751072, -0.000000000000000451312067>,0.010000000000000000208167}cylinder{<-1.214123996413371653346758, 1.051462224238267850751072, -0.000000000000000451312067>,<-2.124716993723400726423733, 1.577193336357401998171213, -0.000000000000000315918447>,0.010000000000000000208167}cylinder{<-2.124716993723400726423733, 1.577193336357401998171213, -0.000000000000000315918447>,<-2.124716993723401614602153, 0.525731112119134813553956, -0.000000000000000383615257>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{3.000000000000000000000000,<-2.124716993723401614602153, 0.525731112119134813553956, -0.000000000000000383615257>,<-1.214123996413371653346758, 1.051462224238267850751072, -0.000000000000000451312067>,<-2.124716993723400726423733, 1.577193336357401998171213, -0.000000000000000315918447>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<-2.124716993723402502780573, 0.525731112119134813553956, -0.000000000000000936472540>,0.029999999999999998889777}sphere{<-3.035309991033431131768339, 1.051462224238269627107911, -0.000000000000000767230514>,0.029999999999999998889777}sphere{<-2.124716993723401170512943, 1.577193336357401998171213, -0.000000000000000868775729>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<-2.124716993723402502780573, 0.525731112119134813553956, -0.000000000000000936472540>,<-3.035309991033431131768339, 1.051462224238269627107911, -0.000000000000000767230514>,0.010000000000000000208167}cylinder{<-3.035309991033431131768339, 1.051462224238269627107911, -0.000000000000000767230514>,<-2.124716993723401170512943, 1.577193336357401998171213, -0.000000000000000868775729>,0.010000000000000000208167}cylinder{<-2.124716993723401170512943, 1.577193336357401998171213, -0.000000000000000868775729>,<-2.124716993723402502780573, 0.525731112119134813553956, -0.000000000000000936472540>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{3.000000000000000000000000,<-2.124716993723402502780573, 0.525731112119134813553956, -0.000000000000000936472540>,<-3.035309991033431131768339, 1.051462224238269627107911, -0.000000000000000767230514>,<-2.124716993723401170512943, 1.577193336357401998171213, -0.000000000000000868775729>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<0.607061998206685715651076, -1.051462224238266962572652, 0.000000000000000022565603>,0.029999999999999998889777}sphere{<-0.303530999103343301914748, -0.525731112119133481286326, -0.000000000000000011282802>,0.029999999999999998889777}sphere{<0.607061998206686048717984, 0.000000000000000000000000, -0.000000000000000011282802>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<0.607061998206685715651076, -1.051462224238266962572652, 0.000000000000000022565603>,<-0.303530999103343301914748, -0.525731112119133481286326, -0.000000000000000011282802>,0.010000000000000000208167}cylinder{<-0.303530999103343301914748, -0.525731112119133481286326, -0.000000000000000011282802>,<0.607061998206686048717984, 0.000000000000000000000000, -0.000000000000000011282802>,0.010000000000000000208167}cylinder{<0.607061998206686048717984, 0.000000000000000000000000, -0.000000000000000011282802>,<0.607061998206685715651076, -1.051462224238266962572652, 0.000000000000000022565603>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{3.000000000000000000000000,<0.607061998206685715651076, -1.051462224238266962572652, 0.000000000000000022565603>,<-0.303530999103343301914748, -0.525731112119133481286326, -0.000000000000000011282802>,<0.607061998206686048717984, 0.000000000000000000000000, -0.000000000000000011282802>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<0.607061998206685604628774, -1.051462224238266962572652, -0.000000000000000022565603>,0.029999999999999998889777}sphere{<1.517654995516715121794959, -0.525731112119134147420141, 0.000000000000000011282802>,0.029999999999999998889777}sphere{<0.607061998206686048717984, -0.000000000000000019859114, -0.000000000000000056414008>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<0.607061998206685604628774, -1.051462224238266962572652, -0.000000000000000022565603>,<1.517654995516715121794959, -0.525731112119134147420141, 0.000000000000000011282802>,0.010000000000000000208167}cylinder{<1.517654995516715121794959, -0.525731112119134147420141, 0.000000000000000011282802>,<0.607061998206686048717984, -0.000000000000000019859114, -0.000000000000000056414008>,0.010000000000000000208167}cylinder{<0.607061998206686048717984, -0.000000000000000019859114, -0.000000000000000056414008>,<0.607061998206685604628774, -1.051462224238266962572652, -0.000000000000000022565603>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{3.000000000000000000000000,<0.607061998206685604628774, -1.051462224238266962572652, -0.000000000000000022565603>,<1.517654995516715121794959, -0.525731112119134147420141, 0.000000000000000011282802>,<0.607061998206686048717984, -0.000000000000000019859114, -0.000000000000000056414008>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<0.607061998206685604628774, -1.051462224238267406661862, -0.000000000000000203090430>,0.029999999999999998889777}sphere{<1.517654995516715121794959, -0.525731112119134369464746, -0.000000000000000169242025>,0.029999999999999998889777}sphere{<1.517654995516714233616540, -1.577193336357401776126608, -0.000000000000000101545215>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<0.607061998206685604628774, -1.051462224238267406661862, -0.000000000000000203090430>,<1.517654995516715121794959, -0.525731112119134369464746, -0.000000000000000169242025>,0.010000000000000000208167}cylinder{<1.517654995516715121794959, -0.525731112119134369464746, -0.000000000000000169242025>,<1.517654995516714233616540, -1.577193336357401776126608, -0.000000000000000101545215>,0.010000000000000000208167}cylinder{<1.517654995516714233616540, -1.577193336357401776126608, -0.000000000000000101545215>,<0.607061998206685604628774, -1.051462224238267406661862, -0.000000000000000203090430>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{3.000000000000000000000000,<0.607061998206685604628774, -1.051462224238267406661862, -0.000000000000000203090430>,<1.517654995516715121794959, -0.525731112119134369464746, -0.000000000000000169242025>,<1.517654995516714233616540, -1.577193336357401776126608, -0.000000000000000101545215>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<2.428247992826744194871935, -1.051462224238268738929492, -0.000000000000000315918447>,0.029999999999999998889777}sphere{<1.517654995516715343839564, -0.525731112119134369464746, -0.000000000000000451312067>,0.029999999999999998889777}sphere{<1.517654995516714677705750, -1.577193336357401776126608, -0.000000000000000383615257>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<2.428247992826744194871935, -1.051462224238268738929492, -0.000000000000000315918447>,<1.517654995516715343839564, -0.525731112119134369464746, -0.000000000000000451312067>,0.010000000000000000208167}cylinder{<1.517654995516715343839564, -0.525731112119134369464746, -0.000000000000000451312067>,<1.517654995516714677705750, -1.577193336357401776126608, -0.000000000000000383615257>,0.010000000000000000208167}cylinder{<1.517654995516714677705750, -1.577193336357401776126608, -0.000000000000000383615257>,<2.428247992826744194871935, -1.051462224238268738929492, -0.000000000000000315918447>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{3.000000000000000000000000,<2.428247992826744194871935, -1.051462224238268738929492, -0.000000000000000315918447>,<1.517654995516715343839564, -0.525731112119134369464746, -0.000000000000000451312067>,<1.517654995516714677705750, -1.577193336357401776126608, -0.000000000000000383615257>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<-0.303530999103342913336689, 0.525731112119133481286326, -0.000000000000000011282802>,0.029999999999999998889777}sphere{<0.607061998206686159740286, 1.051462224238266962572652, 0.000000000000000022565603>,0.029999999999999998889777}sphere{<0.607061998206686048717984, -0.000000000000000111022302, -0.000000000000000011282802>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<-0.303530999103342913336689, 0.525731112119133481286326, -0.000000000000000011282802>,<0.607061998206686159740286, 1.051462224238266962572652, 0.000000000000000022565603>,0.010000000000000000208167}cylinder{<0.607061998206686159740286, 1.051462224238266962572652, 0.000000000000000022565603>,<0.607061998206686048717984, -0.000000000000000111022302, -0.000000000000000011282802>,0.010000000000000000208167}cylinder{<0.607061998206686048717984, -0.000000000000000111022302, -0.000000000000000011282802>,<-0.303530999103342913336689, 0.525731112119133481286326, -0.000000000000000011282802>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{3.000000000000000000000000,<-0.303530999103342913336689, 0.525731112119133481286326, -0.000000000000000011282802>,<0.607061998206686159740286, 1.051462224238266962572652, 0.000000000000000022565603>,<0.607061998206686048717984, -0.000000000000000111022302, -0.000000000000000011282802>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<-1.214123996413372985614387, -1.051462224238266518483442, 0.000000000000000011282802>,0.029999999999999998889777}sphere{<-0.303530999103343357425899, -0.525731112119133703330931, -0.000000000000000056414008>,0.029999999999999998889777}sphere{<-1.214123996413372541525177, 0.000000000000000666133815, -0.000000000000000022565603>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<-1.214123996413372985614387, -1.051462224238266518483442, 0.000000000000000011282802>,<-0.303530999103343357425899, -0.525731112119133703330931, -0.000000000000000056414008>,0.010000000000000000208167}cylinder{<-0.303530999103343357425899, -0.525731112119133703330931, -0.000000000000000056414008>,<-1.214123996413372541525177, 0.000000000000000666133815, -0.000000000000000022565603>,0.010000000000000000208167}cylinder{<-1.214123996413372541525177, 0.000000000000000666133815, -0.000000000000000022565603>,<-1.214123996413372985614387, -1.051462224238266518483442, 0.000000000000000011282802>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{3.000000000000000000000000,<-1.214123996413372985614387, -1.051462224238266518483442, 0.000000000000000011282802>,<-0.303530999103343357425899, -0.525731112119133703330931, -0.000000000000000056414008>,<-1.214123996413372541525177, 0.000000000000000666133815, -0.000000000000000022565603>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<-0.303530999103342691292085, 0.525731112119133703330931, -0.000000000000000067696810>,0.029999999999999998889777}sphere{<-1.214123996413371653346758, 1.051462224238267850751072, 0.000000000000000000000000>,0.029999999999999998889777}sphere{<-0.303530999103341914135967, 1.577193336357400887948188, -0.000000000000000000000000>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<-0.303530999103342691292085, 0.525731112119133703330931, -0.000000000000000067696810>,<-1.214123996413371653346758, 1.051462224238267850751072, 0.000000000000000000000000>,0.010000000000000000208167}cylinder{<-1.214123996413371653346758, 1.051462224238267850751072, 0.000000000000000000000000>,<-0.303530999103341914135967, 1.577193336357400887948188, -0.000000000000000000000000>,0.010000000000000000208167}cylinder{<-0.303530999103341914135967, 1.577193336357400887948188, -0.000000000000000000000000>,<-0.303530999103342691292085, 0.525731112119133703330931, -0.000000000000000067696810>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{3.000000000000000000000000,<-0.303530999103342691292085, 0.525731112119133703330931, -0.000000000000000067696810>,<-1.214123996413371653346758, 1.051462224238267850751072, 0.000000000000000000000000>,<-0.303530999103341914135967, 1.577193336357400887948188, -0.000000000000000000000000>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<-2.124716993723401170512943, 0.525731112119135035598561, -0.000000000000000146676422>,0.029999999999999998889777}sphere{<-2.124716993723402058691363, -0.525731112119132482085604, -0.000000000000000146676422>,0.029999999999999998889777}sphere{<-1.214123996413372541525177, 0.000000000000000666133815, -0.000000000000000248221637>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<-2.124716993723401170512943, 0.525731112119135035598561, -0.000000000000000146676422>,<-2.124716993723402058691363, -0.525731112119132482085604, -0.000000000000000146676422>,0.010000000000000000208167}cylinder{<-2.124716993723402058691363, -0.525731112119132482085604, -0.000000000000000146676422>,<-1.214123996413372541525177, 0.000000000000000666133815, -0.000000000000000248221637>,0.010000000000000000208167}cylinder{<-1.214123996413372541525177, 0.000000000000000666133815, -0.000000000000000248221637>,<-2.124716993723401170512943, 0.525731112119135035598561, -0.000000000000000146676422>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{3.000000000000000000000000,<-2.124716993723401170512943, 0.525731112119135035598561, -0.000000000000000146676422>,<-2.124716993723402058691363, -0.525731112119132482085604, -0.000000000000000146676422>,<-1.214123996413372541525177, 0.000000000000000666133815, -0.000000000000000248221637>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<-1.214123996413371653346758, 2.102924448476535257412934, -0.000000000000000361049654>,0.029999999999999998889777}sphere{<-1.214123996413371653346758, 1.051462224238267850751072, -0.000000000000000462594869>,0.029999999999999998889777}sphere{<-2.124716993723400726423733, 1.577193336357401776126608, -0.000000000000000327201249>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<-1.214123996413371653346758, 2.102924448476535257412934, -0.000000000000000361049654>,<-1.214123996413371653346758, 1.051462224238267850751072, -0.000000000000000462594869>,0.010000000000000000208167}cylinder{<-1.214123996413371653346758, 1.051462224238267850751072, -0.000000000000000462594869>,<-2.124716993723400726423733, 1.577193336357401776126608, -0.000000000000000327201249>,0.010000000000000000208167}cylinder{<-2.124716993723400726423733, 1.577193336357401776126608, -0.000000000000000327201249>,<-1.214123996413371653346758, 2.102924448476535257412934, -0.000000000000000361049654>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{3.000000000000000000000000,<-1.214123996413371653346758, 2.102924448476535257412934, -0.000000000000000361049654>,<-1.214123996413371653346758, 1.051462224238267850751072, -0.000000000000000462594869>,<-2.124716993723400726423733, 1.577193336357401776126608, -0.000000000000000327201249>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<-2.124716993723402502780573, 0.525731112119134813553956, -0.000000000000000981603746>,0.029999999999999998889777}sphere{<-3.035309991033431131768339, 1.051462224238269627107911, -0.000000000000000812361721>,0.029999999999999998889777}sphere{<-3.035309991033432464035968, 0.000000000000002664535259, -0.000000000000000846210126>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<-2.124716993723402502780573, 0.525731112119134813553956, -0.000000000000000981603746>,<-3.035309991033431131768339, 1.051462224238269627107911, -0.000000000000000812361721>,0.010000000000000000208167}cylinder{<-3.035309991033431131768339, 1.051462224238269627107911, -0.000000000000000812361721>,<-3.035309991033432464035968, 0.000000000000002664535259, -0.000000000000000846210126>,0.010000000000000000208167}cylinder{<-3.035309991033432464035968, 0.000000000000002664535259, -0.000000000000000846210126>,<-2.124716993723402502780573, 0.525731112119134813553956, -0.000000000000000981603746>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{3.000000000000000000000000,<-2.124716993723402502780573, 0.525731112119134813553956, -0.000000000000000981603746>,<-3.035309991033431131768339, 1.051462224238269627107911, -0.000000000000000812361721>,<-3.035309991033432464035968, 0.000000000000002664535259, -0.000000000000000846210126>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<0.607061998206685937695681, -1.051462224238266962572652, -0.000000000000000022565603>,0.029999999999999998889777}sphere{<-0.303530999103343357425899, -0.525731112119133370264024, -0.000000000000000056414008>,0.029999999999999998889777}sphere{<-0.303530999103343357425899, -1.577193336357400887948188, 0.000000000000000011282802>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<0.607061998206685937695681, -1.051462224238266962572652, -0.000000000000000022565603>,<-0.303530999103343357425899, -0.525731112119133370264024, -0.000000000000000056414008>,0.010000000000000000208167}cylinder{<-0.303530999103343357425899, -0.525731112119133370264024, -0.000000000000000056414008>,<-0.303530999103343357425899, -1.577193336357400887948188, 0.000000000000000011282802>,0.010000000000000000208167}cylinder{<-0.303530999103343357425899, -1.577193336357400887948188, 0.000000000000000011282802>,<0.607061998206685937695681, -1.051462224238266962572652, -0.000000000000000022565603>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{3.000000000000000000000000,<0.607061998206685937695681, -1.051462224238266962572652, -0.000000000000000022565603>,<-0.303530999103343357425899, -0.525731112119133370264024, -0.000000000000000056414008>,<-0.303530999103343357425899, -1.577193336357400887948188, 0.000000000000000011282802>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<1.517654995516715787928774, 0.525731112119132926174814, 0.000000000000000000000000>,0.029999999999999998889777}sphere{<1.517654995516715343839564, -0.525731112119134369464746, -0.000000000000000000000000>,0.029999999999999998889777}sphere{<0.607061998206686048717984, 0.000000000000000000000000, -0.000000000000000067696810>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<1.517654995516715787928774, 0.525731112119132926174814, 0.000000000000000000000000>,<1.517654995516715343839564, -0.525731112119134369464746, -0.000000000000000000000000>,0.010000000000000000208167}cylinder{<1.517654995516715343839564, -0.525731112119134369464746, -0.000000000000000000000000>,<0.607061998206686048717984, 0.000000000000000000000000, -0.000000000000000067696810>,0.010000000000000000208167}cylinder{<0.607061998206686048717984, 0.000000000000000000000000, -0.000000000000000067696810>,<1.517654995516715787928774, 0.525731112119132926174814, 0.000000000000000000000000>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{3.000000000000000000000000,<1.517654995516715787928774, 0.525731112119132926174814, 0.000000000000000000000000>,<1.517654995516715343839564, -0.525731112119134369464746, -0.000000000000000000000000>,<0.607061998206686048717984, 0.000000000000000000000000, -0.000000000000000067696810>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<0.607061998206685715651076, -1.051462224238267406661862, -0.000000000000000248221637>,0.029999999999999998889777}sphere{<0.607061998206684605428052, -2.102924448476534813323724, -0.000000000000000146676422>,0.029999999999999998889777}sphere{<1.517654995516714677705750, -1.577193336357401332037398, -0.000000000000000146676422>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<0.607061998206685715651076, -1.051462224238267406661862, -0.000000000000000248221637>,<0.607061998206684605428052, -2.102924448476534813323724, -0.000000000000000146676422>,0.010000000000000000208167}cylinder{<0.607061998206684605428052, -2.102924448476534813323724, -0.000000000000000146676422>,<1.517654995516714677705750, -1.577193336357401332037398, -0.000000000000000146676422>,0.010000000000000000208167}cylinder{<1.517654995516714677705750, -1.577193336357401332037398, -0.000000000000000146676422>,<0.607061998206685715651076, -1.051462224238267406661862, -0.000000000000000248221637>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{3.000000000000000000000000,<0.607061998206685715651076, -1.051462224238267406661862, -0.000000000000000248221637>,<0.607061998206684605428052, -2.102924448476534813323724, -0.000000000000000146676422>,<1.517654995516714677705750, -1.577193336357401332037398, -0.000000000000000146676422>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<2.428247992826744194871935, -1.051462224238268516884887, -0.000000000000000327201249>,0.029999999999999998889777}sphere{<1.517654995516715343839564, -0.525731112119134147420141, -0.000000000000000462594869>,0.029999999999999998889777}sphere{<2.428247992826745083050355, -0.000000000000001776356839, -0.000000000000000361049654>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<2.428247992826744194871935, -1.051462224238268516884887, -0.000000000000000327201249>,<1.517654995516715343839564, -0.525731112119134147420141, -0.000000000000000462594869>,0.010000000000000000208167}cylinder{<1.517654995516715343839564, -0.525731112119134147420141, -0.000000000000000462594869>,<2.428247992826745083050355, -0.000000000000001776356839, -0.000000000000000361049654>,0.010000000000000000208167}cylinder{<2.428247992826745083050355, -0.000000000000001776356839, -0.000000000000000361049654>,<2.428247992826744194871935, -1.051462224238268516884887, -0.000000000000000327201249>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{3.000000000000000000000000,<2.428247992826744194871935, -1.051462224238268516884887, -0.000000000000000327201249>,<1.517654995516715343839564, -0.525731112119134147420141, -0.000000000000000462594869>,<2.428247992826745083050355, -0.000000000000001776356839, -0.000000000000000361049654>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}}
