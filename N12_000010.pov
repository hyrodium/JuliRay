#version 3.7;
global_settings{assumed_gamma 1.0}
camera{perspective location <18.750000000000000000000000, 10.825317547305481014063844, 12.499999999999996447286321> right <0.133333333333333331482962, -0.230940107675850353485814, -0.000000000000000000000000> up <-0.064951905283832905846353, -0.037499999999999998612221, 0.129903810567665811692706> direction <0.750000000000000111022302, 0.433012701892219298294151, 0.499999999999999944488849> sky <-0.433012701892219298294151, -0.249999999999999944488849, 0.866025403784438707610605> look_at <18.000000000000000000000000, 10.392304845413262270881205, 11.999999999999996447286321>}
light_source{<18.750000000000000000000000, 10.825317547305481014063844, 12.499999999999996447286321> rgb<1.0,1.0,1.0>}
background{rgb<1,1,1>}
union{union{object{union{sphere{<0.187592474085079868872938, 0.577350269189625620036566, -0.576996638683937024261184>,0.029999999999999998889777}sphere{<-0.491123473188422809965203, 0.356822089773089878850243, -0.576996638683937024261184>,0.029999999999999998889777}sphere{<-0.491123473188422865476355, -0.356822089773089767827940, -0.576996638683937024261184>,0.029999999999999998889777}sphere{<0.187592474085079730095060, -0.577350269189625620036566, -0.576996638683937024261184>,0.029999999999999998889777}sphere{<0.607061998206686048717984, -0.000000000000000148687307, -0.576996638683937024261184>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<0.187592474085079868872938, 0.577350269189625620036566, -0.576996638683937024261184>,<-0.491123473188422809965203, 0.356822089773089878850243, -0.576996638683937024261184>,0.010000000000000000208167}cylinder{<-0.491123473188422809965203, 0.356822089773089878850243, -0.576996638683937024261184>,<-0.491123473188422865476355, -0.356822089773089767827940, -0.576996638683937024261184>,0.010000000000000000208167}cylinder{<-0.491123473188422865476355, -0.356822089773089767827940, -0.576996638683937024261184>,<0.187592474085079730095060, -0.577350269189625620036566, -0.576996638683937024261184>,0.010000000000000000208167}cylinder{<0.187592474085079730095060, -0.577350269189625620036566, -0.576996638683937024261184>,<0.607061998206686048717984, -0.000000000000000148687307, -0.576996638683937024261184>,0.010000000000000000208167}cylinder{<0.607061998206686048717984, -0.000000000000000148687307, -0.576996638683937024261184>,<0.187592474085079868872938, 0.577350269189625620036566, -0.576996638683937024261184>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{5.000000000000000000000000,<0.187592474085079868872938, 0.577350269189625620036566, -0.576996638683937024261184>,<-0.491123473188422809965203, 0.356822089773089878850243, -0.576996638683937024261184>,<-0.491123473188422865476355, -0.356822089773089767827940, -0.576996638683937024261184>,<0.187592474085079730095060, -0.577350269189625620036566, -0.576996638683937024261184>,<0.607061998206686048717984, -0.000000000000000148687307, -0.576996638683937024261184>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<1.003868410778580466669041, -0.015708913405614588487680, 0.015950168358295685466786>,0.029999999999999998889777}sphere{<0.829638736580318170155124, -0.602767825006238600060726, 0.382411448631120587382526>,0.029999999999999998889777}sphere{<0.325152463505077204963811, -0.949881272368329865329883, 0.015950168358295241377576>,0.029999999999999998889777}sphere{<0.187592474085079646828333, -0.577350269189625731058868, -0.576996638683937024261184>,0.029999999999999998889777}sphere{<0.607061998206686048717984, -0.000000000000000222044605, -0.576996638683937024261184>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<1.003868410778580466669041, -0.015708913405614588487680, 0.015950168358295685466786>,<0.829638736580318170155124, -0.602767825006238600060726, 0.382411448631120587382526>,0.010000000000000000208167}cylinder{<0.829638736580318170155124, -0.602767825006238600060726, 0.382411448631120587382526>,<0.325152463505077204963811, -0.949881272368329865329883, 0.015950168358295241377576>,0.010000000000000000208167}cylinder{<0.325152463505077204963811, -0.949881272368329865329883, 0.015950168358295241377576>,<0.187592474085079646828333, -0.577350269189625731058868, -0.576996638683937024261184>,0.010000000000000000208167}cylinder{<0.187592474085079646828333, -0.577350269189625731058868, -0.576996638683937024261184>,<0.607061998206686048717984, -0.000000000000000222044605, -0.576996638683937024261184>,0.010000000000000000208167}cylinder{<0.607061998206686048717984, -0.000000000000000222044605, -0.576996638683937024261184>,<1.003868410778580466669041, -0.015708913405614588487680, 0.015950168358295685466786>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{5.000000000000000000000000,<1.003868410778580466669041, -0.015708913405614588487680, 0.015950168358295685466786>,<0.829638736580318170155124, -0.602767825006238600060726, 0.382411448631120587382526>,<0.325152463505077204963811, -0.949881272368329865329883, 0.015950168358295241377576>,<0.187592474085079646828333, -0.577350269189625731058868, -0.576996638683937024261184>,<0.607061998206686048717984, -0.000000000000000222044605, -0.576996638683937024261184>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<1.003868410778580466669041, -0.015708913405614560732104, 0.015950168358295546688908>,0.029999999999999998889777}sphere{<0.829638736580318170155124, -0.602767825006238600060726, 0.382411448631120531871375>,0.029999999999999998889777}sphere{<0.575994966328297608448850, -0.399067554132528146126901, 1.017595715169939740718519>,0.029999999999999998889777}sphere{<0.593464169476141778503120, 0.313885048385609266574647, 1.043699900737277141971049>,0.029999999999999998889777}sphere{<0.857904501029906563225325, 0.550813718235801830225284, 0.424648908127707502835335>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<1.003868410778580466669041, -0.015708913405614560732104, 0.015950168358295546688908>,<0.829638736580318170155124, -0.602767825006238600060726, 0.382411448631120531871375>,0.010000000000000000208167}cylinder{<0.829638736580318170155124, -0.602767825006238600060726, 0.382411448631120531871375>,<0.575994966328297608448850, -0.399067554132528146126901, 1.017595715169939740718519>,0.010000000000000000208167}cylinder{<0.575994966328297608448850, -0.399067554132528146126901, 1.017595715169939740718519>,<0.593464169476141778503120, 0.313885048385609266574647, 1.043699900737277141971049>,0.010000000000000000208167}cylinder{<0.593464169476141778503120, 0.313885048385609266574647, 1.043699900737277141971049>,<0.857904501029906563225325, 0.550813718235801830225284, 0.424648908127707502835335>,0.010000000000000000208167}cylinder{<0.857904501029906563225325, 0.550813718235801830225284, 0.424648908127707502835335>,<1.003868410778580466669041, -0.015708913405614560732104, 0.015950168358295546688908>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{5.000000000000000000000000,<1.003868410778580466669041, -0.015708913405614560732104, 0.015950168358295546688908>,<0.829638736580318170155124, -0.602767825006238600060726, 0.382411448631120531871375>,<0.575994966328297608448850, -0.399067554132528146126901, 1.017595715169939740718519>,<0.593464169476141778503120, 0.313885048385609266574647, 1.043699900737277141971049>,<0.857904501029906563225325, 0.550813718235801830225284, 0.424648908127707502835335>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<0.557528031069579155243332, -0.424485109949140460017247, 1.017595715169939740718519>,0.029999999999999998889777}sphere{<0.829638736580318392199729, -0.602767825006238711083029, 0.382411448631120698404828>,0.029999999999999998889777}sphere{<0.325152463505077149452660, -0.949881272368329976352186, 0.015950168358295158110849>,0.029999999999999998889777}sphere{<-0.258747905623922191953312, -0.986126465733151769121889, 0.424648908127706614656915>,0.029999999999999998889777}sphere{<-0.115131906714012632875210, -0.661413779799332579578675, 1.043699900737276919926444>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<0.557528031069579155243332, -0.424485109949140460017247, 1.017595715169939740718519>,<0.829638736580318392199729, -0.602767825006238711083029, 0.382411448631120698404828>,0.010000000000000000208167}cylinder{<0.829638736580318392199729, -0.602767825006238711083029, 0.382411448631120698404828>,<0.325152463505077149452660, -0.949881272368329976352186, 0.015950168358295158110849>,0.010000000000000000208167}cylinder{<0.325152463505077149452660, -0.949881272368329976352186, 0.015950168358295158110849>,<-0.258747905623922191953312, -0.986126465733151769121889, 0.424648908127706614656915>,0.010000000000000000208167}cylinder{<-0.258747905623922191953312, -0.986126465733151769121889, 0.424648908127706614656915>,<-0.115131906714012632875210, -0.661413779799332579578675, 1.043699900737276919926444>,0.010000000000000000208167}cylinder{<-0.115131906714012632875210, -0.661413779799332579578675, 1.043699900737276919926444>,<0.557528031069579155243332, -0.424485109949140460017247, 1.017595715169939740718519>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{5.000000000000000000000000,<0.557528031069579155243332, -0.424485109949140460017247, 1.017595715169939740718519>,<0.829638736580318392199729, -0.602767825006238711083029, 0.382411448631120698404828>,<0.325152463505077149452660, -0.949881272368329976352186, 0.015950168358295158110849>,<-0.258747905623922191953312, -0.986126465733151769121889, 0.424648908127706614656915>,<-0.115131906714012632875210, -0.661413779799332579578675, 1.043699900737276919926444>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<-0.785443933658839132938567, -0.634877228845508345678184, 0.042054353925632892430286>,0.029999999999999998889777}sphere{<-0.276217108771766028940675, -0.985434888705109313100650, 0.398544722560369546471293>,0.029999999999999998889777}sphere{<0.325152463505077260474962, -0.949881272368329976352186, 0.015950168358295283010939>,0.029999999999999998889777}sphere{<0.187592474085079674583909, -0.577350269189625731058868, -0.576996638683937135283486>,0.029999999999999998889777}sphere{<-0.498793847145397928333210, -0.382667063698870768551075, -0.560863364754688231705870>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<-0.785443933658839132938567, -0.634877228845508345678184, 0.042054353925632892430286>,<-0.276217108771766028940675, -0.985434888705109313100650, 0.398544722560369546471293>,0.010000000000000000208167}cylinder{<-0.276217108771766028940675, -0.985434888705109313100650, 0.398544722560369546471293>,<0.325152463505077260474962, -0.949881272368329976352186, 0.015950168358295283010939>,0.010000000000000000208167}cylinder{<0.325152463505077260474962, -0.949881272368329976352186, 0.015950168358295283010939>,<0.187592474085079674583909, -0.577350269189625731058868, -0.576996638683937135283486>,0.010000000000000000208167}cylinder{<0.187592474085079674583909, -0.577350269189625731058868, -0.576996638683937135283486>,<-0.498793847145397928333210, -0.382667063698870768551075, -0.560863364754688231705870>,0.010000000000000000208167}cylinder{<-0.498793847145397928333210, -0.382667063698870768551075, -0.560863364754688231705870>,<-0.785443933658839132938567, -0.634877228845508345678184, 0.042054353925632892430286>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{5.000000000000000000000000,<-0.785443933658839132938567, -0.634877228845508345678184, 0.042054353925632892430286>,<-0.276217108771766028940675, -0.985434888705109313100650, 0.398544722560369546471293>,<0.325152463505077260474962, -0.949881272368329976352186, 0.015950168358295283010939>,<0.187592474085079674583909, -0.577350269189625731058868, -0.576996638683937135283486>,<-0.498793847145397928333210, -0.382667063698870768551075, -0.560863364754688231705870>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<1.003868410778580466669041, -0.015708913405614588487680, 0.015950168358295723630702>,0.029999999999999998889777}sphere{<0.851848491539995089638637, 0.567214208669458153089238, 0.398544722560369713004746>,0.029999999999999998889777}sphere{<0.361088601911640938446624, 0.943189424345671656446655, 0.042054353925632712019045>,0.029999999999999998889777}sphere{<0.209802229044756455289544, 0.592631764486071133113398, -0.560863364754688120683568>,0.029999999999999998889777}sphere{<0.607061998206686048717984, -0.000000000000000219442520, -0.576996638683937135283486>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<1.003868410778580466669041, -0.015708913405614588487680, 0.015950168358295723630702>,<0.851848491539995089638637, 0.567214208669458153089238, 0.398544722560369713004746>,0.010000000000000000208167}cylinder{<0.851848491539995089638637, 0.567214208669458153089238, 0.398544722560369713004746>,<0.361088601911640938446624, 0.943189424345671656446655, 0.042054353925632712019045>,0.010000000000000000208167}cylinder{<0.361088601911640938446624, 0.943189424345671656446655, 0.042054353925632712019045>,<0.209802229044756455289544, 0.592631764486071133113398, -0.560863364754688120683568>,0.010000000000000000208167}cylinder{<0.209802229044756455289544, 0.592631764486071133113398, -0.560863364754688120683568>,<0.607061998206686048717984, -0.000000000000000219442520, -0.576996638683937135283486>,0.010000000000000000208167}cylinder{<0.607061998206686048717984, -0.000000000000000219442520, -0.576996638683937135283486>,<1.003868410778580466669041, -0.015708913405614588487680, 0.015950168358295723630702>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{5.000000000000000000000000,<1.003868410778580466669041, -0.015708913405614588487680, 0.015950168358295723630702>,<0.851848491539995089638637, 0.567214208669458153089238, 0.398544722560369713004746>,<0.361088601911640938446624, 0.943189424345671656446655, 0.042054353925632712019045>,<0.209802229044756455289544, 0.592631764486071133113398, -0.560863364754688120683568>,<0.607061998206686048717984, -0.000000000000000219442520, -0.576996638683937135283486>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<-0.821380072065401645176053, 0.577350269189625620036566, 0.015950168358295435666605>,0.029999999999999998889777}sphere{<-0.491123473188422865476355, 0.356822089773089878850243, -0.576996638683937024261184>,0.029999999999999998889777}sphere{<-0.491123473188422865476355, -0.356822089773089767827940, -0.576996638683937024261184>,0.029999999999999998889777}sphere{<-0.821380072065401756198355, -0.577350269189625620036566, 0.015950168358295435666605>,0.029999999999999998889777}sphere{<-1.025489875180315202385373, 0.000000000000000003266362, 0.382411448631120698404828>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<-0.821380072065401645176053, 0.577350269189625620036566, 0.015950168358295435666605>,<-0.491123473188422865476355, 0.356822089773089878850243, -0.576996638683937024261184>,0.010000000000000000208167}cylinder{<-0.491123473188422865476355, 0.356822089773089878850243, -0.576996638683937024261184>,<-0.491123473188422865476355, -0.356822089773089767827940, -0.576996638683937024261184>,0.010000000000000000208167}cylinder{<-0.491123473188422865476355, -0.356822089773089767827940, -0.576996638683937024261184>,<-0.821380072065401756198355, -0.577350269189625620036566, 0.015950168358295435666605>,0.010000000000000000208167}cylinder{<-0.821380072065401756198355, -0.577350269189625620036566, 0.015950168358295435666605>,<-1.025489875180315202385373, 0.000000000000000003266362, 0.382411448631120698404828>,0.010000000000000000208167}cylinder{<-1.025489875180315202385373, 0.000000000000000003266362, 0.382411448631120698404828>,<-0.821380072065401645176053, 0.577350269189625620036566, 0.015950168358295435666605>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{5.000000000000000000000000,<-0.821380072065401645176053, 0.577350269189625620036566, 0.015950168358295435666605>,<-0.491123473188422865476355, 0.356822089773089878850243, -0.576996638683937024261184>,<-0.491123473188422865476355, -0.356822089773089767827940, -0.576996638683937024261184>,<-0.821380072065401756198355, -0.577350269189625620036566, 0.015950168358295435666605>,<-1.025489875180315202385373, 0.000000000000000003266362, 0.382411448631120698404828>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<-0.821380072065401645176053, 0.577350269189625620036566, 0.015950168358295546688908>,0.029999999999999998889777}sphere{<-0.370299140544567695343403, 0.949881272368329643285279, 0.424648908127706947723823>,0.029999999999999998889777}sphere{<-0.295625596302641791002230, 0.602767825006238266993819, 1.043699900737276697881839>,0.029999999999999998889777}sphere{<-0.700555739421546586065404, 0.015708913405614338687499, 1.017595715169939518673914>,0.029999999999999998889777}sphere{<-1.025489875180315202385373, 0.000000000000000055511151, 0.382411448631120753915980>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<-0.821380072065401645176053, 0.577350269189625620036566, 0.015950168358295546688908>,<-0.370299140544567695343403, 0.949881272368329643285279, 0.424648908127706947723823>,0.010000000000000000208167}cylinder{<-0.370299140544567695343403, 0.949881272368329643285279, 0.424648908127706947723823>,<-0.295625596302641791002230, 0.602767825006238266993819, 1.043699900737276697881839>,0.010000000000000000208167}cylinder{<-0.295625596302641791002230, 0.602767825006238266993819, 1.043699900737276697881839>,<-0.700555739421546586065404, 0.015708913405614338687499, 1.017595715169939518673914>,0.010000000000000000208167}cylinder{<-0.700555739421546586065404, 0.015708913405614338687499, 1.017595715169939518673914>,<-1.025489875180315202385373, 0.000000000000000055511151, 0.382411448631120753915980>,0.010000000000000000208167}cylinder{<-1.025489875180315202385373, 0.000000000000000055511151, 0.382411448631120753915980>,<-0.821380072065401645176053, 0.577350269189625620036566, 0.015950168358295546688908>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{5.000000000000000000000000,<-0.821380072065401645176053, 0.577350269189625620036566, 0.015950168358295546688908>,<-0.370299140544567695343403, 0.949881272368329643285279, 0.424648908127706947723823>,<-0.295625596302641791002230, 0.602767825006238266993819, 1.043699900737276697881839>,<-0.700555739421546586065404, 0.015708913405614338687499, 1.017595715169939518673914>,<-1.025489875180315202385373, 0.000000000000000055511151, 0.382411448631120753915980>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<-0.821380072065401423131448, 0.577350269189625509014263, 0.015950168358295546688908>,0.029999999999999998889777}sphere{<-0.370299140544567695343403, 0.949881272368329421240674, 0.424648908127706836701520>,0.029999999999999998889777}sphere{<0.256567651996951806570024, 0.985434888705109202078347, 0.085441037089571403306820>,0.029999999999999998889777}sphere{<0.192911704685406082404597, 0.634877228845508900789696, -0.532899696232898678616152>,0.029999999999999998889777}sphere{<-0.473296626880721271746211, 0.382667063698871379173738, -0.575847415016584651681342>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<-0.821380072065401423131448, 0.577350269189625509014263, 0.015950168358295546688908>,<-0.370299140544567695343403, 0.949881272368329421240674, 0.424648908127706836701520>,0.010000000000000000208167}cylinder{<-0.370299140544567695343403, 0.949881272368329421240674, 0.424648908127706836701520>,<0.256567651996951806570024, 0.985434888705109202078347, 0.085441037089571403306820>,0.010000000000000000208167}cylinder{<0.256567651996951806570024, 0.985434888705109202078347, 0.085441037089571403306820>,<0.192911704685406082404597, 0.634877228845508900789696, -0.532899696232898678616152>,0.010000000000000000208167}cylinder{<0.192911704685406082404597, 0.634877228845508900789696, -0.532899696232898678616152>,<-0.473296626880721271746211, 0.382667063698871379173738, -0.575847415016584651681342>,0.010000000000000000208167}cylinder{<-0.473296626880721271746211, 0.382667063698871379173738, -0.575847415016584651681342>,<-0.821380072065401423131448, 0.577350269189625509014263, 0.015950168358295546688908>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{5.000000000000000000000000,<-0.821380072065401423131448, 0.577350269189625509014263, 0.015950168358295546688908>,<-0.370299140544567695343403, 0.949881272368329421240674, 0.424648908127706836701520>,<0.256567651996951806570024, 0.985434888705109202078347, 0.085441037089571403306820>,<0.192911704685406082404597, 0.634877228845508900789696, -0.532899696232898678616152>,<-0.473296626880721271746211, 0.382667063698871379173738, -0.575847415016584651681342>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<0.270872714232332389983071, 0.986126465733151436054982, 0.113404705611361067418841>,0.029999999999999998889777}sphere{<-0.370299140544567584321101, 0.949881272368329532262976, 0.424648908127707003234974>,0.029999999999999998889777}sphere{<-0.295625596302641846513382, 0.602767825006238155971516, 1.043699900737276475837234>,0.029999999999999998889777}sphere{<0.391697046876187338071418, 0.424485109949140404506096, 1.115050252423004817359242>,0.029999999999999998889777}sphere{<0.741812257356194382218462, 0.661413779799332468556372, 0.540096202264473412846257>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<0.270872714232332389983071, 0.986126465733151436054982, 0.113404705611361067418841>,<-0.370299140544567584321101, 0.949881272368329532262976, 0.424648908127707003234974>,0.010000000000000000208167}cylinder{<-0.370299140544567584321101, 0.949881272368329532262976, 0.424648908127707003234974>,<-0.295625596302641846513382, 0.602767825006238155971516, 1.043699900737276475837234>,0.010000000000000000208167}cylinder{<-0.295625596302641846513382, 0.602767825006238155971516, 1.043699900737276475837234>,<0.391697046876187338071418, 0.424485109949140404506096, 1.115050252423004817359242>,0.010000000000000000208167}cylinder{<0.391697046876187338071418, 0.424485109949140404506096, 1.115050252423004817359242>,<0.741812257356194382218462, 0.661413779799332468556372, 0.540096202264473412846257>,0.010000000000000000208167}cylinder{<0.741812257356194382218462, 0.661413779799332468556372, 0.540096202264473412846257>,<0.270872714232332389983071, 0.986126465733151436054982, 0.113404705611361067418841>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{5.000000000000000000000000,<0.270872714232332389983071, 0.986126465733151436054982, 0.113404705611361067418841>,<-0.370299140544567584321101, 0.949881272368329532262976, 0.424648908127707003234974>,<-0.295625596302641846513382, 0.602767825006238155971516, 1.043699900737276475837234>,<0.391697046876187338071418, 0.424485109949140404506096, 1.115050252423004817359242>,<0.741812257356194382218462, 0.661413779799332468556372, 0.540096202264473412846257>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<0.397016277476513912425560, -0.313885048385609877197311, 1.159147194874043274026576>,0.029999999999999998889777}sphere{<0.382711215241133440034815, 0.399067554132527591015389, 1.131183526352253831959160>,0.029999999999999998889777}sphere{<-0.295625596302641957535684, 0.602767825006238378016121, 1.043699900737276697881839>,0.029999999999999998889777}sphere{<-0.700555739421546919132311, 0.015708913405614338687499, 1.017595715169939518673914>,0.029999999999999998889777}sphere{<-0.272479519394613900473701, -0.550813718235802274314494, 1.088946066855667860195922>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<0.397016277476513912425560, -0.313885048385609877197311, 1.159147194874043274026576>,<0.382711215241133440034815, 0.399067554132527591015389, 1.131183526352253831959160>,0.010000000000000000208167}cylinder{<0.382711215241133440034815, 0.399067554132527591015389, 1.131183526352253831959160>,<-0.295625596302641957535684, 0.602767825006238378016121, 1.043699900737276697881839>,0.010000000000000000208167}cylinder{<-0.295625596302641957535684, 0.602767825006238378016121, 1.043699900737276697881839>,<-0.700555739421546919132311, 0.015708913405614338687499, 1.017595715169939518673914>,0.010000000000000000208167}cylinder{<-0.700555739421546919132311, 0.015708913405614338687499, 1.017595715169939518673914>,<-0.272479519394613900473701, -0.550813718235802274314494, 1.088946066855667860195922>,0.010000000000000000208167}cylinder{<-0.272479519394613900473701, -0.550813718235802274314494, 1.088946066855667860195922>,<0.397016277476513912425560, -0.313885048385609877197311, 1.159147194874043274026576>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{5.000000000000000000000000,<0.397016277476513912425560, -0.313885048385609877197311, 1.159147194874043274026576>,<0.382711215241133440034815, 0.399067554132527591015389, 1.131183526352253831959160>,<-0.295625596302641957535684, 0.602767825006238378016121, 1.043699900737276697881839>,<-0.700555739421546919132311, 0.015708913405614338687499, 1.017595715169939518673914>,<-0.272479519394613900473701, -0.550813718235802274314494, 1.088946066855667860195922>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<-0.818092606760401408472205, -0.592631764486071133113398, 0.043203577592985445421370>,0.029999999999999998889777}sphere{<-0.364979909944241287522715, -0.943189424345671656446655, 0.468745850578745626435762>,0.029999999999999998889777}sphere{<-0.292338130997641554298383, -0.567214208669458153089238, 1.070953309971966804781118>,0.029999999999999998889777}sphere{<-0.700555739421546586065404, 0.015708913405614342156946, 1.017595715169939740718519>,0.029999999999999998889777}sphere{<-1.025489875180315202385373, 0.000000000000000055511151, 0.382411448631120753915980>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<-0.818092606760401408472205, -0.592631764486071133113398, 0.043203577592985445421370>,<-0.364979909944241287522715, -0.943189424345671656446655, 0.468745850578745626435762>,0.010000000000000000208167}cylinder{<-0.364979909944241287522715, -0.943189424345671656446655, 0.468745850578745626435762>,<-0.292338130997641554298383, -0.567214208669458153089238, 1.070953309971966804781118>,0.010000000000000000208167}cylinder{<-0.292338130997641554298383, -0.567214208669458153089238, 1.070953309971966804781118>,<-0.700555739421546586065404, 0.015708913405614342156946, 1.017595715169939740718519>,0.010000000000000000208167}cylinder{<-0.700555739421546586065404, 0.015708913405614342156946, 1.017595715169939740718519>,<-1.025489875180315202385373, 0.000000000000000055511151, 0.382411448631120753915980>,0.010000000000000000208167}cylinder{<-1.025489875180315202385373, 0.000000000000000055511151, 0.382411448631120753915980>,<-0.818092606760401408472205, -0.592631764486071133113398, 0.043203577592985445421370>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{5.000000000000000000000000,<-0.818092606760401408472205, -0.592631764486071133113398, 0.043203577592985445421370>,<-0.364979909944241287522715, -0.943189424345671656446655, 0.468745850578745626435762>,<-0.292338130997641554298383, -0.567214208669458153089238, 1.070953309971966804781118>,<-0.700555739421546586065404, 0.015708913405614342156946, 1.017595715169939740718519>,<-1.025489875180315202385373, 0.000000000000000055511151, 0.382411448631120753915980>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}}
