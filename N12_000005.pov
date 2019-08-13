#version 3.7;
global_settings{assumed_gamma 1.0}
camera{perspective location <18.750000000000000000000000, 10.825317547305481014063844, 12.499999999999996447286321> right <0.133333333333333331482962, -0.230940107675850353485814, -0.000000000000000000000000> up <-0.064951905283832905846353, -0.037499999999999998612221, 0.129903810567665811692706> direction <0.750000000000000111022302, 0.433012701892219298294151, 0.499999999999999944488849> sky <-0.433012701892219298294151, -0.249999999999999944488849, 0.866025403784438707610605> look_at <18.000000000000000000000000, 10.392304845413262270881205, 11.999999999999996447286321>}
light_source{<18.750000000000000000000000, 10.825317547305481014063844, 12.499999999999996447286321> rgb<1.0,1.0,1.0>}
background{rgb<1,1,1>}
union{union{object{union{sphere{<0.187592474085079868872938, 0.577350269189625620036566, -0.788260990513641734978023>,0.029999999999999998889777}sphere{<-0.491123473188422809965203, 0.356822089773089878850243, -0.788260990513641734978023>,0.029999999999999998889777}sphere{<-0.491123473188422865476355, -0.356822089773089767827940, -0.788260990513641734978023>,0.029999999999999998889777}sphere{<0.187592474085079730095060, -0.577350269189625620036566, -0.788260990513641734978023>,0.029999999999999998889777}sphere{<0.607061998206686048717984, -0.000000000000000148687307, -0.788260990513641734978023>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<0.187592474085079868872938, 0.577350269189625620036566, -0.788260990513641734978023>,<-0.491123473188422809965203, 0.356822089773089878850243, -0.788260990513641734978023>,0.010000000000000000208167}cylinder{<-0.491123473188422809965203, 0.356822089773089878850243, -0.788260990513641734978023>,<-0.491123473188422865476355, -0.356822089773089767827940, -0.788260990513641734978023>,0.010000000000000000208167}cylinder{<-0.491123473188422865476355, -0.356822089773089767827940, -0.788260990513641734978023>,<0.187592474085079730095060, -0.577350269189625620036566, -0.788260990513641734978023>,0.010000000000000000208167}cylinder{<0.187592474085079730095060, -0.577350269189625620036566, -0.788260990513641734978023>,<0.607061998206686048717984, -0.000000000000000148687307, -0.788260990513641734978023>,0.010000000000000000208167}cylinder{<0.607061998206686048717984, -0.000000000000000148687307, -0.788260990513641734978023>,<0.187592474085079868872938, 0.577350269189625620036566, -0.788260990513641734978023>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{5.000000000000000000000000,<0.187592474085079868872938, 0.577350269189625620036566, -0.788260990513641734978023>,<-0.491123473188422809965203, 0.356822089773089878850243, -0.788260990513641734978023>,<-0.491123473188422865476355, -0.356822089773089767827940, -0.788260990513641734978023>,<0.187592474085079730095060, -0.577350269189625620036566, -0.788260990513641734978023>,<0.607061998206686048717984, -0.000000000000000148687307, -0.788260990513641734978023>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<0.982261278771691115352382, -0.000010413094383293630685, -0.181207850537238868149359>,0.029999999999999998889777}sphere{<0.794677662593765599119422, -0.577367117930265849601312, 0.193971622945525612458084>,0.029999999999999998889777}sphere{<0.303545331498188020180606, -0.934182772057098764761918, -0.181207850537239034682813>,0.029999999999999998889777}sphere{<0.187592474085079674583909, -0.577350269189625620036566, -0.788260990513641846000326>,0.029999999999999998889777}sphere{<0.607061998206686048717984, 0.000000000000000027755576, -0.788260990513641846000326>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<0.982261278771691115352382, -0.000010413094383293630685, -0.181207850537238868149359>,<0.794677662593765599119422, -0.577367117930265849601312, 0.193971622945525612458084>,0.010000000000000000208167}cylinder{<0.794677662593765599119422, -0.577367117930265849601312, 0.193971622945525612458084>,<0.303545331498188020180606, -0.934182772057098764761918, -0.181207850537239034682813>,0.010000000000000000208167}cylinder{<0.303545331498188020180606, -0.934182772057098764761918, -0.181207850537239034682813>,<0.187592474085079674583909, -0.577350269189625620036566, -0.788260990513641846000326>,0.010000000000000000208167}cylinder{<0.187592474085079674583909, -0.577350269189625620036566, -0.788260990513641846000326>,<0.607061998206686048717984, 0.000000000000000027755576, -0.788260990513641846000326>,0.010000000000000000208167}cylinder{<0.607061998206686048717984, 0.000000000000000027755576, -0.788260990513641846000326>,<0.982261278771691115352382, -0.000010413094383293630685, -0.181207850537238868149359>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{5.000000000000000000000000,<0.982261278771691115352382, -0.000010413094383293630685, -0.181207850537238868149359>,<0.794677662593765599119422, -0.577367117930265849601312, 0.193971622945525612458084>,<0.303545331498188020180606, -0.934182772057098764761918, -0.181207850537239034682813>,<0.187592474085079674583909, -0.577350269189625620036566, -0.788260990513641846000326>,<0.607061998206686048717984, 0.000000000000000027755576, -0.788260990513641846000326>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<0.982261278771691115352382, -0.000010413094383321386260, -0.181207850537238895904935>,0.029999999999999998889777}sphere{<0.794677662593765599119422, -0.577367117930265849601312, 0.193971622945525473680206>,0.029999999999999998889777}sphere{<0.491178712376336423783840, -0.356849352099808059257668, 0.801053427293343056625474>,0.029999999999999998889777}sphere{<0.491189661769978358485389, 0.356794827142487558369055, 0.801071142849143602582274>,0.029999999999999998889777}sphere{<0.794695379084834452321218, 0.577333419957290816526552, 0.194000287316940467352566>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<0.982261278771691115352382, -0.000010413094383321386260, -0.181207850537238895904935>,<0.794677662593765599119422, -0.577367117930265849601312, 0.193971622945525473680206>,0.010000000000000000208167}cylinder{<0.794677662593765599119422, -0.577367117930265849601312, 0.193971622945525473680206>,<0.491178712376336423783840, -0.356849352099808059257668, 0.801053427293343056625474>,0.010000000000000000208167}cylinder{<0.491178712376336423783840, -0.356849352099808059257668, 0.801053427293343056625474>,<0.491189661769978358485389, 0.356794827142487558369055, 0.801071142849143602582274>,0.010000000000000000208167}cylinder{<0.491189661769978358485389, 0.356794827142487558369055, 0.801071142849143602582274>,<0.794695379084834452321218, 0.577333419957290816526552, 0.194000287316940467352566>,0.010000000000000000208167}cylinder{<0.794695379084834452321218, 0.577333419957290816526552, 0.194000287316940467352566>,<0.982261278771691115352382, -0.000010413094383321386260, -0.181207850537238895904935>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{5.000000000000000000000000,<0.982261278771691115352382, -0.000010413094383321386260, -0.181207850537238895904935>,<0.794677662593765599119422, -0.577367117930265849601312, 0.193971622945525473680206>,<0.491178712376336423783840, -0.356849352099808059257668, 0.801053427293343056625474>,<0.491189661769978358485389, 0.356794827142487558369055, 0.801071142849143602582274>,<0.794695379084834452321218, 0.577333419957290816526552, 0.194000287316940467352566>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<0.491166471049717778818433, -0.356866200840448066777810, 0.801053427293343167647777>,0.029999999999999998889777}sphere{<0.794677662593765377074817, -0.577367117930265960623615, 0.193971622945525695724811>,0.029999999999999998889777}sphere{<0.303545331498188075691758, -0.934182772057098764761918, -0.181207850537239006927237>,0.029999999999999998889777}sphere{<-0.303502333636893217860830, -0.934206056935690365428115, 0.194000287316939884485478>,0.029999999999999998889777}sphere{<-0.187546092386060736512121, -0.577404793655250880846097, 0.801071142849143269515366>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<0.491166471049717778818433, -0.356866200840448066777810, 0.801053427293343167647777>,<0.794677662593765377074817, -0.577367117930265960623615, 0.193971622945525695724811>,0.010000000000000000208167}cylinder{<0.794677662593765377074817, -0.577367117930265960623615, 0.193971622945525695724811>,<0.303545331498188075691758, -0.934182772057098764761918, -0.181207850537239006927237>,0.010000000000000000208167}cylinder{<0.303545331498188075691758, -0.934182772057098764761918, -0.181207850537239006927237>,<-0.303502333636893217860830, -0.934206056935690365428115, 0.194000287316939884485478>,0.010000000000000000208167}cylinder{<-0.303502333636893217860830, -0.934206056935690365428115, 0.194000287316939884485478>,<-0.187546092386060736512121, -0.577404793655250880846097, 0.801071142849143269515366>,0.010000000000000000208167}cylinder{<-0.187546092386060736512121, -0.577404793655250880846097, 0.801071142849143269515366>,<0.491166471049717778818433, -0.356866200840448066777810, 0.801053427293343167647777>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{5.000000000000000000000000,<0.491166471049717778818433, -0.356866200840448066777810, 0.801053427293343167647777>,<0.794677662593765377074817, -0.577367117930265960623615, 0.193971622945525695724811>,<0.303545331498188075691758, -0.934182772057098764761918, -0.181207850537239006927237>,<-0.303502333636893217860830, -0.934206056935690365428115, 0.194000287316939884485478>,<-0.187546092386060736512121, -0.577404793655250880846097, 0.801071142849143269515366>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<-0.794648997385814959137917, -0.577387944422915744269176, -0.181190134981438405459286>,0.029999999999999998889777}sphere{<-0.303513283030535041540077, -0.934206056631806114332051, 0.193982571761140004662494>,0.029999999999999998889777}sphere{<0.303545331498187964669455, -0.934182772057098764761918, -0.181207850537239006927237>,0.029999999999999998889777}sphere{<0.187592474085079730095060, -0.577350269189625509014263, -0.788260990513641957022628>,0.029999999999999998889777}sphere{<-0.491128947417614480919212, -0.356838938701540153708436, -0.788250041698027592573794>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<-0.794648997385814959137917, -0.577387944422915744269176, -0.181190134981438405459286>,<-0.303513283030535041540077, -0.934206056631806114332051, 0.193982571761140004662494>,0.010000000000000000208167}cylinder{<-0.303513283030535041540077, -0.934206056631806114332051, 0.193982571761140004662494>,<0.303545331498187964669455, -0.934182772057098764761918, -0.181207850537239006927237>,0.010000000000000000208167}cylinder{<0.303545331498187964669455, -0.934182772057098764761918, -0.181207850537239006927237>,<0.187592474085079730095060, -0.577350269189625509014263, -0.788260990513641957022628>,0.010000000000000000208167}cylinder{<0.187592474085079730095060, -0.577350269189625509014263, -0.788260990513641957022628>,<-0.491128947417614480919212, -0.356838938701540153708436, -0.788250041698027592573794>,0.010000000000000000208167}cylinder{<-0.491128947417614480919212, -0.356838938701540153708436, -0.788250041698027592573794>,<-0.794648997385814959137917, -0.577387944422915744269176, -0.181190134981438405459286>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{5.000000000000000000000000,<-0.794648997385814959137917, -0.577387944422915744269176, -0.181190134981438405459286>,<-0.303513283030535041540077, -0.934206056631806114332051, 0.193982571761140004662494>,<0.303545331498187964669455, -0.934182772057098764761918, -0.181207850537239006927237>,<0.187592474085079730095060, -0.577350269189625509014263, -0.788260990513641957022628>,<-0.491128947417614480919212, -0.356838938701540153708436, -0.788250041698027592573794>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<0.982261278771691115352382, -0.000010413094383293627297, -0.181207850537238840393783>,0.029999999999999998889777}sphere{<0.794691995247110205546903, 0.577343833355557944919667, 0.193982571761140115684796>,0.029999999999999998889777}sphere{<0.303568522218448433314109, 0.934178794305088433524986, -0.181190134981438349948135>,0.029999999999999998889777}sphere{<0.187606806738424197744664, 0.577360682096198285506716, -0.788250041698027592573794>,0.029999999999999998889777}sphere{<0.607061998206686048717984, 0.000000000000000027759811, -0.788260990513641957022628>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<0.982261278771691115352382, -0.000010413094383293627297, -0.181207850537238840393783>,<0.794691995247110205546903, 0.577343833355557944919667, 0.193982571761140115684796>,0.010000000000000000208167}cylinder{<0.794691995247110205546903, 0.577343833355557944919667, 0.193982571761140115684796>,<0.303568522218448433314109, 0.934178794305088433524986, -0.181190134981438349948135>,0.010000000000000000208167}cylinder{<0.303568522218448433314109, 0.934178794305088433524986, -0.181190134981438349948135>,<0.187606806738424197744664, 0.577360682096198285506716, -0.788250041698027592573794>,0.010000000000000000208167}cylinder{<0.187606806738424197744664, 0.577360682096198285506716, -0.788250041698027592573794>,<0.607061998206686048717984, 0.000000000000000027759811, -0.788260990513641957022628>,0.010000000000000000208167}cylinder{<0.607061998206686048717984, 0.000000000000000027759811, -0.788260990513641957022628>,<0.982261278771691115352382, -0.000010413094383293627297, -0.181207850537238840393783>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{5.000000000000000000000000,<0.982261278771691115352382, -0.000010413094383293627297, -0.181207850537238840393783>,<0.794691995247110205546903, 0.577343833355557944919667, 0.193982571761140115684796>,<0.303568522218448433314109, 0.934178794305088433524986, -0.181190134981438349948135>,<0.187606806738424197744664, 0.577360682096198285506716, -0.788250041698027592573794>,<0.607061998206686048717984, 0.000000000000000027759811, -0.788260990513641957022628>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<-0.794672188106074761648756, 0.577350269189625731058868, -0.181207850537239228971842>,0.029999999999999998889777}sphere{<-0.491123473188422643431750, 0.356822089773089878850243, -0.788260990513641846000326>,0.029999999999999998889777}sphere{<-0.491123473188422754454052, -0.356822089773089767827940, -0.788260990513641734978023>,0.029999999999999998889777}sphere{<-0.794672188106074872671059, -0.577350269189625509014263, -0.181207850537239256727418>,0.029999999999999998889777}sphere{<-0.982275611166536166685148, 0.000000000000000008096753, 0.193971622945525112857723>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<-0.794672188106074761648756, 0.577350269189625731058868, -0.181207850537239228971842>,<-0.491123473188422643431750, 0.356822089773089878850243, -0.788260990513641846000326>,0.010000000000000000208167}cylinder{<-0.491123473188422643431750, 0.356822089773089878850243, -0.788260990513641846000326>,<-0.491123473188422754454052, -0.356822089773089767827940, -0.788260990513641734978023>,0.010000000000000000208167}cylinder{<-0.491123473188422754454052, -0.356822089773089767827940, -0.788260990513641734978023>,<-0.794672188106074872671059, -0.577350269189625509014263, -0.181207850537239256727418>,0.010000000000000000208167}cylinder{<-0.794672188106074872671059, -0.577350269189625509014263, -0.181207850537239256727418>,<-0.982275611166536166685148, 0.000000000000000008096753, 0.193971622945525112857723>,0.010000000000000000208167}cylinder{<-0.982275611166536166685148, 0.000000000000000008096753, 0.193971622945525112857723>,<-0.794672188106074761648756, 0.577350269189625731058868, -0.181207850537239228971842>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{5.000000000000000000000000,<-0.794672188106074761648756, 0.577350269189625731058868, -0.181207850537239228971842>,<-0.491123473188422643431750, 0.356822089773089878850243, -0.788260990513641846000326>,<-0.491123473188422754454052, -0.356822089773089767827940, -0.788260990513641734978023>,<-0.794672188106074872671059, -0.577350269189625509014263, -0.181207850537239256727418>,<-0.982275611166536166685148, 0.000000000000000008096753, 0.193971622945525112857723>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<-0.794672188106074761648756, 0.577350269189625620036566, -0.181207850537239228971842>,0.029999999999999998889777}sphere{<-0.303573997124399752589596, 0.934182772057098320672708, 0.194000287316940245307961>,0.029999999999999998889777}sphere{<-0.187662046344598903990786, 0.577367117930265072445195, 0.801071142849143269515366>,0.029999999999999998889777}sphere{<-0.607122712042052037340056, 0.000010413094382960563777, 0.801053427293342390491659>,0.029999999999999998889777}sphere{<-0.982275611166536166685148, 0.000000000000000027755576, 0.193971622945525112857723>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<-0.794672188106074761648756, 0.577350269189625620036566, -0.181207850537239228971842>,<-0.303573997124399752589596, 0.934182772057098320672708, 0.194000287316940245307961>,0.010000000000000000208167}cylinder{<-0.303573997124399752589596, 0.934182772057098320672708, 0.194000287316940245307961>,<-0.187662046344598903990786, 0.577367117930265072445195, 0.801071142849143269515366>,0.010000000000000000208167}cylinder{<-0.187662046344598903990786, 0.577367117930265072445195, 0.801071142849143269515366>,<-0.607122712042052037340056, 0.000010413094382960563777, 0.801053427293342390491659>,0.010000000000000000208167}cylinder{<-0.607122712042052037340056, 0.000010413094382960563777, 0.801053427293342390491659>,<-0.982275611166536166685148, 0.000000000000000027755576, 0.193971622945525112857723>,0.010000000000000000208167}cylinder{<-0.982275611166536166685148, 0.000000000000000027755576, 0.193971622945525112857723>,<-0.794672188106074761648756, 0.577350269189625620036566, -0.181207850537239228971842>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{5.000000000000000000000000,<-0.794672188106074761648756, 0.577350269189625620036566, -0.181207850537239228971842>,<-0.303573997124399752589596, 0.934182772057098320672708, 0.194000287316940245307961>,<-0.187662046344598903990786, 0.577367117930265072445195, 0.801071142849143269515366>,<-0.607122712042052037340056, 0.000010413094382960563777, 0.801053427293342390491659>,<-0.982275611166536166685148, 0.000000000000000027755576, 0.193971622945525112857723>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<-0.794672188106074872671059, 0.577350269189625509014263, -0.181207850537239173460691>,0.029999999999999998889777}sphere{<-0.303573997124399752589596, 0.934182772057098209650405, 0.194000287316940217552386>,0.029999999999999998889777}sphere{<0.303502332701619081944955, 0.934206056631805892287446, -0.181161470093032767048058>,0.029999999999999998889777}sphere{<0.187597947317965085822777, 0.577387944422916299380688, -0.788232325305718073416017>,0.029999999999999998889777}sphere{<-0.491111232120318264016134, 0.356838938701540986375704, -0.788260989996651173505882>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<-0.794672188106074872671059, 0.577350269189625509014263, -0.181207850537239173460691>,<-0.303573997124399752589596, 0.934182772057098209650405, 0.194000287316940217552386>,0.010000000000000000208167}cylinder{<-0.303573997124399752589596, 0.934182772057098209650405, 0.194000287316940217552386>,<0.303502332701619081944955, 0.934206056631805892287446, -0.181161470093032767048058>,0.010000000000000000208167}cylinder{<0.303502332701619081944955, 0.934206056631805892287446, -0.181161470093032767048058>,<0.187597947317965085822777, 0.577387944422916299380688, -0.788232325305718073416017>,0.010000000000000000208167}cylinder{<0.187597947317965085822777, 0.577387944422916299380688, -0.788232325305718073416017>,<-0.491111232120318264016134, 0.356838938701540986375704, -0.788260989996651173505882>,0.010000000000000000208167}cylinder{<-0.491111232120318264016134, 0.356838938701540986375704, -0.788260989996651173505882>,<-0.794672188106074872671059, 0.577350269189625509014263, -0.181207850537239173460691>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{5.000000000000000000000000,<-0.794672188106074872671059, 0.577350269189625509014263, -0.181207850537239173460691>,<-0.303573997124399752589596, 0.934182772057098209650405, 0.194000287316940217552386>,<0.303502332701619081944955, 0.934206056631805892287446, -0.181161470093032767048058>,<0.187597947317965085822777, 0.577387944422916299380688, -0.788232325305718073416017>,<-0.491111232120318264016134, 0.356838938701540986375704, -0.788260989996651173505882>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<0.303513280741717128563550, 0.934206056935690476450418, -0.181143753700723303401432>,0.029999999999999998889777}sphere{<-0.303573997124399808100748, 0.934182772057098542717313, 0.194000287316940245307961>,0.029999999999999998889777}sphere{<-0.187662046344598903990786, 0.577367117930265183467498, 0.801071142849143491559971>,0.029999999999999998889777}sphere{<0.491062756805739519805343, 0.356866200840447622688600, 0.801117524129858482595523>,0.029999999999999998889777}sphere{<0.794625803380429873712387, 0.577404793655251102890702, 0.194075333805579425261456>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<0.303513280741717128563550, 0.934206056935690476450418, -0.181143753700723303401432>,<-0.303573997124399808100748, 0.934182772057098542717313, 0.194000287316940245307961>,0.010000000000000000208167}cylinder{<-0.303573997124399808100748, 0.934182772057098542717313, 0.194000287316940245307961>,<-0.187662046344598903990786, 0.577367117930265183467498, 0.801071142849143491559971>,0.010000000000000000208167}cylinder{<-0.187662046344598903990786, 0.577367117930265183467498, 0.801071142849143491559971>,<0.491062756805739519805343, 0.356866200840447622688600, 0.801117524129858482595523>,0.010000000000000000208167}cylinder{<0.491062756805739519805343, 0.356866200840447622688600, 0.801117524129858482595523>,<0.794625803380429873712387, 0.577404793655251102890702, 0.194075333805579425261456>,0.010000000000000000208167}cylinder{<0.794625803380429873712387, 0.577404793655251102890702, 0.194075333805579425261456>,<0.303513280741717128563550, 0.934206056935690476450418, -0.181143753700723303401432>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{5.000000000000000000000000,<0.303513280741717128563550, 0.934206056935690476450418, -0.181143753700723303401432>,<-0.303573997124399808100748, 0.934182772057098542717313, 0.194000287316940245307961>,<-0.187662046344598903990786, 0.577367117930265183467498, 0.801071142849143491559971>,<0.491062756805739519805343, 0.356866200840447622688600, 0.801117524129858482595523>,<0.794625803380429873712387, 0.577404793655251102890702, 0.194075333805579425261456>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<0.491068230038624431443850, -0.356794827142488335525172, 0.801146189337781589046017>,0.029999999999999998889777}sphere{<0.491057281998526495847557, 0.356849352099807171079249, 0.801128472945472736022055>,0.029999999999999998889777}sphere{<-0.187662046344598987257513, 0.577367117930265183467498, 0.801071142849143269515366>,0.029999999999999998889777}sphere{<-0.607122712042052148362359, 0.000010413094383043830504, 0.801053427293342390491659>,0.029999999999999998889777}sphere{<-0.187644332043610362337205, -0.577333419957290927548854, 0.801099808574057381527211>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<0.491068230038624431443850, -0.356794827142488335525172, 0.801146189337781589046017>,<0.491057281998526495847557, 0.356849352099807171079249, 0.801128472945472736022055>,0.010000000000000000208167}cylinder{<0.491057281998526495847557, 0.356849352099807171079249, 0.801128472945472736022055>,<-0.187662046344598987257513, 0.577367117930265183467498, 0.801071142849143269515366>,0.010000000000000000208167}cylinder{<-0.187662046344598987257513, 0.577367117930265183467498, 0.801071142849143269515366>,<-0.607122712042052148362359, 0.000010413094383043830504, 0.801053427293342390491659>,0.010000000000000000208167}cylinder{<-0.607122712042052148362359, 0.000010413094383043830504, 0.801053427293342390491659>,<-0.187644332043610362337205, -0.577333419957290927548854, 0.801099808574057381527211>,0.010000000000000000208167}cylinder{<-0.187644332043610362337205, -0.577333419957290927548854, 0.801099808574057381527211>,<0.491068230038624431443850, -0.356794827142488335525172, 0.801146189337781589046017>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{5.000000000000000000000000,<0.491068230038624431443850, -0.356794827142488335525172, 0.801146189337781589046017>,<0.491057281998526495847557, 0.356849352099807171079249, 0.801128472945472736022055>,<-0.187662046344598987257513, 0.577367117930265183467498, 0.801071142849143269515366>,<-0.607122712042052148362359, 0.000010413094383043830504, 0.801053427293342390491659>,<-0.187644332043610362337205, -0.577333419957290927548854, 0.801099808574057381527211>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}union{object{union{sphere{<-0.794668805462123706284672, -0.577360682096197730395204, -0.181190134464448399098657>,0.029999999999999998889777}sphere{<-0.303568523891515174017997, -0.934178794305087989435776, 0.194028952524862963180396>,0.029999999999999998889777}sphere{<-0.187658663700647904137853, -0.577343833355557944919667, 0.801088858921934154899702>,0.029999999999999998889777}sphere{<-0.607122712042052037340056, 0.000010413094382960568860, 0.801053427293342390491659>,0.029999999999999998889777}sphere{<-0.982275611166536166685148, 0.000000000000000027754305, 0.193971622945525140613299>,0.029999999999999998889777}} pigment{rgb<0.1,0.1,0.1>}}object{union{cylinder{<-0.794668805462123706284672, -0.577360682096197730395204, -0.181190134464448399098657>,<-0.303568523891515174017997, -0.934178794305087989435776, 0.194028952524862963180396>,0.010000000000000000208167}cylinder{<-0.303568523891515174017997, -0.934178794305087989435776, 0.194028952524862963180396>,<-0.187658663700647904137853, -0.577343833355557944919667, 0.801088858921934154899702>,0.010000000000000000208167}cylinder{<-0.187658663700647904137853, -0.577343833355557944919667, 0.801088858921934154899702>,<-0.607122712042052037340056, 0.000010413094382960568860, 0.801053427293342390491659>,0.010000000000000000208167}cylinder{<-0.607122712042052037340056, 0.000010413094382960568860, 0.801053427293342390491659>,<-0.982275611166536166685148, 0.000000000000000027754305, 0.193971622945525140613299>,0.010000000000000000208167}cylinder{<-0.982275611166536166685148, 0.000000000000000027754305, 0.193971622945525140613299>,<-0.794668805462123706284672, -0.577360682096197730395204, -0.181190134464448399098657>,0.010000000000000000208167}} pigment{rgb<0.1,0.1,0.1>}}object{polygon{5.000000000000000000000000,<-0.794668805462123706284672, -0.577360682096197730395204, -0.181190134464448399098657>,<-0.303568523891515174017997, -0.934178794305087989435776, 0.194028952524862963180396>,<-0.187658663700647904137853, -0.577343833355557944919667, 0.801088858921934154899702>,<-0.607122712042052037340056, 0.000010413094382960568860, 0.801053427293342390491659>,<-0.982275611166536166685148, 0.000000000000000027754305, 0.193971622945525140613299>} pigment{rgbft<0.5,1.0,1.0,0.1,0.1>}}}}
