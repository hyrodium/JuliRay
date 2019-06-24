//Hyrodium's POV-Ray format
//Hy_color.inc

//#include "Hy_color.inc"
#declare Red=rgb <1,0,0>;
#declare Green=rgb <0,1,0>;
#declare Blue=rgb <0,0,1>;

#declare Cyan=rgb <0,1,1>;
#declare Magenta=rgb <1,0,1>;
#declare Yellow=rgb <1,1,0>;

#declare White=rgb 1;
#declare Gray=rgb 0.5;
#declare Black=rgb 0;

#declare RedTr=rgbft <1,0,0,0.3,0.3>;
#declare GreenTr=rgbft <0,1,0,0.3,0.3>;
#declare BlueTr=rgbft <0,0,1,0.3,0.3>;

#declare CyanTr=rgbft <0,1,1,0.3,0.3>;
#declare MagentaTr=rgbft <1,0,1,0.3,0.3>;
#declare YellowTr=rgbft <1,1,0,0.3,0.3>;

#declare WhiteTr=rgbft <1,1,1,0.3,0.3>;
#declare GrayTr=rgbft <0.5,0.5,0.5,0.3,0.3>;
#declare BlackTr=rgbft <0,0,0,0.3,0.3>;


#declare Clear=rgbf 1;


#macro HSVtoRGB(r,g,b)
    <r,g,b>
#end

#macro Color(RGB,F,T)
    rgbft<vdot(RGB,x),vdot(RGB,y),vdot(RGB,z),F,T>
#end
