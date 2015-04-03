function dydt = odes_div(t, y);

global p;

dydt = zeros(4901,1);

%% DIVkP
dydt(1) = -p.kdeg_dkp*y(1) + p.kph1_pc*y(501) - p.kpc_ph1*y(1)*y(401)   + p.kpk0_pk1*y(701) -  p.kpk1_pk0*y(801)*y(1) + p.kpk1_pk2p*y(801) - p.kpk2p_pk1*y(1301)*y(1) + p.kpt2_pk1h*y(1101) - p.kpk1h_pt2*y(1201)*y(1) + p.kjkp_dj*y(2601) - p.kdj_jkp*y(1)*y(201)+ p.kpk2_pk3*y(1001) - p.kpk3_pk2*y(1401)*y(1) + p.kpt3_pk1p*y(1501) - p.kpk1p_pt3*y(1601)*y(1) - p.kpk1p_1h*y(1601)*y(1) + p.k1h_pk1p*y(1201) - p.kph2_p12*y(2001)*y(1) + p.kp12_ph2*y(901)  + p.kpt4_pk3h*y(1901) - p.kpk3h_pt4*y(2101)*y(1) + p.kh1_h2*y(501)*y(101) - p.kh2_h1*y(2001)*y(1) -  p.kph1_p11*y(501)*y(1) + p.kp11_ph1*y(601) + (p.D_dk/y(4901)^2)*(y(2)-y(1))  - p.kbdl_dldk*y(2801)*y(1) + p.kudl_dldk*y(2901)+ p.kdeg_plc*(y(501)+2*y(601)+2*y(701)+y(801)+y(901)+y(1001)+2*y(1101)+y(1201)+y(1501)+y(1901));%+ p.kdeg_dl3*y(3101);
for v1a=2:49
    dydt(v1a) = -p.kdeg_dkp*y(v1a)  + p.kph1_pc*y(v1a+500) - p.kpc_ph1*y(v1a)*y(v1a+400)  + p.kpk0_pk1*y(v1a+700) -  p.kpk1_pk0*y(v1a+800)*y(v1a) + p.kpk1_pk2p*y(v1a+800) - p.kpk2p_pk1*y(v1a+1300)*y(v1a) + p.kpt2_pk1h*y(v1a+1100) - p.kpk1h_pt2*y(v1a+1200)*y(v1a) + p.kjkp_dj*y(v1a+2600) - p.kdj_jkp*y(v1a)*y(v1a+200) + p.kpk2_pk3*y(v1a+1000) - p.kpk3_pk2*y(v1a+1400)*y(v1a) + p.kpt3_pk1p*y(v1a+1500) - p.kpk1p_pt3*y(v1a+1600)*y(v1a) - p.kpk1p_1h*y(v1a+1600)*y(v1a) + p.k1h_pk1p*y(v1a+1200) - p.kph2_p12*y(v1a+2000)*y(v1a) + p.kp12_ph2*y(v1a+900)  + p.kpt4_pk3h*y(v1a+1900) - p.kpk3h_pt4*y(v1a+2100)*y(v1a) + p.kh1_h2*y(v1a+500)*y(v1a+100) - p.kh2_h1*y(v1a+2000)*y(v1a) -  p.kph1_p11*y(v1a+500)*y(v1a) + p.kp11_ph1*y(v1a+600)+ (p.D_dk/y(4901)^2)*(y(v1a-1)-2*y(v1a)+y(v1a+1))  - p.kbdl_dldk*y(v1a+2800)*y(v1a) + p.kudl_dldk*y(v1a+2900)+ p.kdeg_plc*(y(v1a+500)+2*y(v1a+600)+2*y(v1a+700)+y(v1a+800)+y(v1a+900)+y(v1a+1000)+2*y(v1a+1100)+y(v1a+1200)+y(v1a+1500)+y(v1a+1900));%+ p.kdeg_dl3*y(v1+3100);
end

dydt(50) = -p.kdeg_dkp*y(50)  + p.kph1_pc*y(550) - p.kpc_ph1*y(50)*y(450)  + p.kpk0_pk1*y(750) -  p.kpk1_pk0*y(850)*y(50) + p.kpk1_pk2p*y(850) - p.kpk2p_pk1*y(1350)*y(50) + p.kpt2_pk1h*y(1150) - p.kpk1h_pt2*y(1250)*y(50) + p.kjkp_dj*y(2650) - p.kdj_jkp*y(50)*y(250) + p.kpk2_pk3*y(1050) - p.kpk3_pk2*y(1450)*y(50) + p.kpt3_pk1p*y(1550) - p.kpk1p_pt3*y(1650)*y(50) - p.kpk1p_1h*y(1650)*y(50) + p.k1h_pk1p*y(1250) - p.kph2_p12*y(2050)*y(50) + p.kp12_ph2*y(950)  + p.kpt4_pk3h*y(1950) - p.kpk3h_pt4*y(2150)*y(50) + p.kh1_h2*y(550)*y(150) - p.kh2_h1*y(2050)*y(50) -  p.kph1_p11*y(550)*y(50) + p.kp11_ph1*y(650) +(p.D_dk/y(4901)^2)*(y(49)-y(50))  - p.kbdl_dldk*y(2850)*y(50) + p.kudl_dldk*y(2950)+ p.kdeg_plc*(y(550)+2*y(650)+2*y(750)+y(850)+y(950)+y(1050)+2*y(1150)+y(1250)+y(1550)+y(1950));
dydt(51) = -p.kdeg_dkp*y(51)  + p.kph1_pc*y(551) - p.kpc_ph1*y(51)*y(451)  + p.kpk0_pk1*y(751) -  p.kpk1_pk0*y(851)*y(51) + p.kpk1_pk2p*y(851) - p.kpk2p_pk1*y(1351)*y(51) + p.kpt2_pk1h*y(1151) - p.kpk1h_pt2*y(1251)*y(51) + p.kjkp_dj*y(2651) - p.kdj_jkp*y(51)*y(251) + p.kpk2_pk3*y(1051) - p.kpk3_pk2*y(1451)*y(51) + p.kpt3_pk1p*y(1551) - p.kpk1p_pt3*y(1651)*y(51) - p.kpk1p_1h*y(1651)*y(51) + p.k1h_pk1p*y(1251) - p.kph2_p12*y(2051)*y(51) + p.kp12_ph2*y(951)  + p.kpt4_pk3h*y(1951) - p.kpk3h_pt4*y(2151)*y(51) + p.kh1_h2*y(550)*y(151) - p.kh2_h1*y(2050)*y(51) -  p.kph1_p11*y(551)*y(51) + p.kp11_ph1*y(651) +(p.D_dk/y(4901)^2)*(y(52)-y(51))  - p.kbdl_dldk*y(2851)*y(51) + p.kudl_dldk*y(2951)+ p.kdeg_plc*(y(551)+2*y(651)+2*y(751)+y(851)+y(951)+y(1051)+2*y(1151)+y(1251)+y(1551)+y(1951));
for v1b=52:99
    dydt(v1b) = -p.kdeg_dkp*y(v1b)  + p.kph1_pc*y(v1b+500) - p.kpc_ph1*y(v1b)*y(v1b+400)  + p.kpk0_pk1*y(v1b+700) -  p.kpk1_pk0*y(v1b+800)*y(v1b) + p.kpk1_pk2p*y(v1b+800) - p.kpk2p_pk1*y(v1b+1300)*y(v1b) + p.kpt2_pk1h*y(v1b+1100) - p.kpk1h_pt2*y(v1b+1200)*y(v1b) + p.kjkp_dj*y(v1b+2600) - p.kdj_jkp*y(v1b)*y(v1b+200) + p.kpk2_pk3*y(v1b+1000) - p.kpk3_pk2*y(v1b+1400)*y(v1b) + p.kpt3_pk1p*y(v1b+1500) - p.kpk1p_pt3*y(v1b+1600)*y(v1b) - p.kpk1p_1h*y(v1b+1600)*y(v1b) + p.k1h_pk1p*y(v1b+1200) - p.kph2_p12*y(v1b+2000)*y(v1b) + p.kp12_ph2*y(v1b+900)  + p.kpt4_pk3h*y(v1b+1900) - p.kpk3h_pt4*y(v1b+2100)*y(v1b) + p.kh1_h2*y(v1b+500)*y(v1b+100) - p.kh2_h1*y(v1b+2000)*y(v1b) -  p.kph1_p11*y(v1b+500)*y(v1b) + p.kp11_ph1*y(v1b+600)+ (p.D_dk/y(4901)^2)*(y(v1b-1)-2*y(v1b)+y(v1b+1))  - p.kbdl_dldk*y(v1b+2800)*y(v1b) + p.kudl_dldk*y(v1b+2900)+ p.kdeg_plc*(y(v1b+500)+2*y(v1b+600)+2*y(v1b+700)+y(v1b+800)+y(v1b+900)+y(v1b+1000)+2*y(v1b+1100)+y(v1b+1200)+y(v1b+1500)+y(v1b+1900));%+ p.kdeg_dl3*y(v1+3100);
end 

dydt(100) = -p.kdeg_dkp*y(100)  + p.kph1_pc*y(600) - p.kpc_ph1*y(100)*y(500)  + p.kpk0_pk1*y(800) -  p.kpk1_pk0*y(900)*y(100) + p.kpk1_pk2p*y(900) - p.kpk2p_pk1*y(1400)*y(100) + p.kpt2_pk1h*y(1200) - p.kpk1h_pt2*y(1300)*y(100) + p.kjkp_dj*y(2700) - p.kdj_jkp*y(100)*y(300) + p.kpk2_pk3*y(1100) - p.kpk3_pk2*y(1500)*y(100) + p.kpt3_pk1p*y(1600) - p.kpk1p_pt3*y(1700)*y(100) - p.kpk1p_1h*y(1700)*y(100) + p.k1h_pk1p*y(1300) - p.kph2_p12*y(2100)*y(100) + p.kp12_ph2*y(1000)  + p.kpt4_pk3h*y(2000) - p.kpk3h_pt4*y(2200)*y(100) + p.kh1_h2*y(600)*y(200) - p.kh2_h1*y(2100)*y(100) -  p.kph1_p11*y(600)*y(100) + p.kp11_ph1*y(700) +(p.D_dk/y(4901)^2)*(y(99)-y(100))  - p.kbdl_dldk*y(2900)*y(100) + p.kudl_dldk*y(3000)+ p.kdeg_plc*(y(600)+2*y(700)+2*y(800)+y(900)+y(1000)+y(1100)+2*y(1200)+y(1300)+y(1600)+y(2000));%+ p.kdeg_dl3*y(3200);
 
%% DIVk
dydt(101) = p.ksyn_dk - p.kdeg_dk*y(101)  - p.kpk1_pk2*y(801)*y(101) + p.kpk2_pk1*y(1001)  -p.kj_jk*y(201)*y(101) + p.kjk_j*y(301) + p.kpk3_pk2p*y(1401) - p.kpk2p_pk3*y(1301)*y(101) - p.kpk3_pk4*y(1401)*y(101) + p.kpk4_pk3*y(1801) - p.kpc_ph2*y(401)*y(101) + p.kph2_pc*y(2001)  - p.kpk1p_p3h*y(1601)*y(101) + p.kp3h_pk1p*y(2101) - p.kph1_p12*y(501)*y(101) + p.kp12_ph1*y(901) + p.kh2_h1*y(2001)*y(1) - p.kh1_h2*y(501)*y(101) - p.kph2_p22*y(2001)*y(101) + p.kp22_ph2*y(1701)+ (p.D_dkp/y(4901)^2)*(y(102)-y(101));
for v2a=102:149
    dydt(v2a) = p.ksyn_dk - p.kdeg_dk*y(v2a)  - p.kpk1_pk2*y(v2a+700)*y(v2a) + p.kpk2_pk1*y(v2a+900)  -p.kj_jk*y(v2a+100)*y(v2a) + p.kjk_j*y(v2a+200) + p.kpk3_pk2p*y(v2a+1300) - p.kpk2p_pk3*y(v2a+1200)*y(v2a) - p.kpk3_pk4*y(v2a+1300)*y(v2a) + p.kpk4_pk3*y(v2a+1700) - p.kpc_ph2*y(v2a+300)*y(v2a) + p.kph2_pc*y(v2a+1900)  - p.kpk1p_p3h*y(v2a+1500)*y(v2a) + p.kp3h_pk1p*y(v2a+2000) - p.kph1_p12*y(v2a+400)*y(v2a) + p.kp12_ph1*y(v2a+800) + p.kh2_h1*y(v2a+1900)*y(v2a-100) - p.kh1_h2*y(v2a+400)*y(v2a) - p.kph2_p22*y(v2a+1900)*y(v2a) + p.kp22_ph2*y(v2a+1600)+ (p.D_dkp/y(4901)^2)*(y(v2a-1)-2*y(v2a)+y(v2a+1));
end

dydt(150) = p.ksyn_dk - p.kdeg_dk*y(150)  - p.kpk1_pk2*y(850)*y(150) + p.kpk2_pk1*y(1050)  -p.kj_jk*y(250)*y(150) + p.kjk_j*y(350) + p.kpk3_pk2p*y(1450) - p.kpk2p_pk3*y(1350)*y(150) - p.kpk3_pk4*y(1450)*y(150) + p.kpk4_pk3*y(1850) - p.kpc_ph2*y(450)*y(150) + p.kph2_pc*y(2050)  - p.kpk1p_p3h*y(1650)*y(150) + p.kp3h_pk1p*y(2150) - p.kph1_p12*y(550)*y(150) + p.kp12_ph1*y(950) + p.kh2_h1*y(2050)*y(50) - p.kh1_h2*y(550)*y(150) - p.kph2_p22*y(2050)*y(150) + p.kp22_ph2*y(1750)+ (p.D_dkp/y(4901)^2)*(y(149)-y(150));
dydt(151) = p.ksyn_dk - p.kdeg_dk*y(151)  - p.kpk1_pk2*y(851)*y(151) + p.kpk2_pk1*y(1051)  -p.kj_jk*y(251)*y(151) + p.kjk_j*y(351) + p.kpk3_pk2p*y(1451) - p.kpk2p_pk3*y(1351)*y(151) - p.kpk3_pk4*y(1451)*y(151) + p.kpk4_pk3*y(1851) - p.kpc_ph2*y(451)*y(151) + p.kph2_pc*y(2051)  - p.kpk1p_p3h*y(1651)*y(151) + p.kp3h_pk1p*y(2151) - p.kph1_p12*y(551)*y(151) + p.kp12_ph1*y(951) + p.kh2_h1*y(2051)*y(51) - p.kh1_h2*y(551)*y(151) - p.kph2_p22*y(2051)*y(151) + p.kp22_ph2*y(1751)+ (p.D_dkp/y(4901)^2)*(y(152)-y(151));
for v2b=152:199
    dydt(v2b) = p.ksyn_dk - p.kdeg_dk*y(v2b)  - p.kpk1_pk2*y(v2b+700)*y(v2b) + p.kpk2_pk1*y(v2b+900)  -p.kj_jk*y(v2b+100)*y(v2b) + p.kjk_j*y(v2b+200) + p.kpk3_pk2p*y(v2b+1300) - p.kpk2p_pk3*y(v2b+1200)*y(v2b) - p.kpk3_pk4*y(v2b+1300)*y(v2b) + p.kpk4_pk3*y(v2b+1700) - p.kpc_ph2*y(v2b+300)*y(v2b) + p.kph2_pc*y(v2b+1900)  - p.kpk1p_p3h*y(v2b+1500)*y(v2b) + p.kp3h_pk1p*y(v2b+2000) - p.kph1_p12*y(v2b+400)*y(v2b) + p.kp12_ph1*y(v2b+800) + p.kh2_h1*y(v2b+1900)*y(v2b-100) - p.kh1_h2*y(v2b+400)*y(v2b) - p.kph2_p22*y(v2b+1900)*y(v2b) + p.kp22_ph2*y(v2b+1600)+ (p.D_dkp/y(4901)^2)*(y(v2b-1)-2*y(v2b)+y(v2b+1));
end

dydt(200) = p.ksyn_dk - p.kdeg_dk*y(200)  - p.kpk1_pk2*y(900)*y(200) + p.kpk2_pk1*y(1100)  -p.kj_jk*y(300)*y(200) + p.kjk_j*y(400) + p.kpk3_pk2p*y(1500) - p.kpk2p_pk3*y(1400)*y(200) - p.kpk3_pk4*y(1500)*y(200) + p.kpk4_pk3*y(1900) - p.kpc_ph2*y(500)*y(200) + p.kph2_pc*y(2100)  - p.kpk1p_p3h*y(1700)*y(200) + p.kp3h_pk1p*y(2200) - p.kph1_p12*y(600)*y(200) + p.kp12_ph1*y(1000) + p.kh2_h1*y(2100)*y(100) - p.kh1_h2*y(600)*y(200) - p.kph2_p22*y(2100)*y(200) + p.kp22_ph2*y(1800)+ (p.D_dkp/y(4901)^2)*(y(199)-y(200));

%%DIVJ
dydt(201) =  p.kdj_djp*y(2201)*y(2301) - p.kdjp_dj*y(201) - p.kj_jk*y(201)*y(101) + p.kjk_j*y(301) + p.kjkp_dj*y(2601) - p.kdj_jkp*y(1)*y(201) - p.kdegpp2*y(201);
for v3=202:299
    dydt(v3) =  p.kdj_djp*y(v3+2000)*y(v3+2100) - p.kdjp_dj*y(v3) - p.kj_jk*y(v3)*y(v3-100) + p.kjk_j*y(v3+100) + p.kjkp_dj*y(v3+2400) - p.kdj_jkp*y(v3-200)*y(v3)- p.kdegpp2*y(v3);
end
dydt(300) = p.kdj_djp*y(2300)*y(2400) - p.kdjp_dj*y(300) - p.kj_jk*y(300)*y(200) + p.kjk_j*y(400) + p.kjkp_dj*y(2700) - p.kdj_jkp*y(100)*y(300)- p.kdegpp2*y(300);

%%DJk
dydt(301) =   p.kj_jk*y(201)*y(101) - p.kjk_j*y(301) - p.kjk_jkp*y(301) + p.kjkp_jk*y(2601)- p.kdegpp2*y(301);
for v4=302:399
    dydt(v4) =  p.kj_jk*y(v4-100)*y(v4-200) - p.kjk_j*y(v4) - p.kjk_jkp*y(v4) + p.kjkp_jk*y(v4+2300)- p.kdegpp2*y(v4);
end
dydt(400) =  p.kj_jk*y(300)*y(200) - p.kjk_j*y(400) - p.kjk_jkp*y(400) + p.kjkp_jk*y(2700)- p.kdegpp2*y(400);
%% PLEC
dydt(401) =  - p.kdeg_plc*y(401)- p.kpc_ph1*y(401)*y(1) + p.kph1_pc*y(501)  + p.kpk2p_pc*y(1301) +  p.kpk1p_pc*y(1601) - p.kpc_ph2*y(401)*y(101) + p.kph2_pc*y(2001)+ p.kplc_plcb*y(2401)*y(2501) - p.kplcb_plc*y(401);
for v5=402:499
    dydt(v5) =  - p.kdeg_plc*y(v5)- p.kpc_ph1*y(v5)*y(v5-400) + p.kph1_pc*y(v5+100)  + p.kpk2p_pc*y(v5+900) +  p.kpk1p_pc*y(v5+1200) - p.kpc_ph2*y(v5)*y(v5-300) + p.kph2_pc*y(v5+1600)+ p.kplc_plcb*y(v5+2000)*y(v5+2100) - p.kplcb_plc*y(v5);
end
dydt(500) =  - p.kdeg_plc*y(500)- p.kpc_ph1*y(500)*y(100) + p.kph1_pc*y(600)  + p.kpk2p_pc*y(1400) +  p.kpk1p_pc*y(1700) - p.kpc_ph2*y(500)*y(200) + p.kph2_pc*y(2100)+ p.kplc_plcb*y(2500)*y(2600) - p.kplcb_plc*y(500);
%% PLECH1
dydt(501) =  p.kpc_ph1*y(401)*y(1) - p.kph1_pc*y(501) -p.kph1_p11*y(501)*y(1) + p.kp11_ph1*y(601)  - p.kph1_p12*y(501)*y(101) + p.kp12_ph1*y(901) + p.kh2_h1*y(2001)*y(1) - p.kh1_h2*y(501)*y(101) - p.kph1_ph2*y(501) + p.kph2_ph1*y(2001) - p.kdeg_plc*y(501);
for v6=502:599
    dydt(v6) =  p.kpc_ph1*y(v6-100)*y(v6-500) - p.kph1_pc*y(v6) -p.kph1_p11*y(v6)*y(v6-500) + p.kp11_ph1*y(v6+100)  - p.kph1_p12*y(v6)*y(v6-400) + p.kp12_ph1*y(v6+400) + p.kh2_h1*y(v6+1500)*y(v6-500) - p.kh1_h2*y(v6)*y(v6-400) - p.kph1_ph2*y(v6) + p.kph2_ph1*y(v6+1500) - p.kdeg_plc*y(v6);
end
dydt(600) =  p.kpc_ph1*y(500)*y(100) - p.kph1_pc*y(600) -p.kph1_p11*y(600)*y(100) + p.kp11_ph1*y(700)  - p.kph1_p12*y(600)*y(200) + p.kp12_ph1*y(1000) + p.kh2_h1*y(2100)*y(100) - p.kh1_h2*y(600)*y(200) - p.kph1_ph2*y(600) + p.kph2_ph1*y(2100) - p.kdeg_plc*y(600);
%% Pk11
dydt(601) = p.kph1_p11*y(501)*y(1)- p.kp11_ph1*y(601) + p.kpt4_p11*y(1901) - p.kp11_pt4*y(601) - p.kp11_pk0*y(601) + p.kpk0_p11*y(701) - p.kdeg_plc*y(601);
for v7=602:699
    dydt(v7) = p.kph1_p11*y(v7-100)*y(v7-600)- p.kp11_ph1*y(v7) + p.kpt4_p11*y(v7+1300) - p.kp11_pt4*y(v7) - p.kp11_pk0*y(v7) + p.kpk0_p11*y(v7+100) - p.kdeg_plc*y(v7);
end
dydt(700) = p.kph1_p11*y(600)*y(100)- p.kp11_ph1*y(700) + p.kpt4_p11*y(2000) - p.kp11_pt4*y(700) - p.kp11_pk0*y(700) + p.kpk0_p11*y(800) - p.kdeg_plc*y(700);
%% Pk0
dydt(701) = p.kp11_pk0*y(601) - p.kpk0_p11*y(701) - p.kpk0_pk1*y(701) + p.kpk1_pk0*y(801)*y(1) - p.kdeg_plc*y(701);
for v8=702:799
    dydt(v8) = p.kp11_pk0*y(v8-100) - p.kpk0_p11*y(v8) - p.kpk0_pk1*y(v8) + p.kpk1_pk0*y(v8+100)*y(v8-700) - p.kdeg_plc*y(v8);
end
dydt(800) = p.kp11_pk0*y(700) - p.kpk0_p11*y(800) - p.kpk0_pk1*y(800) + p.kpk1_pk0*y(900)*y(100) - p.kdeg_plc*y(800);
%%Pk1
dydt(801) =  p.kpk0_pk1*y(701) -  p.kpk1_pk0*y(801)*y(1)  - p.kpk1_pk2*y(801)*y(101) + p.kpk2_pk1*y(1001) - p.kpk1_pk2p*y(801) + p.kpk2p_pk1*y(1301)*y(1) -  p.kpk1_pk1h*y(801) + p.kpk1h_pk1*y(1201) - p.kdeg_plc*y(801);% - p.kpk1_pk5*y(801)*y(2201) +  p.kpk5_pk1*y(2401)
for v9=802:899
    dydt(v9) =  p.kpk0_pk1*y(v9-100) -  p.kpk1_pk0*y(v9)*y(v9-800)  - p.kpk1_pk2*y(v9)*y(v9-700) + p.kpk2_pk1*y(v9+200) - p.kpk1_pk2p*y(v9) + p.kpk2p_pk1*y(v9+500)*y(v9-800) -  p.kpk1_pk1h*y(v9) + p.kpk1h_pk1*y(v9+400) - p.kdeg_plc*y(v9); % - p.kpk1_pk5*y(v9)*y(v9+1400) +  p.kpk5_pk1*y(v9+1600)
end
dydt(900) =  p.kpk0_pk1*y(800) -  p.kpk1_pk0*y(900)*y(100)  - p.kpk1_pk2*y(900)*y(200) + p.kpk2_pk1*y(1100) - p.kpk1_pk2p*y(900) + p.kpk2p_pk1*y(1400)*y(100) -  p.kpk1_pk1h*y(900) + p.kpk1h_pk1*y(1300)  - p.kdeg_plc*y(900);%- p.kpk1_pk5*y(900)*y(2300) +  p.kpk5_pk1*y(2500)
%%Pk12
dydt(901) = p.kph1_p12*y(501)*y(101)- p.kp12_ph1*y(901) + p.kph2_p12*y(2001)*y(1)- p.kp12_ph2*y(901)  - p.kp12_pk2*y(901) + p.kpk2_p12*y(1001) - p.kdeg_plc*y(901);
for v10=902:999
    dydt(v10) = p.kph1_p12*y(v10-400)*y(v10-800)- p.kp12_ph1*y(v10) + p.kph2_p12*y(v10+1100)*y(v10-900)- p.kp12_ph2*y(v10)  - p.kp12_pk2*y(v10) + p.kpk2_p12*y(v10+100) - p.kdeg_plc*y(v10);
end
dydt(1000) = p.kph1_p12*y(600)*y(200)- p.kp12_ph1*y(1000) + p.kph2_p12*y(2100)*y(100)- p.kp12_ph2*y(1000)  - p.kp12_pk2*y(1000) + p.kpk2_p12*y(1100) - p.kdeg_plc*y(1000);
%% Pk2
dydt(1001) = p.kpk1_pk2*y(801)*y(101) - p.kpk2_pk1*y(1001) + p.kp12_pk2*y(901) - p.kpk2_p12*y(1001) -p.kpk2_pt2*y(1001) + p.kpt2_pk2*y(1101) - p.kpk2_pk3*y(1001) + p.kpk3_pk2*y(1401)*y(1) - p.kdeg_plc*y(1001); 
for v11=1002:1099
    dydt(v11) = p.kpk1_pk2*y(v11-200)*y(v11-900) - p.kpk2_pk1*y(v11) + p.kp12_pk2*y(v11-100) - p.kpk2_p12*y(v11) -p.kpk2_pt2*y(v11) + p.kpt2_pk2*y(v11+100) - p.kpk2_pk3*y(v11) + p.kpk3_pk2*y(v11+400)*y(v11-1000) - p.kdeg_plc*y(v11);
end
dydt(1100) = p.kpk1_pk2*y(900)*y(200) - p.kpk2_pk1*y(1100) + p.kp12_pk2*y(1000) - p.kpk2_p12*y(1100) -p.kpk2_pt2*y(1100) + p.kpt2_pk2*y(1200) - p.kpk2_pk3*y(1100) + p.kpk3_pk2*y(1500)*y(100) - p.kdeg_plc*y(1100);
%% PT2
dydt(1101) = p.kpk2_pt2*y(1001) - p.kpt2_pk2*y(1101) -p.kpt2_pk1h*y(1101) + p.kpk1h_pt2*y(1201)*y(1) - p.kdeg_plc*y(1101);
for v12=1102:1199
    dydt(v12) = p.kpk2_pt2*y(v12-100) - p.kpt2_pk2*y(v12) -p.kpt2_pk1h*y(v12) + p.kpk1h_pt2*y(v12+100)*y(v12-1100) - p.kdeg_plc*y(v12);
end
dydt(1200) = p.kpk2_pt2*y(1100) - p.kpt2_pk2*y(1200) -p.kpt2_pk1h*y(1200) + p.kpk1h_pt2*y(1300)*y(100) - p.kdeg_plc*y(1200);
%% Pk1H
dydt(1201) =  p.kpt2_pk1h*y(1101) - p.kpk1h_pt2*y(1201)*y(1) + p.kpk1_pk1h*y(801) - p.kpk1h_pk1*y(1201) + p.kpk1p_1h*y(1601)*y(1) - p.k1h_pk1p*y(1201) - p.kdeg_plc*y(1201);% + p.kpt5_pk1h*y(2601) - p.kpk1h_pt5*y(1201)*y(2301)
for v13=1202:1299
    dydt(v13) =  p.kpt2_pk1h*y(v13-100) - p.kpk1h_pt2*y(v13)*y(v13-1200) + p.kpk1_pk1h*y(v13-400) - p.kpk1h_pk1*y(v13) + p.kpk1p_1h*y(v13+400)*y(v13-1200) - p.k1h_pk1p*y(v13) - p.kdeg_plc*y(v13);%+ p.kpt5_pk1h*y(v13+1400) - p.kpk1h_pt5*y(v13)*y(v13+1100)
end
dydt(1300) =  p.kpt2_pk1h*y(1200) - p.kpk1h_pt2*y(1300)*y(100) + p.kpk1_pk1h*y(900) - p.kpk1h_pk1*y(1300) + p.kpk1p_1h*y(1700)*y(100) - p.k1h_pk1p*y(1300)  - p.kdeg_plc*y(1300);%+ p.kpt5_pk1h*y(2700) - p.kpk1h_pt5*y(1300)*y(2400)
%% Pk2P
dydt(1301) = p.kpk1_pk2p*y(801) - p.kpk2p_pk1*y(1301)*y(1) - p.kpk2p_pc*y(1301)  + p.kpk3_pk2p*y(1401) - p.kpk2p_pk3*y(1301)*y(101) - p.kdeg_plc*y(1301); 
for v14=1302:1399
    dydt(v14) = p.kpk1_pk2p*y(v14-500) - p.kpk2p_pk1*y(v14)*y(v14-1300) - p.kpk2p_pc*y(v14)  + p.kpk3_pk2p*y(v14+100) - p.kpk2p_pk3*y(v14)*y(v14-1200) - p.kdeg_plc*y(v14);
end
dydt(1400) = p.kpk1_pk2p*y(900) - p.kpk2p_pk1*y(1400)*y(100) - p.kpk2p_pc*y(1400)  + p.kpk3_pk2p*y(1500) - p.kpk2p_pk3*y(1400)*y(200) - p.kdeg_plc*y(1400);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pk3
dydt(1401) = p.kpk2_pk3*y(1001) - p.kpk3_pk2*y(1401)*y(1) -  p.kpk3_pt3*y(1401) + p.kpt3_pk3*y(1501) - p.kpk3_pk2p*y(1401) + p.kpk2p_pk3*y(1301)*y(101) - p.kpk3_pk4*y(1401)*y(101) + p.kpk4_pk3*y(1801) - p.kpk3_pk3h*y(1401) + p.kpk3h_pk3*y(2101)  - p.kdeg_plc*y(1401);%- p.kpk3_pk6*y(1401)*y(2201) +  p.kpk6_pk3*y(2501)
for v15=1402:1499
    dydt(v15) = p.kpk2_pk3*y(v15-400) - p.kpk3_pk2*y(v15)*y(v15-1400) -  p.kpk3_pt3*y(v15) + p.kpt3_pk3*y(v15+100) - p.kpk3_pk2p*y(v15) + p.kpk2p_pk3*y(v15-100)*y(v15-1300) - p.kpk3_pk4*y(v15)*y(v15-1300) + p.kpk4_pk3*y(v15+400) - p.kpk3_pk3h*y(v15) + p.kpk3h_pk3*y(v15+700) - p.kdeg_plc*y(v15);%- p.kpk3_pk6*y(v15)*y(v15+800) +  p.kpk6_pk3*y(v15+1100)
end
dydt(1500) = p.kpk2_pk3*y(1100) - p.kpk3_pk2*y(1500)*y(100) -  p.kpk3_pt3*y(1500) + p.kpt3_pk3*y(1600) - p.kpk3_pk2p*y(1500) + p.kpk2p_pk3*y(1400)*y(200) - p.kpk3_pk4*y(1500)*y(200) + p.kpk4_pk3*y(1900) - p.kpk3_pk3h*y(1500) + p.kpk3h_pk3*y(2200)  - p.kdeg_plc*y(1500);%- p.kpk3_pk6*y(1500)*y(2300) +  p.kpk6_pk3*y(2600)
%% PT3
dydt(1501) = p.kpk3_pt3*y(1401) - p.kpt3_pk3*y(1501) + p.kpk1p_pt3*y(1601)*y(1) - p.kpt3_pk1p*y(1501) - p.kdeg_plc*y(1501);
for v16=1502:1599
    dydt(v16) = p.kpk3_pt3*y(v16-100) - p.kpt3_pk3*y(v16) + p.kpk1p_pt3*y(v16+100)*y(v16-1500) - p.kpt3_pk1p*y(v16) - p.kdeg_plc*y(v16);
end
dydt(1600) = p.kpk3_pt3*y(1500) - p.kpt3_pk3*y(1600) + p.kpk1p_pt3*y(1700)*y(100) - p.kpt3_pk1p*y(1600) - p.kdeg_plc*y(1600);
%% Pk1P
dydt(1601) = p.kpt3_pk1p*y(1501) - p.kpk1p_pt3*y(1601)*y(1) -p.kpk1p_1h*y(1601)*y(1) + p.k1h_pk1p*y(1201) - p.kpk1p_pc*y(1601) - p.kpk1p_p3h*y(1601)*y(101) +p.kp3h_pk1p*y(2101) - p.kdeg_plc*y(1601);
for v17=1602:1699
    dydt(v17) = p.kpt3_pk1p*y(v17-100) - p.kpk1p_pt3*y(v17)*y(v17-1600) -p.kpk1p_1h*y(v17)*y(v17-1600) + p.k1h_pk1p*y(v17-400) - p.kpk1p_pc*y(v17) - p.kpk1p_p3h*y(v17)*y(v17-1500) +p.kp3h_pk1p*y(v17+500) - p.kdeg_plc*y(v17);
end
dydt(1700) = p.kpt3_pk1p*y(1600) - p.kpk1p_pt3*y(1700)*y(100) -p.kpk1p_1h*y(1700)*y(100) + p.k1h_pk1p*y(1300) - p.kpk1p_pc*y(1700) - p.kpk1p_p3h*y(1700)*y(200) +p.kp3h_pk1p*y(2200) - p.kdeg_plc*y(1700);


%%% Pk22 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dydt(1701) = p.kph2_p22*y(2001)*y(101) - p.kp22_ph2*y(1701) - p.kp22_pk4*y(1701) + p.kpk4_p22*y(1801) - p.kdeg_plc*y(1701);
for v18=1702:1799
    dydt(v18) = p.kph2_p22*y(v18+300)*y(v18-1600) - p.kp22_ph2*y(v18) - p.kp22_pk4*y(v18) + p.kpk4_p22*y(v18+100) - p.kdeg_plc*y(v18);
end
dydt(1800) = p.kph2_p22*y(2100)*y(200) - p.kp22_ph2*y(1800) - p.kp22_pk4*y(1800) + p.kpk4_p22*y(1900) - p.kdeg_plc*y(1800);
%% Pk4
dydt(1801) = p.kpk3_pk4*y(1401)*y(101) - p.kpk4_pk3*y(1801) - p.kpk4_pt4*y(1801) + p.kpt4_pk4*y(1901) + p.kp22_pk4*y(1701) - p.kpk4_p22*y(1801) - p.kdeg_plc*y(1801);
for v19=1802:1899
    dydt(v19) = p.kpk3_pk4*y(v19-400)*y(v19-1700) - p.kpk4_pk3*y(v19) - p.kpk4_pt4*y(v19) + p.kpt4_pk4*y(v19+100) + p.kp22_pk4*y(v19-100) - p.kpk4_p22*y(v19) - p.kdeg_plc*y(v19);
end
dydt(1900) = p.kpk3_pk4*y(1500)*y(200) - p.kpk4_pk3*y(1900) - p.kpk4_pt4*y(1900) + p.kpt4_pk4*y(2000) + p.kp22_pk4*y(1800) - p.kpk4_p22*y(1900) - p.kdeg_plc*y(1900);
%% PT4
dydt(1901) = p.kpk4_pt4*y(1801) - p.kpt4_pk4*y(1901) - p.kpt4_p11*y(1901) + p.kp11_pt4*y(601)- p.kpt4_pk3h*y(1901) + p.kpk3h_pt4*y(2101)*y(1) - p.kdeg_plc*y(1901);
for v20=1902:1999
    dydt(v20) = p.kpk4_pt4*y(v20-100) - p.kpt4_pk4*y(v20) - p.kpt4_p11*y(v20) + p.kp11_pt4*y(v20-1300)- p.kpt4_pk3h*y(v20) + p.kpk3h_pt4*y(v20+200)*y(v20-1900) - p.kdeg_plc*y(v20);
end
dydt(2000) = p.kpk4_pt4*y(1900) - p.kpt4_pk4*y(2000) - p.kpt4_p11*y(2000) + p.kp11_pt4*y(700)- p.kpt4_pk3h*y(2000) + p.kpk3h_pt4*y(2200)*y(100) - p.kdeg_plc*y(2000);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLECH2
dydt(2001) = p.kpc_ph2*y(401)*y(101) - p.kph2_pc*y(2001) - p.kph2_p22*y(2001)*y(101) + p.kp22_ph2*y(1701)- p.kph2_p12*y(2001)*y(1) + p.kp12_ph1*y(901) + p.kh1_h2*y(501)*y(101) - p.kh2_h1*y(2001)*y(1)  +  p.kph1_ph2*y(501) - p.kph2_ph1*y(2001) - p.kdeg_plc*y(2001);
for v21=2002:2099
    dydt(v21) = p.kpc_ph2*y(v21-1600)*y(v21-1900) - p.kph2_pc*y(v21) - p.kph2_p22*y(v21)*y(v21-1900) + p.kp22_ph2*y(v21-300)- p.kph2_p12*y(v21)*y(v21-2000) + p.kp12_ph1*y(v21-1100) + p.kh1_h2*y(v21-1500)*y(v21-1900) - p.kh2_h1*y(v21)*y(v21-2000)  +  p.kph1_ph2*y(v21-1500) - p.kph2_ph1*y(v21) - p.kdeg_plc*y(v21);
end
dydt(2100) = p.kpc_ph2*y(500)*y(200) - p.kph2_pc*y(2100) - p.kph2_p22*y(2100)*y(200) + p.kp22_ph2*y(1800)- p.kph2_p12*y(2100)*y(100) + p.kp12_ph1*y(1000) + p.kh1_h2*y(600)*y(200) - p.kh2_h1*y(2100)*y(100)  +  p.kph1_ph2*y(600) - p.kph2_ph1*y(2100) - p.kdeg_plc*y(2100);

%% Pk3H
dydt(2101) = p.kpt4_pk3h*y(1901) - p.kpk3h_pt4*y(2101)*y(1) + p.kpk3_pk3h*y(1401) - p.kpk3h_pk3*y(2101) + p.kpk1p_p3h*y(1601)*y(101) - p.kp3h_pk1p*y(2101)  - p.kdeg_plc*y(2101);%+ p.kpt6_pk3h*y(2701)- p.kpk3h_pt6*y(2101)*y(2301)
for v22=2102:2199
    dydt(v22) = p.kpt4_pk3h*y(v22-200) - p.kpk3h_pt4*y(v22)*y(v22-2100) + p.kpk3_pk3h*y(v22-700) - p.kpk3h_pk3*y(v22) + p.kpk1p_p3h*y(v22-500)*y(v22-2000) - p.kp3h_pk1p*y(v22) - p.kdeg_plc*y(v22);%+ p.kpt6_pk3h*y(v22+600)- p.kpk3h_pt6*y(v22)*y(v22+200)
end
dydt(2200) = p.kpt4_pk3h*y(2000) - p.kpk3h_pt4*y(2200)*y(100) + p.kpk3_pk3h*y(1500) - p.kpk3h_pk3*y(2200) + p.kpk1p_p3h*y(1700)*y(200) - p.kp3h_pk1p*y(2200) - p.kdeg_plc*y(2200);% + p.kpt6_pk3h*y(2800)- p.kpk3h_pt6*y(2200)*y(2400)

% equations for DivJ

dydt(2201)= p.ksyn_divj - p.kdeg_divj*y(2201) + p.D_divj*(y(2202)-y(2201))/(y(4901)^2) - p.growth*y(2201) - p.kdj_djp*y(2201)*y(2301) + p.kdjp_dj*y(201);

for v23a=2202:2249
dydt(v23a) = p.ksyn_divj - p.kdeg_divj*y(v23a) + p.D_divj*(y(v23a+1)-2*y(v23a) + y(v23a-1))/(y(4901)^2) - p.growth*y(v23a) - p.kdj_djp*y(v23a)*y(v23a+100) + p.kdjp_dj*y(v23a-2000);
end
dydt(2250) = p.ksyn_divj - p.kdeg_divj*y(2250) + p.D_divj*(y(2249) - y(2250))/(y(4901)^2) - p.growth*y(2250)- p.kdj_djp*y(2250)*y(2350) + p.kdjp_dj*y(250);
dydt(2251) = p.ksyn_divj - p.kdeg_divj*y(2251) + p.D_divj*(y(2252) - y(2251))/(y(4901)^2) - p.growth*y(2251)- p.kdj_djp*y(2251)*y(2351) + p.kdjp_dj*y(251);
for v23b=2252:2299
dydt(v23b) = p.ksyn_divj - p.kdeg_divj*y(v23b) + p.D_divj*(y(v23b+1)-2*y(v23b) + y(v23b-1))/(y(4901)^2) - p.growth*y(v23b) - p.kdj_djp*y(v23b)*y(v23b+100) + p.kdjp_dj*y(v23b-2000);
end

dydt(2300) = p.ksyn_divj - p.kdeg_divj*y(2300) + p.D_divj*(y(2299) - y(2300))/(y(4901)^2) - p.growth*y(2300)- p.kdj_djp*y(2300)*y(2400) + p.kdjp_dj*y(300);

% DivJ sticky
dydt(2301:2400) = 0;

% PleC
dydt(2401) = p.ksyn_plc2 - p.kdeg_plc2*y(2401) - p.kplc_plcb*y(2401)*y(2501) + p.kplcb_plc*y(401) + p.D_plc*(y(2402)-y(2401))/(y(4901)^2);% + p.kprot2*y(6201)*y(5101);
for v25a=2402:2449
    dydt(v25a) = p.ksyn_plc2 - p.kdeg_plc2*y(v25a) - p.kplc_plcb*y(v25a)*y(v25a+100) + p.kplcb_plc*y(v25a-2000) + p.D_plc*(y(v25a-1)-2*y(v25a)+y(v25a+1))/(y(4901)^2);% +  p.kprot2*y(v25-100)*y(v25-1200);
end
dydt(2450) = p.ksyn_plc2 - p.kdeg_plc2*y(2450) - p.kplc_plcb*y(2450)*y(2550) + p.kplcb_plc*y(450) + p.D_plc*(y(2449)-y(2450))/(y(4901)^2);
dydt(2451) = p.ksyn_plc2 - p.kdeg_plc2*y(2451) - p.kplc_plcb*y(2451)*y(2551) + p.kplcb_plc*y(451) + p.D_plc*(y(2452)-y(2451))/(y(4901)^2);
for v25b=2452:2499
    dydt(v25b) = p.ksyn_plc2 - p.kdeg_plc2*y(v25b) - p.kplc_plcb*y(v25b)*y(v25b+100) + p.kplcb_plc*y(v25b-2000) + p.D_plc*(y(v25b-1)-2*y(v25b)+y(v25b+1))/(y(4901)^2);% +  p.kprot2*y(v25-100)*y(v25-1200);
end

dydt(2500) = p.ksyn_plc2 - p.kdeg_plc2*y(2500) - p.kplc_plcb*y(2500)*y(2600) + p.kplcb_plc*y(500) + p.D_plc*(y(2499)-y(2500))/(y(4901)^2);% + p.kprot2*y(6300)*y(5200);

% PleC sticky
dydt(2501:2600) = 0;

%% future equations
% DIvJ:DivK-P
for v27=2601:2700
dydt(v27) = p.kjk_jkp*y(v27-2300) - p.kjkp_jk*y(v27) - p.kjkp_dj*y(v27) + p.kdj_jkp*y(v27-2600)*y(v27-2400) - p.kdegpp2*y(v27);
end

% DivLfree
dydt(2701) = p.ksyn_dl - p.kdeg_dl*y(2701) - p.kdlf_dlb*y(3001)*y(2701) + p.kdlb_dlf*y(2801)+ p.D_dl*(y(2702)-y(2701))/(y(4901)^2);%-kcp_ck*y(3401)*(y(3001)^2/(y(3001)^2+Kmdl^2)) + kck_cp*y(3501);
for v28a=2702:2749
    dydt(v28a) = p.ksyn_dl - p.kdeg_dl*y(v28a)- p.kdlf_dlb*y(v28a+300)*y(v28a) + p.kdlb_dlf*y(v28a+100)+ p.D_dl*(y(v28a-1)-2*y(v28a)+y(v28a+1))/(y(4901)^2);%-kcp_ck*y(v31+400)*(y(v31)^2/(y(v31)^2+Kmdl^2)) + kck_cp*y(v31+500);
end
dydt(2750) = p.ksyn_dl - p.kdeg_dl*y(2750)- p.kdlf_dlb*y(3050)*y(2750) + p.kdlb_dlf*y(2850)+ p.D_dl*(y(2749)-y(2750))/(y(4901)^2);
dydt(2751) = p.ksyn_dl - p.kdeg_dl*y(2751)- p.kdlf_dlb*y(3051)*y(2751) + p.kdlb_dlf*y(2851)+ p.D_dl*(y(2752)-y(2751))/(y(4901)^2);
for v28b=2752:2799
    dydt(v28b) = p.ksyn_dl - p.kdeg_dl*y(v28b)- p.kdlf_dlb*y(v28b+300)*y(v28b) + p.kdlb_dlf*y(v28b+100)+ p.D_dl*(y(v28b-1)-2*y(v28b)+y(v28b+1))/(y(4901)^2);%-kcp_ck*y(v31+400)*(y(v31)^2/(y(v31)^2+Kmdl^2)) + kck_cp*y(v31+500);
end
dydt(2800) = p.ksyn_dl - p.kdeg_dl*y(2800)- p.kdlf_dlb*y(3100)*y(2800) + p.kdlb_dlf*y(2900)+ p.D_dl*(y(2799)-y(2800))/(y(4901)^2);%-kcp_ck*y(3500)*(y(3100)^2/(y(3100)^2+Kmdl^2)) + kck_cp*y(3600);

% DivL bound
for v29 = 2801:2900
    dydt(v29) = - p.kdeg_dl*y(v29) + p.kdlf_dlb*y(v29+200)*y(v29-100) - p.kdlb_dlf*y(v29) - p.kbdl_dldk*y(v29)*y(v29-2800) + p.kudl_dldk*y(v29+100);
end

% DivL:DivK~P
for v30=2901:3000
    dydt(v30) = - p.kdeg_dl*y(v30) + p.kbdl_dldk*y(v30-100)*y(v30-2900) - p.kudl_dldk*y(v30);
end

% DivL Sticky
dydt(3001:3100) = 0;

% CtrA and CtrA-P
%dydt(3101:3300) = 0;
%% CTRA
dydt(3101) =  p.ksyn_ctr- p.kdeg_ctr*y(3101)+ (p.D_ctr/y(4901)^2)*(y(3102)-y(3101)) + p.kctr_phos*y(3201)*y(3401) - p.kctr_kin*y(3101)*y(3501);
for v32a=3102:3149
    dydt(v32a) =  p.ksyn_ctr - p.kdeg_ctr*y(v32a)+ (p.D_ctr/y(4901)^2)*(y(v32a-1)-2*y(v32a)+y(v32a+1))+ p.kctr_phos*y(v32a+100)*y(v32a+300) - p.kctr_kin*y(v32a)*y(v32a+400);
end
dydt(3150) =  p.ksyn_ctr - p.kdeg_ctr*y(3150) + (p.D_ctr/y(4901)^2)*(y(3149)-y(3150))+ p.kctr_phos*y(3250)*y(3450) - p.kctr_kin*y(3150)*y(3550);
dydt(3151) =  p.ksyn_ctr - p.kdeg_ctr*y(3151)+ (p.D_ctr/y(4901)^2)*(y(3152)-y(3151))+ p.kctr_phos*y(3251)*y(3451) - p.kctr_kin*y(3151)*y(3551);
for v32b=3152:3199
    dydt(v32b) =  p.ksyn_ctr - p.kdeg_ctr*y(v32b) + (p.D_ctr/y(4901)^2)*(y(v32b-1)-2*y(v32b)+y(v32b+1))+ p.kctr_phos*y(v32b+100)*y(v32b+300) - p.kctr_kin*y(v32b)*y(v32b+400);
end
dydt(3200) =  p.ksyn_ctr - p.kdeg_ctr*y(3200) + (p.D_ctr/y(4901)^2)*(y(3199)-y(3200))+ p.kctr_phos*y(3300)*y(3500) - p.kctr_kin*y(3200)*y(3600);

%% CTRAP
dydt(3201) = -p.kdeg_ctr*y(3201) + (p.D_ctrp/y(4901)^2)*(y(3202)-y(3201)) - p.kctr_phos*y(3201)*y(3401) + p.kctr_kin*y(3101)*y(3501);
for v33a=3202:3249
    dydt(v33a) = -p.kdeg_ctr*y(v33a)  + (p.D_ctrp/y(4901)^2)*(y(v33a-1)-2*y(v33a)+y(v33a+1))-  p.kctr_phos*y(v33a)*y(v33a+200) + p.kctr_kin*y(v33a-100)*y(v33a+300);
end
dydt(3250) = -p.kdeg_ctr*y(3250) + (p.D_ctrp/y(4901)^2)*(y(3249)-y(3250)) -  p.kctr_phos*y(3250)*y(3450) + p.kctr_kin*y(3150)*y(3550);
dydt(3251) = -p.kdeg_ctr*y(3251) + (p.D_ctrp/y(4901)^2)*(y(3252)-y(3251)) -  p.kctr_phos*y(3251)*y(3451) + p.kctr_kin*y(3151)*y(3551);
for v33b=3252:3299
    dydt(v33b) = -p.kdeg_ctr*y(v33b) + (p.D_ctrp/y(4901)^2)*(y(v33b-1)-2*y(v33b)+y(v33b+1))-  p.kctr_phos*y(v33b)*y(v33b+200) + p.kctr_kin*y(v33b-100)*y(v33b+300);
end
dydt(3300) = -p.kdeg_ctr*y(3300) + (p.D_ctrp/y(4901)^2)*(y(3299)-y(3300)) -  p.kctr_phos*y(3300)*y(3500) + p.kctr_kin*y(3200)*y(3600);

% CckA free
dydt(3301) = p.ksyn_ccka - p.kdeg_ccka*y(3301) + p.D_cck*(y(3302) - y(3301))/(y(4901)^2) - p.kcf_cb*y(3301)*y(3901) + p.kcb_cf*y(3401);
for v34a = 3302:3349
    dydt(v34a) = p.ksyn_ccka - p.kdeg_ccka*y(v34a) + p.D_cck*(y(v34a-1) - 2*y(v34a) + y(v34a+1))/(y(4901)^2)- p.kcf_cb*y(v34a)*y(v34a+600) + p.kcb_cf*y(v34a+100);
end
dydt(3350) = p.ksyn_ccka - p.kdeg_ccka*y(3350) + p.D_cck*(y(3349) - y(3350))/(y(4901)^2)- p.kcf_cb*y(3350)*y(3950) + p.kcb_cf*y(3450);
dydt(3351) = p.ksyn_ccka - p.kdeg_ccka*y(3351) + p.D_cck*(y(3352) - y(3351))/(y(4901)^2)- p.kcf_cb*y(3351)*y(3951) + p.kcb_cf*y(3451);
for v34b = 3352:3399
    dydt(v34b) = p.ksyn_ccka - p.kdeg_ccka*y(v34b) + p.D_cck*(y(v34b-1) - 2*y(v34b) + y(v34b+1))/(y(4901)^2)- p.kcf_cb*y(v34b)*y(v34b+600) + p.kcb_cf*y(v34b+100);
end
dydt(3400) = p.ksyn_ccka - p.kdeg_ccka*y(3400) + p.D_cck*(y(3399) - y(3400))/(y(4901)^2)- p.kcf_cb*y(3400)*y(4000) + p.kcb_cf*y(3500);

% CCKAphos
%dydt(3401) = p.kcf_cb*y(3301)*y(3901) - p.kcb_cf*y(3401)  -p.kcp_ck*y(3401)*(y(2801)^p.h/(y(2801)^p.h+p.Kmdl^p.h)) + p.kck_cp*y(3501)- p.kdeg_ccka*y(3401);% - p.kbcp_cps2*y(3401)*y(4001) + p.kucp_cps2*y(4301) + p.ke_cps2*y(4301) ;
for v35=3401:3500
    dydt(v35) = p.kcf_cb*y(v35-100)*y(v35+500) - p.kcb_cf*y(v35)- p.kdeg_ccka*y(v35) -p.kcp_ck*y(v35)*(y(v35-600)^p.h/(y(v35-600)^p.h+p.Kmdl^p.h)) + p.kck_cp*y(v35+100);% - p.kbcp_cps2*y(v35)*y(v35+600) + p.kucp_cps2*y(v35+900) + p.ke_cps2*y(v35+900)- p.kdeg_ccka*y(v35);
end
%dydt(3500) =  p.kcf_cb*y(3400)*y(4000) - p.kcb_cf*y(3500)  -p.kcp_ck*y(3500)*(y(2900)^p.h/(y(2900)^p.h+p.Kmdl^p.h)) + p.kck_cp*y(3600)- p.kdeg_ccka*y(3500);% - p.kbcp_cps2*y(3500)*y(4100) + p.kucp_cps2*y(4400) + p.ke_cps2*y(4400)- p.kdeg_ccka*y(3500);
%% Cckakin
dydt(3501) = p.kcp_ck*y(3401)*(y(2801)^p.h/(y(2801)^p.h+p.Kmdl^p.h)) - p.kck_cp*y(3501)- p.kdeg_ccka*y(3501);% - p.kck_ck1*y(3501)*y(3201) + p.kck1_ck*y(3601) + p.kct1_ck1h*y(3701) - p.kck1h_ct1*y(3501)*y(3301)  - p.kck_ck2*y(3501)*y(3901) + p.kck2_ck*y(4101) + p.kct2_ck2h*y(4201) - p.kck2h_ct2*y(3501)*y(4001)- p.kdeg_ccka*y(3501);
    for v36=3502:3599
        dydt(v36) = p.kcp_ck*y(v36-100)*(y(v36-700)^p.h/(y(v36-700)^p.h+p.Kmdl^p.h)) - p.kck_cp*y(v36)- p.kdeg_ccka*y(v36);% - p.kck_cp*y(v36) - p.kck_ck1*y(v36)*y(v36-300) + p.kck1_ck*y(v36+100) + p.kct1_ck1h*y(v36+200) - p.kck1h_ct1*y(v36)*y(v36-200)  - p.kck_ck2*y(v36)*y(v36+400) + p.kck2_ck*y(v36+600) + p.kct2_ck2h*y(v36+700) - p.kck2h_ct2*y(v36)*y(v36+500)- p.kdeg_ccka*y(v36);
    end
dydt(3600) = p.kcp_ck*y(3500)*(y(2900)^p.h/(y(2900)^p.h+p.Kmdl^p.h)) - p.kck_cp*y(3600)- p.kdeg_ccka*y(3600);% - p.kck_ck1*y(3600)*y(3300) + p.kck1_ck*y(3700) + p.kct1_ck1h*y(3800) - p.kck1h_ct1*y(3600)*y(3400)  - p.kck_ck2*y(3600)*y(4000) + p.kck2_ck*y(4200) + p.kct2_ck2h*y(4300) - p.kck2h_ct2*y(3600)*y(4100)- p.kdeg_ccka*y(3600);

% %% CK1
% dydt(3601)= kck_ck1*y(3501)*y(3201) - kck1_ck*y(3601) - kck1_ct1*y(3601) + kct1_ck1*y(3701)- kdeg_ctr*y(3601);
% for v37=3602:3699
%     dydt(v37)= kck_ck1*y(v37-100)*y(v37-400) - kck1_ck*y(v37) - kck1_ct1*y(v37) + kct1_ck1*y(v37+100)- kdeg_ctr*y(v37);
% end
% dydt(3700)= kck_ck1*y(3600)*y(3300) - kck1_ck*y(3700) - kck1_ct1*y(3700) + kct1_ck1*y(3800)- kdeg_ctr*y(3700);
% 
% %% CT1
% 
% dydt(3701) = kck1_ct1*y(3601) - kct1_ck1*y(3701) - kct1_ck1h*y(3701) + kck1h_ct1*y(3501)*y(3301)- kdeg_ctr*y(3701);
% for v38=3702:3799
%     dydt(v38) = kck1_ct1*y(v38-100) - kct1_ck1*y(v38) - kct1_ck1h*y(v38) + kck1h_ct1*y(v38-200)*y(v38-400)- kdeg_ctr*y(v38);
% end
% dydt(3800) = kck1_ct1*y(3700) - kct1_ck1*y(3800) - kct1_ck1h*y(3800) + kck1h_ct1*y(3600)*y(3400)- kdeg_ctr*y(3800);
% 
% %% CPS
% dydt(3801) = kbcp_cps*y(3401)*y(3301) - kucp_cps*y(3801) - ke_cps*y(3801)- kdeg_ctr*y(3801);
%     for v39=3802:3899
%         dydt(v39) = kbcp_cps*y(v39-400)*y(v39-500) - kucp_cps*y(v39) - ke_cps*y(v39)- kdeg_ctr*y(v39);
%     end
% dydt(3900) = kbcp_cps*y(3500)*y(3400) - kucp_cps*y(3900) - ke_cps*y(3900)- kdeg_ctr*y(3900);

dydt(3601:3900) = 0;
%% CckA stcky
for v40 = 3901:4000
dydt(v40) = 0;%- p.kcf_cb*y(v40-600)*y(v40) + p.kcb_cf*y(v40-500);
end

dydt(4001:4900) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dydt(4901)=p.growth*y(4901);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PLED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
