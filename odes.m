function dydt = odes(t, y);

global p;

dydt = zeros(4901,1);

%% DIVkP
dydt(1) = -p.kdeg_dkp*y(1) + p.kph1_pc*y(501) - p.kpc_ph1*y(1)*y(401)   + p.kpk0_pk1*y(701) -  p.kpk1_pk0*y(801)*y(1) + p.kpk1_pk2p*y(801) - p.kpk2p_pk1*y(1301)*y(1) + p.kpt2_pk1h*y(1101) - p.kpk1h_pt2*y(1201)*y(1) + p.kjkp_dj*y(2601) - p.kdj_jkp*y(1)*y(201)+ p.kpk2_pk3*y(1001) - p.kpk3_pk2*y(1401)*y(1) + p.kpt3_pk1p*y(1501) - p.kpk1p_pt3*y(1601)*y(1) - p.kpk1p_1h*y(1601)*y(1) + p.k1h_pk1p*y(1201) - p.kph2_p12*y(2001)*y(1) + p.kp12_ph2*y(901)  + p.kpt4_pk3h*y(1901) - p.kpk3h_pt4*y(2101)*y(1) + p.kh1_h2*y(501)*y(101) - p.kh2_h1*y(2001)*y(1) -  p.kph1_p11*y(501)*y(1) + p.kp11_ph1*y(601) + (p.D_dk/y(4901)^2)*(y(2)-y(1))  - p.kbdl_dldk*y(2801)*y(1) + p.kudl_dldk*y(2901)+ p.kdeg_plc*(y(501)+2*y(601)+2*y(701)+y(801)+y(901)+y(1001)+2*y(1101)+y(1201)+y(1501)+y(1901));%+ p.kdegpp2*y(2601)+ p.kdeg_dl*y(2901);
for v1=2:99
    dydt(v1) = -p.kdeg_dkp*y(v1)  + p.kph1_pc*y(v1+500) - p.kpc_ph1*y(v1)*y(v1+400)  + p.kpk0_pk1*y(v1+700) -  p.kpk1_pk0*y(v1+800)*y(v1) + p.kpk1_pk2p*y(v1+800) - p.kpk2p_pk1*y(v1+1300)*y(v1) + p.kpt2_pk1h*y(v1+1100) - p.kpk1h_pt2*y(v1+1200)*y(v1) + p.kjkp_dj*y(v1+2600) - p.kdj_jkp*y(v1)*y(v1+200) + p.kpk2_pk3*y(v1+1000) - p.kpk3_pk2*y(v1+1400)*y(v1) + p.kpt3_pk1p*y(v1+1500) - p.kpk1p_pt3*y(v1+1600)*y(v1) - p.kpk1p_1h*y(v1+1600)*y(v1) + p.k1h_pk1p*y(v1+1200) - p.kph2_p12*y(v1+2000)*y(v1) + p.kp12_ph2*y(v1+900)  + p.kpt4_pk3h*y(v1+1900) - p.kpk3h_pt4*y(v1+2100)*y(v1) + p.kh1_h2*y(v1+500)*y(v1+100) - p.kh2_h1*y(v1+2000)*y(v1) -  p.kph1_p11*y(v1+500)*y(v1) + p.kp11_ph1*y(v1+600)+ (p.D_dk/y(4901)^2)*(y(v1-1)-2*y(v1)+y(v1+1))  - p.kbdl_dldk*y(v1+2800)*y(v1) + p.kudl_dldk*y(v1+2900)+ p.kdeg_plc*(y(v1+500)+2*y(v1+600)+2*y(v1+700)+y(v1+800)+y(v1+900)+y(v1+1000)+2*y(v1+1100)+y(v1+1200)+y(v1+1500)+y(v1+1900));%+ p.kdegpp2*y(v1+2600)+ p.kdeg_dl*y(v1+2900) ;%+ p.kdeg_dl3*y(v1+3100);
end 
dydt(100) = -p.kdeg_dkp*y(100)  + p.kph1_pc*y(600) - p.kpc_ph1*y(100)*y(500)  + p.kpk0_pk1*y(800) -  p.kpk1_pk0*y(900)*y(100) + p.kpk1_pk2p*y(900) - p.kpk2p_pk1*y(1400)*y(100) + p.kpt2_pk1h*y(1200) - p.kpk1h_pt2*y(1300)*y(100) + p.kjkp_dj*y(2700) - p.kdj_jkp*y(100)*y(300) + p.kpk2_pk3*y(1100) - p.kpk3_pk2*y(1500)*y(100) + p.kpt3_pk1p*y(1600) - p.kpk1p_pt3*y(1700)*y(100) - p.kpk1p_1h*y(1700)*y(100) + p.k1h_pk1p*y(1300) - p.kph2_p12*y(2100)*y(100) + p.kp12_ph2*y(1000)  + p.kpt4_pk3h*y(2000) - p.kpk3h_pt4*y(2200)*y(100) + p.kh1_h2*y(600)*y(200) - p.kh2_h1*y(2100)*y(100) -  p.kph1_p11*y(600)*y(100) + p.kp11_ph1*y(700) +(p.D_dk/y(4901)^2)*(y(99)-y(100))  - p.kbdl_dldk*y(2900)*y(100) + p.kudl_dldk*y(3000) + p.kdeg_plc*(y(600)+2*y(700)+2*y(800)+y(900)+y(1000)+y(1100)+2*y(1200)+y(1300)+y(1600)+y(2000));%+ p.kdegpp2*y(2700) + p.kdeg_dl*y(3000);
 
%% DIVk
dydt(101) = p.ksyn_dk - p.kdeg_dk*y(101)  - p.kpk1_pk2*y(801)*y(101) + p.kpk2_pk1*y(1001)  -p.kj_jk*y(201)*y(101) + p.kjk_j*y(301) + p.kpk3_pk2p*y(1401) - p.kpk2p_pk3*y(1301)*y(101) - p.kpk3_pk4*y(1401)*y(101) + p.kpk4_pk3*y(1801) - p.kpc_ph2*y(401)*y(101) + p.kph2_pc*y(2001)  - p.kpk1p_p3h*y(1601)*y(101) + p.kp3h_pk1p*y(2101) - p.kph1_p12*y(501)*y(101) + p.kp12_ph1*y(901) + p.kh2_h1*y(2001)*y(1) - p.kh1_h2*y(501)*y(101) - p.kph2_p22*y(2001)*y(101) + p.kp22_ph2*y(1701)+ (p.D_dkp/y(4901)^2)*(y(102)-y(101));
for v2=102:199
    dydt(v2) = p.ksyn_dk - p.kdeg_dk*y(v2)  - p.kpk1_pk2*y(v2+700)*y(v2) + p.kpk2_pk1*y(v2+900)  -p.kj_jk*y(v2+100)*y(v2) + p.kjk_j*y(v2+200) + p.kpk3_pk2p*y(v2+1300) - p.kpk2p_pk3*y(v2+1200)*y(v2) - p.kpk3_pk4*y(v2+1300)*y(v2) + p.kpk4_pk3*y(v2+1700) - p.kpc_ph2*y(v2+300)*y(v2) + p.kph2_pc*y(v2+1900)  - p.kpk1p_p3h*y(v2+1500)*y(v2) + p.kp3h_pk1p*y(v2+2000) - p.kph1_p12*y(v2+400)*y(v2) + p.kp12_ph1*y(v2+800) + p.kh2_h1*y(v2+1900)*y(v2-100) - p.kh1_h2*y(v2+400)*y(v2) - p.kph2_p22*y(v2+1900)*y(v2) + p.kp22_ph2*y(v2+1600)+ (p.D_dkp/y(4901)^2)*(y(v2-1)-2*y(v2)+y(v2+1));
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

for v23=2202:2299
dydt(v23) = p.ksyn_divj - p.kdeg_divj*y(v23) + p.D_divj*(y(v23+1)-2*y(v23) + y(v23-1))/(y(4901)^2) - p.growth*y(v23) - p.kdj_djp*y(v23)*y(v23+100) + p.kdjp_dj*y(v23-2000);
end

dydt(2300) = p.ksyn_divj - p.kdeg_divj*y(2300) + p.D_divj*(y(2299) - y(2300))/(y(4901)^2) - p.growth*y(2300)- p.kdj_djp*y(2300)*y(2400) + p.kdjp_dj*y(300);

% DivJ sticky
dydt(2301:2400) = 0;

% PleC
dydt(2401) = p.ksyn_plc2 - p.kdeg_plc2*y(2401) - p.kplc_plcb*y(2401)*y(2501) + p.kplcb_plc*y(401) + p.D_plc*(y(2402)-y(2401))/(y(4901)^2);% + p.kprot2*y(6201)*y(5101);
for v25=2402:2499
    dydt(v25) = p.ksyn_plc2 - p.kdeg_plc2*y(v25) - p.kplc_plcb*y(v25)*y(v25+100) + p.kplcb_plc*y(v25-2000) + p.D_plc*(y(v25-1)-2*y(v25)+y(v25+1))/(y(4901)^2);% +  p.kprot2*y(v25-100)*y(v25-1200);
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
for v28=2702:2799
    dydt(v28) = p.ksyn_dl - p.kdeg_dl*y(v28)- p.kdlf_dlb*y(v28+300)*y(v28) + p.kdlb_dlf*y(v28+100)+ p.D_dl*(y(v28-1)-2*y(v28)+y(v28+1))/(y(4901)^2);%-kcp_ck*y(v31+400)*(y(v31)^2/(y(v31)^2+Kmdl^2)) + kck_cp*y(v31+500);
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
dydt(3101) =  p.ksyn_ctr- p.kdeg_ctr*y(3101) + (p.D_ctr/y(4901)^2)*(y(3102)-y(3101)) + p.kctr_phos*y(3201)*y(3401) - p.kctr_kin*y(3101)*y(3501);
for v32=3102:3199
    dydt(v32) =  p.ksyn_ctr - p.kdeg_ctr*y(v32) + (p.D_ctr/y(4901)^2)*(y(v32-1)-2*y(v32)+y(v32+1))+ p.kctr_phos*y(v32+100)*y(v32+300) - p.kctr_kin*y(v32)*y(v32+400);
end
dydt(3200) =  p.ksyn_ctr - p.kdeg_ctr*y(3200)+ (p.D_ctr/y(4901)^2)*(y(3199)-y(3200))+ p.kctr_phos*y(3300)*y(3500) - p.kctr_kin*y(3200)*y(3600);
%% CTRAP
dydt(3201) = -p.kdeg_ctr*y(3201)+ (p.D_ctrp/y(4901)^2)*(y(3202)-y(3201)) - p.kctr_phos*y(3201)*y(3401) + p.kctr_kin*y(3101)*y(3501);
for v33=3202:3299
    dydt(v33) = -p.kdeg_ctr*y(v33) + (p.D_ctrp/y(4901)^2)*(y(v33-1)-2*y(v33)+y(v33+1))-  p.kctr_phos*y(v33)*y(v33+200) + p.kctr_kin*y(v33-100)*y(v33+300);
end
dydt(3300) = -p.kdeg_ctr*y(3300) + (p.D_ctrp/y(4901)^2)*(y(3299)-y(3300)) -  p.kctr_phos*y(3300)*y(3500) + p.kctr_kin*y(3200)*y(3600);

% CckA free
dydt(3301) = p.ksyn_ccka - p.kdeg_ccka*y(3301) + p.D_cck*(y(3302) - y(3301))/(y(4901)^2) - p.kcf_cb*y(3301)*y(3901) + p.kcb_cf*y(3401);
for v34 = 3302:3399
    dydt(v34) = p.ksyn_ccka - p.kdeg_ccka*y(v34) + p.D_cck*(y(v34-1) - 2*y(v34) + y(v34+1))/(y(4901)^2)- p.kcf_cb*y(v34)*y(v34+600) + p.kcb_cf*y(v34+100);
end
dydt(3400) = p.ksyn_ccka - p.kdeg_ccka*y(3400) + p.D_cck*(y(3399) - y(3400))/(y(4901)^2)- p.kcf_cb*y(3400)*y(4000) + p.kcb_cf*y(3500);

% CCKAphos
for v35=3401:3500
    dydt(v35) = p.kcf_cb*y(v35-100)*y(v35+500) - p.kcb_cf*y(v35)- p.kdeg_ccka*y(v35) -p.kcp_ck*y(v35)*(y(v35-600)^p.h/(y(v35-600)^p.h+p.Kmdl^p.h)) + p.kck_cp*y(v35+100); 
end

%% Cckakin
dydt(3501) = p.kcp_ck*y(3401)*(y(2801)^p.h/(y(2801)^p.h+p.Kmdl^p.h)) - p.kck_cp*y(3501)- p.kdeg_ccka*y(3501); 
    for v36=3502:3599
        dydt(v36) = p.kcp_ck*y(v36-100)*(y(v36-700)^p.h/(y(v36-700)^p.h+p.Kmdl^p.h)) - p.kck_cp*y(v36)- p.kdeg_ccka*y(v36) ;
    end
dydt(3600) = p.kcp_ck*y(3500)*(y(2900)^p.h/(y(2900)^p.h+p.Kmdl^p.h)) - p.kck_cp*y(3600)- p.kdeg_ccka*y(3600); 

dydt(3601:3900) = 0;
%% CckA stcky
for v40 = 3901:4000
dydt(v40) = 0;
end

dydt(4001:4900) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dydt(4901)=p.growth*y(4901);

