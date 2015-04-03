function param(cond, type)

global p phen_type 

if cond == 0
    p.growth = 0;
    p.ksyn_divj = 0;
else p.growth = 0.0055;
    p.ksyn_divj = 0.05;     % 0 for delta_divJ
end

%% Synthesis rate constants [units --> 1/min]
p.ksyn_dk = 0.05;       % 0.1 for overexpression	
p.kdeg_dk = 0.005;	
p.kdeg_dkp = 0.005;

p.kdeg_plc = 0.05;
p.ksyn_plc2 = 0.1;     % 0 for delta_pleC
p.kdeg_plc2 = 0.05;

p.kdeg_divj = 0.05;
p.kdegpp2 = 0.05;

p.ksyn_ccka = 0.025;
p.kdeg_ccka= 0.05;

p.ksyn_ctr = 0.05;
p.kdeg_ctr = 0.05;

p.ksyn_dl = 0.025;
p.kdeg_dl = 0.05;

%% Binding constants

p.kdj_djp = 1;          % DivJ binding
p.kdjp_dj = 0;

p.kplc_plcb = 1;        % pleC binding
p.kplcb_plc = 0.5;

p.kdlf_dlb = 1;         % DivL binding
p.kdlb_dlf = 0.1;

p.kbdl_dldk = 2.5;      % DivL:DivK-P binding
p.kudl_dldk = 0.5;

p.kcf_cb = 1;           % CckA binding
p.kcb_cf = 0.1;
 
p.kcp_ck=10;            % CckA phos-to-kin
p.kck_cp=1;
p.Kmdl=0.5;
p.h = 4;

%% DivJ phosphorylation paramters
p.kj_jk = 5; 
p.kjk_j = 0.0016;	
p.kjk_jkp = 5; 
p.kjkp_jk = 0.16;
p.kjkp_dj = 1;
p.kdj_jkp = 5;

%% PleC phosphorylation paramters
p.kpc_ph1 = 5;	
p.kph1_pc = 5;	

p.kph1_p11 = 5;	
p.kp11_ph1 = 2.5;

p.kp11_pk0 = 2.5;	
p.kpk0_p11 = 5;

p.kpk0_pk1 = 0.16;	
p.kpk1_pk0 = 5;

p.kpk1_pk2 = 5;	
p.kpk2_pk1 = 0.0016;
 
p.kpk1_pk1h = 5;	
p.kpk1h_pk1 = 5;

p.kpk1_pk2p = 0.16;	
p.kpk2p_pk1 = 5;

p.kpk2_pt2 = 5;	
p.kpt2_pk2 = 0.16;

p.kpt2_pk1h = 0.16;	
p.kpk1h_pt2 = 5;

p.kpk2p_pc = 5;	
p.kpk1p_pc = 5;

p.kpk2_pk3 = 0.16;	
p.kpk3_pk2 = 5;

p.kpk3_pt3 = 5;	
p.kpt3_pk3 = 0.16;

p.kpt3_pk1p = 0.16;	
p.kpk1p_pt3 = 5;

p.kpk1p_1h = 5;	
p.k1h_pk1p = 0.16;

p.kpk3_pk2p = 0.0016;	
p.kpk2p_pk3 = 5;

p.kpk3_pk4 = 5;	
p.kpk4_pk3 = 0.0016;

p.kpk4_pt4 = 5;	
p.kpt4_pk4 = 0.16;

p.kpc_ph2 = 0.05;	
p.kph2_pc = 5;

p.kph2_p22 = 0.016;	
p.kp22_ph2 = 1.6e-08;

p.kp22_pk4 = 5;	
p.kpk4_p22 = 5;

p.kpt4_pk3h = 0.16;	
p.kpk3h_pt4 = 5;

p.kpk3_pk3h = 5;	
p.kpk3h_pk3 = 5;

p.kpk1p_p3h = 5;	
p.kp3h_pk1p = 0.0016;

p.kph1_p12 = 1.6e-02;	
p.kp12_ph1 = 1.6e-04;

p.kp12_pk2 = 5;	
p.kpk2_p12 = 5;

p.kph2_p12 = 1.6;	
p.kp12_ph2 = 1.6e-04;

p.kh1_h2 = 1.6e-02;	
p.kh2_h1 = 1.6;

p.kp11_pt4 = 0.0755;	
p.kpt4_p11 = 5;
 
p.kph1_ph2 = 10;	
p.kph2_ph1 = 5e-03;

%% CtrA phosphorylation
p.kctr_kin = 600;
p.kctr_phos = 600;

%% Diffusion parameters [units --> um^2/min]

p.D_dk = 100;
p.D_dkp = 100;
p.D_divj = 100;
p.D_dl = 100;
p.D_plc = 10;
p.D_ctr = 100;
p.D_ctrp = 100;
p.D_cck = 100;
 
%% Mutants

if type == 1
    phen_type = 1;
    
elseif type == 2            %% del_plec
    p.ksyn_plc2 = 0;
    
    phen_type = 2;

elseif type == 3;       %% plec_F778L
    p.kp11_pk0 = 0;	
    p.kpk0_p11 = 0;
    
    p.kp22_pk4 = 0;	
    p.kpk4_p22 = 0;

    p.kp12_pk2 = 0;         	
    p.kpk2_p12 = 0;
    
    p.kpk2_pt2 = 0;         	
    p.kpt2_pk2 = 0;

    p.kpk3_pt3 = 0;        		
    p.kpt3_pk3 = 0;

    p.kpk4_pt4 = 0;         		
    p.kpt4_pk4 = 0;        

    p.kp11_pt4 = 0;         	
    p.kpt4_p11 = 0;         
    
    phen_type = 3;

elseif type == 4;       %% divk_ovex
    p.ksyn_dk = 0.4;
    phen_type = 4;

elseif type == 5;       %% divl_misloc
    phen_type = 5;

elseif type == 6;       %% divl_desloc
    phen_type = 6;
    
elseif type == 7;       %% divk_D90G
   % p.ksyn_plc2 = 0;
    %p.kbdl_dldk = 0;%0.1;          
    %p.kudl_dldk = 0.5;

    p.kph1_p11 = 2.5;       
    p.kp11_ph1 = 50;        
    
    p.kph2_p22 = 0.025;         	
    p.kp22_ph2 = 2.5e-03;   	
    
    p.kph1_p12 = 1.6e-02;       
    p.kp12_ph1 = 16;            

    p.kph2_p12 = 1.6;           	
    p.kp12_ph2 = 16;
    
    phen_type = 7;

elseif type == 8;       %% divk_D53N
    p.kjk_jkp = 0;                  
    p.kjkp_jk = 0;                  

    p.kpk2_pt2 = 0;         	
    p.kpt2_pk2 = 0;

    p.kpk3_pt3 = 0;        		
    p.kpt3_pk3 = 0;

    p.kpk4_pt4 = 0;         		
    p.kpt4_pk4 = 0;        

    p.kp11_pt4 = 0;         	
    p.kpt4_p11 = 0;

    phen_type = 8;

elseif type == 9;       %% del_divJ
    phen_type = 9;
    p.ksyn_divj = 0;

elseif type == 10;      %% divj_h338A
    phen_type = 10;
    p.kjk_jkp = 0; 
    p.kjkp_jk = 0;

elseif type == 11;       %% plec_H610A
    p.kp11_pk0 = 0;	
    p.kpk0_p11 = 0;
    
    p.kp22_pk4 = 0;	
    p.kpk4_p22 = 0;

    p.kp12_pk2 = 0;         	
    p.kpk2_p12 = 0;
    
    p.kpk2_pt2 = 0;         	
    p.kpt2_pk2 = 0;

    p.kpk3_pt3 = 0;        		
    p.kpt3_pk3 = 0;

    p.kpk4_pt4 = 0;         		
    p.kpt4_pk4 = 0;        

    p.kp11_pt4 = 0;         	
    p.kpt4_p11 = 0;
    
    p.kph1_ph2 = 0;        % hydrolysis 
    p.kph2_ph1 = 0;
    
    phen_type = 11;    

end
