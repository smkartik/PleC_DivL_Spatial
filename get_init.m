function initial_cond = get_init(type)

% values assigned to types
wild_type = 1;
del_plec = 2;
plec_f778l = 3;
divk_ovex = 4;
divl_misloc = 5;
divl_deloc = 6;
divk_d90g = 7;
divk_d53n = 8;
del_divj = 9;
divj_h338a = 10;
plec_h610a = 11;

%% declare globals and numerical placeholders for mutants
global p

%% choose mutant
param(0, eval(type));

%% INITIAL CONDITIONS
y0=zeros(4901,1);

y0(2591:2600) = 1;          % PleC sticky
y0(3001:3100) = 1;          % DivL sticky %% cgange for DivL misloc
y0(3901:4000) = 1;          % CckA sticky

y0(4901) = 0.013;           % cell size

%% INTEGRATION PARAMETERS
t0 = 0;
tf = 150;%150;
         
[t,y]=ode15s(@odes,[t0 tf],y0);

%% Define grid M
for n=51:100
    M(:,n)=y(:,4901)*(n/100)-0.5*(y(:,4901)*(.5) + y(:,4901)*.51);  % Each element from column 51 to 100 is assigned the 'n'th fraction of total cell length at a given time step followed by subtracting the mid point value so that the centre of the grid takes the value 0
end

M(:,1:50)=fliplr(M(:,51:100));  % The 1 st 50 colums is the flipped version of the next 50 example: 3 2 1 1 2 3
M(:,1:50)=-M(:,1:50);           %  now its -3 -2 -1 1 2 3
M=100*M;                        % exoressed in micrometres

initial_cond.init = y(end,:);
initial_cond.type = eval(type);
initial_cond.strain = type;

% Run swarmer to early PD spatiotemporal simulation
sw_to_early_pd = main_events(initial_cond)
end
