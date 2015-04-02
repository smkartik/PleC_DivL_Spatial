# PleC_DivL_Spatial
Matlab code for simulating results shown in the manuscript 'Dynamical Localization of DivL and PleC in the Asymmetric Division Cycle  of Caulobacter crescentus: A Theoretical Investigation of Alternative Models'

The master file to run the code is get_init.m. To simulate wild type or a particular mutant, enter the appropriate string from the list below
>>['wild_type', 'del_plec', 'plec_f778l', 'divk_ovex', 'divl_misloc', 'divl_deloc', 'divk_d90g', 'plec_h610a']

For example, to simulate PleC deletion, enter the following line of code in Matlab

>>initial_cond = get_init('del_plec')

Or to simulate DivK overexpression,

>>initial_cond = get_init('divk_ovex')

The space-time simulation will be plotted at then end. The complete simulation takes approximately 20 mins (2.8 Ghz intel core I5)

##############################################################################################################################

THE DETAILS OF EACH .m FILE ARE EXPLAINED BELOW

PARAM.m contains the full list of parameters (rate and diffusion constants). It takes 2 arguments. the first argument specifies whether we are simulating a cell of fixed size (0) or a growing cell (1). The second argument specifies which cell type is to be simulated.

ODES.m contains the complete list of differential equations.

GET_INIT.m is the main file. As explained above, it takes a single argument in the form of a string, i.e, the type of cell to be simulated. It calls the odesolver (ode15s) and simulates a non-growing swarmer cell. The endpoints of the simulation serve as initial conditions to simulate a growing cell. get_initcalls the file main_events.m (see below)

MAIN_EVENTS.m simulates a Cauloacter cell growing from swarmer up to predivisional stage prior to compartmentalization (early PD).  

EVENTS.m contains information on when an event takes place. In our model, an event is defined as an enforced change in the location of proteins CckA, DivL, DivJ or PleC.  

ODES_DIV.m is similar to odes.m except that proteins are not allowed to diffuse across mid-cell. 

DIV_MAIN.m It is used to simulate the Caulobacter cell after compartmentalization (early_to_late pd). It calls the ode_solver using odes_div.m as argument. It uses the endpoints of main_events.m simulation as initial conditions.

