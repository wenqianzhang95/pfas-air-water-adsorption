%% Full-transport forward solve used by inverse PFOS/PFOA SSL calculations
% Flow:      Richards equation is solved by Celia et al (1990) mixed form 
%            method with Newton iteration
% Transport: Advection-diffusion-adsorption equation is solved by fully
%            implicit in time; upwind for advection, and central
%            difference for diffusion in space
%            PFAS adsorption on solid, air-water interface are kinetic based on a two-domain model     
% Coupling:  Flow and transport are fully coupled implicitly
% Unit system: cgs with micro mol for the concentration of species

% Transport backbone: Schaefer_Fwd_Model_SMP1_SMP3_LONG_TERM.m
% Kralchevsky equations and parameters: make_fig3_surface_tension_2x3_matlab.m
%
% SSL design:
%   * Every trial is a complete transient Richards/transport calculation.
%     Only intermediate file output can be suppressed; no transport term is
%     frozen, averaged, or replaced by the accelerated SSL approximation.
%   * The water flow, solid sorption, Aaw model, forcing, and initial mass
%     are identical between the two interfacial-isotherm cases.
%   * isotherm_flg = 1 uses the 10-mM Langmuir fits paired with the user's
%     Kralchevsky comparison; the Ryan legacy parameters are not selectable.
%   * isotherm_flg = 2 uses the Kralchevsky/Stern model.
%   * A prescribed NaCl concentration is required by the Kralchevsky model.
%     The default is 10 mM, a common condition for both PFOS and PFOA data.
clc; clear;
scriptDir = fileparts(mfilename('fullpath'));
startingDir = pwd;
cleanupDir = onCleanup(@() cd(startingDir));
cd(scriptDir);
setpath;

% Turn on gravity
gravity reset on,

% Shared comparison settings. Environment variables are used so that the
% same immutable model file can run each case on Slurm.
Aawflg = 1;                 % fixed for both isotherms: thermo Aaw with scaling factor
NaCl_run_mM = 10;           % prescribed NaCl concentration (mM = mol/m^3 numerically)
simulation_cycles = 5;      % repeat the nominal 10-year daily forcing five times = 50 years
interface_temperature_K = 298.0; % temperature used by the Fig. 3 fits
output_interval_days = 30;  % save full profiles and diagnostic time series every 30 days
afff_active_years = 30;      % explicit source-window setting; current input source column is all zero
progress_interval_days = 30; % concise live/Slurm heartbeat interval
verbose_timesteps = false;   % per-step output is available for debugging but is very large
ryan_validation = false; % SSL package intentionally excludes legacy regression mode
langmuir_parameter_set = "matched_10mM";

initial_scale = str2double(getenv('PFAS_INITIAL_SCALE'));
if ~isfinite(initial_scale) || initial_scale <= 0
    initial_scale = 1;
end
ssl_response_only = strcmpi(getenv('PFAS_SSL_RESPONSE_ONLY'), 'true');
receptor_mcl_ug_L = str2double(getenv('PFAS_RECEPTOR_MCL_UG_L'));
if ~isfinite(receptor_mcl_ug_L) || receptor_mcl_ug_L <= 0
    receptor_mcl_ug_L = 0.004;
end

cycles_env = str2double(getenv('PFAS_SIMULATION_CYCLES'));
if isfinite(cycles_env) && cycles_env >= 1 && cycles_env == floor(cycles_env)
    simulation_cycles = cycles_env;
end
output_days_env = str2double(getenv('PFAS_OUTPUT_INTERVAL_DAYS'));
if isfinite(output_days_env) && output_days_env >= 1 && output_days_env == floor(output_days_env)
    output_interval_days = output_days_env;
end
afff_years_env = str2double(getenv('PFAS_AFFF_ACTIVE_YEARS'));
if isfinite(afff_years_env) && afff_years_env >= 0
    afff_active_years = afff_years_env;
end
progress_days_env = str2double(getenv('PFAS_PROGRESS_INTERVAL_DAYS'));
if isfinite(progress_days_env) && progress_days_env > 0
    progress_interval_days = progress_days_env;
end
verbose_timesteps = strcmpi(getenv('PFAS_VERBOSE_TIMESTEPS'), 'true');
case_id = str2double(getenv('PFAS_CASE_ID'));
if ~ismember(case_id, 1:4)
    error(['PFAS_CASE_ID must be one of 1, 2, 3, or 4: ', ...
        '1=PFOS-Langmuir, 2=PFOS-Kralchevsky, ', ...
        '3=PFOA-Langmuir, 4=PFOA-Kralchevsky.']);
end
if ryan_validation && ~ismember(case_id, [1 3])
    error('Ryan validation is defined only for Langmuir cases 1 (PFOS) and 3 (PFOA).');
end
case_pfas = [1 1 2 2];
case_isotherm = [1 2 1 2];
PFAS_run = case_pfas(case_id);
SMP_run = 1;               % four requested cases all use SMP1
isotherm_run = case_isotherm(case_id);

transport_smoke_test = strcmpi(getenv('PFAS_TRANSPORT_SMOKE_TEST'), 'true');
short_test_days = str2double(getenv('PFAS_SHORT_TEST_DAYS'));
short_term_test = isfinite(short_test_days) && short_test_days > 0;
if transport_smoke_test
    simulation_cycles = 1;
    defaultOutputRoot = fullfile(scriptDir, 'Results_smoke_test');
elseif short_term_test
    short_test_days = floor(short_test_days);
    simulation_cycles = 1;
    defaultOutputRoot = fullfile(scriptDir, sprintf( ...
        'Results_short_%dday_%ddayOutput', ...
        short_test_days, output_interval_days));
else
    defaultOutputRoot = fullfile(scriptDir, 'Results_50yr');
end
if ryan_validation
    defaultOutputRoot = fullfile(scriptDir, 'Results_Ryan_validation');
end
outputRoot = getenv('PFAS_OUTPUT_ROOT');
if isempty(outputRoot)
    outputRoot = defaultOutputRoot;
end
if ~exist(outputRoot, 'dir')
    mkdir(outputRoot);
end

fprintf(['\nHPC case %d: PFAS flag %d, isotherm flag %d, SMP1, ', ...
    '%d forcing cycles, %d-day output, %.3g-year AFFF window\n'], ...
    case_id, PFAS_run, isotherm_run, simulation_cycles, ...
    output_interval_days, afff_active_years);
fprintf('Live progress every %.3g simulated days; output root: %s\n', ...
    progress_interval_days, outputRoot);

for PFAS_loop = PFAS_run
    for SMP_loop = SMP_run
        for isotherm_loop = isotherm_run
scenario_wall_timer = tic;
%% Specify Modeling Scenario
SMPflg = SMP_loop;     % 1 -> SMP-1 |  3 -> SMP-3
interpflg = 2;  % 1 -> linear interpolation |  2 -> exponential interpolation | 3 -> "step" interpolation
PFASflg = PFAS_loop;    % 1 -> PFOS | 2 -> PFOA | 3 -> PFHxS | 4 -> PFBS
paramflg = 2;   % 1 -> parameters calibrated using flushing data | 2 -> calibrated using ambient data
cs2flg = 3;       % 1 -> interpolated conc. @ points | 2 -> averaged interp. conc. for 3 regions | 3 -> avg interp. conc. for test cell
labile_flg = 0;         % 1 -> account for non-labile PFAS

if isotherm_loop == 1
    isotherm_flg = 1;
    isotherm_name = "LangmuirSzyszkowski";
elseif isotherm_loop == 2
    isotherm_flg = 2;
    isotherm_name = "Kralchevsky";
else
    error('Unsupported interfacial isotherm flag.');
end

% Optional conservative re-parameterization for the nonlinear interface.
% In this mode the Newton primary variable is the bounded logit of surface
% excess (Gamma), while aqueous concentration and surface pressure are
% obtained from the same Kralchevsky lookup table by monotone inversion.
gamma_primary_mode = strcmpi(getenv('PFAS_INTERFACE_PRIMARY'), 'gamma');
if gamma_primary_mode && isotherm_flg ~= 2
    error('PFAS_INTERFACE_PRIMARY=gamma is defined only for Kralchevsky cases.');
end
conservative_inventory_mode = strcmpi(getenv('PFAS_CONSERVATIVE_INVENTORY'), 'true');
closed_box_mode = strcmpi(getenv('PFAS_CLOSED_BOX_TEST'), 'true');
conservative_step_mass_tolerance = str2double(getenv('PFAS_CONSERVATIVE_STEP_MASS_TOL'));
if ~isfinite(conservative_step_mass_tolerance) || conservative_step_mass_tolerance <= 0
    conservative_step_mass_tolerance = 1e-8;
end
step_mass_tolerance = str2double(getenv('PFAS_STEP_MASS_TOL'));
if ~isfinite(step_mass_tolerance) || step_mass_tolerance <= 0
    step_mass_tolerance = inf; % disabled unless explicitly requested
end
step_mass_rejections = 0;
time_step_cuts = 0;
linear_solve_rejections = 0;
state_step_rejections = 0;
minimum_dt_seconds = str2double(getenv('PFAS_MINIMUM_DT_SECONDS'));
if ~isfinite(minimum_dt_seconds) || minimum_dt_seconds <= 0
    minimum_dt_seconds = 1e-6;
end

% if interpflg == 1
%     interpname = "cs2_lininterp";
% elseif interpflg == 2
%     interpname = "cs2_expinterp";
% elseif interpflg == 3
%     interpname = "cs2_stepinterp";
% else
%     fprintf('Please specify interpflg corresponding to one of the appropriate scenarios\n');
%     return 
% end

rhob = 1.627;                                       % g/cm^3            -> soil bulk density

if PFASflg == 1     % PFOS Parameters
    COC = "PFOS";
    z = 5;                          % specifies row of excel table to pull input data from
    M = 500.13;                                     % g/mol             -> molecular weight of PFOS
    sigma0 = 71;                                    % dyn/cm            -> surface tension
    Kf = 6.9;   % used Kd for PFOS in Schaefer 2022 Desorption Paper (units L/kg)
    non_labile = (1000000/1000)*((rhob/1000)/M); % units are micro mol/cm^3
    % Independent Langmuir fit to the same 10-mM surface-tension dataset
    % used in the user's Kralchevsky comparison.
    langmuir_KL = 123.794092031724;                 % m^3/mol
    langmuir_gamma_m = 2.49648324768916e-6;         % mol/m^2
    langmuir_fit_temperature_K = 298.0;
    a = 1/langmuir_KL;                              % mol/m^3 (= umol/cm^3)
    b = 1000*8.314*langmuir_fit_temperature_K*langmuir_gamma_m/sigma0;
    % Fig. 3 Kralchevsky parameters [K1 (m^3/mol), K2 (m^6/mol^2), Gamma_m (mol/m^2)]
    kralchevsky_pars = [1553, 3.8825, 3.0695e-6];
    use_activity_correction = false;
elseif PFASflg == 2     % PFOA Parameters
    COC = "PFOA";
    z = 10;                          % specifies row of excel table to pull input data from
    M = 414.07;                                     % g/mol             -> molecular weight of PFOA
    sigma0 = 71.85;                                    % dyn/cm            -> surface tension
    Kf = 2.2;   % used Kd for PFOS in Schaefer 2022 Desorption Paper (units L/kg)
    non_labile = (3000/1000)*((rhob/1000)/M);  % units are micro mol/cm^3
    % Independent Langmuir fit paired with the 10-mM PFOA Kralchevsky case.
    langmuir_KL = 6.227660735;                      % m^3/mol
    langmuir_gamma_m = 2.988939298e-6;              % mol/m^2
    langmuir_fit_temperature_K = 298.0;
    a = 1/langmuir_KL;                              % mol/m^3 (= umol/cm^3)
    b = 1000*8.314*langmuir_fit_temperature_K*langmuir_gamma_m/sigma0;
    kralchevsky_pars = [49.2, 0.13, 3.36e-6];
    use_activity_correction = false;
else
    error('This comparison version supports only PFOS and PFOA.');
end
% Hard preflight for the sole approved physical comparison. This package
% cannot silently run an activity-corrected or legacy Fig. 3 parameter set.
if use_activity_correction
    error('Preflight failure: use_activity_correction must be false.');
end
if PFASflg == 1
    expected_kral = [1553, 3.8825, 3.0695e-6];
else
    expected_kral = [49.2, 0.13, 3.36e-6];
end
if any(abs(kralchevsky_pars-expected_kral) > 1e-14*max(1,abs(expected_kral)))
    error('Preflight failure: mass-action parameters differ from SI Table S1.');
end
if NaCl_run_mM ~= 10 || abs(interface_temperature_K-298) > 1e-12
    error('Preflight failure: require NaCl=10 mM and T=298 K.');
end
if ~conservative_inventory_mode
    error('Preflight failure: require single_total_inventory_BE formulation.');
end
%% Read datasets
if ryan_validation
    atmbc = readtable(fullfile(scriptDir, '..', 'Rainfall_Inputs_Schaefer.xlsx'), ...
        'Sheet', 'Pert Retention Analysis (30min)', ...
        'VariableNamingRule','preserve');
    forcing_interval_seconds = 30*60;
else
    atmbc = readtable('Rainfall_Inputs_10year.xlsx', 'Sheet', 'Model Input', 'VariableNamingRule','preserve');
    forcing_interval_seconds = day;
end
receptor_forcing = atmbc; % retain the complete forcing record for mean DF
if transport_smoke_test && ~ryan_validation
    atmbc = atmbc(1:min(2,height(atmbc)), :);
elseif short_term_test && ~ryan_validation
    if short_test_days > height(atmbc)
        % Diagnostic runs may need to cross a failure in a later repeated
        % forcing cycle. Tile the same forcing record exactly as production
        % does, then trim to the requested diagnostic horizon.
        forcingRepeats = ceil(short_test_days/height(atmbc));
        atmbc = repmat(atmbc,forcingRepeats,1);
    end
    atmbc = atmbc(1:short_test_days, :);
    atmbc.("Time Step") = (1:short_test_days)';
end
% 1st column of atmbc: index of data
% 2nd column of atmbc: precipitation (input as cm/day and it will get
%                      converted to cm3/sec)
% 3rd column of atmbc: evapotranspiration
% 4th column of atmbc: infiltration due to AFFF release

%% Set up model geometry
[nx,ny,nz] = deal(1, 1, 150);               % # of cells
[Lx,Ly,Lz] = deal(1, 1, 150);               % cm
G = cartGrid([nx, ny, nz], [Lx, Ly, Lz]);   % generate grid
G = computeGeometry(G);                     % compute the data structure for the grid
dz = Lz/nz;                                 % grid size along z 

%% Parameters that Change with SMP Location
% if SMPflg == 1 && paramflg == 1
%     Ks = 1478/day*(Lx/nx)*(Ly/ny);                   %Input as cm/day, converts to cm/s -> Saturated hydraulic conductivity
%     poro = 0.3444;                        % [-]       -> porosity
%     thetar = 0.06263;                     % [-]       -> residual water content
%     thetas = 0.3444;                     % [-]       -> saturated water content
%     n1 = 2.774;                         % [-]       -> VG model
%     m1 = 1 - 1./n1;                     % [-]       -> VG model
%     alpha1 = 0.02456;
%     L = 1.488;
%     SMP = "SMP1";
%     hSMP = "SMP1";      % used for importing initial h and soil concentration conditions later in code
%     param = "flushing_params";
% elseif SMPflg == 3 && paramflg == 1
%     Ks = 411.3/day*(Lx/nx)*(Ly/ny);                   %Input as cm/day, converts to cm/s -> Saturated hydraulic conductivity
%     poro = 0.4217;                        % [-]       -> porosity
%     thetar = 0.03003;                     % [-]       -> residual water content
%     thetas = 0.4217;                     % [-]       -> saturated water content
%     n1 = 2.549;                         % [-]       -> VG model
%     m1 = 1 - 1./n1;                     % [-]       -> VG model
%     alpha1 = 0.01;
%     L = 0.4465;
%     SMP = "SMP3";
%     hSMP = "SMP3";      % used for importing initial h and soil concentration conditions later in code
%     param = "flushing_params";
if SMPflg == 1 && paramflg == 2
    Ks = 1478/day*(Lx/nx)*(Ly/ny);                   %Input as cm/day, converts to cm/s -> Saturated hydraulic conductivity
    poro = 0.3444;                        % [-]       -> porosity
    thetar = 0.06263;                     % [-]       -> residual water content
    thetas = 0.3444;                     % [-]       -> saturated water content
    n1 = 2.774;                         % [-]       -> VG model
    m1 = 1 - 1./n1;                     % [-]       -> VG model
    alpha1 = 0.02456;
    L = 1.488;
    SMP = "SMP1";
    hSMP = "SMP1";      % used for importing initial h and soil concentration conditions later in code
    param = "ambient_params";
elseif SMPflg == 3 && paramflg == 2
    Ks = 919.7/day*(Lx/nx)*(Ly/ny);                   %Input as cm/day, converts to cm/s -> Saturated hydraulic conductivity
    poro = 0.3792;                        % [-]       -> porosity
    thetar = 0.04889;                     % [-]       -> residual water content
    thetas = 0.3792;                     % [-]       -> saturated water content
    n1 = 2.773;                         % [-]       -> VG model
    m1 = 1 - 1./n1;                     % [-]       -> VG model
    alpha1 = 0.02497;
    L = 1.164;
    SMP = "SMP3";
    hSMP = "SMP3";      % used for importing initial h and soil concentration conditions later in code
    param = "ambient_params";
elseif SMPflg == 4
    Ks = 869.4/day*(Lx/nx)*(Ly/ny);                   %Input as cm/day, converts to cm/s -> Saturated hydraulic conductivity
    poro = 0.3658;                        % [-]       -> porosity
    thetar = 0.0482;                     % [-]       -> residual water content
    thetas = 0.3658;                     % [-]       -> saturated water content
    n1 = 3.2788;                         % [-]       -> VG model
    m1 = 1 - 1./n1;                     % [-]       -> VG model
    alpha1 = 0.0342;
    L = 0.5;
    SMP = "Rosetta3";
    hSMP = "SMP1";      % used for importing initial h and soil concentration conditions later in code
    param = "ambient_params";
elseif SMPflg == 5
    Ks = 511.545/day*(Lx/nx)*(Ly/ny);                   %Input as cm/day, converts to cm/s -> Saturated hydraulic conductivity
    poro = 0.34443;                        % [-]       -> porosity
    thetar = 0.0336817;                     % [-]       -> residual water content
    thetas = 0.34443;                     % [-]       -> saturated water content
    n1 = 2.7733;                         % [-]       -> VG model
    m1 = 1 - 1./n1;                     % [-]       -> VG model
    alpha1 = 0.02363;
    L = 0.5;
    SMP = "Rosetta3_lowSTDEV";
    hSMP = "SMP1";      % used for importing initial h and soil concentration conditions later in code
    param = "ambient_params";
elseif SMPflg == 6
    Ks = 1477.6/day*(Lx/nx)*(Ly/ny);                   %Input as cm/day, converts to cm/s -> Saturated hydraulic conductivity
    poro = 0.387129;                        % [-]       -> porosity
    thetar = 0.0626286;                     % [-]       -> residual water content
    thetas = 0.387129;                     % [-]       -> saturated water content
    n1 = 3.87634;                         % [-]       -> VG model
    m1 = 1 - 1./n1;                     % [-]       -> VG model
    alpha1 = 0.049414;
    L = 0.5;
    SMP = "Rosetta3_upSTDEV";
    hSMP = "SMP1";      % used for importing initial h and soil concentration conditions later in code
    param = "ambient_params";
else
    fprintf('Please specify SMPflg and paramflg corresponding to one of the SMP Locations and parameter sets\n');
    return 
end 

% Ryan's archived 60-day result used an earlier SMP1 ambient hydraulic
% parameter set. Restore it only for the regression comparison.
if ryan_validation
    thetar = 0.06219;
    alpha1 = 0.02645;
    L = 1.491;
end



%% Soil parameters
rock = makeRock(G, Ks, poro);       % Note: perm is replaced by Ks
rock.Ks = rock.perm;                % cm/s      -> copy the replaced perm (by Ks) to Ks

%% Fluid parameters
mu    = 0.01;                       % g/(cm s)  -> viscosity of water
rhow  = 1;                          % g/cm^3    -> density of water
thetacritl = 0.002;                 % [-]       -> critical theta to switch to pressure head bc under dry condition
hcritl = -1./alpha1 * (thetacritl^(-1/m1) - 1).^(1/n1);
%% Transport parameters
% Solute
Dm = 5.4e-6;                                % cm^2/s            -> molecular diffusion coeff. [ref: Schaefer et al. 2019 JEE]
alphaL = 100*0.83*(log10(Lz/100))^2.414;    % cm                -> Longitudinal dispersivity. 
% Solid adsorption 
% NB1: It can be either linear (n=1) or nonlinear (n!=1)
% NB2: unit for Kf is mg/kg / (mg/L)^n3
n3 = 1;         % = 1 for linear adsorption
%sf = 4.412;                                          % scaling factor for thermodynamic-based Aaw
% Air-water adsorption
% if Aawflg == 1
%     m2 = 1305;                                      % [-]               -> Aaw fitting parameter
%     n2 = -2848.6;                                   % [-]               -> Aaw fitting parameter
%     alpha2 = 1543.6;                                % [-]               -> Aaw fitting parameter
% elseif Aawflg == 0
%     m2 = 313.39;                                    % [-]               -> Aaw fitting parameter
%     n2 = -674.38;                                   % [-]               -> Aaw fitting parameter
%     alpha2 = 360.99;                                % [-]               -> Aaw fitting parameter
% end
R = 8.314;                                      % J/K/mol           -> gas constant
Tp = interface_temperature_K;                   % K                 -> shared interface-model temperature
if isotherm_flg == 1 && abs(Tp-langmuir_fit_temperature_K) > 1e-12
    error('Matched 10 mM Langmuir parameters are defined at 298 K.');
end
Fs = 1;                                      % fraction of solid adsorption, 1-all equilibrium, 0-all kinetic.                            
alphas = 0.9/hour;                              % kinetic coefficient for solid adsorption
Faw = 1;                                        % fraction of AW adsorption, 1-all equilibrium, 0-all kinetic. 
alphaaw = 5/hour;                               % kinetic coefficient for AW adsorption
if conservative_inventory_mode && (Fs ~= 1 || Faw ~= 1)
    error(['Conservative inventory mode requires fully equilibrium solid and ', ...
        'air-water adsorption (Fs=Faw=1).']);
end

% Precompute the nonlinear Kralchevsky relation once per run. The lookup
% table avoids solving a Stern charge-balance problem in every cell at every
% transport Newton iteration. Concentration is numerically identical in
% mol/m^3 and micro mol/cm^3.
if isotherm_flg == 2
    interface_lookup = build_kralchevsky_lookup(NaCl_run_mM, ...
        kralchevsky_pars, use_activity_correction, Tp);
else
    interface_lookup = struct();
end
    
%% Extract grid information
nf = G.faces.num;                                 % Number of faces
nc = G.cells.num;                                 % Number of cells
N  = double(G.faces.neighbors);                   % Map: face -> cell
F  = G.cells.faces(:,1);                          % Map: cell -> face
intInx = all(N ~= 0, 2);                          % Interior faces
N  = N(intInx, :);                                % Interior neighbors
AF = G.faces.areas;
FC = G.faces.neighbors;
LF = zeros(length(FC),1);
LF(intInx) = sqrt(sum((G.cells.centroids(FC(intInx,1),:)-G.cells.centroids(FC(intInx,2),:)).^2,2));
AF=AF(intInx);
LF=LF(intInx);
vC  = G.cells.volumes;                            % Volume of grid cells

%% Compute transmissibilities
hT = computeTrans(G, rock);                       % Half-transmissibilities. Note: perm is replaced by Ks
T  = 1 ./ accumarray(F, 1 ./ hT, [nf, 1]);        % Harmonic average. Note: perm is replaced by Ks
T  = T(intInx);                                   % Restricted to interior

%% Define discrete operators
nf   = size(N,1);
M_g  = sparse([(1:nf)'; (1:nf)'], N, ones(nf,1)*[-1 1], nf, nc);
grad  = @(x)  M_g*x;
div   = @(x)-M_g'*x;
avg   = @(x) 0.5 * (x(N(:,1)) + x(N(:,2)));
upw   = @(flag, x) flag.*x(N(:, 1)) + ~flag.*x(N(:, 2));
gradz = grad(G.cells.centroids(:,3));

%% Schedule and infiltration/evaporation rates
halfh = 60*30;

% Output is monthly-scale, while Ti below still forces the nonlinear solver
% to hit every daily forcing change. Thus lowering output frequency does not
% average or otherwise alter the daily rainfall boundary condition.
Tf = simulation_cycles*length(atmbc.("Time Step"))*forcing_interval_seconds;
Ts = (0:output_interval_days*day:Tf)';                  % regular profile-output times
if Ts(end) < Tf
    Ts(end+1,1) = Tf;                                   % always save the final state
end
nstep = numel(Ts)-1;
Ti = (1:1:simulation_cycles*atmbc.("Time Step")(end,1))*forcing_interval_seconds;
inIx1 = 1;                                              % index for the top cell
inIx2 = nc;                                             % index for the bottom cell    
total_in_w = 0;                                         % total water entering the domain    
total_out_w = 0;                                        % total water exiting the domain
total_in_c = 0;                                         % total PFAS entering the domain
total_out_c = 0;                                        % total PFAS exiting the domain

% The transport model has no receptor-well geometry.  These are the only
% inputs retained from Smith et al. (2024), Table 4 AFFF scenario.
receptor.W_cm = 3000;
receptor.Ugw_cm_yr = 100*365.25;
receptor.bsat_cm = 305;
receptor.alpha_v_cm = 0.0056*receptor.W_cm;
% ET is stored as cm/s in the transport workbook; rainfall and the AFFF
% water column are cm/day.  The current workbook has zero ET and zero AFFF.
net_infiltration_cm_yr = mean(0.47*receptor_forcing.("Rainfall (cm/day)") + ...
    receptor_forcing.("PFAS infiltration") - receptor_forcing.ET*day)*365.25;
if ~(isfinite(net_infiltration_cm_yr) && net_infiltration_cm_yr > 0)
    error('Mean net infiltration must be positive for the receptor DF calculation.');
end
groundwater_DF_value = groundwater_DF(net_infiltration_cm_yr, receptor);
% Report both common conventions explicitly. DAF > 1 is the divisor used
% below; AF is the transmitted concentration fraction and equals 1/DAF.
groundwater_AF_value = 1/groundwater_DF_value;
if ~ryan_validation
    atmbc.("Rainfall (cm/day)") = atmbc.("Rainfall (cm/day)").*0.47; % fraction reaching WT (Schaefer, 2023 SI)
end
atmbc.("Rainfall (cm/day)") = atmbc.("Rainfall (cm/day)")/day*(Lx/nx)*(Ly/ny); % convert the unit to cm^3/s -> Infiltration or evaporation rate
% atmbc.("Rainfall (cm/day)") = (atmbc.("Rainfall (cm/day)")/day*(Lx/nx)*(Ly/ny)); % convert the unit to cm^3/s -> Infiltration or evaporation rate
atmbc.ET = atmbc.ET*(Lx/nx)*(Ly/ny);
atmbc.("PFAS infiltration") = atmbc.("PFAS infiltration")/day*(Lx/nx)*(Ly/ny);

%% Initialize for flow and Transport and set BC
ci = 100/M;
% NB: c0 needs to be set as a very small value for nonlinear adsorption. It cannot
%     be zero.

if cs2flg == 1
    initial_soil = hSMP + '_' + COC;
    cs2_cond = 'cs2_at_Location';
elseif cs2flg ==2
    initial_soil = hSMP + '_' + COC + '_AreaAvg';
    cs2_cond = 'cs2_AreaAvg';
elseif cs2flg == 3
    initial_soil = COC + '_Test_Plot_Avg';
    cs2_cond = 'cs2_TestPlotAvg';
else
    fprintf('Please specify cs2flg corresponding to an appropriate set of initial conditions\n');
    return 
end    

soil_data = readtable('Schaefer2022_Soil_Data_1.xlsx','Sheet', initial_soil, 'VariableNamingRule','preserve');     % 1st column = deep soil data, 2nd = shallow soil data
% rows coincide with PFAS => (1)PFBS (2)PFPeS (3)PFHxS (4) PFHpS (5)PFOS
% (6)PFBA (7)PFPeA (8)PFHxA (9)PFHpA (10)PFOA (11)FOSA (12)FHxSA (13)8:2 FtS (14)6:2
% FtS ...
if PFASflg == 1
    initial_profile_ug_kg = soil_data.PFOS;
    ctot = initial_scale*initial_profile_ug_kg*((rhob/1000)/M);
elseif PFASflg == 2
    initial_profile_ug_kg = soil_data.PFOA;
    ctot = initial_scale*initial_profile_ug_kg*((rhob/1000)/M);
elseif PFASflg == 3
    ctot = soil_data.PFHxS*((rhob/1000)/M);    
elseif PFASflg == 4
    ctot = soil_data.PFBS*((rhob/1000)/M); 
else
    fprintf('Please specify PFASflg corresponding to an appropriate PFAS \n');
    return 
end

if labile_flg == 1
    ctot = ctot - non_labile;
    COC = COC + '_nonlabile';
end

if Aawflg == 1
    SMP = SMP + "_thermoAaw";
elseif Aawflg == 2
    SMP = SMP + "_nosfAaw";
elseif Aawflg == 3
    SMP = SMP + "_noAaw";
else
    fprintf('Please specify Aawflg corresponding to an appropriate Aaw fxn scenario \n');
    return   
end

%c0 = 18000.*ones(G.cells.num,1);
c0 = 1E-40*ones(G.cells.num,1);
cs2 = 0*ones(G.cells.num,1);

%Use Linear Interpolation to Intialize Soil Concentration Values
% cs2 = zeros(150,1);
% if interpflg == 1           % linear interpolation
%     cs2_90 = soil_data(z,2)*((rhob/1000)/M)*ones(90,1);   % measured values for 1st 90cm of profile
%     cs2_60 = ((soil_data(z,1)*(rhob/1000))/M)*ones(60,1);   % measured values for remaining 60cm of profile
%     cs2(1) = cs2_90(1);
%     cs2(91:end) = cs2_60;
%     cs2 = interp1([1, 91:150]', [cs2(1); cs2(91:150)], [1:150]','linear');
% elseif interpflg == 2       % exp interplotation
%     intery = [soil_data(z,2); soil_data(z,1)];
%     interx = [1; 91];
%     cs2 = interp1([1, 91:150]', [cs2(1); cs2(91:150)], [1:150]','linear');
%     cs_exp = fit(interx, intery, 'exp1');
%     %plot(cs_exp,interx,intery)
%     cs2(1:91) = cs_exp(1:91).*(rhob/1000)./M;
%     cs2(92:end) = intery(2).*(rhob/1000)./M;
% elseif interpflg == 3       % step-fxn interpolation
%     cs2 = [((soil_data(z,2)*(rhob/1000))/M)*ones(90,1);((soil_data(z,1)*(rhob/1000))/M)*ones(60,1)] ;        % micro mol / cm^3 -> initial c (1300(ug/kg) = shallow PFOS soil conc in Schaefer 2022b)
% else
%     fprintf('Please specify interpflg corresponding to interpolation scheme that is applicable\n');
%     return 
% end
caw2 = 0*ones(G.cells.num,1);                           % micro mol / cm^3 -> initial caw2
c = c0;
[~, surface_pressure] = interface_state_numeric(c, isotherm_flg, ...
    sigma0, a, b, R, Tp, interface_lookup);
sigma = sigma0 - surface_pressure;                       % dyn/cm (= mN/m numerically)

init_cond = readtable('Initial_Conditions.xlsx', 'Sheet', 'Retention_' + hSMP,  'VariableNamingRule','preserve');
h0 = zeros(150,1);
% Use Linear Interpolation to Intialize h Values
% Convert theta to h before interpolation to avoid discontinuities in h
% profiles
theta_to_h = @(h0) -(((((alpha1*sigma0)./sigma).^(-n1)).*(((thetar-thetas)./(thetar-h0)).^(n1/(n1-1))-1)).^(1/n1)); % convert theta to h (pressure head)
h01 = theta_to_h(init_cond.("Water Content")(1));
h20 = theta_to_h(init_cond.("Water Content")(20));
h60 = theta_to_h(init_cond.("Water Content")(60));
h100 = theta_to_h(init_cond.("Water Content")(100));
h150 = theta_to_h(init_cond.("Water Content")(150));
h0(1) = h01(1);
h0(20) = h20(20);
h0(60) = h60(60);
h0(100) = h100(100);
h0(150) = h150(150);

h0 = interp1([1, 60, 100, 150]', [h0(1); h0(60); h0(100); h0(150)], (1:150)','linear');
%h0 = init_cond.WaterContent;
%h0 = -(((((alpha1*sigma0)./sigma).^(-n1)).*(((thetar-thetas)./(thetar-h0)).^(n1/(n1-1))-1)).^(1/n1)); % convert theta to h (pressure head)
init_cond.h = h0;

% Constant h for initial conditions
%h0 = ones(150,1).*-300;

%% Initialize for solution loop
[h, c, cs2, caw2]    = initVariablesADI(h0, c0, cs2, caw2);                                     % create h, c, cs2, caw2 as AD variables
[hIx, cIx, cs2Ix, caw2Ix] = deal(1:nc, (nc+1):(2*nc), (2*nc+1):(3*nc), (3*nc+1):(4*nc));        % Index for h and c
[tol,htol,ttol,ctol,cs2tol, caw2tol, Aawtol] = deal(1e-3,5e-4,1e-3,...
                                      1e-4,1e-2,1e-2,1e-1);             % Tolerances for local res, h, theta, c
[maxits,maxchop,it_L,it_R] = deal(10,20,6,8);                           % Maximum iteration, maximu chops, range of iterations 
                                                                        % beyond which to cut or increase time step sizes
[dt_init,dtmax] = deal(2e-5*day,1*day);                               	% Initial and maximum time step sizes
pargs = {};% {'EdgeColor','k','EdgeAlpha',.1};                          % Custom settings for plotCellData

% Initialize parameters that depend on h, c, and caw2
[Gamma_aw, surface_pressure] = interface_state_ad(c, isotherm_flg, ...
    sigma0, a, b, R, Tp, interface_lookup);
sigma = sigma0 - surface_pressure;
theta = swc(h, thetas, thetar, m1, n1, alpha1, sigma0, sigma);          % water content
kr = relperm(theta, thetas, thetar, m1, L);                             % relative permeability
cap = smc(h, thetas, thetar, m1, n1, alpha1, sigma0, sigma);            % soil moisture capacity
Aaw = Aaw_func_thermo_rfr(sigma0,poro,alpha1,n1,theta,thetar,thetas, Aawflg);% air-water area [cm^2/cm^3]

% Recover aqueous concentration from measured total soil concentration.
% Iterate because surface tension changes theta and Aaw through capillarity.
for init_iter = 1:12
    c_previous = c.val;
    c.val = solve_initial_aqueous_concentration(ctot, theta.val, Aaw.val, ...
        rhob, Kf, M, n3, isotherm_flg, sigma0, a, b, R, Tp, interface_lookup);
    [Gamma_aw, surface_pressure] = interface_state_ad(c, isotherm_flg, ...
        sigma0, a, b, R, Tp, interface_lookup);
    sigma = sigma0 - surface_pressure;
    theta = swc(h, thetas, thetar, m1, n1, alpha1, sigma0, sigma);
    kr = relperm(theta, thetas, thetar, m1, L);
    cap = smc(h, thetas, thetar, m1, n1, alpha1, sigma0, sigma);
    Aaw = Aaw_func_thermo_rfr(sigma0,poro,alpha1,n1,theta,thetar,thetas, Aawflg);
    if norm(c.val-c_previous, inf) <= 1e-8*max(1, norm(c.val, inf))
        break
    end
end

if gamma_primary_mode
    gamma_max = 100*kralchevsky_pars(3); % mol/m2 -> micromol/cm2
    gamma_initial = min(max(Gamma_aw.val, gamma_max*1e-12), gamma_max*(1-1e-12));
    gamma_logit_initial = log(gamma_initial./(gamma_max-gamma_initial));
    [h, gamma_logit, cs2, caw2] = initVariablesADI(h.val, gamma_logit_initial, ...
        cs2.val, caw2.val);
    [hIx, cIx, cs2Ix, caw2Ix] = deal(1:nc, (nc+1):(2*nc), ...
        (2*nc+1):(3*nc), (3*nc+1):(4*nc));
    [c, Gamma_aw, surface_pressure] = interface_state_gamma_primary_ad( ...
        gamma_logit, interface_lookup, gamma_max);
else
    gamma_logit = [];
end
c0 = c.val;
Kaw = Gamma_aw./(c + 1e-30);                                            % effective Kaw for diagnostics [cm]
caw1 = Faw * Gamma_aw.*Aaw;                                             % equilibrium AWI PFAS [micro mol/cm^3]
cs1 = Fs * rhob*Kf*(c.*M).^n3/M;                                        % Solid-phase adsorption [micromol/cm^3]
ctot = theta.*c + caw1 + caw2 + cs1 + rhob*cs2;                         % total mass per numerical grid cell volume [micromol/cm^3]
pfas_vd = sum(ctot.*vC);                                                % total PFAS in the vadose zone (assuming the entire domain is vadose zone)    
pfas_vd_c = sum(theta.*c.*vC);                                          % total aqueous PFAS in the vadose zone (assuming the entire domain is vadose zone)    
pfas_vd_caw = sum((caw1 + caw2).*vC);                                   % total air-water interfacial PFAS in the vadose zone
pfas_vd_cs = sum((cs1 + rhob*cs2).*vC);                                 % total solid-phase PFAS in the vadose zone (assuming the entire domain is vadose zone)    
qbottom = kr(inIx2)*rock.Ks(inIx2)*(Lx/nx)*(Ly/ny);                     % water flux at the bottom (assuming unit gradient)   
qcbottom = c(inIx2)*kr(inIx2)*rock.Ks(inIx2)*(Lx/nx)*(Ly/ny);           % PFAS flux at the bottom (assuming unit gradient for water and zero gradient for concentration)   
cum_discharge = 0*qcbottom;                                             % initialize cumulative discharge. *qcbottom is to make an AD variable
pore_v_mean = 0;                                                        % initialize mean porewater velocity
v = 0.*kr;
peak_leachate_ug_L = max(c.val(end),0)*M*1000;
peak_receptor_ug_L = peak_leachate_ug_L/groundwater_DF_value;
peak_time_years = 0;
forcing_cycle_seconds = height(receptor_forcing)*forcing_interval_seconds;
last_cycle_peak_receptor_ug_L = 0;
previous_cycle_peak_receptor_ug_L = 0;

sol = repmat(struct('time',[],'h',[], 'theta', [], 'c', [], 'Aaw', [], ...
                                    'caw1', [], 'caw2', [], 'cs1', [], ...
                                    'cs2', [], 'ctot', [], 'v', [], 'cap', [], ...
                                    'Kaw', [], 'Gamma_aw', [], ...
                                    'surface_pressure', [], 'sigma', []), [nstep+1,1]);
sol(1) = struct('time', 0, 'h', h.val,'theta', theta.val, 'c', c.val, ...   % stores initial values in solution structure
                                    'Aaw', Aaw.val, 'caw1', caw1.val, ...
                                    'caw2', caw2.val, 'cs1', cs1.val, ...
                                    'cs2', cs2.val, 'ctot', ctot.val, 'v', v.val, ...
                                    'cap', cap.val, 'Kaw', Kaw.val, ...
                                    'Gamma_aw', Gamma_aw.val, ...
                                    'surface_pressure', surface_pressure.val, ...
                                    'sigma', sigma.val);
sol_obs = repmat(struct('time',[],'pfas_vd',[],...
                                  'pfas_vd_c',[],...
                                  'pfas_vd_caw',[],...
                                  'pfas_vd_cs',[],...
                                  'qbottom',[],...
                                  'qcbottom',[],...
                                  'cum_discharge',[],...
                                  'pore_v_mean',[],...
                                  'cbottom',[],...
                                  'leachate_ug_L',[],...
                                  'receptor_ug_L',[],...
                                  'groundwater_DAF',[],...
                                  'groundwater_AF',[]),[length(Ts),1]);
sol_obs(1) = struct('time', 0, 'pfas_vd',pfas_vd.val,...                    % stores initial observation values 
                                  'pfas_vd_c',pfas_vd_c.val,...
                                  'pfas_vd_caw',pfas_vd_caw.val,...
                                  'pfas_vd_cs',pfas_vd_cs.val,...
                                  'qbottom',qbottom.val,...
                                  'qcbottom',qcbottom.val,...
                                  'cum_discharge',cum_discharge.val,...
                                  'pore_v_mean',pore_v_mean,...
                                  'cbottom',c.val(end)/ci,...
                                  'leachate_ug_L',peak_leachate_ug_L,...
                                  'receptor_ug_L',peak_receptor_ug_L,...
                                  'groundwater_DAF',groundwater_DF_value,...
                                  'groundwater_AF',groundwater_AF_value);
                                         
%% Main loop
t   = 0;                        % Initial time
dt  = dt_init;                  % Initial time step size
%t_step(1) = dt;
htop0 = 0;                      % Initial ponding depth (will be used to save temporary values in the loop)
htop = htop0;                   % Initial ponding depth
flag_ponding = 0;               % flag to denote whether ponding occurs
n = 1;                        % Time step                    
next_progress_day = progress_interval_days;
progressFile = fullfile(outputRoot, sprintf('progress_case_%d.csv', case_id));
fidProgress = fopen(progressFile,'w');
if fidProgress < 0
    error('Could not create live progress file: %s', progressFile);
end
fprintf(fidProgress, ['status,case_id,pfas,isotherm,simulated_days,', ...
    'total_days,percent,peak_receptor_ug_L,current_dt_days,time_step_cuts,', ...
    'linear_solve_rejections,state_step_rejections,step_mass_rejections,wall_seconds\n']);
fprintf(fidProgress, 'running,%d,%s,%s,0,%.10g,0,0,%.10g,0,0,0,0,0\n', ...
    case_id,char(COC),char(isotherm_name),Tf/day,dt/day);
fclose(fidProgress);
while t<Tf
    t = t + dt;
    if verbose_timesteps
        fprintf('\nTime step %d: Time %.2f -> %.2f days %.3f day \n', ...
            n, convertTo(t - dt, day), convertTo(t, day), convertTo(dt, day));
    end
    
    % Newton loop
    reshNorm   = 1e99;          % Initialize Norm with a big number
    rescNorm   = 1e99;          % Initialize Norm with a big number
    h0 = h.val;                 % save temporary values for the Newton loop
    theta0 = theta.val;
    c0 = c.val;
    if gamma_primary_mode
        gamma_logit0 = gamma_logit.val;
    end
    cs10 = cs1.val;
    caw10 = caw1.val;
    cs20 = cs2.val;
    caw20 = caw2.val;
    cap0 = cap.val;
    Gamma_aw0 = Gamma_aw.val;
    Aaw0 = Aaw.val;             % save temporary values for the Newton loop
    pfas_inventory0 = pfas_vd.val;
    nit = 0;                    % # of iterations
    nchop=0;                    % count the # of times dt is chopped
    flag_ponding = 0;           % set ponding flag to 0
    mass_step_reject = false;
    while (reshNorm > tol || rescNorm > tol || hNorm > htol || tNorm > ttol || cs2Norm > cs2tol || caw2Norm > caw2tol ...
            || cNorm > ctol || AawNorm > Aawtol || (~isreal(h.val)) || (~isreal(c.val)) || (~isreal(cs2.val)) || (~isreal(caw2.val)) ...
            || (gamma_primary_mode && (~isreal(gamma_logit.val) || sum(isnan(gamma_logit.val)))) ...
            || sum(isnan(h.val)) || sum(isnan(c.val)) || sum(isnan(cs2.val)) || sum(isnan(caw2.val)))
               
        % Relative permeability
        kr = relperm(theta, thetas, thetar, m1, L);
        % Pressure head differences across each interface
        dh = grad(h)-gradz;
        % Phase Darcy fluxes
        v  = -avg(kr).*T.*dh;               % use arithmetic average for kr
        %v  = -upw(dh.val<=0, kr).*T.*dh;   % use upstream weighting for kr
        
        % Compute phase flux for transport
        D = abs(v) ./AF ./ avg(theta) .* alphaL + avg(theta).^(7/3) / thetas^2 *Dm;     % compute dispersion coefficient
        c_upw = upw(v.val >= 0, c);                                                     % use upstream weighting for c
        c_flux = v .*c_upw - avg(theta).*D.*AF./LF.*grad(c);                            % solute flux
        
        % Conservation of water
        water = vC./dt.*(theta.val - theta0) + vC./dt.* cap.val .*(h - h.val) + div(v);
        qtop = water(inIx1);                                                            % water flux from top cell to the cell below
        % Conservation of PFAS. In conservative-inventory mode, assemble
        % the entire equilibrium inventory as one backward-Euler difference.
        % This prevents independent component linearizations from leaving a
        % persistent mass residual in the fully coupled Gamma formulation.
        if conservative_inventory_mode
            pfas_inventory = theta.*c + rhob*Kf*M^(n3-1).*c.^n3 + Aaw.*Gamma_aw;
            pfas_inventory_previous = theta0.*c0 + rhob*Kf*M^(n3-1).*c0.^n3 + ...
                Aaw0.*Gamma_aw0;
            mc = vC./dt.*(pfas_inventory-pfas_inventory_previous) + div(c_flux);
        else
            mc = vC./dt.*(theta.*c - theta0.*c0) + div(c_flux);
            % Note: n3 < 1, Jacobian will be INF when c=0
            mc = mc + Fs * vC./dt * rhob*Kf*M^(n3-1).*(c.^n3-c0.^n3);
        end
        % Include kinetic solid-phase adsorption
        mc = mc + vC .* alphas *rhob .* ((1-Fs) *Kf*M^(n3-1).*c.^n3- cs2);
        
        % Include equilibrium air-water interfacial adsorption only in the
        % component-wise legacy assembly.  The conservative assembly above
        % already contains Aaw*Gamma_aw in pfas_inventory; adding it again
        % would double-count the dominant PFOS storage term.
        if ~conservative_inventory_mode
            mc = mc + Faw * vC./dt .* (Gamma_aw .* Aaw - Gamma_aw0 .* Aaw0);
        end
        % Include kinetic air-water interfacial adsorption
        mc = mc + vC .* alphaaw .* ((1-Faw).*Aaw.*Gamma_aw - caw2);
        
        % Conservation of immobile PFAS at solid surfaces
        imcs = (cs2 - cs20)./dt - alphas * ((1-Fs)*Kf*M^(n3-1).*c.^n3 - cs2);
        
        % Conservation of immobile PFAS at air-water interfaces
        imcaw = (caw2 - caw20)./dt - alphaaw * ((1-Faw).*Aaw.*Gamma_aw - caw2);
       
        
        % Bottom boundary: free drainage for transport runs, no-flow for
        % the closed-box conservation verification.
        if closed_box_mode
            qbottom_eq = 0*kr(inIx2)*rock.Ks(inIx2)*(Lx/nx)*(Ly/ny);
            qcbottom_eq = 0*c(inIx2);
        else
            qbottom_eq = kr(inIx2)*rock.Ks(inIx2)*(Lx/nx)*(Ly/ny);
            water(inIx2) = water(inIx2) + qbottom_eq;
            qcbottom_eq = c(inIx2)*qbottom_eq;
            mc(inIx2) = mc(inIx2) + qcbottom_eq;
        end
        
        % Pressure head boundary
        %h_bottom = Lz - interp1(gw(:,1)*day,gw(:,2),t,'linear');
        %water(inIx2) = water(inIx2) - kr(inIx2)*rock.Ks(inIx2)*(2*(h_bottom - h(inIx2))/dz - 1)*(Lx/nx)*(Ly/ny);
        %mc(inIx2) = mc(inIx2) - c(inIx2)*kr(inIx2)*rock.Ks(inIx2)*(2*(h_bottom - h(inIx2))/dz - 1)*(Lx/nx)*(Ly/ny);
        
        % Infiltration rate: cm^3 /s
        eps1 = 1e-3;
        % number of time intervals with rainfall input
        nd = floor((t-eps1)/forcing_interval_seconds) + 1;
        % Repeat the 10 year data for longer simulation times
        nd = mod(nd-1, length(atmbc.ET)) + 1;
        atmbc2 = atmbc.("Rainfall (cm/day)")(nd,1);           % precipitation
        atmbc3 = atmbc.ET(nd,1);                     % evapotranspiration    
        atmbc4 = atmbc.("PFAS infiltration")(nd,1);       % water infiltration due to AFFF application
% net infiltration rate (precipitation - evapotranspiration)
        inRate = atmbc2 - atmbc3;
        %inRate_save(n) = inRate.*(24*60*60);
       % t_step(n+1) = t_step(n) + dt*(n>1);
        
        %Include influx of water due to AFFF application
        if t <= min(Tf, afff_active_years*365.25*day)
           % Apply a nonzero AFFF source only during the configured window.
           inRate = inRate + atmbc4;
           %inRate_save(n) = inRate;
        end
        
        % Set top boundary condition for water infiltration
        if closed_box_mode
            % No top water or PFAS flux in the conservation test.
            flag_ponding = 0;
        elseif h.val(inIx1) <= (hcritl(inIx1)+eps) && inRate < 0
            % switch to pressure head BC when it is dry
            water(inIx1) = h(inIx1) - hcritl(inIx1);
            % no ponding (flag = 0)
            flag_ponding = 0;
        elseif h.val(inIx1) >= 0 && htop >= 0
            % switch to pressure head when there is water ponding
            htop = htop0 + (inRate - qtop.val)*dt/(Lx/nx)/(Ly/ny);
            water(inIx1) = h(inIx1) - htop;
            % ponding arises (flag = 1)
            flag_ponding = 1;
        else
            % If not dry or ponding, then use flux BC
            water(inIx1) = water(inIx1) - inRate;
            flag_ponding = 0;
        end
        
        %Influx of PFAS: micro mol / s
        if ~closed_box_mode && t <= min(Tf, afff_active_years*365.25*day)
           mc(inIx1) = mc(inIx1) - atmbc4*ci;
        end
        
        % Collect and concatenate all equations (i.e., assemble and
        % linearize system)
        eqs = {water,mc,imcs,imcaw};
        eq  = cat(eqs{:});
        
        % Compute Newton update and update variable
        res = eq.val;
        reshFresh = norm(res(hIx).*dt./(vC.*rock.poro),inf);
        if conservative_inventory_mode
            % Scale the PFAS residual by actual local total PFAS storage,
            % rather than the legacy 100-ug/L reference concentration. This
            % makes the nonlinear solve resolve the mass that is present.
            pfas_storage_scale = max(abs(pfas_inventory_previous), 1e-12);
            rescFresh = norm(res(cIx).*dt./(vC.*pfas_storage_scale),inf);
        else
            rescFresh = norm(res(cIx).*dt./(vC.*max(theta.val, 1e-20).*ci),inf);
        end
        conservative_mass_ok = true;
        if conservative_inventory_mode
            source_rate_fresh = 0;
            if ~closed_box_mode && t <= min(Tf, afff_active_years*365.25*day)
                source_rate_fresh = ci*atmbc4;
            end
            pfas_inventory_fresh = sum(pfas_inventory.*vC);
            conservative_step_mass_error = abs(source_rate_fresh*dt - ...
                qcbottom_eq*dt - (pfas_inventory_fresh-pfas_inventory0))/ ...
                max(abs(pfas_inventory0), eps);
            conservative_mass_ok = conservative_step_mass_error <= ...
                conservative_step_mass_tolerance;
        end
        % Do not accept an iterate merely because the residual from the
        % preceding state was small. Reassembly at the current state is
        % required before an accepted time step.
        if conservative_inventory_mode && nit > 0 && reshFresh <= tol && rescFresh <= tol && ...
                hNorm <= htol && tNorm <= ttol && cNorm <= ctol && ...
                cs2Norm <= cs2tol && caw2Norm <= caw2tol && AawNorm <= Aawtol && ...
                conservative_mass_ok
            qbottom = kr(inIx2)*rock.Ks(inIx2)*(Lx/nx)*(Ly/ny);
            qcbottom = qcbottom_eq;
            reshNorm = reshFresh;
            rescNorm = rescFresh;
            break
        end
        % The coupled equation blocks use very different physical units
        % (head, dissolved PFAS, and immobile PFAS stores).  Row/column
        % equilibration preserves every Newton root while avoiding false
        % singular-matrix decisions caused by those units.  Apply it to both
        % Langmuir and Kralchevsky formulations; Langmuir's raw Jacobian can
        % otherwise become badly scaled during wet/dry forcing transitions.
        [upd, linear_solve_failed] = scaled_newton_update(eq.jac{1}, res);
        if gamma_primary_mode && conservative_inventory_mode
            % A full Newton step can temporarily drive the logit many orders
            % of magnitude away from the admissible Kralchevsky branch.  A
            % coupled damping factor keeps Gamma changes within a factor of
            % exp(2) per iteration and prevents the next Jacobian from being
            % evaluated on a numerically flat tail of the logistic map.
            maxGammaLogitStep = norm(upd(cIx), inf);
            if maxGammaLogitStep > 2
                upd = upd .* (2/maxGammaLogitStep);
            end
        elseif conservative_inventory_mode
            % Langmuir retains aqueous concentration as the primary unknown.
            % Use a fraction-to-boundary rule to preserve C > 0.  Do not
            % impose a lower bound on the damping factor: cells initialized
            % near zero may legitimately require alpha << 1e-8.  The old
            % 1e-8 floor could cross C through zero, create fractional-power
            % complex values, and enter an unlimited time-step chop loop.
            cUpdate = upd(cIx);
            negativeUpdate = cUpdate < 0;
            if any(negativeUpdate)
                alphaC = 0.99*min(c.val(negativeUpdate) ./ ...
                    (-cUpdate(negativeUpdate)));
                if isfinite(alphaC) && alphaC < 1
                    upd = upd .* max(alphaC, 0);
                end
            end
        end
        h.val  = h.val   + upd(hIx);
        if gamma_primary_mode
            gamma_logit.val = gamma_logit.val + upd(cIx);
        else
            c.val = c.val + upd(cIx);
        end
        cs2.val = cs2.val + upd(cs2Ix);
        caw2.val = caw2.val + upd(caw2Ix);
        
        % Update parameters that depend on h and c
        if gamma_primary_mode
            [c, Gamma_aw, surface_pressure] = interface_state_gamma_primary_ad( ...
                gamma_logit, interface_lookup, gamma_max);
        else
            [Gamma_aw, surface_pressure] = interface_state_ad(c, isotherm_flg, ...
                sigma0, a, b, R, Tp, interface_lookup);
        end
        sigma = sigma0 - surface_pressure;
        theta = swc(h, thetas, thetar, m1, n1, alpha1, sigma0, sigma);
        % swc is already bounded by construction.  Do not clip only the ADI
        % values here: changing .val without changing the Jacobian makes the
        % residual and its derivative describe different states.
        cap = smc(h, thetas, thetar, m1, n1, alpha1, sigma0, sigma);
        kr = relperm(theta, thetas, thetar, m1, L);
        %v  = -avg(kr).*T.*dh;
        Kaw = Gamma_aw./(c + 1e-30);
        Aawm = Aaw.val;
		%Aaw = iac(theta, thetas, m2, n2, alpha2);
        Aaw = Aaw_func_thermo_rfr(sigma0,poro,alpha1,n1,theta,thetar,thetas, Aawflg);
        caw1 = Faw * Gamma_aw.*Aaw;
        cs1 = Fs * rhob*Kf*(c.*M).^n3/M; 
        ctot = theta.*c + caw1 + caw2 + cs1 + rhob*cs2;
        pfas_vd = sum(ctot.*vC);
        pfas_vd_c = sum(theta.*c.*vC);
        pfas_vd_caw = sum((caw1 + caw2).*vC);
        pfas_vd_cs = sum((cs1 + rhob*cs2).*vC);
        % Free drainage boundary (bottom)
        if closed_box_mode
            qbottom = 0*kr(inIx2)*rock.Ks(inIx2)*(Lx/nx)*(Ly/ny);
            qcbottom = 0*c(inIx2);
        else
            qbottom = kr(inIx2)*rock.Ks(inIx2)*(Lx/nx)*(Ly/ny);
            qcbottom = c(inIx2)*qbottom;
        end
        % Pressure head boundary (bottom)
        %qbottom = -kr(inIx2)*rock.Ks(inIx2)*(2*(h_bottom - h(inIx2))/dz - 1)*(Lx/nx)*(Ly/ny);
        %qcbottom = -c(inIx2)*kr(inIx2)*rock.Ks(inIx2)*(2*(h_bottom - h(inIx2))/dz - 1)*(Lx/nx)*(Ly/ny);
        
        % Compute mean pore velocity in the vadose zone
        pore_v = abs(v)./ avg(theta);
        pore_v_mean = mean(pore_v.val);
        
        % Compute norms
        reshNorm = norm(res(hIx).*dt./(vC.*rock.poro),inf);
        rescNorm = norm(res(cIx).*dt./(vC.*theta.val.*ci),inf);
        hNorm = norm(upd(hIx),inf);
        tNorm = norm(cap.val.*upd(hIx),inf);
        cNorm = norm(upd(cIx),inf);
        cs2Norm = norm(upd(cs2Ix),inf);
        caw2Norm = norm(upd(caw2Ix),inf);
        AawNorm = norm(Aaw.val-Aawm,inf);
        nit     = nit + 1;
        converged_this_iteration = reshNorm <= tol && rescNorm <= tol && ...
            hNorm <= htol && tNorm <= ttol && cNorm <= ctol && ...
            cs2Norm <= cs2tol && caw2Norm <= caw2tol && AawNorm <= Aawtol;
        if converged_this_iteration && isfinite(step_mass_tolerance)
            source_rate = 0;
            if ~closed_box_mode && t <= min(Tf, afff_active_years*365.25*day)
                source_rate = ci*atmbc4;
            end
            step_mass_residual = source_rate*dt - qcbottom.val*dt - ...
                (pfas_vd.val-pfas_inventory0);
            step_mass_error = abs(step_mass_residual)/max(abs(pfas_inventory0), eps);
            mass_step_reject = step_mass_error > step_mass_tolerance;
            if mass_step_reject
                fprintf('reject mass-imbalanced step: rel error %.3e > %.3e\n', ...
                    step_mass_error, step_mass_tolerance);
            end
        end
        %fprintf('  Iteration %3d:  Res = %.4e, hRes = %.4e, tRes = %.4e, cRes = %.4e, dt = %.4e\n', nit, reshNorm, hNorm, tNorm, cNorm, dt);
        
        % Reject any inadmissible state through one bounded retry path.
        % Previously, NaN/complex and mass-rejection branches bypassed
        % maxchop and could halve dt forever without another progress line.
        state_invalid = (~isreal(h.val)) || (~isreal(c.val)) || ...
            (~isreal(cs2.val)) || (~isreal(caw2.val)) || ...
            any(~isfinite(h.val)) || any(~isfinite(c.val)) || ...
            any(~isfinite(cs2.val)) || any(~isfinite(caw2.val)) || ...
            any(c.val <= 0) || any(~isfinite(theta.val)) || ...
            any(~isreal(theta.val)) || any(~isfinite(sigma.val)) || ...
            any(~isreal(sigma.val)) || any(sigma.val <= 0) || ...
            (gamma_primary_mode && ((~isreal(gamma_logit.val)) || ...
            any(~isfinite(gamma_logit.val)))) || ...
            (~isreal(pfas_vd.val)) || any(~isfinite(pfas_vd.val));
        retry_for_iterations = nit >= maxits-2;
        retry_required = retry_for_iterations || linear_solve_failed || ...
            mass_step_reject || state_invalid;
        if retry_required
            if linear_solve_failed
                retry_reason = 'linear_solve';
            elseif state_invalid
                retry_reason = 'inadmissible_state';
            elseif mass_step_reject
                retry_reason = 'step_mass_balance';
            else
                retry_reason = 'newton_iterations';
            end
            proposed_dt = dt/2.0;
            step_start_time = t-dt;
            cannot_advance = (step_start_time + proposed_dt == step_start_time);
            if nchop >= maxchop || proposed_dt < minimum_dt_seconds || cannot_advance
                diagnostic = struct( ...
                    'case_id',case_id,'pfas',char(COC),'model',char(isotherm_name), ...
                    'reason',retry_reason,'step_start_day',step_start_time/day, ...
                    'attempted_end_day',t/day,'dt_seconds',dt, ...
                    'proposed_dt_seconds',proposed_dt,'newton_iteration',nit, ...
                    'step_chops',nchop,'total_time_step_cuts',time_step_cuts, ...
                    'linear_solve_rejections',linear_solve_rejections, ...
                    'state_step_rejections',state_step_rejections, ...
                    'h',h.val,'c',c.val,'theta',theta.val,'sigma',sigma.val, ...
                    'residual',res,'h_update_norm',hNorm,'c_update_norm',cNorm, ...
                    'water_residual_norm',reshFresh,'pfas_residual_norm',rescFresh);
                diagnosticFile = fullfile(outputRoot, sprintf( ...
                    'FAILED_STEP_case_%d_day_%.8f.mat',case_id,step_start_time/day));
                save(diagnosticFile,'diagnostic','-v7');
                error(['Adaptive time-step retry limit reached: case=%d, ', ...
                    'day=%.10g, reason=%s, chops=%d, dt=%.6g s. ', ...
                    'Diagnostic saved to %s'],case_id,step_start_time/day, ...
                    retry_reason,nchop,dt,diagnosticFile);
            end
            if verbose_timesteps
                fprintf(['cut time step: reason=%s chop=%d/%d ', ...
                    'dt = %.4e -> %.4e s\n'],retry_reason,nchop+1,maxchop,dt,proposed_dt);
            end
            time_step_cuts = time_step_cuts + 1;
            nchop=nchop+1;
            t = t - dt;
            dt = proposed_dt;
            t = t + dt;
            nit=0;
            h.val=h0;
            if gamma_primary_mode
                gamma_logit.val = gamma_logit0;
                [c, Gamma_aw, surface_pressure] = interface_state_gamma_primary_ad( ...
                    gamma_logit, interface_lookup, gamma_max);
            else
                c.val = c0;
                [Gamma_aw, surface_pressure] = interface_state_ad(c, isotherm_flg, ...
                    sigma0, a, b, R, Tp, interface_lookup);
            end
            caw1.val = caw10;
            cs1.val = cs10;
            caw2.val = caw20;
            cs2.val = cs20;
            theta.val = theta0;
            cap.val = cap0;
            Aaw.val = Aaw0;
            sigma = sigma0 - surface_pressure;
            theta = swc(h, thetas, thetar, m1, n1, alpha1, sigma0, sigma);
            cap = smc(h, thetas, thetar, m1, n1, alpha1, sigma0, sigma);
            Kaw = Gamma_aw./(c + 1e-30);
			%Aaw = iac(theta, thetas, m2, n2, alpha2);
            Aaw = Aaw_func_thermo_rfr(sigma0,poro,alpha1,n1,theta,thetar,thetas, Aawflg);
            caw1 = Faw * Gamma_aw.*Aaw;
            htop = htop0;
            reshNorm = 1e99;
            if mass_step_reject
                step_mass_rejections = step_mass_rejections + 1;
                mass_step_reject = false;
            end
            if linear_solve_failed
                linear_solve_rejections = linear_solve_rejections + 1;
            end
            if state_invalid
                state_step_rejections = state_step_rejections + 1;
            end
            continue
        end
        % Conservative-inventory verification requires a fresh residual
        % assembly after each update; legacy production mode retains its
        % original convergence sequence.
        if conservative_inventory_mode
            reshNorm = 1e99;
            rescNorm = 1e99;
        end
        if nit>(maxits + conservative_inventory_mode)
             break
        end
    end
     
    if nit > maxits
        error('Newton solves did not converge')
    else 
        if closed_box_mode
            q = 0;
        else
            q = -(kr.val(inIx1)+kr.val(inIx1+1))/2*rock.Ks(inIx1)*((h.val(inIx1+1)-h.val(inIx1))/dz-1)*(Lx/nx)*(Ly/ny);
            total_in_w = total_in_w + q*dt;
            total_out_w = total_out_w + qbottom.val*dt;
            if t <= min(Tf, afff_active_years*365.25*day)
                total_in_c = total_in_c + ci*atmbc4*dt;
            end
            total_out_c = total_out_c + qcbottom.val*dt;
            cum_discharge = cum_discharge + qcbottom.val * dt;
        end
        current_leachate_ug_L = max(c.val(end),0)*M*1000;
        current_receptor_ug_L = current_leachate_ug_L/groundwater_DF_value;
        if current_receptor_ug_L > peak_receptor_ug_L
            peak_leachate_ug_L = current_leachate_ug_L;
            peak_receptor_ug_L = current_receptor_ug_L;
            peak_time_years = t/(365.25*day);
        end
        if t > max(0,Tf-forcing_cycle_seconds)
            last_cycle_peak_receptor_ug_L = max( ...
                last_cycle_peak_receptor_ug_L,current_receptor_ug_L);
        elseif t > max(0,Tf-2*forcing_cycle_seconds)
            previous_cycle_peak_receptor_ug_L = max( ...
                previous_cycle_peak_receptor_ug_L,current_receptor_ug_L);
        end
        simulated_days_now = t/day;
        if simulated_days_now + 1e-9 >= next_progress_day || t >= Tf
            progress_percent = min(100,100*t/Tf);
            progress_wall_seconds = toc(scenario_wall_timer);
            fprintf(['PROGRESS case=%d pfas=%s model=%s day=%.2f/%.2f ', ...
                '(%.1f%%) peak_receptor=%.6g ug/L dt=%.4g d cuts=%d ', ...
                'linear_retry=%d state_retry=%d mass_reject=%d wall=%.1f s\n'], ...
                case_id,char(COC),char(isotherm_name),simulated_days_now,Tf/day, ...
                progress_percent,peak_receptor_ug_L,dt/day,time_step_cuts, ...
                linear_solve_rejections,state_step_rejections,step_mass_rejections, ...
                progress_wall_seconds);
            fidProgress = fopen(progressFile,'w');
            if fidProgress >= 0
                fprintf(fidProgress, ['status,case_id,pfas,isotherm,simulated_days,', ...
                    'total_days,percent,peak_receptor_ug_L,current_dt_days,time_step_cuts,', ...
                    'linear_solve_rejections,state_step_rejections,step_mass_rejections,wall_seconds\n']);
                fprintf(fidProgress, 'running,%d,%s,%s,%.10g,%.10g,%.8g,%.10g,%.10g,%d,%d,%d,%d,%.10g\n', ...
                    case_id,char(COC),char(isotherm_name),simulated_days_now,Tf/day, ...
                    progress_percent,peak_receptor_ug_L,dt/day,time_step_cuts, ...
                    linear_solve_rejections,state_step_rejections,step_mass_rejections, ...
                    progress_wall_seconds);
                fclose(fidProgress);
            end
            next_progress_day = (floor(simulated_days_now/progress_interval_days)+1) ...
                *progress_interval_days;
        end
        % Save profile results for specified time points
        if ismember(t,Ts)
            istep=find(Ts==t);
            sol(istep) = struct('time', t, ...
                          'h', h.val, ...
                          'theta', theta.val, ...
                          'c', c.val, ...
                          'Aaw', Aaw.val,...
                          'caw1', caw1.val,...
                          'caw2', caw2.val,...
                          'cs1', cs1.val,...
                          'cs2', cs2.val,...
                          'ctot', ctot.val,...
                          'v', v.val,...
                          'cap', cap.val,...
                          'Kaw', Kaw.val,...
                          'Gamma_aw', Gamma_aw.val,...
                          'surface_pressure', surface_pressure.val,...
                          'sigma', sigma.val);
            sol_obs(istep) = struct('time', t, ...
                          'pfas_vd',pfas_vd.val,...
                          'pfas_vd_c',pfas_vd_c.val,...
                          'pfas_vd_caw',pfas_vd_caw.val,...
                          'pfas_vd_cs',pfas_vd_cs.val,...
                          'qbottom',qbottom.val,...
                          'qcbottom',qcbottom.val,...
                          'cum_discharge',cum_discharge.val,...
                          'pore_v_mean',pore_v_mean,...
                          'cbottom',c.val(end)/ci,...
                          'leachate_ug_L',current_leachate_ug_L,...
                          'receptor_ug_L',current_receptor_ug_L,...
                          'groundwater_DAF',groundwater_DF_value,...
                          'groundwater_AF',groundwater_AF_value);
        end
       
    end
    
    % Increase or cut time step size based on it_L or it_R
    if(nit<=it_L)
        % increase time step size when iteration # is small
        dt = dt*1.8;
        dt = min(dt,dtmax);
    elseif(nit>it_R)
        % reduce time step size when iteration # is big
        dt = dt/2.;
    end
    
    % Adjust time step size to hit the time points where results are saved
    DT1 = Ti-t;
    DT1 = DT1(DT1>0);
    DT1 = min(DT1);
    DT2 = Ts-t;
    DT2 = DT2(DT2>0);
    DT2 = min(DT2);
    dt = min(dt,DT1);
    dt = min(dt,DT2);
    
    % prevent ponding depth htop from getting negative
    if flag_ponding == 1
       htop0 = max(htop,0);
    end
    
    n=n+1;
end

%inRate_save = inRate_save'.*(24*60*60);
%t_step = t_step'./(24*60*60);

% Check global mass balance
total_w = sum((sol(end).theta - thetar) .* vC);
total_w_added = sum((sol(end).theta - sol(1).theta).*vC);
if ~closed_box_mode
    % q is evaluated below the top control volume, so its storage change is
    % included as a boundary correction only for the open transport model.
    total_in_w = total_in_w + ...
        (sol(end).theta(inIx1)-sol(1).theta(inIx1))*vC(inIx1);
end
water_balance_residual = total_in_w - total_out_w - total_w_added;
if closed_box_mode
    initial_water_inventory = sum(sol(1).theta.*vC);
    merrw = abs(water_balance_residual)/max(abs(initial_water_inventory), eps);
else
    merrw = abs(water_balance_residual)/max(abs(total_in_w-total_out_w), eps);
end

total_c = sum((sol(end).ctot) .* vC);
total_c_added = sum((sol(end).ctot - sol(1).ctot).*vC);
pfas_balance_residual = total_in_c - total_out_c - total_c_added;
merrc = abs(pfas_balance_residual)/max(abs(total_in_c-total_out_c), eps);
initial_pfas_inventory = sum(sol(1).ctot.*vC);
merrc_inventory = abs(pfas_balance_residual)/max(abs(initial_pfas_inventory), eps);

%% Save results
run_wall_seconds = toc(scenario_wall_timer);
if isotherm_flg == 1
    interface_parameter_tag = char(langmuir_parameter_set);
else
    interface_parameter_tag = 'not_applicable';
end
if gamma_primary_mode
    interface_primary_tag = 'Gamma_logit';
else
    interface_primary_tag = 'aqueous_concentration';
end
if conservative_inventory_mode
    inventory_assembly_tag = 'single_total_inventory_BE';
else
    inventory_assembly_tag = 'componentwise_legacy';
end
solver_formulation_version = 'PFAS_transport_SSL_2026-08_SIparams_concentration_v5_bounded_retry';
if ryan_validation
    runTag = sprintf('%s_%s_%s_%s_RyanValidation60day_%s_%s', char(SMP), ...
        char(COC), char(isotherm_name), interface_parameter_tag, char(param), char(cs2_cond));
else
    runTag = sprintf('%s_%s_%s_%s_NaCl%gmM_%dcycles_%s_%s', char(SMP), ...
        char(COC), char(isotherm_name), interface_parameter_tag, NaCl_run_mM, simulation_cycles, ...
        char(param), char(cs2_cond));
end
simulated_years = Tf/(365.25*day);
ssl_profile_peak_ug_kg = initial_scale*max(initial_profile_ug_kg);
ssl_profile_mean_ug_kg = initial_scale*mean(initial_profile_ug_kg);
peak_in_last_forcing_cycle = peak_time_years >= max(0, simulated_years-10.0);
horizon_peak_covered = simulated_years > 19.9 && ...
    ~peak_in_last_forcing_cycle && previous_cycle_peak_receptor_ug_L > 0 && ...
    last_cycle_peak_receptor_ug_L <= 1.02*previous_cycle_peak_receptor_ug_L;
responseTable = table(case_id,string(COC),string(isotherm_name),initial_scale, ...
    ssl_profile_peak_ug_kg,ssl_profile_mean_ug_kg,simulated_years, ...
    peak_leachate_ug_L,groundwater_DF_value,groundwater_AF_value,peak_receptor_ug_L, ...
    peak_time_years,receptor_mcl_ug_L,peak_receptor_ug_L/receptor_mcl_ug_L, ...
    peak_in_last_forcing_cycle,last_cycle_peak_receptor_ug_L, ...
    previous_cycle_peak_receptor_ug_L, ...
    last_cycle_peak_receptor_ug_L/max(previous_cycle_peak_receptor_ug_L,eps), ...
    horizon_peak_covered,net_infiltration_cm_yr,receptor.W_cm, ...
    receptor.Ugw_cm_yr,receptor.bsat_cm,receptor.alpha_v_cm,merrc_inventory,merrw, ...
    run_wall_seconds,langmuir_KL,langmuir_gamma_m,a,b, ...
    string(interface_primary_tag),string(inventory_assembly_tag), ...
    string(solver_formulation_version),conservative_step_mass_tolerance, ...
    step_mass_rejections,time_step_cuts,linear_solve_rejections, ...
    state_step_rejections,minimum_dt_seconds,NaCl_run_mM,interface_temperature_K, ...
    use_activity_correction,kralchevsky_pars(1),kralchevsky_pars(2),kralchevsky_pars(3), ...
    'VariableNames',{'Case_ID','PFAS','Isotherm','Initial_Profile_Scale', ...
    'SSL_Profile_Peak_ug_kg','SSL_Profile_Mean_ug_kg','Simulated_Years', ...
    'Peak_Leachate_ug_L','Groundwater_DAF','Groundwater_AF','Peak_Receptor_ug_L', ...
    'Peak_Time_Years','Target_MCL_ug_L','Peak_to_Target_Ratio', ...
    'Peak_In_Last_10yr','Last_10yr_Max_Receptor_ug_L', ...
    'Previous_10yr_Max_Receptor_ug_L','Last_to_Previous_10yr_Peak_Ratio', ...
    'Horizon_Peak_Covered','Net_Infiltration_cm_yr','Source_Width_cm', ...
    'Groundwater_Velocity_cm_yr','Saturated_Thickness_cm','Vertical_Dispersivity_cm', ...
    'PFAS_Mass_Error_Inventory_Normalized','Water_Mass_Error_Flux_Normalized', ...
    'Wall_Seconds','Langmuir_KL_m3_mol','Langmuir_Gamma_m_mol_m2', ...
    'Langmuir_a_mol_m3','Langmuir_b','Interface_Primary', ...
    'PFAS_Inventory_Assembly','Solver_Formulation_Version', ...
    'Conservative_Step_Mass_Tolerance','Step_Mass_Rejections','Time_Step_Cuts', ...
    'Linear_Solve_Rejections','State_Step_Rejections','Minimum_dt_seconds', ...
    'NaCl_mM','Temperature_K','Use_Activity_Correction','MassAction_K1_m3_mol', ...
    'MassAction_K2_m6_mol2','MassAction_Gamma_m_mol_m2'});
writetable(responseTable,fullfile(outputRoot,'ssl_forward_response.csv'));

if ~ssl_response_only
matFile = fullfile(outputRoot, "results_" + string(runTag) + ".mat");
save(matFile, '-v7.3', 'sol', 'G', 'total_in_c', 'total_c_added', ...
    'total_out_c', 'total_in_w', 'total_out_w', 'total_w_added', ...
    'merrc', 'merrc_inventory', 'merrw', 'pfas_balance_residual', ...
    'water_balance_residual', 'isotherm_name', 'isotherm_flg', 'NaCl_run_mM', ...
    'interface_temperature_K', 'kralchevsky_pars', 'use_activity_correction', 'a', 'b', 'sigma0', ...
    'langmuir_parameter_set', 'langmuir_KL', 'langmuir_gamma_m', ...
    'langmuir_fit_temperature_K', 'initial_scale', 'responseTable', ...
    'groundwater_DF_value', 'groundwater_AF_value', 'receptor', 'net_infiltration_cm_yr', ...
    'peak_leachate_ug_L', 'peak_receptor_ug_L', 'peak_time_years', ...
    'output_interval_days', 'simulation_cycles', 'afff_active_years', ...
    'case_id', 'run_wall_seconds', 'ryan_validation', ...
    'gamma_primary_mode', 'conservative_inventory_mode', ...
    'conservative_step_mass_tolerance', 'step_mass_rejections', ...
    'time_step_cuts', 'linear_solve_rejections', 'state_step_rejections', ...
    'minimum_dt_seconds', ...
    'solver_formulation_version', ...
    'forcing_interval_seconds');

runtimeFile = fullfile(outputRoot, "results_runtime_" + string(runTag) + ".csv");
runtimeTable = table(string(COC),string(SMP),string(isotherm_name), ...
    string(interface_parameter_tag),Tf/day,output_interval_days, ...
    run_wall_seconds,run_wall_seconds/(Tf/day), ...
    'VariableNames',{'PFAS','SMP','Isotherm','InterfaceParameterSet','SimulatedDays', ...
    'OutputIntervalDays','WallSeconds','SecondsPerSimulatedDay'});
writetable(runtimeTable,runtimeFile);
fprintf('\nCompleted %s / %s / %s: %.1f simulated days in %.2f s\n', ...
    char(COC),char(SMP),char(isotherm_name),Tf/day,run_wall_seconds);

% Write time series
temp_time = [sol_obs.time]';
temp_vd = [sol_obs.pfas_vd]';
temp_vdc = [sol_obs.pfas_vd_c]';
temp_vdcaw = [sol_obs.pfas_vd_caw]';
temp_vdcs = [sol_obs.pfas_vd_cs]';
temp_qbottom = [sol_obs.qbottom]';
temp_qcbottom = [sol_obs.qcbottom]';
temp_cum_discharge = [sol_obs.cum_discharge]';
temp_pore_v_mean = [sol_obs.pore_v_mean]';
temp_leachate_ug_L = [sol_obs.leachate_ug_L]';
temp_receptor_ug_L = [sol_obs.receptor_ug_L]';
temp_groundwater_DAF = [sol_obs.groundwater_DAF]';
temp_groundwater_AF = [sol_obs.groundwater_AF]';
% write header
cHeader = {'time' 'vd' 'vd_c' 'vd_caw' 'vd_cs' ...
            'qbottom' 'qcbottom' 'cum_discharge' 'pore_v_mean' ...
            'leachate_ug_L' 'receptor_ug_L' 'groundwater_DAF' 'groundwater_AF'}; %header
textHeader = strjoin(cHeader, ',');
%write header to file
timeSeriesFile = fullfile(outputRoot, "results_time_series_" + string(runTag) + ".csv");
fid = fopen(timeSeriesFile,"w"); 
fprintf(fid,'%s\n',textHeader);
fclose(fid);
writematrix([temp_time temp_vd temp_vdc temp_vdcaw temp_vdcs ...
    temp_qbottom temp_qcbottom temp_cum_discharge temp_pore_v_mean ...
    temp_leachate_ug_L temp_receptor_ug_L temp_groundwater_DAF temp_groundwater_AF], ...
    timeSeriesFile, 'WriteMode', 'append');

% Write mass balance information
cHeader = {'total_in_c' 'total_c_added' 'total_out_c' ...
    'pfas_balance_residual' 'pfas_error_flux_normalized' ...
    'pfas_error_inventory_normalized' 'total_in_w' 'total_out_w' ...
    'total_w_added' 'water_balance_residual' 'water_error_flux_normalized'};
textHeader = strjoin(cHeader, ',');
%write header to file
massBalanceFile = fullfile(outputRoot, "results_mass_balance_" + string(runTag) + ".csv");
fid = fopen(massBalanceFile,'w'); 
fprintf(fid,'%s\n',textHeader);
fclose(fid);
writematrix([total_in_c total_c_added total_out_c pfas_balance_residual ...
    merrc merrc_inventory total_in_w total_out_w total_w_added ...
    water_balance_residual merrw], massBalanceFile, 'WriteMode', 'append');
end

fprintf(['\nSSL response: scale=%.8g, profile peak=%.8g ug/kg, ', ...
    'receptor peak=%.8g ug/L at %.3f years, target=%.8g ug/L, DAF=%.6g, AF=%.6g\n'], ...
    initial_scale,ssl_profile_peak_ug_kg,peak_receptor_ug_L, ...
    peak_time_years,receptor_mcl_ug_L,groundwater_DF_value,groundwater_AF_value);

fidProgress = fopen(progressFile,'w');
if fidProgress >= 0
    fprintf(fidProgress, ['status,case_id,pfas,isotherm,simulated_days,', ...
        'total_days,percent,peak_receptor_ug_L,current_dt_days,time_step_cuts,', ...
        'linear_solve_rejections,state_step_rejections,step_mass_rejections,wall_seconds\n']);
    fprintf(fidProgress, 'complete,%d,%s,%s,%.10g,%.10g,100,%.10g,%.10g,%d,%d,%d,%d,%.10g\n', ...
        case_id,char(COC),char(isotherm_name),Tf/day,Tf/day, ...
        peak_receptor_ug_L,dt/day,time_step_cuts,linear_solve_rejections, ...
        state_step_rejections,step_mass_rejections,run_wall_seconds);
    fclose(fidProgress);
end

        end
    end
end
%% Functions used in the code
function kr = relperm(theta, thetas, thetar, m1, L)
% relative permeability
    se = (theta-thetar)./(thetas-thetar);
    temp = se.^L .* (1 - (1 - se.^(1/m1)).^m1).^2;        %add L from hydrus as input parameter (where 0.5 is)
    kr = temp.*(theta.val<thetas) + 1.*(theta.val==thetas);
end

function theta = swc(h, thetas, thetar, m1, n1, alpha1, sigma0, sigma)
% SWC function to link h and theta    
    % turn on surfactant-induced flow    
    temp = thetar + (thetas-thetar)./(1+(sigma0./sigma.*alpha1.*abs(h)).^n1).^m1;
    % turn off surfactant-induced flow
    %temp = thetar + (thetas-thetar)./(1+(alpha1.*abs(h)).^n1).^m1;
    theta = temp.*(h.val<0) + thetas.*(h.val>=0);
end

function dthetadh = smc(h, thetas, thetar, m1, n1, alpha1, sigma0, sigma)
% soil moisture capacity
    % turn on surfactant-induced flow
    temp = m1*n1*(sigma0./sigma.*alpha1).^n1.*(thetas-thetar).*(1+(sigma0./sigma.*alpha1.*abs(h)).^n1).^(-1-m1).*abs(h).^(n1-1);
    % turn off surfactant-induced flow
    %temp = m1*n1*(alpha1).^n1.*(thetas-thetar).*(1+(alpha1.*abs(h)).^n1).^(-1-m1).*abs(h).^(n1-1);
    dthetadh = temp.*(h.val<0) + 0.*(h.val>=0);
end

function Aaw = Aaw_func_thermo_rfr(sigma0,poro,alpha1,n1,theta,thetar,thetas, Aawflg)
%Computing air-water interfacial area using the thermodynamic approach

% Aaw       air-water interfacial area (cm^2/cm^3)
% sigma0    surface tension (dyn/cm)
% poro      porosity (-)
% alphal    V-G parameter (cm^-1)
% n1         V-G parameter (-)
% theta        water content (-)
% thetar       residual water content (-)
% thetas       saturated water content (-)
% rhow      water density (kg/m^3)
% g         gravity acceleration (m/s^2)
% Pc        capillary pressure
% sf        scaling factor to correct the thermodyanmic-based Aaw (see Guo et al., 2022) (-)

%{
% Copyright 2021-2022 Bo Guo (University of Arizona, Email: boguo@arizona.edu or guobo07@gmail.com).

% This file is part of the implementation for the PFAS-LEACH-Screening analytical 
% solver presented in the article Guo et al. (2022) AWR

% Guo, B., Zeng, J., Brusseau, M.L. and Zhang, Y., 2022. 
% A screening model for quantifying PFAS leaching in the vadose zone and 
% mass discharge to groundwater. Advances in Water Resources, 160, p.104102.

% Development of PFAS-LEACH-Screening is supported by the ESTCP Project ER21-5041. 
% Project Webpage: https://www.serdp-estcp.org/Program-Areas/Environmental-Restoration/ER21-5041

% The PFAS-LEACH-Screening analytical solver is distributed in the hope that 
% it will be useful, but WITHOUT ANY WARRANTY without even the implied warranty 
% of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details, <http://www.gnu.org/licenses/>.
%}

%sf = 4.412; %4.15;              % value set based on what was used for Vinton soil in Guo et al., 2022
d50 = 0.3;
rhow = 1000;
g = 9.81;
m = 1 - 1/n1;
Sr = thetar/thetas;
N = 200;
dSw = (1-theta/thetas)/N;
Sw = zeros(N+1,length(theta.val));
Sw(1,:) = (theta.val/thetas)' ;
Sw(N+1,:) = ones(1,length(theta));
Sw(2:N,:) = Sw(1,:) + (1:N-1)' .* dSw.val';
Pc = @(Sw) (((1-Sr)./(Sw-Sr)).^(1/m) - 1).^(1/n1)/alpha1/100*rhow*g;
if Aawflg == 1
    sf = (-0.65.*Sw(1,:) + 1.33).*(-0.45*d50 + 5);
elseif Aawflg == 2
    sf = 1;
elseif Aawflg == 3
    sf = 0;
else
    fprintf('Please specify Aawflg corresponding to an appropriate Aaw fxn scenario \n');
    return   
end

Aaw = 0.5 * (Pc(Sw(1,:)) +  sum(2*Pc(Sw(2:N,:)),1)  + Pc(Sw(N+1,:)));
Aaw = 10 * dSw .* Aaw' * poro/sigma0;
Aaw = sf'.*Aaw;

end

function validate_interface_models()
% Fast numeric self-test. Run from PowerShell with:
%   $env:PFAS_INTERFACE_VALIDATION_ONLY='true'; matlab -batch ...
    R = 8.314;
    temperature_K = 298.0;
    salt_mM = 10;
    testConcentration = [1e-6; 1e-3; 1; 100];
    compounds = { ...
        struct('name','PFOS','sigma0',71,'a',1/123.794092, ...
            'b',1000*R*temperature_K*2.496483e-6/71, ...
            'pars',[1653, 0.002*1653, 3.0695e-6]); ...
        struct('name','PFOA','sigma0',71.85,'a',1/6.227660735, ...
            'b',1000*R*temperature_K*2.988939298e-6/71.85, ...
            'pars',[97.3971891355371, 0.1927241590803115, ...
            3.4168952485221763e-6])};
    for j = 1:numel(compounds)
        compound = compounds{j};
        lookup = build_kralchevsky_lookup(salt_mM, compound.pars, true, temperature_K);
        [gammaLookup, pressureLookup] = evaluate_kralchevsky_lookup( ...
            testConcentration, lookup);
        gammaDirect = zeros(size(testConcentration));
        pressureDirect = zeros(size(testConcentration));
        for i = 1:numel(testConcentration)
            [gammaMolM2, pressureDirect(i)] = kralchevsky_state_numeric( ...
                testConcentration(i), salt_mM, compound.pars, true, temperature_K);
            gammaDirect(i) = 100*gammaMolM2;
        end
        relativeGammaError = max(abs(gammaLookup-gammaDirect)./max(gammaDirect,1e-30));
        relativePressureError = max(abs(pressureLookup-pressureDirect)./max(pressureDirect,1e-30));
        assert(relativeGammaError < 2e-3, 'Kralchevsky Gamma lookup validation failed.');
        assert(relativePressureError < 2e-3, 'Kralchevsky pressure lookup validation failed.');
        [gammaLangmuir, pressureLangmuir] = interface_state_numeric( ...
            testConcentration, 1, compound.sigma0, compound.a, compound.b, ...
            R, temperature_K, struct());
        assert(all(isfinite(gammaLangmuir)) && all(diff(gammaLangmuir)>0));
        assert(all(isfinite(pressureLangmuir)) && all(diff(pressureLangmuir)>0));
        fprintf('%s validation passed: max lookup errors Gamma=%.3g, pressure=%.3g\n', ...
            compound.name, relativeGammaError, relativePressureError);
    end
end

function lookup = build_kralchevsky_lookup(salt_mM, pars, useActivity, temperature_K)
% Precompute Kralchevsky surface excess and surface pressure on a log grid.
% c is mol/m^3, numerically equal to micro mol/cm^3. Gamma is converted
% from mol/m^2 to micro mol/cm^2 by multiplying by 100.
    cGrid = logspace(-12, 3, 700)';
    gammaGrid = zeros(size(cGrid));
    pressureGrid = zeros(size(cGrid));
    for i = 1:numel(cGrid)
        [gammaMolM2, pressureGrid(i)] = kralchevsky_state_numeric( ...
            cGrid(i), salt_mM, pars, useActivity, temperature_K);
        gammaGrid(i) = 100*gammaMolM2;
    end
    if any(~isfinite(gammaGrid)) || any(~isfinite(pressureGrid))
        error('Kralchevsky lookup contains non-finite values.');
    end
    if any(diff(gammaGrid) < -1e-12)
        warning('Kralchevsky surface excess is not monotone on the lookup grid.');
    end
    lookup.cMin = cGrid(1);
    lookup.cMax = cGrid(end);
    lookup.gammaMax = max(gammaGrid);
    lookup.gammaLowSlope = gammaGrid(1)/cGrid(1);
    lookup.pressureLowSlope = pressureGrid(1)/cGrid(1);
    lookup.logC = log10(cGrid);
    lookup.gammaPP = pchip(lookup.logC, gammaGrid);
    lookup.pressurePP = pchip(lookup.logC, pressureGrid);
    % Inverse maps are used only by the optional Gamma-primary solver.
    % unique protects pchip from machine-precision plateaus near saturation.
    [gammaUnique, uniqueIndex] = unique(gammaGrid, 'stable');
    cUnique = cGrid(uniqueIndex);
    pressureUnique = pressureGrid(uniqueIndex);
    if numel(gammaUnique) < 4 || any(diff(gammaUnique) <= 0)
        error('Kralchevsky Gamma lookup is not strictly monotone; cannot invert it.');
    end
    lookup.logcFromGammaPP = pchip(gammaUnique, log(cUnique));
    lookup.pressureFromGammaPP = pchip(gammaUnique, pressureUnique);
    lookup.gammaInverseMin = gammaUnique(1);
    lookup.gammaInverseMax = gammaUnique(end);
end

function [Gamma, surfacePressure] = interface_state_ad(c, isothermFlg, ...
        sigma0, a, b, R, temperature_K, lookup)
% Return ADI-compatible surface excess [micro mol/cm^2] and surface
% pressure [mN/m]. For Kralchevsky, a local linearization of the lookup
% table supplies dGamma/dc and dPi/dc to the transport Jacobian.
    if isothermFlg == 1
        Gamma = 0.1/R/temperature_K*sigma0*b .* c./(a+c);
        surfacePressure = sigma0*b.*log(1+c./a);
    elseif isothermFlg == 2
        [gammaValue, pressureValue, dGammaDc, dPressureDc] = ...
            evaluate_kralchevsky_lookup(c.val, lookup);
        Gamma = gammaValue + dGammaDc.*(c-c.val);
        surfacePressure = pressureValue + dPressureDc.*(c-c.val);
    else
        error('Unknown interfacial-isotherm flag.');
    end
end

function [Gamma, surfacePressure] = interface_state_numeric(c, isothermFlg, ...
        sigma0, a, b, R, temperature_K, lookup)
% Numeric counterpart of interface_state_ad.
    c = max(c, 0);
    if isothermFlg == 1
        Gamma = 0.1/R/temperature_K*sigma0*b .* c./(a+c);
        surfacePressure = sigma0*b.*log(1+c./a);
    elseif isothermFlg == 2
        [Gamma, surfacePressure] = evaluate_kralchevsky_lookup(c, lookup);
    else
        error('Unknown interfacial-isotherm flag.');
    end
end

function [gammaValue, pressureValue, dGammaDc, dPressureDc] = ...
        evaluate_kralchevsky_lookup(cValue, lookup)
% Evaluate pchip tables in log concentration and calculate local slopes.
    originalSize = size(cValue);
    cValue = max(cValue(:), 0);
    gammaValue = zeros(size(cValue));
    pressureValue = zeros(size(cValue));
    dGammaDc = zeros(size(cValue));
    dPressureDc = zeros(size(cValue));

    low = cValue <= lookup.cMin;
    high = cValue >= lookup.cMax;
    middle = ~(low | high);

    gammaValue(low) = lookup.gammaLowSlope*cValue(low);
    pressureValue(low) = lookup.pressureLowSlope*cValue(low);
    dGammaDc(low) = lookup.gammaLowSlope;
    dPressureDc(low) = lookup.pressureLowSlope;

    if any(middle)
        cMid = cValue(middle);
        x = log10(cMid);
        gammaValue(middle) = ppval(lookup.gammaPP, x);
        pressureValue(middle) = ppval(lookup.pressurePP, x);
        dx = 1e-5;
        dGammaDx = (ppval(lookup.gammaPP, x+dx) - ...
            ppval(lookup.gammaPP, x-dx))/(2*dx);
        dPressureDx = (ppval(lookup.pressurePP, x+dx) - ...
            ppval(lookup.pressurePP, x-dx))/(2*dx);
        dGammaDc(middle) = dGammaDx./(log(10)*cMid);
        dPressureDc(middle) = dPressureDx./(log(10)*cMid);
    end

    if any(high)
        xMax = log10(lookup.cMax);
        gammaValue(high) = ppval(lookup.gammaPP, xMax);
        pressureValue(high) = ppval(lookup.pressurePP, xMax);
        % Saturated extrapolation above the calibrated lookup range.
        dGammaDc(high) = 0;
        dPressureDc(high) = 0;
    end

    gammaValue = reshape(gammaValue, originalSize);
    pressureValue = reshape(pressureValue, originalSize);
    dGammaDc = reshape(dGammaDc, originalSize);
    dPressureDc = reshape(dPressureDc, originalSize);
end

function [c, Gamma, surfacePressure] = interface_state_gamma_primary_ad( ...
        gammaLogit, lookup, gammaMax)
% Gamma-primary state: a logit variable guarantees 0 < Gamma < Gamma_max.
    Gamma = gammaMax./(1 + exp(-gammaLogit));
    [cValue, pressureValue, dCdGamma, dPressureDgamma] = ...
        evaluate_kralchevsky_inverse_lookup(Gamma.val, lookup);
    c = cValue + dCdGamma.*(Gamma-Gamma.val);
    surfacePressure = pressureValue + dPressureDgamma.*(Gamma-Gamma.val);
end

function [cValue, pressureValue, dCdGamma, dPressureDgamma] = ...
        evaluate_kralchevsky_inverse_lookup(gammaValue, lookup)
% Monotone inverse of the precomputed Kralchevsky state table.
    originalSize = size(gammaValue);
    gammaValue = max(gammaValue(:), 0);
    cValue = zeros(size(gammaValue));
    pressureValue = zeros(size(gammaValue));
    dCdGamma = zeros(size(gammaValue));
    dPressureDgamma = zeros(size(gammaValue));
    low = gammaValue <= lookup.gammaInverseMin;
    high = gammaValue >= lookup.gammaInverseMax;
    middle = ~(low | high);
    cValue(low) = gammaValue(low)/lookup.gammaLowSlope;
    pressureValue(low) = gammaValue(low)*lookup.pressureLowSlope/lookup.gammaLowSlope;
    dCdGamma(low) = 1/lookup.gammaLowSlope;
    dPressureDgamma(low) = lookup.pressureLowSlope/lookup.gammaLowSlope;
    if any(middle)
        g = gammaValue(middle);
        logc = ppval(lookup.logcFromGammaPP, g);
        cValue(middle) = exp(logc);
        pressureValue(middle) = ppval(lookup.pressureFromGammaPP, g);
        dg = max(lookup.gammaInverseMax*1e-6, 1e-14);
        dLogCdGamma = (ppval(lookup.logcFromGammaPP, g+dg) - ...
            ppval(lookup.logcFromGammaPP, g-dg))/(2*dg);
        dCdGamma(middle) = cValue(middle).*dLogCdGamma;
        dPressureDgamma(middle) = (ppval(lookup.pressureFromGammaPP, g+dg) - ...
            ppval(lookup.pressureFromGammaPP, g-dg))/(2*dg);
    end
    if any(high)
        g = lookup.gammaInverseMax;
        cValue(high) = exp(ppval(lookup.logcFromGammaPP, g));
        pressureValue(high) = ppval(lookup.pressureFromGammaPP, g);
        dCdGamma(high) = 0;
        dPressureDgamma(high) = 0;
    end
    cValue = reshape(max(cValue, 1e-40), originalSize);
    pressureValue = reshape(pressureValue, originalSize);
    dCdGamma = reshape(dCdGamma, originalSize);
    dPressureDgamma = reshape(dPressureDgamma, originalSize);
end

function aqueous = solve_initial_aqueous_concentration(total, theta, Aaw, ...
        rhob, Kf, M, n3, isothermFlg, sigma0, a, b, R, temperature_K, lookup)
% Solve total = aqueous + solid + AWI mass independently in each cell.
    aqueous = zeros(size(total));
    for i = 1:numel(total)
        target = max(total(i), 0);
        if target == 0
            aqueous(i) = 1e-30;
            continue
        end
        residual = @(ci) initial_mass_residual(ci, target, theta(i), ...
            Aaw(i), rhob, Kf, M, n3, isothermFlg, sigma0, a, b, R, ...
            temperature_K, lookup);
        upper = max(target/max(theta(i), 1e-12), 1e-12);
        aqueous(i) = max(fzero(residual, [0, upper]), 1e-30);
    end
end

function residual = initial_mass_residual(c, target, theta, Aaw, rhob, ...
        Kf, M, n3, isothermFlg, sigma0, a, b, R, temperature_K, lookup)
    Gamma = interface_state_numeric(c, isothermFlg, sigma0, a, b, R, ...
        temperature_K, lookup);
    solidMass = rhob*Kf*(c*M).^n3/M;
    residual = theta*c + solidMass + Aaw*Gamma - target;
end

function [gamma, surfacePressure] = kralchevsky_state_numeric( ...
        cRaw, salt_mM, pars, useActivity, temperature_K)
% Equations and units follow make_fig3_surface_tension_2x3_matlab.m.
    if cRaw <= 0
        gamma = 0;
        surfacePressure = 0;
        return
    end
    k1 = pars(1);
    k2 = pars(2);
    gammaM = pars(3);
    F = 96485.0;
    R = 8.314;
    epsilon = 8.85e-12*78.4;

    cTotal = salt_mM + cRaw;
    if useActivity
        ionicStrengthMolar = cTotal/1000;
        rootI = sqrt(ionicStrengthMolar);
        activity = 10^(-0.5*(rootI/(1+rootI) - 0.3*ionicStrengthMolar));
    else
        activity = 1;
    end
    cSurfactant = activity*cRaw;
    cElectrolyte = activity*cTotal;
    chargeBalance = @(potential) stern_charge_balance(potential, k1, k2, ...
        gammaM, cSurfactant, cElectrolyte, F, R, temperature_K, epsilon);
    potential = solve_stern_potential(chargeBalance);
    expTerm = exp(F*potential/(R*temperature_K));
    adsorption = k1*cSurfactant*expTerm + k2*cElectrolyte*cSurfactant;
    gamma = gammaM*adsorption/(1+adsorption);
    halfPotential = F*potential/(2*R*temperature_K);
    surfacePressure = 1000*(-R*temperature_K*gammaM*log(1-gamma/gammaM) + ...
        4*R*temperature_K*sqrt(2*epsilon*R*temperature_K*cElectrolyte)/F* ...
        (cosh(halfPotential)-1));
end

function potential = solve_stern_potential(chargeBalance)
% Preserve the source model's [-0.3, 0] V bracket and expand only if needed.
    bracket = [-0.3, 0];
    fBracket = [chargeBalance(bracket(1)), chargeBalance(bracket(2))];
    if all(isfinite(fBracket)) && fBracket(1)*fBracket(2) <= 0
        potential = fzero(chargeBalance, bracket);
        return
    end
    candidates = linspace(-1, 0, 401);
    values = arrayfun(chargeBalance, candidates);
    crossing = find(isfinite(values(1:end-1)) & isfinite(values(2:end)) & ...
        values(1:end-1).*values(2:end) <= 0, 1, 'last');
    if isempty(crossing)
        error('Unable to bracket the Kralchevsky Stern potential.');
    end
    potential = fzero(chargeBalance, candidates(crossing:crossing+1));
end

function residual = stern_charge_balance(potential, k1, k2, gammaM, ...
        cSurfactant, cElectrolyte, F, R, temperature_K, epsilon)
    expTerm = exp(F*potential/(R*temperature_K));
    adsorption = k1*cSurfactant*expTerm + k2*cElectrolyte*cSurfactant;
    gamma = gammaM*adsorption/(1+adsorption);
    counterionExcess = gamma*k2*cElectrolyte/(k1*expTerm+k2*cElectrolyte);
    residual = F*(-gamma+counterionExcess) - ...
        sqrt(8*epsilon*R*temperature_K*cElectrolyte)* ...
        sinh(F*potential/(2*R*temperature_K));
end

function DF = groundwater_DF(infiltration_cm_yr, receptor)
% Smith et al. (2024) finite-width source dilution/attenuation factor.
    mixing_depth = sqrt(2*receptor.alpha_v_cm*receptor.W_cm) + ...
        receptor.bsat_cm*(1-exp(-infiltration_cm_yr*receptor.W_cm/ ...
        (receptor.Ugw_cm_yr*receptor.bsat_cm)));
    mixing_depth = min(mixing_depth,receptor.bsat_cm);
    DF = 1 + receptor.Ugw_cm_yr*mixing_depth/ ...
        (infiltration_cm_yr*receptor.W_cm);
end

function [update, nearSingular] = scaled_newton_update(jacobian, residual)
% Solve J*update=-residual after max-norm row/column equilibration.  The
% diagonal transformations are algebraically reversible and therefore do
% not regularize or otherwise alter the Newton equation. A near-singular
% linear system is deliberately rejected so the outer adaptive integrator
% can restore the last accepted state and retry with a smaller time step.
    nearSingular = false;
    n = size(jacobian, 1);
    rowMagnitude = full(max(abs(jacobian), [], 2));
    rowMagnitude(~isfinite(rowMagnitude) | rowMagnitude == 0) = 1;
    rowScale = 1./rowMagnitude;
    rowScaledJacobian = spdiags(rowScale, 0, n, n)*jacobian;

    columnMagnitude = full(max(abs(rowScaledJacobian), [], 1))';
    columnMagnitude(~isfinite(columnMagnitude) | columnMagnitude == 0) = 1;
    columnScale = 1./columnMagnitude;
    scaledJacobian = rowScaledJacobian*spdiags(columnScale, 0, n, n);
    scaledRightHandSide = -rowScale.*residual;

    nearWarningState = warning('query', 'MATLAB:nearlySingularMatrix');
    singularWarningState = warning('query', 'MATLAB:singularMatrix');
    warning('error', 'MATLAB:nearlySingularMatrix');
    warning('error', 'MATLAB:singularMatrix');
    restoreWarnings = onCleanup(@() restore_linear_warning_states( ...
        nearWarningState, singularWarningState));
    try
        scaledUpdate = scaledJacobian\scaledRightHandSide;
    catch ME
        if strcmp(ME.identifier, 'MATLAB:nearlySingularMatrix') || ...
                strcmp(ME.identifier, 'MATLAB:singularMatrix')
            nearSingular = true;
            scaledUpdate = zeros(size(scaledRightHandSide));
        else
            rethrow(ME);
        end
    end
    update = columnScale.*scaledUpdate;
end

function restore_linear_warning_states(nearWarningState, singularWarningState)
    warning(nearWarningState.state, nearWarningState.identifier);
    warning(singularWarningState.state, singularWarningState.identifier);
end
