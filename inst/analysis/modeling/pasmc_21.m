function m = normoxia()
 
% setup
oxygen = '21%';
treatment = 'None';
in_path = '/Users/will/Dropbox (Partners HealthCare)/Copeland.2021.hypoxia.flux/inst/analysis/modeling/matlab-input/';
out_path = '/Users/will/Dropbox (Partners HealthCare)/Copeland.2021.hypoxia.flux/inst/analysis/modeling/';
cell_type = 'pasmc';
addpath(genpath(in_path));
output = [out_path cell_type '_21/' cell_type '_21_' datestr(now,'yyyy-mm-dd-HH-MM-SS') '.mat'];

% biomass equation
biomass = readtable([ cell_type '_biomass.csv']);
biomasseqn = '';
for i=1:height(biomass)-1
    biomasseqn = [biomasseqn num2str(biomass{i, 2}) ' ' char(biomass{i, 1}) ' + '];
end
biomasseqn = [biomasseqn num2str(biomass{height(biomass), 2}) ' ' char(biomass{height(biomass), 1}) ' -> biomass'];
clear biomass;

% reaction equations
reactions = readtable('reactions.csv');
% reactions = rmmissing(reactions);
reaction_list = [biomasseqn; reactions.equation]; 
reaction_ids = ['BIOMASS'; reactions.name];
r = reaction(reaction_list); 
r.id = reaction_ids;
clear biomasseqn reaction_ids reaction_list;

% define tracers
t{1} = tracer({'[1,2-13C2]GLC: GLC.x @ 1 2'});
t{1}.atoms.it = [0.01; 0.99];

t{2} = tracer({'[U-13C6]GLC: GLC.x @ 1 2 3 4 5 6'});
t{2}.atoms.it = [0.01; 0.99];

t{3} = tracer({'[U-13C5]GLN: GLN.x @ 1 2 3 4 5'});
t{3}.atoms.it = [0.006; 0.994];

% define metabolites with mids
frags = {...          
    'PYR: PYR.ms @ 1 2 3'; ...
    'LAC: LAC @ 1 2 3'; ...
    'ALA: ALA @ 1 2 3'; ...
    'AKG: AKG @ 1 2 3 4 5'; ...
    'MAL: MAL @ 1 2 3 4'; ...
    'ASP: ASP @ 1 2 3 4'; ...
    'GLU: GLU @ 1 2 3 4 5'; ...
    'GLN: GLN @ 1 2 3 4 5'; ...
    'CIT: CIT @ 1 2 3 4 5 6'; ...
    'SER: SER @ 1 2 3';...
    '3PG: 3PG @ 1 2 3';...
    'FBP: FBP @ 1 2 3 4 5 6';...
    };

fragments{1} = frags;

fragments{2} = frags;

fragments{3} = frags;

% initialize the msdata object for each set of fragments
for i=1:length(fragments)
    d{i} = msdata(fragments{i});
end
clear frags fragments;

% import mid data
mids = readtable([cell_type '_mids.csv'], 'TreatAsEmpty', {'NA'});
mids.se(isnan(mids.se)) = 0.01;
mids = mids(strcmp(mids.treatment, treatment), :);
mids = mids(strcmp(mids.oxygen, oxygen), :);


n = 1; 
for label = {'glc2' 'glc6' 'q5'}
    a = mids(strcmp(mids.tracer, char(label)), :); 
    for metabolite = unique(mids.metabolite)'
        b = a(strcmp(a.metabolite, char(metabolite)), :);
        y = [];
        h = [];
        for time = unique(b.time)'
            c = b(b.time == time, 'mean');
            c = table2array(c);
            y = [y c];
            h = [h time];
        end
        d{n}{char(metabolite)}.mdvs = mdv(y);
        d{n}{char(metabolite)}.mdvs.time = h;
    end
    n = 1 + n;
end

% % assign sem from mid measurements
% n = 1; 
% for label = {'[1,2-13C2] glucose' '[U-13C6] glucose' '[U-13C5] glutamine'}; 
%     a = mids(strcmp(mids.label, char(label)), :); 
%     for compound = unique(mids.compound)';
%         b = a(strcmp(a.compound, char(compound)), :);
%         y = [];
%         for time = unique(b.time)';
%             c = b(b.time == time, 'se');
%             c = table2array(c);
%             y = [y c];
%         end
%         d{n}{char(compound)}.mdvs.std = y;
%     end
%     n = 1 + n;
% end

clear a b c y n x z time compound label;

% assign sem for mid (can repeat loop above to assign actual sem values)
sem_value = 0.013;
for i=1:length(d)
    for j=1:length(d{i})
        d{i}{d{i}.id{j}}.mdvs.std = sem_value * ones(length(d{i}{d{i}.id{j}}.mdvs.val),1);
    end
end

% import flux data
fluxes = readtable([cell_type '_fluxes.csv']);
fluxes = fluxes(strcmp(fluxes.treatment, treatment), :);
fluxes = fluxes(strcmp(fluxes.oxygen, oxygen), :);

f = data(fluxes.name);
f.val = abs(fluxes.mean)';
f.std = fluxes.se';
% f.std = 0.1 * f.val;

% set fixed rates
fixed_flux_equations = {'sPYR'}';

for i=1:length(fixed_flux_equations)
    for j=1:length(r)
        if strcmp(fixed_flux_equations{i},r.id{j})
            r(j).flx.val = 1; 
            r(j).flx.fix = true;
            break;
        end
    end
end

% assemble data from individual experiments
for i=1:length(d)
    struct(i).tracers = t{i};
    struct(i).data_ms = d{i};
    struct(i).data_flx = f;
    struct(i).id = t{i}.id;
end

e = experiment(struct);

m = model(r,'expts',e);

% % label washout in glc6 experiments after 48 h
% for metabolite = {'AKG' 'GLU' 'CIT' 'ASP' 'MAL'}
%     m.expts{'[U-13C6]GLC'}.data_ms{char(metabolite)}.mdvs.on = [true true false];
% end

% % remove low labeled 3-carbon metabolites
% for metabolite = {'LAC' 'PYR' 'ALA'}
%     m.expts{'[1,2-13C2]GLC'}.data_ms{char(metabolite)}.on = [false];
% end

m.mets{'SUC'}.sym = list('rot',atommap('1:4 2:3 3:2 4:1'));
m.mets{'FUM'}.sym = list('rot',atommap('1:4 2:3 3:2 4:1'));

m.states{'CO2','c'}.bal = false;
for metabolite = {'LAC.x' 'PYR.x' 'ALA.x' 'GLN.x' 'GLU.x' 'ASP.x'}
    m.states{char(metabolite)}.bal = false;
end
for pool = {'PYR.ms'}
    m.states{char(pool)}.val = 0;
    m.states{char(pool)}.fix = true;
end

% options
m.options.cont_alpha = 0.05;
m.options.sim_ss = false;
m.options.sim_sens = false;
m.options.sim_more = false;
m.options.sim_na = false;
m.options.fit_reinit = true;
m.options.fit_starts = 10;
m.options.fit_reltol = 0.001;
m.options.fit_tau = 1e-6;
m.options.fit_nudge = 6;
m.options.int_tspan = [0 12 24 36 48];
m.options.int_maxstep = Inf;
m.options.hpc_on = 1;

% % options to focus on extracellular flux balance
% m.options.sim_ss = true;
% m.expts.on = [true false false];
% 
% for expt = {'[1,2-13C2]GLC'}
%     m.expts{char(expt)}.data_ms.on = false;
% end

save(output, 'm');

return
