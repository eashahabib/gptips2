function gp = gpdemo1_config(gp)
%GPDEMO1_CONFIG Config file demonstrating simple (naive) symbolic regression.
%
%   The simple quartic polynomial (y=x+x^2+x^3+x^4) from John Koza's 1992
%   Genetic Programming book is used. It is very easy to solve.
%
%   GP = GPDEMO1_CONFIG(GP) returns the user specified parameter structure
%   GP for the quartic polynomial problem.
%   
%   Example:
%
%   GP = GPTIPS(@GPDEMO1_CONFIG) performs a GPTIPS run using this
%   configuration file and returns the results in a structure called GP.
%
%   Copyright (c) 2009-2015 Dominic Searson 
%
%   GPTIPS 2
% 
%   See also QUARTIC_FITFUN, GPDEMO1

%run control
gp.runcontrol.pop_size = 300;			
gp.runcontrol.num_gen = 100;
gp.runcontrol.runs = 1;
gp.runcontrol.verbose = 25;  
gp.runcontrol.timeout = 10*60*60;            %seconds
gp.fitness.terminate = true;
gp.fitness.terminate_value = 0.1;

gp.runcontrol.parallel.auto = true;

%selection
gp.selection.tournament.size = 20;
gp.selection.tournament.p_pareto = 0.3;
gp.selection.elite_fraction = 0.2;

%fitness function
gp.fitness.fitfun = @quartic_fitfun; 

%Salustowicz
x = [0.05:0.1:9.95]';
y = exp(-x) .* x.^3 .* cos(x) .* sin(x) .* (cos(x) .* sin(x).^2 - 1 );
gp.userdata.ytrain = y;
gp.userdata.xtrain = x;

%test grid
% [x1, x2] = meshgrid(-5:0.2:5, -5:0.2:5);
% x1 = x1(:); x2 = x2(:);
% x = meshgridequal(-1:0.4:1, 4);
% x1 = x(:,1); x2 = x(:,2); x3 = x(:,3); x4 = x(:,4);
% y = 4.*x1.*x2.*cos((x3)) + exp(sin(x4));
% gp.userdata.ytest = y;
% gp.userdata.xtest = [x1 x2 x3 x4];

% x=linspace(-1,1,20)'; 
% gp.userdata.x = x;
% gp.userdata.y = 1./x + 22*x.^2 + x.^3 + x.^4; 
% gp.userdata.name = 'Quartic Polynomial';

% y = 1./(1+x1.^-4) + 1./(1+x2.^-4);
% y = 8./(1 + x1.^2 + x2.^1);
% y = 22 - 4.2.*(cos(x1) - tan(x2)).*(tanh(x3)./sin(x4));

%input configuration 
gp.nodes.inputs.num_inp = size(gp.userdata.xtrain,2); 		         

%quartic example doesn't need constants
gp.nodes.const.num_dec_places = 3;
gp.nodes.const.p_ERC = 0.3;		
gp.nodes.const.p_int = 0.1;

%maximum depth of trees 
gp.treedef.max_depth = 10;

%maximum depth of sub-trees created by mutation operator
gp.treedef.max_mutate_depth = 2;

gp.treedef.max_nodes = inf;
gp.operators.mutation.p_mutate = 0.40;
gp.operators.crossover.p_cross = 0.40;
gp.operators.directrepro.p_direct = 0.20;

%genes
gp.genes.multigene = false;

%define function nodes
gp.nodes.functions.name = {'times','minus','plus','rdivide', 'sin', 'cos', 'tan', 'exp'}; % 'sin','cos','exp', 'power'