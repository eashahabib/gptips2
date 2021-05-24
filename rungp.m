function gpout = rungp(configFile)
%RUNGP Runs GPTIPS 2 using the specified configuration file.
%
%   GP = RUNGP(@CONFIGFILE) runs GPTIPS using the user parameters contained
%   in the configuration file CONFIGFILE.M and returns the results in the
%   GP data structure. Post run analysis commands may then be run on this
%   structure, e.g. POPBROWSER(GP), SUMMARY(GP), RUNTREE(GP,'best') etc.
%
%   Demos:
%
%   Use GPDEMO1 at the commmand line for a demonstration of naive symbolic
%   regression (not multigene).
%
%   For demos involving multigene symbolic regression see GPDEMO2, GPDEMO3
%   and GPDEMO4.
%
%   Addionally, example config files for muligene regression on data
%   generated by several non-linear mathematical functions are
%   CUBIC_CONFIG, RIPPLE_CONFIG, UBALL_CONFIG and SALUSTOWICZ1D_CONFIG.
%
%   Remarks:
%
%   For further information see the GPTIPS websites at:
%
%   https://sites.google.com/site/gptips4matlab/
%
%   Also, see the following paper for a summary of GPTIPS 2 capabilities:
%
%   Searson DP, GPTIPS 2: an open-source software platform for symbolic
%   data mining, Chapter 22 in Gandomi, AH, Alavi, AH, Ryan, C (eds).
%   Springer Handbook of Genetic Programming Applications, Springer, 2015.
%   (author proof freely available at http://arxiv.org/abs/1412.4690)
%   
%   ----------------------------------------------------------------------
%
%   This program is free software: you can redistribute it and/or modify it
%   under the terms of the GNU General Public License as published by the
%   Free Software Foundation, either version 3 of the License, or (at your
%   option) any later version.
%
%   This program is distributed in the hope that it will be useful, but
%   WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   General Public License for more details.
%
%   You should have received a copy of the GNU General Public License along
%   with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%   ----------------------------------------------------------------------
%
%   Copyright (c) 2009-2015 Dominic Searson
%
%   GPTIPS 2
%
%   See also GPDEMO1, GPDEMO2, GPDEMO3, GPDEMO4, CUBIC_CONFIG,
%   UBALL_CONFIG, SALUSTOWICZ1D_CONFIG, RIPPLE_CONFIG

%if no config supplied then report error
if nargin < 1 || isempty(configFile)
    disp('To run GPTIPS a configuration m file function handle must be specified');
    disp('e.g. GP = RUNGP(@CUBIC_CONFIG)');
    return;
end

%generate prepared data structure with some default run parameter values
gp = gpdefaults();

%setup random number generator
gp = gprandom(gp);

%run user configuration file
gp = feval(configFile,gp);

gp.improv = zeros(1, gp.runcontrol.num_gen);

%multiple runs (are merged after each run)
for run=1:gp.runcontrol.runs;
    
    gp.state.run = run;
    
    %run user configuration file each time for subsequent runs, unless
    %suppressed (default)
    if run > 1 && ~gp.runcontrol.suppressConfig;
        gp = feval(configFile,gp);
    end
    
    %attach name of config file for reference
    gp.info.configFile = configFile;
    
    %perform error checks
    gp = gpcheck(gp);
    
    %perform initialisation
    gp = gpinit(gp);
    
    %main generation loop
    for count=1:gp.runcontrol.num_gen
        
        %start gen timer
        gp = gptic(gp);
        
        if count == 1;
            
            %generate the initial population
            gp = initbuild(gp);
            
        else
            
            %use crossover, mutation etc. to generate a new population
            gp = popbuild(gp);
            
        end
        
        %calculate fitnesses of population members
        gp = evalfitness(gp);
        
        %call user defined function 
        gp = gp_userfcn_hybrid(gp); %multicomplex newton
        %gp = gp_userfcn(gp); %multicomplex GD
        %gp = gp_userfcn2(gp); %complex GD
        
        %update run statistics
        gp = updatestats(gp);
      
        %display current stats on screen
        displaystats(gp);
        
        %end gen timer
        gp = gptoc(gp);
        
        %check termination criteria
        gp = gpterminate(gp);
        
        %save gp structure
        if  gp.runcontrol.savefreq && ...
                ~mod(gp.state.count-1,gp.runcontrol.savefreq)
            save gptips_checkpoint gp
        end
        
        %break out of generation loop if termination required
        if gp.state.terminate
            if ~gp.runcontrol.quiet
                disp(['Termination criterion met (' ... 
                    gp.state.terminationReason '). Terminating run.']);
            end
            break;
        end
        
    end  %end generation loop
    
    %finalise the run
    gp = gpfinalise(gp);
    
    %merge independent runs
    if gp.state.run > 1
        gpout = mergegp(gpout,gp,true);
    else
        gpout = gp;
    end
    
end %end independent run loop