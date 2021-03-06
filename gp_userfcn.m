function gp = gp_userfcn(gp)
%% GRADIENT DESCENT METHOD with multicomplex method
%GP_USERFCN Calls a user defined function once per generation if one has been 
%specified in the field GP.USERDATA.USER_FCN.
% 
%   Remarks:
%
%   The user defined function must accept (as 1st) and return (as 2nd) the
%   GP structure as arguments.
%
%   Example:
%
%   In multigene symbolic regression the function
%   regressmulti_fitfun_validate can be called once per generation to
%   evaluate the best individual so far on a validation ('holdout') data
%   set.
%
%   Copyright (c) 2009-2015 Dominic Searson 
%
%   GPTIPS 2
% 
% C = gp.fitness.returnvalues{gp.state.current_individual}; %current coeffs
% evalstr = tree2evalstr(gp.pop{250},gp);
% 
% y = gp.userdata.ytrain;
% numData = gp.userdata.numytrain; % num of data
% numGenes = numel(evalstr);
% 
% pat = 'x(\d+)';
% evalstr = regexprep(evalstr,pat,'gp.userdata.xtrain(:,$1)+0.000001i')
% 
% 
% %eval each gene in the current individual
% for i = 1:numGenes
%     ind = i + 1;
%     eval(['geneOutputs(:,ind)=' evalstr{i} ';']);
%     geneOutputs(:,ind) = imag(geneOutputs(:,ind))/0.000001;
%     
%     %check for nonsensical answers and break out early with an 'inf' if so
%     if  any(~isfinite(geneOutputs(:,ind))) || any(~isreal(geneOutputs(:,ind)))
%         fitness = Inf;
%         gp.fitness.returnvalues{gp.state.current_individual} = [];
%         return
%     end
% end
% 
% geneOutputs
% 
% %only calc. weighting coeffs during an actual run or if forced
% 
%     %set gp.userdata.bootSample to true to resample data for weights computation
%     
%     %prepare LS matrix
%     
%         goptrans = geneOutputs';
%         prj = goptrans * geneOutputs;
%     
%     %calculate tree weight coeffs using SVD based least squares
%     %normal equation
%     
%             theta = pinv(prj) * goptrans * y;
%         
%         %theta = [];
%         %fitness = Inf;
%         %gp.fitness.returnvalues{gp.state.current_individual} = [];
%         %return;
%     
%  


C = gp.pop;
string2Beval = string(C); %string to be evaluated
pat = 'x(\d+)';
string2Beval = regexprep(string2Beval,pat,'x(:,$1)'); %$1 makes x1 and x2 things, multiple x variables

coeffies = zeros(7, 1); % empty array to store the coefficients
coeff_new2 = coeffies;
%FD_re = string2Beval;

h = 1e-5; %complex step size
fitness = gp.fitness.values;
a = min( [min(fitness),0.2] )/2; %step length
xtrain = gp.userdata.xtrain;
y1 = gp.userdata.ytrain;

fitness_new = fitness;

n = length(y1);

idx_new=[0];

%% array of multicomplex numbers upto i7
num2improv = 6;
for i=1:num2improv
    z_temp = zeros(1,2^i);
    z_temp(2^(i-1)+1) = h;
    z(i) = multi(z_temp);
end

%% indices to be improved
frac = gp.selection.elite_fraction+0.2; % 20 percent of the remaining population is also included
if gp.state.count > 1
    if gp.improv(gp.state.count-1)<0.1
        % random indices,
        idx_chosen = randi([1 gp.runcontrol.pop_size], 1, floor(frac*gp.runcontrol.pop_size));
        % how to not include pop with nan fitness?
        idx_chosen = idx_chosen(~isnan( fitness(idx_chosen) ));
    else
        %indices of elite population only
        [~, idx_chosen] = mink(fitness, floor( frac*gp.runcontrol.pop_size));
        idx_chosen = idx_chosen';
        %note "mink" ignores nan values
    end
else
    [~, idx_chosen] = mink(fitness, floor( frac*gp.runcontrol.pop_size));
    idx_chosen = idx_chosen';
end
[~, ia, ~] = unique(string2Beval(idx_chosen));
idx_chosen = sort(idx_chosen(ia));

%% run through all the possible equations to find the ones with constants
for j= idx_chosen %1:length(string2Beval) 
    
    char_temp = convertStringsToChars((string2Beval(j)));
    str_idx = strfind(char_temp, '['); %start index positions of constants in the function
    
    %if equation has constants, this loop is initiated
    if ~isempty(str_idx) %&& ( ~isinf(fitness(j)) && ~isnan(fitness(j)) )
        
        m = min(length(str_idx),num2improv); %only first 7 coefficients are improved
        
        idx_new(end+1) = j; %contains the index of all the changed constants
          
        end_idx = strfind(char_temp, ']'); %end index positions of constants in the function
        %breakdown = ""; % will contain everything except constants
        %breakdown(1) = string(char_temp(1:str_idx(1)));
        
        new_eqn = string(char_temp(1:end_idx(1)-1));
        
        
        %for each coefficient
        for k = 1:m-1
            %%
            %breakdown(k+1) = string(char_temp(end_idx(k):str_idx(k+1)));
            coeffies(k) = str2double(char_temp(str_idx(k)+1:end_idx(k)-1)); %extracting coefficients
            
            %creating equation with each constants c_j is replaced with c_j+i_j
            new_eqn = new_eqn + "+z(" + num2str(k) + ")" + string(char_temp(end_idx(k):end_idx(k+1)-1));
        end
        
        %adding multi to last constant along with the remaining string
        k=m;
        coeffies(k) = str2double(char_temp(str_idx(k)+1:end_idx(k)-1));
        new_eqn = new_eqn+ "+z(" + num2str(k) + ")" + string(char_temp(end_idx(k):end));
        
        %breakdown(m+1) = string(char_temp(end_idx(m):end ));
        
        %regexprep(C{j},exp,'x(:,$1)');
        evalstr_orig = tree2evalstr(cellstr(string2Beval{j}),gp);
        
        x = xtrain;
        eval(['out_orig=' evalstr_orig{1} ';']);
        
        if length(out_orig)==1
            out_orig = ones(length(x),1)*out_orig;
        end
        
        temp = cellstr(convertStringsToChars(new_eqn));
        evalstr = tree2evalstr(temp,gp);
        
        p_num = zeros(m, 1);
        %p_den = p_num;
        
        %for each pt
        for i=1:n
            %%
            x = xtrain(i,:); y=y1(i);
            eval(['out=' evalstr{1} ';']);
            
            if any(isnan(out.zn))
                break;
            end
            
            %for each constant, adding their relevant contrib.
            for k = 1:m
                deriv1 = out.zn(2^(k-1)+1)/h;
                p_num(k) = p_num(k) + (out_orig(i)-y)*deriv1;
                %p_den(k) = p_den(k) + (out_orig(i)-y)^2;
            end
            
        end
        
        %p_den = (p_den./n).^0.5;
        add = a/(n*fitness(j)) * p_num;
        if any(isnan(add))
            add(isnan(add)) = 0;
        end
        if any(isinf(add))
            add(isinf(add)) = 0;
        end
       
        coeff_new2(1:m) = coeffies(1:m) - add;
        
        %% %updating the coefficient
        C{j} = regexprep(C{j}, compose('%.3f',coeffies(1:m)), compose('%.3f',coeff_new2(1:m)));
        
        %preprocess cell array of string expressions into a form that
        %Matlab can evaluate
        evalstr = tree2evalstr(C{j},gp);
        
        
        [fitness_new(j),gp] = feval(gp.fitness.fitfun,evalstr,gp);
        %gp2.fitness.values(i) = fitness;
    end
    
end

%evaluate fitness with the new coeffs
% keep new constants only for those whose  fitness_new < fitness_old

%idx_new = idx_new(2:end);


%sum(gp2.fitness.values < gp.fitness.values)

temp_comparison = fitness_new < fitness;
gp.pop(temp_comparison) = C(temp_comparison);
gp.fitness.values(temp_comparison) = fitness(temp_comparison);

gp.improv(gp.state.count) = ( sum(temp_comparison) )/(length(idx_new)-1);

%gp.improv(gp.state.count)


% if ~isempty(gp.userdata.user_fcn)
%     [~,gp] = feval(gp.userdata.user_fcn,gp);
% end

%% notes for presentation
% tree depth reduction helps speed up runs but decrease in fitness is
% lower though, between generations