function gp = gp_userfcn_hybrid(gp)
%% NEWTON'S METHOD with complex and multicomplex numbers
%GP_USERFCN Calls a user defined function once per generation if one has been
%specified in the field GP.USERDATA.USER_FCN.
%

%% initial equations
C = gp.pop;

string2Beval = string(C); %string to be evaluated

pat = 'x(\d+)';
string2Beval = regexprep(string2Beval,pat,'x(:,$1)'); %$1 makes x1 and x2 things, multiple x variables

%% arrays to store coefficients
coeffies = zeros(7,1); % empty array to store the coefficients
coeff_new2 = coeffies;

%% Initial parameters
h = 1e-5; %multi-complex step size
fitness = gp.fitness.values;
xtrain = gp.userdata.xtrain;
ytrain = gp.userdata.ytrain;

n = length(ytrain);

idx_new=[0]; % array to hold the indices of changed eqns

%% array of multicomplex numbers upto i7
num2improv = 10;
multi_eqn(num2improv) = string(" ");
for i=1:num2improv
    z_temp = zeros(1,2^i);
    z_temp(2^(i-1)+1) = h;
    z(i) = multi(z_temp);
end

%% indices to be improved
frac = gp.selection.elite_fraction+0.2; % 20 percent of the remaining population is also included
if gp.state.count == 1
    return
else
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
end
[~, ia, ~] = unique(string2Beval(idx_chosen));
idx_chosen = sort(idx_chosen(ia));

%%
%run through all the possible equations to find the ones with constants
for j= idx_chosen %1:length(string2Beval)
    
    char_temp = convertStringsToChars((string2Beval(j)));
    str_idx = strfind(char_temp, '['); %start index positions of constants in the function
    
    %if equation has constants, this loop is initiated
    if ~isempty(str_idx) %&& ( ~isinf(fitness(j)) && ~isnan(fitness(j)) )
        
        m = min(length(str_idx),num2improv); %only first 7 coefficients are improved
        
        idx_new(end+1) = j; %contains the index of all the changed constants
        
        end_idx = strfind(char_temp, ']'); %end index positions of constants in the function
        
        %         new_eqn = string(char_temp(1:end_idx(1)-1)); %new multicomplex eqn for 1st order initiated
        
        
        %% evaluating original string
        %regexprep(C{j},exp,'x(:,$1)');
        sum1 = zeros(m, 1);
        sum2 = sum1;
        sum3 = sum1;
        
        %for each coefficient
        for k = 1:m
            %%
            coeffies(k) = str2double(char_temp(str_idx(k)+1:end_idx(k)-1)); %extracting coefficients
            
            multi_eqn(k) = string(char_temp(1:end_idx(k)-1)) + "+ z(1) + z(2) " + string(char_temp(end_idx(k):end));
            
            %array to hold the multicomplex equations for 2nd derivative
        end
        
        %%
        
        temp = cellstr(convertStringsToChars(multi_eqn));
        evalstr_2nd = tree2evalstr(temp,gp);
        
        %solving the eqn for each x point
        for i=1:n
            %%
            x = xtrain(i,:); y=ytrain(i);
            %             eval(['out=' evalstr{1} ';'])
            %eval(['out2(' num2str(k) ')=' evalstr_2nd{k} ';']);
            
            %             if any(isnan(out.zn))% || any(isnan(out2.zn))
            %                 break;
            %             end
            
            %for each constant, adding their relevant contrib.
            for k = 1:m
                eval(['out_multi=' evalstr_2nd{k} ';']);
                out=out_multi.zn;
                deriv1 = out(2)/h;
                deriv2 = out(end)/h^2;
                
                sum1(k) = sum1(k) + (out(1)-y)*deriv1 ;
                sum2(k) = sum2(k) + deriv1^2 + (out(1)-y)*deriv2;
                sum3(k) = sum3(k) + (out(1)-y)^2;
            end
            
        end
        
        add = sum1.*sum3 ./ (sum3.*sum2 - sum1.^2);
        if any(isnan(add))
            add(isnan(add)) = 0;
        end
        if any(isinf(add))
            add(isinf(add)) = 0;
        end
        
        coeff_new2(1:m) = coeffies(1:m) - add;
        
        %% %updating the coefficient
        C{j} = regexprep(C{j}, compose('%.4f',coeffies(1:m)), compose('%.4f',coeff_new2(1:m)));
        
    end
    
end

%evaluate fitness with the new coeffs
% keep new constants only for those whose  fitness_new < fitness_old

idx_new = idx_new(2:end);

fitness_new = fitness;
for i = idx_new %1:gp2.runcontrol.pop_size
    
    %preprocess cell array of string expressions into a form that
    %Matlab can evaluate
    evalstr = tree2evalstr(C{i},gp);
    
    
    [fitness_new(i),gp] = feval(gp.fitness.fitfun,evalstr,gp);
    %gp2.fitness.values(i) = fitness;
    
end

%sum(gp2.fitness.values < gp.fitness.values)

temp_comparison = fitness_new < fitness;
gp.pop(temp_comparison) = C(temp_comparison);
gp.fitness.values(temp_comparison) = fitness(temp_comparison);

gp.improv(gp.state.count) = ( sum(temp_comparison) )/length(idx_new);

%gp.improv(gp.state.count)




% if ~isempty(gp.userdata.user_fcn)
%     [~,gp] = feval(gp.userdata.user_fcn,gp);
% end

%% notes for presentation
% tree depth reduction helps speed up runs but decrease in fitness is
% lower though, between generations
%check for duplicate equations?