%%
%note the coeffs inside a branch is inside the []
clear A idx_new;


gp = gp20;


for j=1:length(gp.pop)
    A{j} = tree2evalstr(gp.pop{j}, gp);
end























































C = gp.pop;
string2Beval = string(C); %string to be evaluated
exp = 'x(\d+)';
string2Beval = regexprep(string2Beval,exp,'x(:,$1)'); %$1 makes x1 and x2 things, multiple x variables

coeffies = zeros(length(string2Beval), 10); % empty array to store the coefficients
coeff_new2 = coeffies;
FD_re = string2Beval; %string to be evaluated for FD

h = 1e-5; %complex step size
fitness = gp.fitness.values;
a = min( [min(fitness),0.2] )/2; %step length
xtrain = gp.userdata.xtrain;
y1 = gp.userdata.ytrain;

n = length(y1);

idx_new=[0];

%indices of top 10% of the population 
[~, idx_minFit] = mink(fitness, floor(gp.selection.elite_fraction*gp.runcontrol.pop_size));
%note "mink" ignores nan values

%run through all the possible equations to find the ones with constants
for j= 63%length(string2Beval) %idx_minFit' %1:length(string2Beval) 
    
    char_temp = convertStringsToChars((string2Beval(j)));
    str_idx = strfind(char_temp, '['); %start index positions of constants in the function
    
    %if equation has constants, this loop is initiated
    if ~isempty(str_idx) %&& ( ~isinf(fitness(j)) && ~isnan(fitness(j)) )
        
        m = min(length(str_idx),7); %only first 10 coefficients are improved
        
        idx_new(end+1) = j; %contains the index of all the changed constants
          
        end_idx = strfind(char_temp, ']'); %end index positions of constants in the function
        %breakdown = ""; % will contain everything except constants
        %breakdown(1) = string(char_temp(1:str_idx(1)));
        
        new_eqn = string(char_temp(1:end_idx(1)-1));
        
        
        %for each coefficient
        for k = 1:m-1
            %%
            %breakdown(k+1) = string(char_temp(end_idx(k):str_idx(k+1)));
            coeffies(j,k) = str2double(char_temp(str_idx(k)+1:end_idx(k)-1)); %extracting coefficients
            
            %creating the respective multicomplex array for the constant
            eval(['z' num2str(k) '(' num2str(2^k) ') = 0 ;']);
            eval(['z' num2str(k) '(' num2str(2^(k-1)+1) ')= h ;']);
            
            %creating equation with each constants c_j is replaced with c_j+i_j
            new_eqn = new_eqn + "+multi(z" + num2str(k) + ")" + string(char_temp(end_idx(k):end_idx(k+1)-1));
        end
        
        %adding multi to last constant along with the remaining string
        k=m;
        coeffies(j,k) = str2double(char_temp(str_idx(k)+1:end_idx(k)-1));
        eval(['z' num2str(k) '(' num2str(2^k) ') = 0 ;']);
        eval(['z' num2str(k) '(' num2str(2^(k-1)+1) ')= h ;']);
        new_eqn = new_eqn+ "+multi(z" + num2str(k) + ")" + string(char_temp(end_idx(k):end));
        
        %breakdown(m+1) = string(char_temp(end_idx(m):end ));
        
        %regexprep(C{j},exp,'x(:,$1)');
        evalstr_orig = tree2evalstr(cellstr(string2Beval{j}),gp);
        
        x = xtrain;
        eval(['out_orig=' evalstr_orig{1} ';']);
        
        if length(out_orig)~=length(x)
            out_orig = ones(1,length(x))*out_orig;
        end
        
        temp = cellstr(convertStringsToChars(new_eqn));
        evalstr = tree2evalstr(temp,gp);
        
        p_num = zeros(1, m);
        p_den = p_num;
        deriv2 = p_num;
        
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
                deriv2(k) = out.zn(2^(k-1)+1)/h;
                p_num(k) = p_num(k) + (out_orig(i)-y)*deriv2(k);
                p_den(k) = p_den(k) + (out_orig(i)-y)^2;
            end
            
        end
        
        p_den = (p_den./n).^0.5;
        
        %a = (0.5 * abs(-fitness(j)./deriv2)) ;
        
        coeff_new2(j, 1:m) = coeffies(j, 1:m) - a/n .* p_num./p_den;
       
        %% %updating the coefficient
        C{j} = regexprep(C{j}, compose('%.4f',coeffies(j, 1:m)), compose('%.4f',coeff_new2(j, 1:m)));
        
    end
    
end

%evaluate fitness with the new coeffs
% keep new constants only for those whose  fitness_new < fitness_old

idx_new = idx_new(2:end);

% B = convertStringsToChars(string2Beval);
% % B2 = convertStringsToChars(FD_re);
% for j=1:length(B)
%     C{j} = cellstr(B{j});
% %     C2{j} = cellstr(B2{j});
% end


gp2 = gp;

gp2.pop = C;

gp2.state.run_completed = false;
for i = idx_new %1:gp2.runcontrol.pop_size
    
        %preprocess cell array of string expressions into a form that
        %Matlab can evaluate
        evalstr = tree2evalstr(gp2.pop{i},gp2);
        
        
        [fitness,gp2] = feval(gp2.fitness.fitfun,evalstr,gp2);
        gp2.fitness.values(i) = fitness;
        
end
gp2.state.run_completed = true;

sum(gp2.fitness.values < gp.fitness.values)
length(idx_new)
min(gp.fitness.values)
min(gp2.fitness.values)

%find(gp2.fitness.values < gp.fitness.values)

% disp(['          old constant     new constants'])
%    for h = 1:length(gp.pop)
%       disp( [gp.fitness.values(h), gp2.fitness.values(h)]);
%    end
