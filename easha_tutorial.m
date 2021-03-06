%clc; clear; 
clc;
%close all;
% 
% %gp = rungp(@spatial_evol_config);
%figure; 
for i=1:5
    i
    
%     gp = rungp(@gpTEST1_config);
%     eval(['GPTIPS1_gp' num2str(i) '= gp;']);
%     gppretty(gp, 'best')
%     
    
%     gp = rungp(@gpTEST2_1_config);
%     eval(['JINAY1_gp_newton' num2str(i) '= gp;']);
%     gppretty(gp, 'best')
%     
%     gp = rungp(@gpTEST2_2_config);
%     eval(['JINAY2_gp_newton' num2str(i) '= gp;']);
%     gppretty(gp, 'best')
%     
    gp = rungp(@gpTEST3_1_config);
    eval(['ARMANI1_gp_newton' num2str(i) '= gp;']);
    gppretty(gp, 'best')
    
    figure(1); 
plot(1:gp.state.count, gp.results.history.bestfitness); hold on;
    
    gp = rungp(@gpTEST3_2_config);
    eval(['ARMANI2_gp_newton' num2str(i) '= gp;']);
    gppretty(gp, 'best')
    
    figure(2); 
plot(1:gp.state.count, gp.results.history.bestfitness); hold on;
    
    gp = rungp(@gpTEST3_3_config);
    eval(['ARMANI3_gp_newton' num2str(i) '= gp;']);
    gppretty(gp, 'best')
    
    figure(3); 
plot(1:gp.state.count, gp.results.history.bestfitness); hold on;
    
    gp = rungp(@gpTEST3_4_config);
    eval(['ARMANI4_gp_newton' num2str(i) '= gp;']);
    gppretty(gp, 'best')
    
    figure(4); 
plot(1:gp.state.count, gp.results.history.bestfitness); hold on;
    
    gp = rungp(@gpTEST3_5_config);
    eval(['ARMANI5_gp_newton' num2str(i) '= gp;']);
    gppretty(gp, 'best')
    
    figure(5); 
plot(1:gp.state.count, gp.results.history.bestfitness); hold on;

    
%     A = gp.results.history;
%     B = gp.state;
%      plot(1:gp.state.count, gp.results.history.bestfitness','-', 'LineWidth', 1.5); 
%      hold on;
%      myStats2.bestfitness(1:length(A.bestfitness),i) = A.bestfitness(1:end);
%      myStats.FITiters(i) = length(A.bestfitness);
%      myStats.FITtime(i) = B.runTimeElapsed;
%      myStats.FITfitness(i) = gp.results.best.fitness;
     
%     
%     myStats_pre.FITiters(i) = length(A.bestfitness);
%     myStats_pre.FITtime(i) = B.runTimeElapsed;

figure(1); 
plot(1:gp.state.count, gp.results.history.bestfitness); hold on;
end

save('filewithresultsTEST');

% grid minor;
% legend;
% xlabel('generation'); ylabel('fitness')
% title('Model 2 - Original')
% 
% popbrowser(gp, 'train')
% 
% runtree(gp, 'best')

% % shows the equation in a readbale way
% gppretty(gp, 'best')
% 
% 
% 
% % modelStruct = gpmodel2struct(gp, 'best')
% % 
% % modelStruct.train
% 
% %visualing
% % modelSym = gpmodel2sym(gp, 'best')
% % vpa(modelSym, 2)
% %modelSym = '0.4952*exp(cos(x1))*exp(cos(x2))'
% % modelSym = 'cos(0.3086*x1)'
% for i = 1:10
%     figure;
%     eval(['gp = gp1' num2str(i) ';']);
% eval(['modelSym = gp1' num2str(i) '.results.best.eval_individual{1};']);
% h = ezsurf(modelSym, -5:5, 300);
% h.LineStyle = 'none';
% colormap winter;
% light;
% material shiny;
% h.FaceAlpha = 0.5;
% 
% hold; plot3(gp.userdata.xtest(:,1), gp.userdata.xtest(:,2), gp.userdata.ytest, 'mx');
% 
% legend('GP eqn','Dataset')
% 
% end
% % figure;
% % spatialEvoFunc = sym('1/(1+x1^-4) + 1/(1+x2^-4)') % sym('8/(1 + x1^2 + x2^1)')
% % modelDiff = spatialEvoFunc - modelSym;
% % h2 = ezsurf(modelDiff, -5:5, 300);
% % h2.LineStyle = 'none';
% % colormap winter;
% % light;
% % material shiny;
