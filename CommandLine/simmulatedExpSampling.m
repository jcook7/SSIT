%% Using the SSIT to fit and Design Single-Cell Experiments 
% In this script, we show how the SSIT can be used to identify a
% time-inhomogeneous model for the activation of Dusp1 mRNA expression
% under Dexamethasome stimulation of Glucocorticoid Receptors.
clear 
close all
clc
addpath('../CommandLine');
addpath('../EricModel/DUSP1_GR_dataframes');
addpath(genpath('../src'));
addpath('../tensor_toolbox-v3.2.1')

%% Define SSIT Model
% Initial Set up for the model
ModelSetUp = SSIT; % Rename SSIT code to make changes with out changing the original code
ModelSetUp.species = {'x1';'x2'}; % x1: number of on alleles,  x2: mRNA 
ModelSetUp.initialCondition = [0;0]; % Set initial conditions
ModelSetUp.propensityFunctions = {'kon*IGR*(2-x1)';'koff*x1';'kr*x1';'gr*x2'}; % Set the propensity functions of the 4 reactions

% Associated stoichiometry of the four reactions
stoich = [1,0;   % Gene allele switching to on state
         -1,0;   % Gene allele switching to off state
          0,1;   % mRNA production
          0,-1]; % mRNA degredation

ModelSetUp.stoichiometry = transpose(stoich);
% Input expression for time varying signals
ModelSetUp.inputExpressions = {'IGR','kcn0/knc+(t>=0)*kcn1/(r1-knc)*(exp(-knc*t)-exp(-r1*t))'};

% Defining the values of each parameter used
ModelSetUp.parameters = ({'koff',0.14;'kon',0.01;'kr',1;'gr',0.01;...
                   'kcn0',0.01;'kcn1',0.1;'knc',1;'r1',0.01});

ModelSetUp.fspOptions.initApproxSS = true;
timeMatrix = [0 10 20 30 40 50 60 75 90 120 150 180];
ModelSetUp.tSpan = timeMatrix;
NCells = [100 0 0 0 0 0 0 0 100 0 0 100];

for j = 1:5
    Model{j} = ModelSetUp; % Start up model experiment iteration with the same conditions 
    %% Run simulated experiment sampling times
    dataFileName =  'DUSP1_3hr_Dex_100nM_total.csv'; 
    [simData,csvFile] = sampleExperimentSim(dataFileName,timeMatrix,NCells);
    
    %% Solve the model using the FSP
    Model{j}.fittingOptions.modelVarsToFit = 1:8; % Picking parameters 1-8 for fitting
    
    % Priors for fitting
    muLog10Prior = [-1 -1 1 -2 -1 -2 -1 -1];
    sigLog10Prior = [2 2 1 1 2 2 2,2];
    
    muLog10Prior = muLog10Prior(Model{j}.fittingOptions.modelVarsToFit);
    sigLog10Prior = sigLog10Prior(Model{j}.fittingOptions.modelVarsToFit);
    Model{j}.fittingOptions.logPrior = @(p)-(log10(p(:))-muLog10Prior').^2./(2*sigLog10Prior'.^2);
    
    Model{j} = Model{j}.loadData(csvFile,{'x2','RNA_nuc'}); % Load experimental data set
    Model{j}.tSpan = timeMatrix;
    Model{j}.dataSet.nCells = NCells';
    
    for i=1:3
        Model{j}.solutionScheme = 'FSP';
        Model{j}.fspOptions.fspTol = 1e-4;
        Model{j}.fspOptions.verbose = 0;
        Model{j}.fspOptions.bounds=[];
        [fspSoln,Model{j}.fspOptions.bounds] = Model{j}.solve;
        Model{j}.fspOptions.bounds
        
        % Load and Fit smFISH Data
        Model{j}.fspOptions.fspTol = inf;
        fitOptions = optimset('Display','iter','MaxIter',500); 
        Model{j}.parameters(Model{j}.fittingOptions.modelVarsToFit,2) =...
          num2cell(Model{j}.maximizeLikelihood(...
          [Model{j}.parameters{Model{j}.fittingOptions.modelVarsToFit,2}],...
          fitOptions));
    end
    Model{j}.makeFitPlot
    %% Metropolis Hastings to Quantify Parameter Uncertainty
    Model{j}.fittingOptions.modelVarsToFit = 1:4; % Picking parameters 1-4 for fitting
    muLog10Prior = [-1 -1 1 -2];
    sigLog10Prior = [2 2 1 1];
    
    muLog10Prior = muLog10Prior(Model{j}.fittingOptions.modelVarsToFit);
    sigLog10Prior = sigLog10Prior(Model{j}.fittingOptions.modelVarsToFit);
    Model{j}.fittingOptions.logPrior = @(p)-(log10(p(:))-muLog10Prior').^2./(2*sigLog10Prior'.^2);
    % Model{j}.fittingOptions.logPrior = []; % Remove prior before fitting
    
    allFitOptions.CovFIMscale = 0.1;% make smaller for higher acceptance rate
    MHOptions = struct('numberOfSamples',5000,'burnin',100,'thin',1,...
      'useFIMforMetHast',true,'suppressFSPExpansion',true);
    [bestParsFound,~,mhResults] = Model{j}.maximizeLikelihood([Model{j}.parameters{Model{j}.fittingOptions.modelVarsToFit,2}]',...
      MHOptions,'MetropolisHastings');
    Model{j}.parameters(Model{j}.fittingOptions.modelVarsToFit,2) = num2cell(bestParsFound);
    
    %% Pick Samples
    Model{j}.tSpan = timeMatrix; % Define the time points
    mhSamplePar = exp(mhResults.mhSamples);
    pickN = 5; % Number of samples to pick
    halfStack = length(mhSamplePar(:,1))/2;
    k1 = randperm(halfStack,pickN);
    k2 = halfStack + k1; % Picks random samples from second half of stack
    par_pick = mhSamplePar(k2,:);
    fims_sim = zeros(pickN,length(timeMatrix),length(Model{j}.parameters),length(Model{j}.parameters));
    


    %% MH samples
    J = floor(linspace(500,1000,10));
    MHSamplesForFIM = exp(mhResults.mhSamples(J,:));
    fimResults = Model{j}.computeFIM([],'lin',MHSamplesForFIM);
    
    % Find and plot total FIM for each parameter sample
    Nc = Model{j}.dataSet.nCells;
    
    % Find optimal experiment design given parameters sets
    NcOptExperiment = Model{j}.optimizeCellCounts(fimResults,300,'Determinant',[],NCells);
%     FIMOptExpt = Model{j}.totalFim(fimResults,NcOptExperiment);
%     Model{j}.plotMHResults(MHResults,[FIM,FIMOptExpt]) 

    NCells = NCells + NcOptExperiment; % New cell measurement design for next experiment 
        
end
    

% %% det(FIM) of each sample
% for n = 1:pickN
%     for m = 1:length(timeMatrix)
%         fimSample(n,m) = det(squeeze(fims_sim(n,m,:,:)));
%     end
% end
% 
% figure()
% time=[];
% fimPlot=[];
% 
% for z = 1:pickN
%     time = [time timeMatrix];
%     fimPlot = [fimPlot fimSample(z,:)];
% end
% 
% plot(time,fimPlot,'--b',timeMatrix,mean(fimSample,1),'k','LineWidth',2)
% legend('Samples','Mean')
% xlabel('time [min]')
% ylabel('|FIM|')
% 
% %% Plot FIM with standard Deviation
% std_dev = std(fimSample,1);
% curve1 = mean(fimSample,1) + std_dev;
% curve2 = mean(fimSample,1) - std_dev;
% x2 = [timeMatrix, fliplr(timeMatrix)];
% inBetween = [curve1, fliplr(curve2)];
% figure()
% fill(x2, inBetween, 'g');
% hold on;
% plot(timeMatrix,mean(fimSample,1),'r','LineWidth',3)
% legend('Standard Dev','Average |FIM|')
% xlabel('time [min]')
% ylabel('|FIM|')
% hold off
