clear all
close all
addpath(genpath('/home/ale/Projects/headModel/'))

%%
hm = headModel.loadFromFile('/home/ale/Projects/rsc/resources/head_modelColin27_2003_xyz_Standard-10-5-Cap346.mat');
hm.plot();
%%
snr = 10;
Nx = size(hm.cortex.vertices,1);
Nt = 100;
x0 = hm.cortex.vertices(unidrnd(Nx, Nt,1),:);

% Open the surface by the Corpus Callosum (unlabeled vertices)
rmIndices = find(hm.atlas.colorTable==0);

% Get SVD decomposition of Kstd/L. 
% Kstd is the LF standardized by source to reduce depth bias.
[Ut, s2,iLV, Kstd,ind,K,L] = hm.svd4sourceLoc({},{rmIndices});

Jall = [];
err = zeros(Nt,1);
I = zeros(Nx,3);
I(ind,:) = 1;
I = find(I(:));
for k=1:Nt
    
    % Generate Gaussian patch centered on a random location
    J = geometricTools.simulateGaussianSource(hm.cortex.vertices,x0(k,:),0.01);
    
    % Simulate xyz components
    J = [J;J;J];
    
    % Simulate EEG
    y = Kstd*J(I);
    
    % Simulate noisy trials to be able to compute a t-score as in sLoreta
    % for ERPs
    y = y*ones(1,100);
    y = awgn(y,snr,'measured');
    
    % Calculate inverse solution
    [tmp,lambdaOpt] = ridgeSVD(y,Ut, s2,iLV,100,0);
    
    % Compute t-score and store
    Jhat = 0*J;
    Jhat(I) = median(tmp,2)./(std(tmp,[],2)+eps);
       
    % Compute localization error
    [~,loc1] = max(sum(reshape(J,[],3).^2,2));
    [~,loc2] = max(sum(reshape(Jhat,[],3).^2,2));
    err(k) = norm(hm.cortex.vertices(loc1,:)-hm.cortex.vertices(loc2,:));
    
    % Collect solutions for visualization
    Jall = [Jall J Jhat];
end
%%
yall = Kstd*Jall(I,:);
hm.plotOnModel(Jall,yall,'',true);
figure;subplot(121);plot(err);subplot(122);hist(err,20)

% 100x to convert to cm
100*[median(err) mean(err) std(err)]

