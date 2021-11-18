
close all
clear all


pathn='.\'; 
% modeling a 3 layer SNN
Weight = [1 1 1 1 ];
Weight2=[1 1];
SEED=1000;
V_reset = -0.075;
V_e = -0.07; 
V_th = -0.04; 

% model to generate simulated data
ModelSubRF_3layers_LIF_snn(pathn,Weight,Weight2,SEED,V_reset,V_e,V_th)

% next time, no need to run the model 
load('Data.mat');

% get STE (spike triggered stimulus ensemble)
STE=get_STE_snn(SS,spklist,tmSVD,nt,nx,ny);
Nx=nx;Ny=ny;
save('STE.mat','STE','Nx','Ny');

% inferring subunits
k=8;
STNMFanalysis_snn('STE.mat','subunit.mat',k,20,1);
close all

% inferring weights
get_module_weight_snn('STE.mat','subunit.mat','module_weight.mat');

% inferring spikes of each subunit
load('module_weight.mat','W_matrix');
load('subunit.mat','unit','Ngood');
spklist_sub=get_subSTA_min_snn(Ngood,spklist,W_matrix,STE,Nx,Ny);
LGN_subunit_sp=get_subunit_spike_train_snn(Ngood,spike,spklist_sub);
save('subunit_sp.mat','LGN_subunit_sp');

% L3-STA and L1-subunit
figure;
datarange = [-0.5 0.5];
imagesc(spSVD/norm(spSVD),datarange);
colormap(colormap_RF(datarange));
title('L3-1 STA');
colorbar;

%%
%%%% adjusting iglist to manually order the inferred subunits same as modeled subunits 
iglist=[ 2 5 3 4 6 1]; 
%%%%
numModel=[1,4,2,5,3,6];
figure;
fcol=Ngood;
frow=2;
for ii=1:Ngood
    % inferred subunits
    subplot(frow,fcol,ii);
    %spsvd=unit{ii};
    spsvd=unit{iglist(ii)};
	imagesc(spsvd/norm(spsvd),datarange);
	colormap(colormap_RF(datarange));
    % modeled subuints
    subplot(frow,fcol,fcol+ii);
    %spsvd=reshape(spM(ii,:),Nx,Ny);
    spsvd=reshape(spM(numModel(ii),:),Nx,Ny);
	imagesc(spsvd/norm(spsvd),datarange);
	colormap(colormap_RF(datarange));
end
sgtitle('Inferred L1 cells')

% spike train and spike correlation
A=[];
B=[];
for i=1:size(iglist,2)
    B(i,:)=sub_SP(numModel(i),:);
    A(i,:)=LGN_subunit_sp(iglist(i),:);
end
for i=1:size(A,1)
    for j=1:size(B,1)
        R=corrcoef(A(i,:),B(j,:));
        CC(i,j)=R(1,2);
        R=corrcoef(B(j,:),spike);
        CC_spike2(j)=R(1,2);
    end
     R=corrcoef(A(i,:),spike);
     CC_spike1(i)=R(1,2);
end


figure;
imagesc(CC);
colorbar;
colormap(jet)
xlabel('modeled subunits','fontsize',12);
ylabel('Inferred subunits','fontsize',12);


figure;
hold on;
beg=2000;
en=6000;
colo={'r','g','b','k','c','y'};
for i=1:size(iglist,2)
    dummy=[A(i,beg:en)',B(i,beg:en)'];
    dummy(dummy>0)=1;
    NT=size(dummy,2);
    maxX = size(dummy,1);
    [xx,yy]=rasterplot_snn(find(dummy),1:NT,maxX);
    subplot(size(iglist,2),1,i);
    plot(xx, yy,colo{i},'linewidth',1);
    title(['subunit - ',num2str(i)],'fontsize',12);
    axis tight;
    xrange = [0 maxX];
    xlim(xrange);
    axis off;
    set(gca,'Yticklabel','');  
end