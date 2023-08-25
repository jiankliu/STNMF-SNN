function calculatTF_subSP_snn(inname1,inname2,outname,K,Nmodules,nt,nx,ny)
%clear all 
close all
load(inname1);  %data from ModelSubRF_simpleNL.mat  cell_simpleNL.mat
load(inname2);  %data from STNMFanalysis.m          subunit_from_STNMF.mat
%K=20;          %            
%Nmodules       number od sub TF SP to display

Nspike=size(SS,1);
hSubunitSubtemporalFig = figure('Position', [50, 50, 1600, 800]);
% get STE effective spike-triggered stimulus ensemble
SS=SS.*repmat(spklist,1,size(SS,2)) ; % mutiplus the spike count

subTemporal=zeros(K,nt);
for m=1:K
    for n=1:Nspike
        line=reshape(SS(n,:),nt,nx*ny);
        for t=1:nt
            subTemporal(m,t)=subTemporal(m,t)+sum(sum(unit{m}.*reshape(line(t,:),nx,ny),1));
        end
    end
end

%normalize
subTemporal=subTemporal/Nspike;
for i=1:K
    subTemporal(i,:)=subTemporal(i,:)./norm(subTemporal(i,:));
end
save(outname,'subTemporal')


figure(hSubunitSubtemporalFig);
color=['r','b','g','k','c','y','r','b','g','k','c','y','r','b','g','k','c','y','r','b','g','k','c','y','r','b','g','k','c','y','r','b','g','k','c','y'];
datarange = [-0.5, 0.5];
for j=1:Nmodules
    subplot(2,Nmodules,(2*(j-1)+1));
    imagesc(unit{j}/norm(unit{j}),datarange);
    colormap(colormap_RF(datarange));
    subplot(2,Nmodules,2*j);
    plot(subTemporal(j,:),color(j),'linewidth',2);
    axis([0,20,-0.6,0.6]);
end

print('-dpng',[outname,'hSubunitSubtemporalFig.png'])

end