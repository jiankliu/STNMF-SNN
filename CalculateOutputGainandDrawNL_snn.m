function CalculateOutputGainandDrawNL_snn(inname1,inname2,inname3,inname4,outname)
close all

               
load(inname1); 
load(inname2);%mat from the STNMF   e.g. subunit_from_STNMF.mat
load(inname3); %mat from the calculatTF e.g. subTemporalKernel.mat
load(inname4);

Nspike=size(SS,1);
nt=20;
K=size(unit,2);
Nstim=size(CB,2);
iisp=[1:1:Nstim];
fullS = makeStimRows(CB', nt, iisp);


%SS=SS.*repmat(spklist,1,size(SS,2)) ; % mutiplus the spike count
%N=Ngood; % test the first and second 
weight=zeros(1,K);
GS=zeros(K,Nstim);
line=[];
for m=1:K
    for n=1:Nstim
        line=reshape(fullS(n,:),nt,nx*ny);
        spConvSti=zeros(1,nt);
        for t=1:nt
           spConvSti(t)=sum(sum((unit{m}/norm(unit{m})).*reshape(line(t,:),nx,ny),1));
             %spConvSti(t)=sum(sum(unit{m}.*reshape(line(t,:),nx,ny),1));
        end
        GS(m,n)=spConvSti*subTemporal(m,:)';
    end
end
GS=GS';
bin=cell(1,K);


ina=ceil(Nstim/40);
for i=1:K
    [sX,index]=sort(GS(:,i));
    for j=1:40
        if j<40
            bin{i}(j,1)=mean(sX(1+ina*(j-1):ina*j));                       %x
            bin{i}(j,2)=sum(spike(index(1+ina*(j-1):ina*j)))/(Nstim/40);  %y
        elseif j==40
            bin{i}(j,1)=mean(sX(1+ina*(j-1):end));
            bin{i}(j,2)=sum(spike(index(1+ina*(j-1):end)))/(Nstim/40);
        end
    end
    weight(i)=max(bin{i}(:,2))-min(bin{i}(:,2));
end

save(outname,'weight','bin','GS');
end