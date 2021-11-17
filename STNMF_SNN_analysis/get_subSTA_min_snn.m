function spklist_sub=get_subSTA_min_snn(Ngood,spklist,W_matrix,STE,nx,ny)
STE_sub=cell(1,Ngood);
temp=zeros(2,1);
spklist_sub=cell(1,Ngood); %each cell storage the spike number and the position of this spike in the whole spklist
for i=1:Ngood
    spklist_sub{i}=temp;
end
Nspike=length(spklist);

% for each subunit ,we get the spike-triggered simulus and spike that are corresponding to the max weight 
for i=1:Nspike
    [C,I]=min(W_matrix(i,1:Ngood));
    if I <=Ngood
		loca=length(spklist_sub{I}(1,:));
		STE_sub{I}(loca+1,:)=STE(i,:);
        spklist_sub{I}(1,loca+1)=spklist(i);%put each subunit spike to each cell in the first raw
        spklist_sub{I}(2,loca+1)=i;         %put the sub spike position in the second row
    end
end

% 
% subSTA=cell(1,Ngood);
% %calculate sta for subunit and plot them
% for j=1:Ngood
% 
%     ste = STE_sub{j}/std(reshape(STE_sub{j},[],1));
%     STA = reshape( mean( ste ), [ny,nx] );
%     subSTA{j}=STA/norm(STA);
% 
% end

disp( 'subSTA finished!' );
end