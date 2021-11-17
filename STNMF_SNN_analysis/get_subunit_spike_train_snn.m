% build the subunit spike train
function GC_subunit_sp=get_subunit_spike_train_snn(Ngood,sp,spklist_sub)

%% get the order of spklist from sp
GC_subunit_sp=[];
iid=[];
j=1;
for i=1:length(sp)
    if sp(i)>0
        iid(1,j)=j;
        iid(2,j)=i;
        j=j+1;
    end
end

%% build the subunit spike train
sp_sub_all=[];
for L=1:Ngood
    sp_sub=zeros(1,length(sp));
    for i=1:length(spklist_sub{L}(1,:))

        loc_spklist=spklist_sub{L}(2,i); %the spike location in the spklist
        val_spike=spklist_sub{L}(1,i); % the value of the spike
      	if loc_spklist>0
            loc_sp=iid(2,loc_spklist); %the spike location in the sp
            sp_sub(loc_sp)=val_spike;  %
        end

    end
    
    sp_sub_all(L,:)=sp_sub;

end

GC_subunit_sp=sp_sub_all;
disp( 'subunit spike train finished!' );
end