function pretreatment_snn(CB,spike,nt,nx,ny)
% calculate STE

% calculate STA
[sta, ~, SS]= simpleSTC(CB',spike,nt);
sta = sta./norm(sta);
sta = reshape(sta,nt,nx,ny);
% singular value decomposition 
temp=[];
[temp,sigma,sp] = svd(reshape(sta,size(sta,1),[]),'econ');
sigma = diag(sigma);
svdsta = cell(size(sigma));
for m=1:length(svdsta)
	svdsta{m} = reshape(temp(:,m)*sigma(m)*sp(:,m)',size(sta));
end
tmSVD = temp(:,1);
spSVD = reshape(sp(:,1)',size(svdsta{1},2),size(svdsta{1},3));

spklist = spike(spike>0);


save('Datapre.mat','spklist','SS','sta','tmSVD','spSVD');
end
