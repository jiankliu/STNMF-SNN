function STE=get_STE_snn(ss,spklist,tmSVD,nt,nx,ny)
% calculate STE
%this ss is SS
Nspike=size(ss,1);

% get STE effective spike-triggered stimulus ensemble
% SS is a stimulus set inwhich stimulus can stimulate spike
% First step : mutiplus the spike count
ss=ss.*repmat(spklist,1,size(ss,2)) ; 
STE=zeros(Nspike,nx*ny);
% Second step :let SS convolute temporal filter got from SVD

for i=1:Nspike
%    for j=1:nx*ny 
        STE(i,:)=tmSVD'*reshape(ss(i,:),nt,nx*ny); 
        
%    end
end
Nx=nx;Ny=ny;
save('STE.mat','STE','Nx','Ny');
end
