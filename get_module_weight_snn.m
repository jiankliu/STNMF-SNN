function get_module_weight_snn(inname1,inname2,outname)
%2018 1 25 Jss
%calculate the weights for each modules W=S*pseudoinverse(M)
%inname1 : STE data set from the STE
%inname2 : subunit from STNMF
%outname : module_weight_.mat
load(inname1);
load(inname2);
Nspike=size(STE,1);
K=size(unit,2);
W_matrix=zeros(Nspike,K);
module=zeros(K,Nx*Ny);
for i=1:K
    module(i,:)=reshape(unit{1,i},1,Nx*Ny);
end

M_inv=pseudoinverse(module,2^32);

W_matrix=STE*M_inv;
W_matrix=normalize_columns(W_matrix);

avarage_W_M=zeros(1,20);
for i=1:K
    avarage_W_M(i)=sum(W_matrix(:,i))/size(W_matrix,1);

end

save(outname,'W_matrix','avarage_W_M')
disp( 'Module Weight finished!' );
end