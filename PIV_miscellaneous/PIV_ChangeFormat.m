function [X,Y,U,V]= PIV_ChangeFormat(x,y,u,v)
%%% from PIVlab weird format to array format



N= size(u,1);
X=x{1};
Y=y{1};
U= zeros([size(X),N]);
V= zeros([size(X),N]);


for i=1:N
     U(:,:,i)=u{i};
     V(:,:,i)=v{i};

end

end