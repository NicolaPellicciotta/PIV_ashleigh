x=resultslist{1,1};
y=resultslist{2,1};
u= zeros([size(x),size(resultslist,2)]);
v= zeros([size(x),size(resultslist,2)]);
u1= zeros([size(x),size(resultslist,2)]);
v1= zeros([size(x),size(resultslist,2)]);

for i=1:size(resultslist,2)
     u(:,:,i)=resultslist{3,i};
     v(:,:,i)=resultslist{4,i};
     u1(:,:,i)=resultslist{7,i};
     v1(:,:,i)=resultslist{8,i};
     
end

%um=mean(u,3)
     