%% -----This is a cript to calculate CBF of cilia from movies given by Ashleight --%
% -- the main difference is how to load the movie file. ----------% 

clear all;
filenames=dir('*.cxd');

for dd=1:numel(filenames)   
    filename=filenames(dd).name
    data = bfopen(filename);
    Nfs= size(data{1,1});Nfs=Nfs(1);
    fs=zeros([size(data{1,1}{1,1}),Nfs]);
    for t=1:Nfs
        fs(:,:,t)= data{1,1}{t,1};
    end;
    clear data
    %%% normalise  fs to be int
    minfs= min(fs(:));
    maxfs= max(fs(:));
    for t=1:Nfs
    fs(:,:,t)= uint8(255*(fs(:,:,t)-minfs )/(maxfs-minfs)) ;
    end
    fs=uint8(fs);
    fps=322; 
    px2mu=0.08; %%% mu
    
    
    f_lim=[4,25];prob_cilia=0.8;box_size=4;area_min=10;
    [F_temp,s,BW_temp] = find_cilia_nomovie(fs,fps,f_lim,prob_cilia,box_size);
    [F4,BW] = remove_debris(F_temp,area_min,f_lim);
    ss=imresize(BW,box_size,'nearest');    
    IF= imresize(F4,box_size,'nearest');
    imagesc(IF);colorbar();axis equal;
    save(strcat(filename(1:end-4),'.mat'),'IF','box_size','prob_cilia','f_lim',...
        'area_min','ss','F4','BW',...
        'F_temp','s','BW_temp',...
        'fps','px2mu');
    clear fs;
    
end

%% load data and make plot
clear all;
number=[17,37,28,27,19,22,15,16,8,4,25,31,41,35,7,34,23,...
    14,2,42,32,11,18,1,26,6,40,21,29,30,33,38,5,39,36,20,45,24,44,10,12,43,13,9,3];

group= [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3];



filenames_mat=dir('*.mat');
g=[];f=[];

for dd=1:numel(filenames_mat) 

load(filenames_mat(dd).name);
filename_mat=(filenames_mat(dd).name);
% take filen number
cc=1;
while ~strcmp(filename_mat(cc),'.')
cc=cc+1;
end

file_num= str2num(filename_mat(1:cc-1));

f_temp= IF(~isnan(IF(:)));
g_temp=ones(size(f_temp))*group(file_num==number) ;
f=cat(1,f,f_temp);
g=cat(1,g,g_temp);

end
%boxplot(f,g);
boxplot(f,g,'Labels',{'DNAI1','g1','NT'});
 x0=0;y0=0;width=800;height=800;
 set(gcf,'position',[x0,y0,width,height]);set(gca,'FontSize',15);
 ylabel('CBF [Hz]');
 saveas(gcf,'CBF_results.pdf');
 
 %% plot single CBF map for example
 load('1.mat');
 imagesc(IF); axis equal; colorbar()
 x0=0;y0=0;width=800;height=800;
 set(gcf,'position',[x0,y0,width,height]);set(gca,'FontSize',15);
 ylabel('CBF [Hz]');
 saveas(gcf,'CBF_colormap.pdf');
 
 
 
 


