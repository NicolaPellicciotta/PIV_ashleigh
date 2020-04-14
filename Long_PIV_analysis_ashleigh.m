%% -----This script is to run PIV analysis on data provided by Ashleigh ---%
% Step 1: you run PIVlab on all the video files
% Step 2: post analysis: average vectors over time and remove useless data
% Step 3: gather results from all the files and plots

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 1------- PIV ANALYSIS ----------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% set Standard PIV Settings
s = cell(10,2); % To make it more readable, let's create a "settings table"
%Parameter                       %Setting           %Options
s{1,1}= 'Int. area 1';           s{1,2}=256;         % window size of first pass
s{2,1}= 'Step size 1';           s{2,2}=64;         % step of first pass
s{3,1}= 'Subpix. finder';        s{3,2}=1;          % 1 = 3point Gauss, 2 = 2D Gauss
s{4,1}= 'Mask';                  s{4,2}=[];         % If needed, generate via: imagesc(image); [temp,Mask{1,1},Mask{1,2}]=roipoly;
s{5,1}= 'ROI';                   s{5,2}=[];         % Region of interest: [x,y,width,height] in pixels, may be left empty
s{6,1}= 'Nr. of passes';         s{6,2}=3;          % 1-4 nr. of passes
s{7,1}= 'Int. area 2';           s{7,2}=128;         % second pass window size
s{8,1}= 'Int. area 3';           s{8,2}=64;         % third pass window size
s{9,1}= 'Int. area 4';           s{9,2}=32;         % fourth pass window size
s{10,1}='Window deformation';    s{10,2}='*linear'; % '*spline' is more accurate, but slower

%%% Standard image preprocessing settings
p = cell(8,1);
%Parameter                       %Setting           %Options
p{1,1}= 'ROI';                   p{1,2}=s{5,2};     % same as in PIV settings
p{2,1}= 'CLAHE';                 p{2,2}=0;          % 1 = enable CLAHE (contrast enhancement), 0 = disable
p{3,1}= 'CLAHE size';            p{3,2}=50;         % CLAHE window size
p{4,1}= 'Highpass';              p{4,2}=0;          % 1 = enable highpass, 0 = disable
p{5,1}= 'Highpass size';         p{5,2}=15;         % highpass size
p{6,1}= 'Clipping';              p{6,2}=0;          % 1 = enable clipping, 0 = disable
p{7,1}= 'Wiener';                p{7,2}=0;          % 1 = enable Wiener2 adaptive denaoise filter, 0 = disable
p{8,1}= 'Wiener size';           p{8,2}=3;          % Wiener2 window size

%% ------run the PIV analysis on all video files-----------

% directory name with all the video files
long={...
   'beads_assay2'}

% long can be more than one directory
%example    long={'dirname1','dirname2', etc...}
   

%loop on each of the directory, that I call insert
for insert=long;

% load files
insert=insert{1};

% data_dir is the path were the directories are stored
%data_dir = strcat('/run/user/10704/gvfs/smb-share:server=sf3.bss.phy.private.cam.ac.uk,share=space/np451/ependymalJune/4.7.18/',insert);
data_dir= strcat('/home/np451/Desktop/ashleigh/',insert)

% results are going to be stored here
a_folder=strcat('/home/np451/Desktop/ashleigh/PIV/',insert);
mkdir(a_folder); 

cd(data_dir) 

  
% store the name of all the video ending with suffix
suffix='.tif';
direc = dir('*.tif');
N_files= size(direc,1);

% loop PIV analysis on each video in the directory insert

for i=1:N_files
    
        cd(data_dir)

% setting savename and loading the video. Ashleigh use 
% a different format of video that can be opend with bfopen.
% I commented the code for loading with moviereader.

        exp_name = direc(i).name;
        exp_name=exp_name(1:end-4);
        disp(exp_name)
%    if exist(strcat(a_folder,'/',exp_name)) == 0
        
%        mkdir(strcat(a_folder,'/',exp_name));

%------- loading with moviereader ---------%

        %mo=moviereader(direc(i).name);
        %movie_path=strcat(data_dir,mo.Filename)
        %frame_stack=mo.read();
        %frame_stack=(frame_stack)-movmin(frame_stack,40,3);
        %max_int= double(max(frame_stack(:)));
        %if max_int>256        
        %frame_stack= uint8((double(frame_stack)/max_int)*255);
        %else frame_stack= uint8(frame_stack);
        %end

%------- loading with bfopen ---------%

        fps=322; % set frame per second
        px2mu=0.08; %%% pixel to micron

        data = bfopen(direc(i).name);
        Nfs= size(data{1,1});Nfs=Nfs(1);
        fs=zeros([size(data{1,1}{1,1}),Nfs]);
        for t=1:Nfs
            fs(:,:,t)= data{1,1}{t,1};
        end;
        clear datawhich 

%------ select max intensity and remove high intensity background --%   
%------ only for video from Ashleigh ------------------------------%
%------- PIV needs uint8 pixel intensity!!!! so normalise ---------%
        
        f1=(fs(:,:,1));
        thresh = multithresh(f1,2);
        maxfs= mean(f1(:))+3*std(f1(:));%thresh(2)+thresh(2)
        minfs= min(fs(:));
        
        for t=1:Nfs
        fs(:,:,t)= uint8(255*(fs(:,:,t)-minfs )/(maxfs-minfs)) ;
        end
        fs=uint8(fs);

 
%-------for convention I use frame_stack-----------------%
        frame_stack=fs;
        clear fs;
        
%------- remove high frequency spatial noise---------------------%
        %        frame_stack=imadjustn(frame_stack);
        for kk=1:size(frame_stack,3); 
            frame_stack(:,:,kk)= wiener2(frame_stack(:,:,kk),[5,5]);
        end
        
        
%        frame_stack= uint8(double(frame_stack)/2^(8));  %%%% converting images from 16 to 8 bit
%        frame_stack= uint8((frame_stack));  %%%% for images at 8 bit


% -----------  PIV Anlaysis with PIVLab ---------%

        [X,Y,U,V] = PIV_GetData(frame_stack,s,p);
        [x,y,u,v] =  PIV_ChangeFormat(X,Y,U,V);

        close('all');
        
        cd(a_folder);

% ------------- PostProcessing ------------------------------%
%--- mainly removing vectors that are too large and then -----%
%--- interpolating the missing vectors------------------------%

% parameter to play with, it is important that you don't cut reasonable vectors
        ulim= [-10,10]; 
        vlim= [-10,10];

        [U1,V1]= PIV_Validation(X,Y,U,V,ulim,vlim);
        [x,y,u1,v1] =  PIV_ChangeFormat(X,Y,U1,V1);

%        save('PIV_post.mat','U1','V1','ulim','vlim','u1','v1');

%%%----------------- Figures ----------------------%

%------ save a figure of standard deviation over time
        ss= std(double(frame_stack),[],3);
        figure(1);
        subplot(2,2,1)
        title(strcat(exp_name,'frame 1'));
        imagesc(ss);
      
%--------save the PIV results before processing
        subplot(2,2,2)
        title(strcat(exp_name,' PreProcessing'));
        quiver(x,-y,nanmean(u,3),-nanmean(v,3),3);
  
%------- save the PIV results after processing
        subplot(2,2,3)
        title(strcat(exp_name,' PostProcessing'));
        quiver(x,-y,mean(u1,3),-mean(v1,3),3);
        fig=figure(1);
        saveas(fig,strcat(exp_name,'_PostPro.png'));
        close(1); 
  
        clear frame_stack; 
%-------save all the variables with mat extension        
        save(strcat(exp_name,'.mat'));

end

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 2------------ post analysis  ------------------ %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this was very dependent on what you need. 
% here I wanted to have, for each video, only the average velocity over 
% time and to remove vectors outside a circle (because of the PIV was on a 
% circular culture plate ). After save everything in a _post.mat file  

clear all
cd '/home/np451/Desktop/ashleigh/PIV/beads_assay2';

filenames=dir('*ome.mat');
for dd=1:numel(filenames)
   filename=filenames(dd).name;
   load(filename)
   u1m=mean(u1,3);
   v1m=mean(v1,3);
   quiver(x,y,u1m,v1m,'r')
   hold on; axis equal
    [px,py] = getpts ;   % click at the center and approximate Radius
    r = sqrt(diff(px).^2+diff(py).^2) ;
    th = linspace(0,2*pi) ;
    xc = px(1)+r*cos(th) ; 
    yc = py(1)+r*sin(th) ; 
    plot(xc(1,:),yc(1,:),'b.') ;
% Keep only points lying inside circle
    idx = inpolygon(x(:),y(:),xc(1,:)',yc(1,:)) ;
    v1m(~idx)=nan; 
    u1m(~idx)=nan;
% normalise vector for the orientation    
    M= sqrt(u1m.^2 +v1m.^2);Mm=nanmedian(M(:));
    nu= u1m./M; nv= v1m./M;
    px2mu=3.25;
%make nice figure
    figure()
    quiver(x*px2mu*1e-3,y*px2mu*1e-3,u1m,v1m,'k','LineWidth',1.4);axis equal
    xlabel('[mm]');ylabel('[mm]');
    x0=0;y0=0;width=800;height=800;
    set(gcf,'position',[x0,y0,width,height]);set(gca,'FontSize',15);
    saveas(gcf,strcat(filename(1:end-4),'_PIV_result.pdf'));
    close all

%---- spatial correlation as a function of distance and angle
%    bin_res= 32; 
%    [idr,idth,cc,ecc,n_cc] = corr_orientation_theta_func(x(:),y(:),nu(:),nv(:),bin_res);
%    plot(idr*3.25,cc)

%---- save all the variables in a post_mat file. 

    save(strcat(filename(1:end-4),'_post.mat'))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 3------- gather results and plot average results-------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


cd '/home/np451/Desktop/ashleigh/PIV/beads_assay2';

clear all
close all

%number=[3,6,9,1,4,2,7,8,5];

%%% define groups knockouts
group={'DNA gAA','DNA gAB','g1','NT','gAA'};
%%% define groups viscosity of the medium 
viscosity = {'4% ','8% ','ON','PBS'};


col_cc={'ko','r.','b>'};

filenames_post=dir('*post.mat');
g=[];MM=[];v=[];dds=[];
av_g=[];av_MM=[];av_v=[];
fps=10;

hh=1;

for dd=1:numel(filenames_post) 
    
filename_post=filenames_post(dd).name;
file_num=str2num(filename_post(1));

%---- load results of this video in class ciao

ciao= load(filename_post);


% ---- look for which knockout group belongs to
cc=1;
while isempty(strfind(filename_post,group{cc}))
cc=cc+1;
end

% ---- look for which viscosity group belongs to
vv=1;
while isempty(strfind(filename_post,viscosity{vv}))
vv=vv+1;
end


% ---- accumolate all the magnitude of velocities 

M=ciao.M; px2mu =ciao.px2mu; fps= ciao.fps; 
M_temp=M(~isnan(M))*px2mu*fps*1e-3;
%g_temp=ones(size(M_temp(:)))*cc;
%g_temp=ones(size(M_temp))*cc;
%g_temp=ones(size(M_temp))*group(file_num==number);
g_temp= repmat({group{cc}},1,numel(M_temp(:)));
v_temp= repmat({viscosity{vv}},1,numel(M_temp(:)));
dd_temp= repmat(dd,1, numel(M_temp(:)));

MM=cat(1,MM,M_temp(:));  % all the velocity magnitude
%g=cat(1,g,g_temp);
g= [g,g_temp];          % all the knock out
v= [v,v_temp];          % all the viscosity
dds=cat(2,dds,dd_temp); % all the file number


% -----accomulate all the  average quantities

av_M(hh) = median(M_temp(:)); % average magnitude velocity
av_g =[av_g, g_temp(end)];    % the knock out group
av_v = [av_v,v_temp(end)];    % the viscosity of the medium  


hh=hh+1;

end

% ---- find the groups for the average and total quantities

[G,idg,idv] = findgroups(g,v);
[av_G,av_idg,av_idv] = findgroups(av_g,av_v)

%% plot data using avergae quantities and threshold on the velocity

close all

str_array= {'4% ','8% ','PBS','ON'};
for ff=1:numel(str_array)

subplot(2,2,ff);
str= str_array{ff};

 %%% average velocity for wildtype in this medium viscosity

 threshold_ind = strcmp(av_v,str) & strcmp(av_g,'NT');
thre_M= median(av_M(threshold_ind));%-std(av_M(threshold_ind));
cc=1;av_perc=[];
for dd = unique(dds)
    av_perc(cc)= sum(MM(dds==dd)>thre_M)/numel(MM(dds==dd))
    cc=cc+1
end

boxplot(100*av_perc(strcmp(av_v,str)),av_G(strcmp(av_v,str)),'Symbol','o','Labels',av_idg(strcmp(av_idv,str)));
hold on;
cc_plot=1;
for bau=unique(av_G(strcmp(av_v,str)));    
    plot(cc_plot,100*av_perc(av_G==bau),'ko','MarkerFace','k');hold on;
    cc_plot=cc_plot+1;
end

%boxplot(MM,g);
x0=0;y0=0;width=800;height=1000;
 set(gcf,'position',[x0,y0,width,height]);set(gca,'FontSize',8);
 ylabel('percentage of active cilia [%]');
 title(strcat('viscosity  ',str,' NT med flow=',num2str(thre_M,2),' mm/s'),'FontSize',8);
 
end
 saveas(gcf,strcat('Percentage of active cilia for different viscosity.pdf'))


str='NT';
figure(5);
boxplot(av_M(strcmp(av_g,'NT')),av_G(strcmp(av_g,'NT')),'Symbol','o','Labels',av_idv(strcmp(av_idg,'NT')));
N_scatter = numel(idv(strcmp(idg,str)));
groups=idv(strcmp(idg,str));
hold on;
x0=0;y0=0;width=800;height=800;
 set(gcf,'position',[x0,y0,width,height]);set(gca,'FontSize',15);
 ylabel('flow velocity [mm/s]');
 title('NT vs viscosity')

 
%% ---- plot data using all the data from the inserts (not average quantities)

close all
str_array= {'4% ','8% ','PBS','ON'};
for ff=1:numel(str_array)

figure(ff);
str= str_array{ff};
boxplot(MM(strcmp(v,str)),G(strcmp(v,str)),'Symbol','o','Labels',idg(strcmp(idv,str)));
hold on;

%boxplot(MM,g);
x0=0;y0=0;width=800;height=800;
 set(gcf,'position',[x0,y0,width,height]);set(gca,'FontSize',15);
 ylabel('flow velocity [mm/s]');
  title(strcat('viscosity ',str))
 saveas(gcf,strcat('Method2 viscosity ',str,'.pdf'))

 
end

str='NT';
figure(5);
boxplot(MM(strcmp(g,'NT')),G(strcmp(g,'NT')),'Symbol','o','Labels',idv(strcmp(idg,'NT')));
N_scatter = numel(idv(strcmp(idg,str)));
groups=idv(strcmp(idg,str));
hold on;
x0=0;y0=0;width=800;height=800;
 set(gcf,'position',[x0,y0,width,height]);set(gca,'FontSize',15);
 ylabel('flow velocity [mm/s]');
 title('NT vs viscosity')
% saveas(gcf,'Method2 NT vs viscosity.pdf')


%% --- end plot staff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --- after this there is a miscelaneous of analysis that I have no memeory of ---%




%%
clear all
close all

number=[3,6,9,1,4,2,7,8,5];
group=[1,1,1,2,2,2,3,3,3];
col_cc={'ko','r.','b>'};

filenames_post=dir('*post.mat');
g=[];MM=[];
fps=10;

for dd=1:numel(filenames_post) 
    
filename_post=filenames_post(dd).name;
file_num=str2num(filename_post(1));
load(filename_post);



M_temp=M(~isnan(M))*px2mu*fps*1e-3;
g_temp=ones(size(M_temp))*group(file_num==number);

MM=cat(1,MM,M_temp);
g=cat(1,g,g_temp);
figure(1)
hold on;
plot(idr*px2mu*1e-3,cc,col_cc{group(file_num==number)});
end
figure(1);
x0=0;y0=0;width=800;height=800;
 set(gcf,'position',[x0,y0,width,height]);set(gca,'FontSize',15);
 ylabel('correlation function'); xlabel('mutual distance [mm]');
 saveas(gcf,'orientation_correlation_results.pdf');


figure(2);
boxplot(MM,g,'Labels',{'DNAI1','g1','NT'});
%boxplot(MM,g);
x0=0;y0=0;width=800;height=800;
 set(gcf,'position',[x0,y0,width,height]);set(gca,'FontSize',15);
 ylabel('flow velocity [mm/s]');
 saveas(gcf,'flow_velocity_results.pdf');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% analysis and graph2%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ----2-----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
cd '/home/np451/Desktop/ashleigh/PIV/beads_assay2';

filenames=dir('*ome.mat');
for dd=1:numel(filenames)
   filename=filenames(dd).name;
   load(filename)
   u1m=mean(u1,3);
   v1m=mean(v1,3);
   quiver(x,y,u1m,v1m,'r')
   hold on; axis equal
    [px,py] = getpts ;   % click at the center and approximate Radius
    r = sqrt(diff(px).^2+diff(py).^2) ;
    th = linspace(0,2*pi) ;
    xc = px(1)+r*cos(th) ; 
    yc = py(1)+r*sin(th) ; 
    plot(xc(1,:),yc(1,:),'b.') ;
% Keep only points lying inside circle
    idx = inpolygon(x(:),y(:),xc(1,:)',yc(1,:)) ;
    v1m(~idx)=nan; 
    u1m(~idx)=nan;
% normalise vector for the orientation    
    M= sqrt(u1m.^2 +v1m.^2);Mm=nanmedian(M(:));
    nu= u1m./M; nv= v1m./M;
    px2mu=3.25;
%make nice figure
    figure()
    quiver(x*px2mu*1e-3,y*px2mu*1e-3,u1m,v1m,'k','LineWidth',1.4);axis equal
    xlabel('[mm]');ylabel('[mm]');
    x0=0;y0=0;width=800;height=800;
    set(gcf,'position',[x0,y0,width,height]);set(gca,'FontSize',15);
    saveas(gcf,strcat(filename(1:end-4),'_PIV_result.pdf'));
close all
 
    bin_res= 32; 
%    [idr,idth,cc,ecc,n_cc] = corr_orientation_theta_func(x(:),y(:),nu(:),nv(:),bin_res);
%    plot(idr*3.25,cc)
    save(strcat(filename(1:end-4),'_post.mat'))
end


%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%% average
 %%%%%%%%%%%%%%%%%%%%%%%%%%%% plots%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
clear all
close all



number=[3,6,9,1,4,2,7,8,5];
group={'DNA gAA','DNA gAB','g1','NT','gAA'};
viscosity = {'4% ','8% ','ON','PBS'};


col_cc={'ko','r.','b>'};

filenames_post=dir('*post.mat');
g=[];MM=[];v=[];MMe=[];tot_files=[];
fps=10;

for dd=1:numel(filenames_post) 
    
filename_post=filenames_post(dd).name;
file_num=str2num(filename_post(1));
ciao= load(filename_post);

cc=1;
while isempty(strfind(filename_post,group{cc}))
cc=cc+1;
end

vv=1;
while isempty(strfind(filename_post,viscosity{vv}))
vv=vv+1;
end

M=ciao.M; px2mu =ciao.px2mu; fps= ciao.fps; 
M_temp=M(~isnan(M))*px2mu*fps*1e-3;
Me_temp= std(M_temp(:));
M_temp= median(M_temp(:));

%g_temp=ones(size(M_temp(:)))*cc;
%g_temp=ones(size(M_temp))*cc;
%g_temp=ones(size(M_temp))*group(file_num==number);
g_temp= repmat({group{cc}},1,numel(M_temp(:)));
v_temp= repmat({viscosity{vv}},1,numel(M_temp(:)));


MM=cat(1,MM,M_temp(:));
MMe=cat(1,MMe,Me_temp(:));
%g=cat(1,g,g_temp);
g= [g,g_temp];
v= [v,v_temp];
tot_files= [tot_files,{filename_post}];

%figure(1)
%hold on;
%plot(idr*px2mu*1e-3,cc,col_cc{group(file_num==number)});
end
%figure(1);
%x0=0;y0=0;width=800;height=800;
% set(gcf,'position',[x0,y0,width,height]);set(gca,'FontSize',15);
% ylabel('correlation function'); xlabel('mutual distance [mm]');
% saveas(gcf,'orientation_correlation_results.pdf');

[G,idg,idv] = findgroups(g,v);

group_name= strcat(v,'~ ',g);
T=table(tot_files', MM,MMe);
 writetable(T,'Method1 median and std for all exp.txt');
  



%%% box plot for viscosity 4%
%%
close all
str_array= {'4% ','8% ','PBS','ON'};
for ff=1:numel(str_array)

figure(ff);
str= str_array{ff};
boxplot(MM(strcmp(v,str)),G(strcmp(v,str)),'Symbol','o','Labels',idg(strcmp(idv,str)));
hold on;
N_scatter = numel(idg(strcmp(idv,str)));
groups=idg(strcmp(idv,str));
for mm=1:N_scatter
    ind= strcmp(v,str) & strcmp(g,groups{mm});
    mm_array=mm*ones([1,numel(MM(ind))]);
plot(mm_array,MM(ind),'ko','MarkerSize',7,'MarkerFaceColor','k','MarkerEdgeColor','g');
end
%boxplot(MM,g);
x0=0;y0=0;width=800;height=800;
 set(gcf,'position',[x0,y0,width,height]);set(gca,'FontSize',15);
 ylabel('flow velocity [mm/s]');
  title(strcat('viscosity ',str))
 saveas(gcf,strcat('Method1 viscosity ',str,'.pdf'))

 
end

str='NT';
figure(5);
boxplot(MM(strcmp(g,'NT')),G(strcmp(g,'NT')),'Symbol','o','Labels',idv(strcmp(idg,'NT')));
N_scatter = numel(idv(strcmp(idg,str)));
groups=idv(strcmp(idg,str));
hold on;
for mm=1:N_scatter
    ind= strcmp(g,str) & strcmp(v,groups{mm});
    mm_array=mm*ones([1,numel(MM(ind))]);
plot(mm_array,MM(ind),'ko','MarkerSize',7,'MarkerFaceColor','k','MarkerEdgeColor','g');
end
x0=0;y0=0;width=800;height=800;
 set(gcf,'position',[x0,y0,width,height]);set(gca,'FontSize',15);
 ylabel('flow velocity [mm/s]');
 title('NT vs viscosity')
 saveas(gcf,'Method1 NT vs viscosity.pdf')

 
 %% statistical anova analysis
 close all
 group_name= strcat(v,'~ ',g);
 [p,tbl,stats]=anova1(MM,group_name);
 x0=0;y0=0;width=800;height=400;
 figure(2);
 set(gcf,'position',[x0,y0,width,height]);set(gca,'FontSize',5);
 ylabel('flow velocity [mm/s]');
 title('Total boxplot from Anova');
 saveas(gcf,'Method2 total boxplot.pdf');
 
 
 %%%% make tables for anova results and average values
 [c,m,h,nms] = multcompare(stats);
 T=table(nms, num2cell(m(:,1)),num2cell(m(:,2)));
 writetable(T,'median and std all.txt');
  
 %%%% make total tables for anova
 T=table(nms(c(:,1)),nms(c(:,2)),c(:,6))
 writetable(T,'anova results all.txt');
 
 %%%% make only significant for anova
 ind = c(:,6)< 0.05
 T=table(nms(c(ind,1)),nms(c(ind,2)),c(ind,6))
 writetable(T,'significant anova results all.txt');
 