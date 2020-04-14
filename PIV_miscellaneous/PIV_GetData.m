function [x,y,u,v] = PIV_GetData(frame_stack,s,p)
%ciao

amount= size(frame_stack,3);


%% PIV analysis loop
if mod(amount,2) == 1 %Uneven number of images?
    disp('Image folder should contain an even number of images.')
    %remove last image from list
    amount=amount-1;
%    filenames(size(filenames,1))=[];
end
disp(['Found ' num2str(amount) ' images (' num2str(amount/2) ' image pairs).'])
x=cell(amount/2,1);
y=x;
u=x;
v=x;
typevector=x; %typevector will be 1 for regular vectors, 0 for masked areas
counter=0;
%% PIV analysis loop:
for i=1:2:amount
    counter=counter+1;
 %   image1=imread(fullfile(directory, filenames{i})); % read images
 %   image2=imread(fullfile(directory, filenames{i+1}));
 %%%%%image preprocessing
 %   image1 = PIVlab_preproc (frame_stack(:,:,i),p{1,2},p{2,2},p{3,2},p{4,2},p{5,2},p{6,2},p{7,2},p{8,2}); %preprocess images
 %   image2 = PIVlab_preproc (frmae_Stack(:,:,i+1),p{1,2},p{2,2},p{3,2},p{4,2},p{5,2},p{6,2},p{7,2},p{8,2});

 
    [x{counter} y{counter} u{counter} v{counter} typevector{counter}] = piv_FFTmulti (frame_stack(:,:,i),frame_stack(:,:,i+1),s{1,2},s{2,2},s{3,2},s{4,2},s{5,2},s{6,2},s{7,2},s{8,2},s{9,2},s{10,2});
    clc
    disp([int2str((i+1)/amount*100) ' %']);
    
    % Graphical output (disable to improve speed)
    %%{
    imagesc(double(frame_stack(:,:,i))+double(frame_stack(:,:,i+1)));colormap('gray');
    hold on
    quiver(x{counter},-y{counter},u{counter},-v{counter},'g','AutoScaleFactor', 1.5);
    hold off;
    axis image;
    title('PIV','interpreter','none')
    set(gca,'xtick',[],'ytick',[])
    drawnow;
    %%}
end



end