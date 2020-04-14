function [F,s,BW] = find_cilia_nomovie(fs,fps,f_lim,prob_cilia,box_size)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin < 3 || isempty(prob_cilia)
    prob_cilia=0.8;
end

if nargin < 4 || isempty(box_size)
    box_size=4;
end

        mo.FrameRate = fps;
       row_dim=size(fs,2);col_dim=size(fs,1);
        frameload=round(2*mo.FrameRate/20);
        s=std(double(fs(:,:,30:30+frameload)),[],3);s=mat2gray(s);
        sm=medfilt2(s,[5,5]);sm=wiener2(sm,[5,5]);
        [pf,edg]=histcounts(sm(:),'Normalization','probability');
        cpf=cumsum(pf); 
        edg=edg(2:end); T=edg(cpf>prob_cilia ); T=T(1);
        BW= imbinarize(sm,T);  %%%% mask whit standard deviation
        
        Nbox=floor((row_dim*col_dim)/box_size^2);  %%%% number of boxes based on the total area
        nboxes_row=floor(row_dim/box_size);
        nboxes_col=floor(col_dim/box_size);       
        [X,Y]=meshgrid(1:box_size:(nboxes_row)*box_size,1: box_size: (nboxes_col)*box_size);

        F= nan([size(X)]);

for xx=1:(size(X,2)-1);
    for yy=1:(size(X,1)-1);    
        
       mask_value= BW(Y(yy,xx):(Y(yy+1,xx+1)-1),X(yy,xx):(X(yy+1,xx+1)-1));
       good= (mean(mask_value(:))==1) ; 
        
        if good==1;
            froi=fs(Y(yy,xx):(Y(yy+1,xx+1)-1),X(yy,xx):(X(yy+1,xx+1)-1),:);
            temp_ss=ones([size(froi,1),size(froi,2)]);
            [fp,fq,m_pxx] = average_fft_meanpk_func(froi,temp_ss,mo.FrameRate,f_lim);
            F(yy,xx)=fp;
         end
        
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fp,fq,m_pxx] = average_fft_meanpk_func(fs,mask,FR,fq_lim,res)
    if nargin < 5 | isempty(res)
        res=2;
    end
    
    if nargin < 4 | isempty(fq_lim)
        fq_lim(1)= 10;
        fq_lim(2)= 30;
    end
    N_frames= size(fs,3);
    fs_roi= fs(repmat(logical(mask),[1,1,N_frames]));
    fs_roi=reshape(fs_roi,[sum(mask(:)),N_frames]);
    %roi=double(fs_roi)- mean(fs_roi,2);
    roi=double(fs_roi)- movmean(fs_roi,60,2);

%%%% periodogram needs time on the row and pixel on the coloun; 
    window = hann(floor(N_frames));
    window= repmat(window,[1,size(roi,1)])';
    n= floor(N_frames/res);
    if mod(n,2)==0; n= n-1;end
    pxx= abs(fft(double(roi).*window,n,2)).^2;
    m_pxx= mean(pxx(:,1:floor(n/2)),1);
    fq= (0:(FR./n):(FR./2-FR./n));
    f_range=fq> (fq_lim(1)) & fq<(fq_lim(2));
    [pks,locs,w,p] = findpeaks(m_pxx(f_range),fq(f_range));%%%% find the peaks frequency in the selected freq range
    [~,ind_sort]= sort(pks);                                           %%%% sort peaks and get an index 
    pks=pks(ind_sort); locs=locs(ind_sort);w=w(ind_sort);p=p(ind_sort);  %%% order all the variables with the same index
    [~,where] = max(pks);                                       
      
    
    %%%%% asign the frequency as the one with the highest peak, with safety
    %%%%% for harmonics
    if numel(locs)==0; fp=nan;    
    elseif numel(locs)==1; fp= locs(end);
    elseif numel(locs)>1; 
        if abs((locs(end)/2)-locs(end-1))< 0.3 & pks(end-1)> 0.66*pks(end);
            fp= locs(end-1); 
        else fp=locs(end);
        end
    end
end



end

