function [u_filt,v_filt]= PIV_Validation(x,y,u,v,ulim,vlim)

% PIV postprocessing loop
% Settings
umin = ulim(1); % minimum allowed u velocity
umax = ulim(2); % maximum allowed u velocity
vmin = vlim(1); % minimum allowed v velocity
vmax = vlim(2); % maximum allowed v velocity
stdthresh=6; % threshold for standard deviation check
epsilon=0.15; % epsilon for normalized median test
thresh=3; % threshold for normalized median test



amount= size(u,1);

typevector=x; %typevector will be 1 for regular vectors, 0 for masked areas
u_filt=cell(floor(amount/2),1);
v_filt=u_filt;
typevector_filt=u_filt;

for PIVresult=1:size(x,1)
    u_filtered=u{PIVresult,1};
    v_filtered=v{PIVresult,1};
    typevector_filtered=typevector{PIVresult,1};
    %vellimit check
    u_filtered(u_filtered<umin)=NaN;
    u_filtered(u_filtered>umax)=NaN;
    v_filtered(v_filtered<vmin)=NaN;
    v_filtered(v_filtered>vmax)=NaN;
    % stddev check
    meanu=nanmean(nanmean(u_filtered));
    meanv=nanmean(nanmean(v_filtered));
    std2u=nanstd(reshape(u_filtered,size(u_filtered,1)*size(u_filtered,2),1));
    std2v=nanstd(reshape(v_filtered,size(v_filtered,1)*size(v_filtered,2),1));
    minvalu=meanu-stdthresh*std2u;
    maxvalu=meanu+stdthresh*std2u;
    minvalv=meanv-stdthresh*std2v;
    maxvalv=meanv+stdthresh*std2v;
    u_filtered(u_filtered<minvalu)=NaN;
    u_filtered(u_filtered>maxvalu)=NaN;
    v_filtered(v_filtered<minvalv)=NaN;
    v_filtered(v_filtered>maxvalv)=NaN;
    % normalized median check
    %Westerweel & Scarano (2005): Universal Outlier detection for PIV data
    [J,I]=size(u_filtered);
    medianres=zeros(J,I);
    normfluct=zeros(J,I,2);
    b=1;
    for c=1:2
        if c==1; velcomp=u_filtered;else;velcomp=v_filtered;end %#ok<*NOSEM>
        for i=1+b:I-b
            for j=1+b:J-b
                neigh=velcomp(j-b:j+b,i-b:i+b);
                neighcol=neigh(:);
                neighcol2=[neighcol(1:(2*b+1)*b+b);neighcol((2*b+1)*b+b+2:end)];
                med=median(neighcol2);
                fluct=velcomp(j,i)-med;
                res=neighcol2-med;
                medianres=median(abs(res));
                normfluct(j,i,c)=abs(fluct/(medianres+epsilon));
            end
        end
    end
    info1=(sqrt(normfluct(:,:,1).^2+normfluct(:,:,2).^2)>thresh);
    u_filtered(info1==1)=NaN;
    v_filtered(info1==1)=NaN;

    typevector_filtered(isnan(u_filtered))=2;
    typevector_filtered(isnan(v_filtered))=2;
    typevector_filtered(typevector{PIVresult,1}==0)=0; %restores typevector for mask
    
    %Interpolate missing data
    u_filtered=inpaint_nans(u_filtered,4);
    v_filtered=inpaint_nans(v_filtered,4);
    
    u_filt{PIVresult,1}=u_filtered;
    v_filt{PIVresult,1}=v_filtered;
    typevector_filt{PIVresult,1}=typevector_filtered;
    
    
    
    
end



end