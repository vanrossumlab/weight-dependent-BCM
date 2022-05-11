% TODO: 
% - learn to convergence
1;
pkg load geometry % needed for centroid in utils

%pkg load parallel
function x_ar=create_stim(N, K, stimtype, stim_params)
global stimLEFT stimRIGHT 

% Typically assume that N=K 
% Indexing such that x_ar( istim, :) is input vector (row vector)
    switch (stimtype)
    
    case ('N2K2')
        printf('create_stim: forcing K = N=2\n')
        phi = stim_params(1);
         % ugede phi=0.3926
        x_ar(1,:)=[cos(phi),sin(phi)];
        x_ar(2,:)=[sin(phi),cos(phi)];
        return
    
    case ('N2K2b')
        printf('create_stim: forcing K = N=2\n')
        phi1 = stim_params(1);
        phi2 = stim_params(2);
         % ugede phi=0.3926
        x_ar(1,:)=[cos(phi1),sin(phi1)];
        x_ar(2,:)=[sin(phi2),cos(phi2)];
        return
    
    case ('N2quadrant')
        % K angles between 0 .. pi/2, K >=2
        phi=(0:K-1)/(K-1)*pi/2;
        x_ar=[cos(phi); sin(phi)]';
    
    case ('unitvec')
        printf('create_stim: forcing K = N\n')
        bg = stim_params(1)
        x0 = stim_params(2)
        x_ar=x0*eye(n)+bg;
        
    case ('vonmises')
        % should normalize 
        A= stim_params(1);
        wid = stim_params(2); %large wid >wide
        offset = stim_params(3);
      %  printf('create_stim: A= %g, wid= %g \n', A ,wid)
       % K=N; printf('forcing K = N')
        phi0_ar=2*pi*(1:K)/K; % centers
        phi_ar=2*pi*(1:N)/N;
        x_ar= A*exp((cos(phi_ar-phi0_ar'*ones(1,N))-1)/wid);
        x_ar += offset;     
        
    case ('triangular')
        A= stim_params(1);
        wid = stim_params(2); %large wid >wide
        offset = stim_params(3);
             
       % K=N; printf('forcing K = N')
        phi0_ar=2*pi*(1:K)/K; % centers
        phi_ar=2*pi*(1:N)/N;
        x_ar= acos(cos(phi_ar-phi0_ar'*ones(1,N))); % x_ar=0..pi..0
        
        % if wid <1, then lowest point hits zero. wid can be >1 though.
        % first non-zero x_ar has value 2*Pi/N, so single non-zero when wid<2/N
        x_ar= max(1-x_ar/wid/pi,0);
        x_ar += offset;
   
    case ('gauss') % without wrap-around!
        A= stim_params(1);
        wid = stim_params(2); %large wid >wide
        offset = stim_params(3);
        printf('create_stim: A= %g, wid= %g \n', A ,wid)
        
       % K=N; printf('forcing K = N')
        phi0_ar=2*pi*(1:K)/K; % centers
        phi_ar=2*pi*(1:N)/N;
        x_ar= A*exp(-(phi_ar-phi0_ar'*ones(1,N)).^2/2/wid^2);  
        
    case ('gausswrap') % with wrap-around
        A= stim_params(1);
        wid = stim_params(2); %large wid >wide
        offset = stim_params(3);
        printf('create_stim: A= %g, wid= %g \n', A ,wid)
        
       % K=N; printf('forcing K = N')
        phi0_ar=2*pi*(1:K)/K; % centers
        phi_ar=2*pi*(1:N)/N;
        dphi=acos(cos(phi_ar-phi0_ar'*ones(1,N)));
        %wid*= 2*pi;
        x_ar= A*exp(-dphi.^2/2/wid^2);      
    
    case ('vonmisesLR') % return K*2N stimulus
        A= stim_params(1);
        wid = stim_params(2); %large wid >wide
        LR= stim_params(3); % =0 -> none; =1 left; =10 right; =11 both
        bg= stim_params(4); % back ground level when turned off.
        
        leftQ=mod(LR,2); % odd
        rightQ=round(LR/10); % hardcoded. Don't change
        phi0_ar=2*pi*(1:K)/K; % centers
        phi_ar=2*pi*(1:N)/N;
        
        x_arM = A*exp((cos(phi_ar-phi0_ar'*ones(1,N))-1)/wid); 
        if (leftQ && rightQ)
            x_ar=[x_arM     x_arM];
            return;
        elseif (leftQ)
            x_ar=[x_arM     0*x_arM+bg];
            return
        elseif (rightQ)
            x_ar=[0*x_arM+bg x_arM];
            return
        else
            x_ar=0*[x_arM   x_arM]+bg;
            printf('bg only !')
        end     
        
    case ('triangularLR') % return K*2N stimulus
        A= stim_params(1);
        wid = stim_params(2); %large wid >wide
        LR= stim_params(3); % =0 -> none; =1 left; =10 right; =11 both
        bg= stim_params(4); % back ground level when turned off.
        
        leftQ=mod(LR,2); % odd
        rightQ=round(LR/10); % hardcoded. Don't change
        strabQ=(LR>=100)
        phi0_ar=2*pi*(1:K)/K; % centers
        phi_ar=2*pi*(1:N)/N;
        
        #x_arM = A*exp((cos(phi_ar-phi0_ar'*ones(1,N))-1)/wid); 
        
        x_arM= acos(cos(phi_ar-phi0_ar'*ones(1,N)));
        x_arM= max(1-x_arM/wid/pi,0); %wid should be <1 
        
        if (strabQ) #strabismus, for every peak in R, peak in L  can be anywhere.
            # in this case we return N^2 * 2N stimulus
            x_ar=[];
            for k=1:K
                x_ar = [x_ar ; [x_arM shift(x_arM,k)]];
            end
            return;
        end
        
        if (leftQ && rightQ)
            x_ar=[x_arM     x_arM];
            return;
        elseif (leftQ)
            x_ar=[x_arM     0*x_arM+bg];
            return;
        elseif (rightQ)
            x_ar=[0*x_arM+bg x_arM];
            return;
        else
            x_ar=0*[x_arM   x_arM]+bg;
            printf('bg only !')
        end         
           
    % rectangular blocks TO TEST
    case ('rects')
        A= stim_params(1)
        wid = stim_params(2) 
        
       % K=N;
      %  printf('forcing K = N')
        phi0_ar=2*pi*(1:K)/K; % centers
        phi_ar=2*pi*(1:N)/N;
        dphi = acos(cos(phi_ar-phi0_ar'*ones(1,N))); % deal with BCs, out=0..pi
        x_ar=(dphi <= wid/N*pi);
        
    otherwise  
        error('no such stimtype')
    end
end

function [dwmat,th] = get_dwmat(w ,x_ar)  
    % return what weight update for each stim would be.
    global tau_w rule u_inh
    K =size(x_ar)(1);
    N =size(x_ar)(2);
    
    % x_ar and dwmat have dimensions (K,N)
    th=0;
    for k=1:K
        y   = w'*x_ar(k,:)';
        th += y*y/K;
    end
    
    dwmat=zeros(K,N);
    for k=1:K
        x   = x_ar(k,:);
        y   = w'*x_ar(k,:)';
        if (rule==1)                        % Standard BCM
            dwmat(k,:) = x'*y*(y-th);
        elseif (rule==2)                    % variant of wBCM
            dwmat(k,:) = x'*y.*(y-abs(w*th));
        elseif (rule==3)                    % wBCM
            dwmat(k,:) = x'*y*(y-th);
            dwmat(k,:) .*= abs(w)'.^(dwmat(k,:)<0);
        elseif (rule==10)                    % wBCM
            dwmat(k,:) = x'*y*(y-th);
            dwmat(k,:) .*= (w+u_inh)'.^(dwmat(k,:)<0);    
        else
            error('No such rule. stop')
        end 
    end    
    dwmat/=tau_w;
end

function y = gety(w,x)  
    global rectyQ si_y si_x N
    y  = (x+si_x*randn(1,N))*w+si_y*randn();   
    y -= rectyQ*y*(y<0);
    % BUG: NOTE x noise needs to be added to pattern in dw as well
end


%%%%%%%%%%% MAIN ROUTINE %%%%%%%%%%%%%
function [w_snaps,rf_snaps,th_ar]= run(winit, x_ar,ndt, nsnaps=1, th0=0)
    global randomstimseqQ tau_theta tau_w rule si_x  si_y rectwQ snaptimes u_inh
   
    N=length(winit); K=size(x_ar)(1);
    
    if (th0==0)
        th0=N;    % N/K??? 
    end    
    th_ar=th0*ones(1,ndt);
    y_ar=zeros(1,ndt);
    
    % C version only, for older nonC versions see ../old
    % note, there is also another version that takes much less memory (see slowdown code)
    
    w_snaps=zeros(N, nsnaps); rf_snaps=zeros(K,nsnaps);
    
    [w_snaps, rf_snaps, th_ar, snaptimes]=runc(winit,x_ar,y_ar,ndt,nsnaps, tau_w, tau_theta, si_x, si_y, rectwQ, rule, th0, u_inh, w_snaps, rf_snaps, th_ar);
    
end    

function  [w_final,rf_final, successQ] = run_nomem_nosample(x_ar,ndt,th0=0,sd_winit=0.1)
% runs bare bcm that only returns final values, but also inits with random.

    global randomstimseqQ tau_theta tau_w rule si_x  si_y rectwQ rectyQ u_inh
    persistent w_finalfound = sd_winit*randn(1,size(x_ar)(2));
    
    K=size(x_ar)(1);  N=size(x_ar)(2);
    winit=sd_winit*randn(1,N);
    
   % if(rectyQ) winit += 100; end
    winithackQ=0;
    if (winithackQ)
        winit(1)=10+u_inh;
        printf("winit hack in wbcm_utils\n")
    else    
        winit=w_finalfound; % initialize with last found.
        % could also init with sBCM solution
    end
    
    use_sBCMQ=1;
    if (use_sBCMQ && rule != 1) % skip over when called from this function.
        printf("wbcm_utils: using sBCM as winit \n")
        saverule=rule;
        rule=1;
        [winit,rfinit,dum] = run_nomem_nosample(x_ar,ndt,th0=0,sd_winit=0.1);
        figure(1001)
        plot(winit); hold on ; plot(rfinit)
        rule=saverule;
    end
    
    if (th0==0)
        th0=K;
    end    
    
    w_final=zeros(N,1); rf_final=zeros(K,1);
    [w_final, rf_final,th_final]=runc_nomem_nosnap(winit,x_ar,ndt, tau_w, tau_theta, si_x, si_y, rectwQ, rule, th0, u_inh, w_final, rf_final);
    
    successQ=(max(th_final)<1e9);
end  

function [wfound, rffound, wfoundmultpl] = get_uniq_samples(ns,x_ar,ndt,eps)
    nproc=3;
    
    % note UniformOutput=false, means that pararrayfun returns cells. access with w_final_ar{is}
  %  [w_final_cell, rf_final_cell, successQ_cell]=pararrayfun(nproc,@() run_nomem_nosample(x_ar,ndt),1:ns,"UniformOutput", false, "VerboseLevel",1);
    
    [w_final_cell, rf_final_cell, successQ_cell]=arrayfun(@() run_nomem_nosample(x_ar,ndt),1:ns,"UniformOutput", false);
    
    wfound=[]; rffound=[]; 
    wfoundmultpl=zeros(1,ns);
    
    %TODO, make uniqing an option (eps=0)
    for is=1:ns
        if (successQ_cell{is}==0) 
            continue;
        end
        w_final= w_final_cell{is};
        rf_final= rf_final_cell{is};
        
        eps_w=0.01;
        if (max(abs(w_final))<eps_w) 
            continue;
        end
       
        reflectionQ=1;
        [w_final, rf_final] = recenter(w_final, rf_final,reflectionQ);
                
        % how should we set eps?, depends on N and max(w)/ var(w)
         % alternative: calculate full matrix of norms (i<j), then draw histogram, pick reasonable threshold
        [wfoundnew, idx]= addto_uniqlistND(wfound, w_final ,eps, reflectionQ);
        wfoundmultpl(idx)++; % multiplicity
        if (columns(wfoundnew) != columns(wfound))
            rffound=[rffound rf_final];
        end
        wfound = wfoundnew;        
    end
    
    wfoundmultpl= wfoundmultpl(find(wfoundmultpl>0));
    % sort 
    [wfoundmultpl,i] = sort(wfoundmultpl,'descend');
    wfound  = wfound(:,i);
    rffound = rffound(:,i);
       
    % nsout=sum(cell2mat(successQ_cell)); % not needed, follows from size wfound
end    

function [winit,wfin]=phaseplot(wmin,wmax,M=25)
    winit=linspace(wmin,wmax,M);
    [winit1, winit2]  = meshgrid(winit,winit);
    [wfinal1,wfinal2] = pararrayfun(4, @(w1,w2) sim_finalw_only(xar, ndt,[w1,w2]),winit1,winit2);

    wfiN=[reshape(wfinal1,M*M,1), reshape(wfinal2,M*M,1)];
    hist3(wfin)
    figure;
    hist(reshape(wfinal1,M*M,1))

    figure;
    dw1=wfinal1-winit1;
    dw2=wfinal2-winit2;
    quiver(winit1, winit2, dw1, dw2);
end

function [r1,r2]=sort2(r1,r2,eps=1e-2)
% to get fixed points from numerics. Removes duplicates (within r< eps), and sort.
% r1 (r2) are lists of x (y) coord. 
    for i=length(r1):-1:2
        for j=i-1:-1:1
            if (hypot(r1(i)-r1(j), r2(i)-r2(j)) <eps)
                r1(i)=[]; r2(i)=[];
                break; % continue to i-1
            end
        end
    end 
    [r1,idx]=sort(r1); r2=r2(idx);
end    

function [outlist,idx]=addto_uniqlist(inlist,x,y,eps=5e-2) %2D only
% returns updated list, and position of tested element (x,y)
    n=columns(inlist);    
    
    for i=1:n
        if (hypot(x-inlist(1,i),y-inlist(2,i)) < eps)   % norm, 2D only
            outlist=inlist;
            idx=i;
            return
        end    
    end    
    % new uniq point found, append
    outlist=[inlist [x,y]'];
    idx=n+1;
    return
end

function [outlist,idx]=addto_uniqlistND(inlist, x, eps=5e-2, reflectionQ=0)
% returns updated list, and position of tested element x; here for N- dimensions
    n=columns(inlist);  
    
    [dum,maxlocx]= max(x);
    if (reflectionQ)
        flipx=flip(x);
        [dum,maxlocflip]= max(flipx);
        flipx=shift(flipx,maxlocx-maxlocflip); 
    else
        flipx=x;
    end
    
    for i=1:n
        if (norm(x-inlist(:,i)) < eps || norm(flipx-inlist(:,i))<eps )
            outlist=inlist;
            idx=i;
            return
        end    
    end    
    outlist=[inlist x];  % new uniq point found, append
    idx=n+1;
    return
end

function com=comfun(ar) % center of mass of 1D array.
% there are cases where this is not stable.
    com =dot(ar,[1:length(ar)])/sum(ar);
end


function com=comfun2(ar) % center of mass  of ABS 1D array.
    com =dot(abs(ar),[1:length(ar)])/sum(abs(ar));
end



function [wnew,rfnew]=recenter(wsnaps, rfsnaps, reflectionQ=0)
    global x_ar

    useRF_Q=1;
    if (useRF_Q)
        % recenter weights and RF based on RF 
        % does not work when N!=K
        ksnaps=size(wsnaps)(2);
        nhalf=floor(size(wsnaps)(1)/2);
        [dum, maxloc]=max(rfsnaps(:,ksnaps));
        rfnew=shift(rfsnaps,nhalf-maxloc); % need to check when ksnaps >1
        wnew=shift(wsnaps,nhalf-maxloc);
        
        if (reflectionQ)
            [dum, maxloc] = max(rfnew);
            if (comfun(rfnew)>maxloc)
                rfnew=flip(rfnew);
                wnew=flip(wnew);
                % now recenter again....
                [dum, maxloc]=max(rfnew);
                rfnew=shift(rfnew,nhalf-maxloc); % need to check when ksnaps >1
                wnew=shift(wnew,nhalf-maxloc);
            end
        end
    else % use max of weight profile to align.
        %ksnaps=size(wsnaps)(2);
        ksnaps=1;
        nhalf=floor(size(wsnaps)(1)/2);
        [dum, maxloc]=max(wsnaps(:,ksnaps));
        wnew=shift(wsnaps,nhalf-maxloc);
        comfun(wnew);
        
        if (reflectionQ)
            [dum, maxloc] = max(wnew);
            if (comfun2(wnew)>maxloc)
                wnew=flip(wnew);
                % now recenter again....
                [dum, maxloc]=max(wnew);
                wnew=shift(wnew,nhalf-maxloc);
            end
        end
        % when K!=N this needs to be recalculated
        rfnew = x_ar*wnew;
    end
end

function subplotindexed(nrows,ncols,ix,iy)
% plot at subplot coordinates given by ix=1..ncols iy= 1..nrows,
% rather than linear index.
% This need to be followed by a call to 'plot()' to do actual plotting
% no initialization is needed
    subplot(nrows,ncols,ix+(iy-1)*ncols);
end

function [rfL,rfR]= get_rf(wsnap,x_arL,x_arR)
    rfL = x_arL*wsnap;
    rfR = x_arR*wsnap;
end

function [selL,selR, maxL, maxR, OD]=getLRstats(rfL,rfR)
    maxL= max(rfL);
    maxR= max(rfR);
    selL = max(0,min(1,1-mean(rfL)./maxL)); % BCM '82, rectify as mean can be <0
    selR = max(0,min(1,1-mean(rfR)./maxR));
    OD= (maxL-maxR)./(maxL+maxR);
end

function out= int_remap(in,maxo,mino=1)
% maps 1d array to integers between mini and maxi. Useful for plotting etc.
    maxin=max(in); minin=min(in); 
    a=(maxo-mino)/(maxin-minin);
    b=mino-a*minin;
    out= round(in*a+b);
end

function asciiwaitbar(x,maxx,dx=1)
# keep track of dx
    printf("\r %g/%g", x, maxx);
    if (x==maxx)
        printf("\n\n");
    end
end

function nz=zerocrossings(x,eps=0.01)
    n=length(x);
    eps=eps*max(abs(x));
    
    perBCQ=1;
    if (perBCQ)
        x(n+1)=x(1);
    end
    nz=0;
    for i=1: n-1+perBCQ
        %nz+= ( sign(x(i)) != sign(x(i+1)) );
        if (x(i)< -eps && x(i+1)> eps) 
            nz++;
        end    
        if (x(i)> eps && x(i+1)< -eps) 
            nz++;
        end    
    end    
end

