clear all
close all
global x_ar tau_w tau_theta rule randomstimseqQ rectwQ rectyQ si_y si_x N K stimLEFT stimRIGHT snaptimes

# fixed number to encode stimtype
stimLEFT=1;
stimRIGHT=10;
stimSTRAB=100;

source('../utils/wbcm_utils.m')
source('../utils/tightfig.m')

source('./nis/sampleimages_mvr.m')
ret=system('make');
if (ret!=0)
    error('Stopping: compilation error')
end

######################################
# rule =1, standard BCM.
# rule =2, theta*|w|, weakly competitive
# rule =3, weight dependent, dw=..*w^(dw>0)

% Rectification and noise options
rectyQ  = 0;
rectwQ  = 0;
si_y    = 0.0% post noise
si_x    = 0.0 % pre noise


% from nis.
% seed random generator
rand('state',0);
% note there is still randomness in .cc stimulus sequence. If set...

simcase=  1
switch simcase

case{0} 
% do nothing, but re-reads functions
    dum=1;

case{1}% images

    set (0, "defaultlinelinewidth",2)
    set (0, "defaultaxeslinewidth",3)

    sk= 16; % patch size
    N= sk^2
    K= 5*N; % # patches
    zeromeanQ=0;
    unitvarQ=0;
    
    ntrials = 4; %12;
    w_final = zeros(ntrials,N);
    tau_w   =  K*1e3
    tau_theta= -1; % should be >>1, and >> n. Not used when  randomstimseqQ=0
    rule    = 10 %1 or 10
    ndt     = 5e6; % for sBCM 5e7 is sufficient
    th0     = 0.1 % N;    % N/K?
    sd_winit= 0.1;

    u_inh=2 %20  % diverge when >~ 100
    
    plotlist=[10,11,12];
    for itrial =1:ntrials
        itrial
        x_ar=sampleimages_mvr(K, sk,zeromeanQ,unitvarQ)';
        % patch is reshape(x_ar(patch_id,:),[sk,sk])
       
        mean_x    = mean(mean(x_ar))
        std_x     = std(reshape(x_ar,1,[]))
       
        nsnaps=40;     
        w_snaps=zeros(N, nsnaps); rf_snaps=zeros(K,nsnaps);
        winit=sd_winit*randn(1,N);

        [w_snaps, rf_snaps, snaptimes]= runc_nomem(winit, x_ar, ndt,nsnaps, tau_w, tau_theta, si_x, si_y, rectwQ, rule, th0, u_inh, w_snaps, rf_snaps);
        
        minw=min(min(w_snaps));
        maxw=max(max(w_snaps));
        
        w_final(itrial,:)=w_snaps(:,nsnaps);
        rf_final(itrial,:) = rf_snaps(:,nsnaps);
        
        if (find(1==plotlist))
            figure(1) % current trial, w snapshots (z normalized)
            pln=ceil(sqrt(nsnaps));
            for isnap=1:nsnaps
                w2=reshape(w_snaps(:,isnap),[sk,sk]);
                subplot(pln,pln,isnap,'align');
                % align supposedly reduces spacing
                %dum=[isnap,min(w_snaps(:,isnap)),max(w_snaps(:,isnap))];        
                colormap gray; 
                imagesc(w2,[minw,maxw])
                axis('off')
            end
            tightfig(1);
            print(1,'fig1.pdf')
        end
        
        
        if (find(22==plotlist))
            %close(22)
            figure(22) % RF snapshots, shows convergence
            for isnap=1:nsnaps
                rf = rf_snaps(:,isnap);
                subplot(pln,pln,isnap,'align');
                plot(rf)
                axis('off')
            % axis('nolabel')
            end
        % tightfig(22);
            print(22,'fig22.pdf')
        end
        
%         figure(2) % response snapshots. See if it converges... Fig 22 works better...
%         for isnap=1:nsnaps
%             rf = rf_snaps(:,isnap);
%             plot(rf)
%             xlabel('input nr.'); ylabel('response')
%             hold on
%         end
        

       % figure(3)
       % plot(rf_snaps(:,nsnaps)) 
       % title('final RF')
        
    end % itrial
    % save / plot
    
    pln2 =ceil(sqrt(ntrials));
    
    minw=min(min(w_final))
    maxw=max(max(w_final))
    
    if (find(10==plotlist))
        figure(10) % final W samples
        for itrial=1:ntrials
            w2=reshape(w_final(itrial,:),[sk,sk]);
            subplot(pln2,pln2,itrial,'align');
            colormap gray;  
            imagesc(w2,[minw,maxw])
            axis( "square");
            axis('off')
        end
        colorbar;
        tightfig(10);
        print(10,'fig10.pdf')
    end
    
    minrf=min(min(rf_final))
    maxfr=max(max(rf_final))
    
    if (find(11==plotlist))
        figure(11) % final RF samples
        for itrial=1:ntrials
            subplot(pln2,pln2,itrial,'align');
            axis('off')
            plot(rf_final(itrial,:))
        end
        %tightfig(11);
        colorbar;
        print(11,'fig11.pdf')
    end
    w_final_flat= reshape(w_final,ntrials*N,1);
    
    
    if (find(12==plotlist))
        figure(12)
        hist(w_final_flat,40);
        print(12,'fig12.pdf')
    end
    save w_final.dat  w_final % dimensions ntrials x (16*16)

    %save rf_final.dat rf_final
    rft=[(1:1280)' , rf_final'];
    save rf_final.dat rft
endswitch
