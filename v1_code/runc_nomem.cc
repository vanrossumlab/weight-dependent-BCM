#include <octave/oct.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>
#include <time.h>

#define max(a, b) ((a) > (b)) ? (a) : (b)
#define min(a, b) ((a) < (b)) ? (a) : (b)

// As runc.cc but without the need to store all y(idt) and th(idt)
// Uses ring-buffer; only returns the snapshots

DEFUN_DLD(runc_nomem,args,nargout,"help? no way"){
    int nargin = args.length (); 
   
    int nargs_required=14;                
    if (nargin != nargs_required){
        octave_stdout << "Incorrect number of arguments \n ";
        return octave_value_list ();
    }
    
    bool debugQ=false;
    
    NDArray winit       = args(0).array_value();
    NDArray x_ar        = args(1).array_value(); // dim: K X N ! FLIP BELOW!
    long long int ndt   = args(2).int64_value();
    int nsnaps          = args(3).int_value();
    
    double tau_w        = args(4).double_value();
    double tau_th       = args(5).double_value();
    const double si_x   = args(6).double_value();
    const double si_y   = args(7).double_value();
    
    const int rectwQ    = args(8).int_value(); // unused ATM
    const int rule      = args(9).int_value();
    const double th0    = args(10).int_value();
    const double u      = args(11).int_value(); // FF inihibition strength
    
    bool rectyQ = true;
    
    NDArray w_snaps = args(12).array_value(); // allocate in octave? N  x nsnaps
    NDArray rf_snaps= args(13).array_value();
       
    int N = winit.numel();
    int K = x_ar.numel()/N;
    if (N!=K){
        printf("Runc.cc: K= %i\n",K);
    }
    
    RowVector w(N);
    w=winit;
    RowVector x_current(N); // current stim, with noise
    x_current=0*winit;       // quick way to initialize
    
    if (debugQ){
    octave_stdout << " Runc_ nomem \n";
    octave_stdout << winit.ndims() << " " << x_ar.ndims() << "\n";
    octave_stdout  << ndt << " "<< nsnaps <<" " << tau_w << " "<< tau_th << " " << si_y << " "<< si_x << "\n";
    octave_stdout << w_snaps.ndims() << " " << rf_snaps.ndims() << "\n";
    octave_stdout << th0 ;
    octave_stdout << " \n Runc_ nomem \n";
    }
    
    int snap_logscaleQ=0;
    
    // Setup snapshot times
    // edited to include t=1, note last snap should be at ndt-1
    RowVector snaptimes(nsnaps);
    if (nsnaps==1){
        snaptimes(0)=ndt-1;
    }else{ 
        if(snap_logscaleQ){
            for (int is=0; is <nsnaps; is++){
                double logsnaptime=is/((double)nsnaps-1)*log(ndt-1);
                snaptimes(is)=round(exp(logsnaptime));
                if (debugQ) printf("runc.cc: snapshot %i at time %g \n",is, snaptimes(is));
            }
        }else{
            for (int is=0; is <nsnaps; is++){
                double snaptime=round((is)/((double)nsnaps-1)*ndt)-1;
                snaptimes(is)=max(1.0,round(snaptime));         
                if (debugQ) printf("runc.cc: snapshot %i at time %g \n",is, snaptimes(is));
            }
        }
    }
    
    // setup random permutation. 
    // All K stimuli are presented
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_permutation * p = gsl_permutation_alloc (K);
    
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    gsl_rng_set(r, time(NULL)); // Seed with time
    //gsl_rng_set(r, 1); // fix seed

    gsl_permutation_init (p);
    gsl_ran_shuffle (r, p->data, K, sizeof(size_t));
    
    /////////////////////////////////////////////

    int nhist = K;      // 10*K; // K;
    int skip_amount = 10; // only update theta once in a while. Should be related to tau_w 
    int iskip=0;
    double y_current;
    RowVector y2hist(nhist);
    for (int i=0; i< nhist; i++) y2hist(i)=0.0;
    int ring_counter=0; // runs over past activity for theta-update
    
    double th_running = th0;  // running estimate of th 
    
    int counter=0; // runs over the stimuli
    for (long long int idt=1; idt < ndt ; idt ++){
           
        int js =p->data[counter];
        counter ++;
        counter = counter % K;
                            
        y_current = 0.0;
        if (si_y>0){
            double u = gsl_ran_gaussian (r, si_y);
            y_current += u;
        }
        
        for (int i = 0; i < N; i++){
            x_current(i) = x_ar(js,i);
            if (si_x>0){
                double u = gsl_ran_gaussian (r, si_x);
                x_current(i)  += u;
            }
            y_current += w(i)*x_current(i);
        }
        
        if (rectyQ) y_current= max(y_current, 0.0);
        
        if (rule==1){                        // Standard BCM
            double psi = y_current*(y_current-th_running);
            for (int i=0; i<N; i++){ 
                double dw= x_current(i)*psi;
                w(i) += dw/tau_w;
            }
        }    
        else if (rule==2){   // variant of wBCM, dw = x'*y.*(y-abs(w*th));
            for (int i=0; i<N; i++){ 
                double dw= x_current(i)*y_current*(y_current-abs(w(i)*th_running ));
                w(i) += dw/tau_w;
            } 
        }
        else if (rule==3){                  // wBCM
            double psi= y_current*(y_current-th_running);
            for (int i=0; i<N; i++){ 
                double dw=x_current(i)*psi;
                if (dw<0) dw *= abs(w(i));
                w(i) += dw/tau_w;
            }
        }
        else if (rule==10){     // wiBCM,  (w inihibition)
            double psi= y_current*(y_current-th_running);
            for (int i=0; i<N; i++){ 
                double dw = x_current(i)*psi;
                if (dw<0) dw *= w(i)+u; // if x>=0, this can go outside 'i'loop
                w(i) += dw/tau_w;
                w(i) = max(w(i),-u);
            }
        }
        
        else{
            printf("No such rule. Exiting \n");
            exit(1);
        }
        
        // theta update, only correct when random perm and fixed length
        // not for original decaying update
       
        ring_counter ++;
        ring_counter = ring_counter % nhist; // ring_counter=0...nhist-1
        y2hist(ring_counter) = y_current*y_current;
        
        iskip = (iskip+1) % skip_amount;
        if (iskip ==0){
            double new_th=0.0;
            for (int ihist=0; ihist < nhist ; ihist++){
                if (idt-ihist <0) // deal with start of simulation
                    new_th += th0;
                else
                    new_th += y2hist(ihist);
            }
            th_running = new_th/nhist;
            
            if (th_running> 1e10){
                printf("runc.cc: Diverging theta. Exiting.\n");
                octave_stdout <<  w <<"\n";
                break;
            }
        }
        
        // take snapshots of RF and W
        for (int isnap=0; isnap < nsnaps; isnap++){
            if (idt==(long long int)snaptimes(isnap)){
                 printf("Runc_nomem.cc: snapshot at time %lli (%lli) = %g. Theta= %g \n",idt,ndt,idt/(float)ndt, th_running);
                for (int i=0; i< N; i++)
                    w_snaps(i,isnap)= w(i);
                
                for (int js=0; js<K; js++){
                     double y=0.0; 
                     for (int i = 0; i < N; i++) {
                         y+= w(i)*x_ar(js,i);    // no noise in x here!
                     }  
                     if (rectyQ) y = max(y, 0.0);
                     rf_snaps(js,isnap)=y;
                }
            }
        }
    }
    
    gsl_permutation_free (p);
    gsl_rng_free (r);
    octave_value_list retval;
    retval(0)= octave_value (w_snaps);
    retval(1)= octave_value (rf_snaps);
    retval(2)= octave_value (snaptimes);

    return retval;
}
