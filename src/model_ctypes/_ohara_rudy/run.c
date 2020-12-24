#include "math.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../liblsoda/src/common.h"
#include "../../liblsoda/src/lsoda.h"
#include "./ohara_rudy.h"


#define FLAG_PRINT 0
#define FLAG_EULER_CHAIN 1


int rhs(double t, double *y, double *ydot, void *data) {
    
    double *C = (double *)data, *A = ((double *)data) + C_SIZE;
    computeRates(t, C, ydot, y, A);
    return 0;   
}


int run(double *S, double *C, int n_beats, double t_sampling, double tol,
        double *output) {

    double data[C_SIZE + A_SIZE];
    for (int i = 0; i < C_SIZE; ++i) {
        data[i] = C[i];
    }
    
    double  atol[S_SIZE], rtol[S_SIZE];
	int     neq = S_SIZE;
	
	struct lsoda_opt_t opt = {0};
    opt.ixpr    = 0;
    opt.rtol    = rtol;
    opt.atol    = atol;
    opt.itask   = 1;
    

    for (int i = 0; i < S_SIZE; ++i) {
        rtol[i] = tol;
        atol[i] = tol;
    }
    
    atol[1]     *= 1e-2;  // CaMKt
    atol[2]     *= 1e-5;  // cass
    atol[9]     *= 1e-5;  // cai

    double  i_Stim_Period = data[15];
    int     n_samples   = i_Stim_Period / t_sampling;

    double  t               = 0;
    double  t_out           = 0;
    int     i_out_global    = 1;
    int     ctx_state       = 0;
    
    memcpy(output, S, S_SIZE * sizeof(double));
    	  
    for (int i_beat = 0; i_beat < n_beats; ++i_beat) {

        struct lsoda_context_t ctx = {
		    .function = rhs,
		    .neq = neq,
		    .data = data,
		    .state = 1,
	    };
        lsoda_prepare(&ctx, &opt);
    
        for (int i_out_local = 1; i_out_local <= n_samples; i_out_local++, i_out_global++) {
        
            t_out = i_out_global * t_sampling;
            
            lsoda(&ctx, S, &t, t_out);
            
	        memcpy(output + i_out_global * S_SIZE, S, S_SIZE * sizeof(double));
	        /*
	        if (i_beat == (n_beats - 1)) {
	            output[i_out_local - 1] = S[0];
	        }
	        */
	        
		    if (ctx.state <= 0) {
		        return ctx.state;
		    }
		    
	    }
	    
	    ctx_state = ctx.state;
	    lsoda_free(&ctx);
	}
	
	return ctx_state;
}


/* - - - - - - - - - - - - - - - - - */



int rhs_chain(double t, double *y, double *ydot, void *data) {
    
    int chain_length = *((double *)data);
    
    for (int i = 0; i < chain_length; ++i) {
    
        double *C = ((double *)data) + 1 + i *  (C_SIZE + A_SIZE);
        double *A = ((double *)data) + 1 + i *  (C_SIZE + A_SIZE) + C_SIZE;
        
        double *y_i    = y    + i * S_SIZE;
        double *ydot_i = ydot + i * S_SIZE;
        
        computeRates(t, C, ydot_i, y_i, A);
        
        double I_gap_junc = 0;
        double g_gap_junc = 5.0;
        
        if (i < chain_length - 1) {
            I_gap_junc += -g_gap_junc * (y[(i + 1) * S_SIZE] - y[i * S_SIZE]);
        }    
        if (i > 0) {
            I_gap_junc += -g_gap_junc * (y[(i - 1) * S_SIZE] - y[i * S_SIZE]);
        }    
        ydot_i[0] -= I_gap_junc;
    }
    
    return 0;   
}


int event_chain_switch(double t, double *S_chain, int chain_length,
                       double v_threshold, double t_safe) {

    double mean_abs_diff = 0;
    for (int i = 1; i < chain_length; ++i) {
        mean_abs_diff = fabs(S_chain[i * S_SIZE] - S_chain[(i - 1) * S_SIZE]);
    }
    mean_abs_diff /= chain_length;
    
    return (mean_abs_diff < v_threshold) && (t > t_safe);
}



int euler(double *t, double *y, void *data,
          double dt, double t_out) {
    
    int chain_length = *((double *)data);
    double ydot[chain_length * S_SIZE];
    
    while (*t < t_out) {
    
        rhs_chain(*t, y, ydot, data);
        
        if (*t + dt <= t_out) {
            *t += dt;
        } else {
            dt = t_out - *t;
            *t = t_out;
        }
            
        for (int i = 0; i < chain_length * S_SIZE; ++i) {
            y[i] += dt * ydot[i];
        }
    }
    
    return 0;    
}


int run_chain(double *S, double *C,
              int chain_length, double v_threshold, double t_safe,
              int n_beats, double t_sampling,
              double tol,
              double *output) {
              
    /* CHAIN */
    double S_chain[chain_length * S_SIZE];
    double data_chain[1 + chain_length * (C_SIZE + A_SIZE)]; // first element is chain_length
    
    for (int i = 0; i < chain_length; ++i) {
        memcpy(S_chain + i * S_SIZE, S, S_SIZE * sizeof(double));
        memcpy(data_chain + 1 + i * (C_SIZE + A_SIZE), C, C_SIZE * sizeof(double));
        if (i > 0) {
            data_chain[1 + i * (C_SIZE + A_SIZE) + 14] = 0; // stim_amplitude
        }
    }
    data_chain[0] = chain_length;

    double  atol_chain[chain_length * S_SIZE], rtol_chain[chain_length * S_SIZE];
	int     neq_chain = chain_length * S_SIZE;

	struct lsoda_opt_t opt_chain = {0};
    opt_chain.ixpr    = 0;
    opt_chain.rtol    = rtol_chain;
    opt_chain.atol    = atol_chain;
    opt_chain.itask   = 1;

    for (int j = 0; j < chain_length; ++j) {
        for (int i = 0; i < S_SIZE; ++i) {
            rtol_chain[j * S_SIZE + i] = tol;
            atol_chain[j * S_SIZE + i] = tol;
        }
        
        atol_chain[j * S_SIZE + 1]     *=  1e-5;  // CaMKt
        atol_chain[j * S_SIZE + 2]     *=  1e-4;  // cass
        atol_chain[j * S_SIZE + 9]     *=  1e-4;  // cai
    }
    /* CHAIN -- DONE */
    
    
    /* SINGLE CELL */
    int    index_target = chain_length / 2;
    double *data        = data_chain + 1 + index_target * (C_SIZE + A_SIZE);
    S                   = S_chain + index_target * S_SIZE; // TODO

    double  atol[S_SIZE], rtol[S_SIZE];
	int     neq = S_SIZE;

	struct lsoda_opt_t opt = {0};
    opt.ixpr    = 0;
    opt.rtol    = rtol;
    opt.atol    = atol;
    opt.itask   = 1;

    for (int i = 0; i < S_SIZE; ++i) {
        rtol[i] = tol;
        atol[i] = tol;
    }
    atol[1]     *= 1e-2;  // CaMKt
    atol[2]     *= 1e-5;  // cass
    atol[9]     *= 1e-5;  // cai
    /* SINGLE CELL -- DONE*/
    
    
    double  i_Stim_Period = data[15];
    int     n_samples   = i_Stim_Period / t_sampling;

    double  t               = 0;
    double  t_out           = 0;
    int     i_out_global    = 1;
    int     ctx_state       = 0;
    
    memcpy(output, S, S_SIZE * sizeof(double));
    
    if (FLAG_PRINT) {
        for (int i = 0; i < chain_length; ++i) {
            printf("%1.3f ", S_chain[i * S_SIZE]);
        }
        printf("\n");
        fflush(stdout);
    }
    
    for (int i_beat = 0; i_beat < n_beats; ++i_beat) {
    
        int i_out_local = 1;

        /* CHAIN */
        struct lsoda_context_t ctx_chain = {
		    .function = rhs_chain,
		    .neq = neq_chain,
		    .data = data_chain,
		    .state = 1,
	    };
        lsoda_prepare(&ctx_chain, &opt_chain);
    
        for (; i_out_local <= n_samples; i_out_local++, i_out_global++) {
        
            double t_local = (i_out_local - 1) * t_sampling;
            if (event_chain_switch(t_local, S_chain, chain_length, v_threshold, t_safe) == 1) {
                if (FLAG_PRINT) {
		            printf("# chain -> single cell\n");
		        }
		        break;
		    }
        
            t_out = i_out_global * t_sampling;
            
            if (FLAG_EULER_CHAIN) {
                euler(&t, S_chain, data_chain, 5e-3, t_out);
            } else {
                lsoda(&ctx_chain, S_chain, &t, t_out);
            }
            
	        memcpy(output + i_out_global * S_SIZE, S, S_SIZE * sizeof(double));
	        
            if (FLAG_PRINT) {
                for (int i = 0; i < chain_length; ++i) {
                    printf("%1.3f ", S_chain[i * S_SIZE]);
                }
                printf("\n");
                fflush(stdout);
            }
	        	        
		    if (ctx_chain.state <= 0) {
		        return ctx_chain.state;
		    }
		   
	    }
	    
	    ctx_state = ctx_chain.state;
	    lsoda_free(&ctx_chain);
	    /* CHAIN -- DONE */
	    
	    
	    /* SINGLE CELL */
        struct lsoda_context_t ctx = {
		    .function = rhs,
		    .neq = neq,
		    .data = data,
		    .state = 1,
	    };
        lsoda_prepare(&ctx, &opt);
    
        for (; i_out_local <= n_samples; i_out_local++, i_out_global++) {
        
            t_out = i_out_global * t_sampling;
            
            lsoda(&ctx, S, &t, t_out);
            
	        memcpy(output + i_out_global * S_SIZE, S, S_SIZE * sizeof(double));
	        
            if (FLAG_PRINT) {
                for (int i = 0; i < chain_length; ++i) {
                    printf("%1.3f ", S_chain[i * S_SIZE]);
                }
                printf("\n");
                fflush(stdout);
            }
	                
		    if (ctx.state <= 0) {
		        return ctx.state;
		    }
	    }
	    
	    ctx_state = ctx.state;
	    lsoda_free(&ctx);
	    
	    if (FLAG_PRINT) {
	        printf("# single cell -> DONE\n");
	        fflush(stdout);
	    }
	    
	    /* SINGLE CELL -- DONE */
	    
	    for (int i = 0; i < chain_length; ++i) {
            memcpy(S_chain + i * S_SIZE, S, S_SIZE * sizeof(double));
        }
	    
	}
	
	return ctx_state;
}
            
