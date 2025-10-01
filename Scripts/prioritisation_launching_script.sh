#!/bin/bash

# Spread only:
gams prioritsation_model.gms --t_max=12 --t_set=12 --NDspr=0.001 --tstep=12 --inpath=Input --outpath=Results --mNo=3 --popdet=0.008 --eff=0.5 --sprfile=spread_2024.txt
 
scp Results/inf_t.txt Results/inf_t0.txt
scp Results/inf_t_det.txt Results/inf_t_det0.txt
scp Results/popul_nt.txt Results/popul_nt0.txt
scp Results/u.txt Results/u0.txt
scp Results/w.txt Results/w0.txt
 
# Heuristic policy - a sequence of 1-period solutions projecting the spread and treatments for 1 period only
gams prioritsation_model.gms --t_max=12 --t_set=2 --tstep=2 --inpath=Input --outpath=Results --mNo=4 --popdet=0.008 --eff=0.5 --budg=150 --sprfile=spread_2024.txt
gams prioritsation_model.gms --t_max=12 --t_set=3 --tstep=3 --inpath=Input --outpath=Results --mNo=4 --popdet=0.008 --eff=0.5 --budg=150 --sprfile=spread_2024.txt
gams prioritsation_model.gms --t_max=12 --t_set=4 --tstep=4 --inpath=Input --outpath=Results --mNo=4 --popdet=0.008 --eff=0.5 --budg=150 --sprfile=spread_2024.txt
gams prioritsation_model.gms --t_max=12 --t_set=5 --tstep=5 --inpath=Input --outpath=Results --mNo=4 --popdet=0.008 --eff=0.5 --budg=150 --sprfile=spread_2024.txt
gams prioritsation_model.gms --t_max=12 --t_set=6 --tstep=6 --inpath=Input --outpath=Results --mNo=4 --popdet=0.008 --eff=0.5 --budg=150 --sprfile=spread_2024.txt
gams prioritsation_model.gms --t_max=12 --t_set=7 --tstep=7 --inpath=Input --outpath=Results --mNo=4 --popdet=0.008 --eff=0.5 --budg=150 --sprfile=spread_2024.txt
gams prioritsation_model.gms --t_max=12 --t_set=8 --tstep=8 --inpath=Input --outpath=Results --mNo=4 --popdet=0.008 --eff=0.5 --budg=150 --sprfile=spread_2024.txt
gams prioritsation_model.gms --t_max=12 --t_set=9 --tstep=9 --inpath=Input --outpath=Results --mNo=4 --popdet=0.008 --eff=0.5 --budg=150 --sprfile=spread_2024.txt
gams prioritsation_model.gms --t_max=12 --t_set=10 --tstep=10 --inpath=Input --outpath=Results --mNo=4 --popdet=0.008 --eff=0.5 --budg=150 --sprfile=spread_2024.txt
gams prioritsation_model.gms --t_max=12 --t_set=11 --tstep=11 --inpath=Input --outpath=Results --mNo=4 --popdet=0.008 --eff=0.5 --budg=150 --sprfile=spread_2024.txt
# Fixing the undetectable infestation pattern:
gams prioritsation_model.gms --t_max=12 --t_set=12 --tstep=11 --NDspr=0.001 --inpath=Input --outpath=Results --mNo=1 --popdet=0.008 --eff=0.5 --budg=150 --sprfile=spread_2024.txt
 
scp Results/inf_t.txt Results/inf_t_1.txt
scp Results/inf_t_det.txt Results/inf_t_det_1.txt
scp Results/popul_nt.txt Results/popul_nt_1.txt
scp Results/treat_dens.txt Results/treat_dens_1.txt
scp Results/u.txt Results/u_1.txt
scp Results/w.txt Results/w_1.txt
 
# Spread with treatments: initialization (a sequence of 1-period solves projecting treatments for one period t and spread over T periods)
gams prioritsation_model.gms --t_max=12 --t_set=12 --tstep=2 --inpath=Input --outpath=Results --mNo=1 --popdet=0.008 --eff=0.5 --budg=150 --sprfile=spread_2024.txt
gams prioritsation_model.gms --t_max=12 --t_set=12 --tstep=3 --inpath=Input --outpath=Results --mNo=1 --popdet=0.008 --eff=0.5 --budg=150 --sprfile=spread_2024.txt
gams prioritsation_model.gms --t_max=12 --t_set=12 --tstep=4 --inpath=Input --outpath=Results --mNo=1 --popdet=0.008 --eff=0.5 --budg=150 --sprfile=spread_2024.txt
gams prioritsation_model.gms --t_max=12 --t_set=12 --tstep=5 --inpath=Input --outpath=Results --mNo=1 --popdet=0.008 --eff=0.5 --budg=150 --sprfile=spread_2024.txt
gams prioritsation_model.gms --t_max=12 --t_set=12 --tstep=6 --inpath=Input --outpath=Results --mNo=1 --popdet=0.008 --eff=0.5 --budg=150 --sprfile=spread_2024.txt
gams prioritsation_model.gms --t_max=12 --t_set=12 --tstep=7 --inpath=Input --outpath=Results --mNo=1 --popdet=0.008 --eff=0.5 --budg=150 --sprfile=spread_2024.txt
gams prioritsation_model.gms --t_max=12 --t_set=12 --tstep=8 --inpath=Input --outpath=Results --mNo=1 --popdet=0.008 --eff=0.5 --budg=150 --sprfile=spread_2024.txt
gams prioritsation_model.gms --t_max=12 --t_set=12 --tstep=9 --inpath=Input --outpath=Results --mNo=1 --popdet=0.008 --eff=0.5 --budg=150 --sprfile=spread_2024.txt
gams prioritsation_model.gms --t_max=12 --t_set=12 --tstep=10 --inpath=Input --outpath=Results --mNo=1 --popdet=0.008 --eff=0.5 --budg=150 --sprfile=spread_2024.txt
# Fixing the undetectable infestation pattern
gams prioritsation_model.gms --t_max=12 --t_set=12 --tstep=11 --NDspr=0.001 --inpath=Input --outpath=Results --mNo=1 --popdet=0.008 --eff=0.5 --budg=150 --sprfile=spread_2024.txt
 
scp Results/inf_t_det.txt Results/inf_t_ini.txt
scp Results/inf_t_det.txt Results/inf_t_det_ini.txt
scp Results/popul_nt.txt Results/popul_nt_ini.txt
scp Results/treat_dens.txt Results/treat_dens_ini.txt
scp Results/u.txt Results/u_.txt
scp Results/w.txt Results/w_.txt
 
# Full problem - spread with treatments
gams prioritsation_model.gms --t_max=12 --t_set=12 --tstep=11 --inpath=Input --outpath=Results --mNo=2 --popdet=0.008 --eff=0.5 --budg=150 --sprfile=spread_2024.txt
# Fixing the undetectable infestation pattern
gams prioritsation_model.gms --t_max=12 --t_set=12 --tstep=11 --NDspr=0.001 --inpath=Input --outpath=Results --mNo=1 --popdet=0.008 --eff=0.5 --budg=150 --sprfile=spread_2024.txt

