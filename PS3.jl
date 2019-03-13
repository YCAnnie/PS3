include("Flux.jl")
using LinearAlgebra

#stoichiometric_matrix
stoichiometric_matrix=[[-1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0];
[1 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0];
[1 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0];
[0 0 0 1 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0];
[0 0 -1 0 2 -2 0 0 0 0 0 1 0 0 0 0 0 0 0 0];
[0 0 0 -1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0];
[-1 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0];
[0 1 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0];
[0 0 1 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0];
[-1 0 0 1 1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
[1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
[0 1 -1 0 -1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
[0 0 1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
[0 0 0 0 -1.5 1.5 0 0 0 0 0 0 0 0 0 1 0 0 0 0];
[0 0 0 0 -1.5 1.5 0 0 0 0 0 0 0 0 0 0 1 0 0 0];
[0 0 0 0 -2 2 0 0 0 0 0 0 0 0 0 0 0 1 0 0];
[0 0 0 0 1 -1 0 0 0 0 0 0 0 0 0 0 0 0 -1 0];
[0 0 0 0 1.5 -1.5 0 0 0 0 0 0 0 0 0 0 0 0 0 -1]];

#define variables
E=0.01;#miumol/gDW
theta=1;
k_cat=[203;34.5;249;88.1;13.7;13.7];

#metabolites concentration(miuM)
con_metabolite_miuM=[[0.0000392333154234199*1000000];#ATP
[0.0000423016355975951*1000000];#AMP
[0];#diphos
[0];#phos
[0];#water
[0];#carbamoyl
[0.0149223204061994*1000000];#aspar
[0.000485075251893972*1000000];#fumar
[0];#urea
[0];#citru
[0];#succin
[0.000255461883107884*1000000];#argin
[0.00448942660807892*1000000];#ornith
[0.0000653952192061409*1000000];#NADPH
[0];#H+
[0];#O2
[0];#NO
[0.000028417783438096*1000000]];#NADP+

#convert concentration to miumol/gDW
cell_mass=9.4e-10;
cell_vol=3.1e-12;
waterfration=0.7;
drycell_mass=(1-water_fraction)*cell_mass;
conv=cell_vol/drycell_mass;
con_metabolite_gDW=con_metabolite_miuM*conv;


#Km(miuM)
Km_miuM=[[0.000392333154234199*1000000];#ATP_6345
[0.000154276895439494*1000000];#aspar_6345
[0];#succi_4321
[0.00154608643830104*1000000];#argin_3531
[0];#carbamoyl_2133
[0.0016*1000000];#ornith_2133
[3.49694159840394e-6*1000000];#argin_1141339 v51
[8.94427190999917e-7*1000000];#NADPH_1141339 v51
[0];#citru v52
[0]];#NADP+ v52
Km=Km_miuM*conv; #Km(miumol/gDW)

#default_bounds_array
#for fluxes
default_bounds_array=[[0 k_cat[1]*E*theta*(con_metabolite_gDW[1]/(con_metabolite_gDW[1]+Km[1]))*(con_metabolite_gDW[7]/(con_metabolite_gDW[7]+Km[2]))];
[0 k_cat[2]*E*theta*1];
[0 k_cat[3]*E*theta*(con_metabolite_gDW[12]/(con_metabolite_gDW[12]+Km[4]))];
[0 k_cat[4]*E*theta*(con_metabolite_gDW[13]/(con_metabolite_gDW[13]+Km[6]))];
[0 k_cat[5]*E*theta*(con_metabolite_gDW[12]/(con_metabolite_gDW[12]+Km[7]))*(con_metabolite_gDW[14]/(con_metabolite_gDW[14]+Km[8]))];
[0 k_cat[6]*E*theta*(con_metabolite_gDW[18]/(con_metabolite_gDW[18]+Km[10]))];
#for bs
[0 10*1000/3600];
[0 10*1000/3600];
[0 10*1000/3600];
[0 10*1000/3600];
[0 10*1000/3600];
[0 10*1000/3600];
[0 10*1000/3600];
[0 10*1000/3600];
[0 10*1000/3600];
[0 10*1000/3600];
[0 10*1000/3600];
[0 10*1000/3600];
[0 10*1000/3600];
[0 10*1000/3600]];

#species_bounds_array
species_bounds_array=zeros(Float64,18,2);

#objective_coefficient_array
objective_coefficient_array=[0;0;0;0;0;0;0;0;0;-1;0;0;0;0;0;0;0;0;0;0];

min_flag = true;

ans=calculate_optimal_flux_distribution(stoichiometric_matrix,default_bounds_array,species_bounds_array,objective_coefficient_array,min_flag);
vec=ans[2]*3600;
flux_urea=vec[10];
println(flux_urea)
