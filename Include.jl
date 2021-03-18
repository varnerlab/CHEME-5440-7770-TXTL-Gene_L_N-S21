# setup some paths -
_PATH_TO_ROOT = pwd()
_PATH_TO_SRC = joinpath(_PATH_TO_ROOT, "src")

# initialize me -
import Pkg
Pkg.activate(_PATH_TO_ROOT)
Pkg.instantiate()

# load the packages that are required for this calculation -
using GLPK
using JSON
using DelimitedFiles
using Plots
using PrettyTables
using LinearAlgebra

# load my codes -
# Flux.jl -
path_to_code = joinpath(_PATH_TO_SRC,"Flux.jl")
include("$(path_to_code)")

# Utility.jl -
path_to_code = joinpath(_PATH_TO_SRC,"Problem.jl")
include("$(path_to_code)")

# setup reaction and species string arrays -
reaction_string_array = [
    "1,transcription_initiation_closed,[],G+RX,GO,0,1000.0"                                             ;
    "2,transcription_elongation,[],GO+924*NTP,mRNA+G+RX+1848*Pi,0,1000.0"                               ;
    "3,mRNA_degradation,[],mRNA,924*NMP,0,1000.0"                                                       ;
    "4,translation_initiation,[],mRNA+RL,RLO,0,1000.0"                                                  ;
    "5,translation_elongation,[],RLO+308*AAtRNA+616*GTP,308*tRNA+616*GDP+616*Pi+RL+mRNA+PROT,0,1000.0"  ;
    "6,tRNA_charging,[],AA+tRNA+ATP,AMP+2*Pi+AAtRNA,0,1000.0"                                           ;
    "7,AA_exchange,[],AA[e],AA,-1000.0,1000.0"                                                          ;
    "8,NTP_exchange,[],NTP[e],NTP,-1000.0,1000.0"                                                       ;
    "9,NMP_exchange,[],NMP[e],NMP,-1000.0,1000.0"                                                       ;
    "10,ATP_exchange,[],ATP[e],ATP,-1000.0,1000.0"                                                      ;
    "11,AMP_exchange,[],AMP[e],AMP,-1000.0,1000.0"                                                      ;
    "12,GTP_exchange,[],GTP[e],GTP,-1000.0,1000.0"                                                      ;
    "13,GMP_exchange,[],GDP[e],GDP,-1000.0,1000.0"                                                      ;
    "14,Pi_exchange,[],Pi[e],Pi,-1000.0,1000.0"                                                         ;
    "15,P_exchange,[],PROT[e],PROT,-1000.0,1000.0"                                                      ;
]

# setup the species string array -
species_string_array = [
    "1,AA"  ;
    "2,AAtRNA"  ;
    "3,AMP"     ;
    "4,ATP"     ;
    "5,G"       ;
    "6,GDP"     ;
    "7,GO"      ;
    "8,GTP"     ;
    "9,NMP"     ;
    "10,NTP"    ;
    "11,PROT"   ;
    "12,Pi"     ;
    "13,RL"     ;
    "14,RLO"    ;
    "15,RX"     ;
    "16,mRNA"   ;
    "17,tRNA"   ;
]
