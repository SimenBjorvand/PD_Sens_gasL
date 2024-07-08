# Author: Simen Bjorvand
# Date: 05.10.22
# Description:
# This file will build create the model parameter struct 
# require Statistics
struct gasLiftParameters
    Nw::Int64               # number of wells
    Mw::Float64             # Moleweight of gas
    Ta::Array{Float64}      # Temperature in Annulus
    R::Float64              # Gas constant
    ρo::Array{Float64}      # Density of oil
    Tw::Array{Float64}      # Temperature in well tube
    Lr::Array{Float64}      # Length of well above injection
    Dr::Array{Float64}      # Diameter of well above injection 
    Lw::Array{Float64}      # Length of well below injection
    Dw::Array{Float64}      # Diameter of well below injection
    g::Float64              # Gravitational constant
    Ha::Array{Float64}      # Height of Annulus
    La::Array{Float64}      # Length of Annulus
    Da::Array{Float64}      # Diameter of Annulus
    Hw::Array{Float64}      # Height of well tube above injection point
    Hr::Array{Float64}      # Height of well tube above injection point
    μ0::Float64             # 1 cp Oil viscosity
    Civ::Array{Float64}     # Injection valve coefficient
    Cpc::Array{Float64}     # Production valve coefficient
    pr::Array{Float64}      # Reservoir pressure
    PI::Array{Float64}      # Reservoir productivity index
    GOR::Array{Float64}     # Gas-oil ratio
    ρmo::Float64            # Mean oil Density
    TR::Float64             # Riser Temperature
    LR::Float64             # Riser Length
    DR::Float64             # Riser Diameter
    HR::Float64             # Riser Height
    CRH::Float64            # Choke valve coefficient 
    ps::Float64             # Separator pressure
end



function gasLiftPar()

    Nw  = 2                 # [#]
    Mw  = 20.               # [kg/kmol]
    Ta  = [28,28.].+273.15  # [K]
    R   = 8314.463          # [J/K/kmol]
    ρo  = [800,800.]        # [kg/m³]
    Tw  = [32,32.].+273.15  # [K]
    Lr  = [500,500.]        # [m]
    Dr  = [0.121,0.121]     # [m]
    Lw  = [1500,1500.]      # [m]
    Dw  = [0.121,0.121]     # [m]
    g   = 9.81              # [m/s²]
    Ha  = [1000,1000.]      # [m]
    La  = [1500,1500.]      # [m]
    Da  = [0.189,0.189]     # [m]
    Hw  = [1000,1000.]      # [m]
    Hr  = [500,500.]        # [m]
    μ0  = 128*0.001         # [kg/m/s]
    Civ = [1e-4,1e-4]       # [m²]
    Cpc = [2e-3,2e-3]       # [m²]
    pr  = [150,155.]        # [bar]
    PI  = [0.7,0.7]         # [kg/s/bar]
    GOR = [0.1,0.12]        # [kg/kg]
    ρmo = mean(ρo)          # [kg/m³]
    TR  = 30 +273.15        # [K]
    LR  = 500.              # [m]
    DR  = 0.121             # [m]
    HR  = 500.              # [m]
    CRH = 1e-2              # [m²]
    ps  = 20.               # [Bar]


    par = gasLiftParameters(Nw,Mw,Ta,R,ρo,Tw,Lr,Dr,Lw,Dw,g,
                Ha,La,Da,Hw,Hr,μ0,Civ,Cpc,pr,PI,GOR,
                ρmo,TR,LR,DR,HR,CRH,ps)

    return par
end
