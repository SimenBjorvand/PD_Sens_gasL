# Date: 21.10.2022
# Changed. 17.09.2023
# Author: Simen Bjorvand
# Description: A file containing usefull functions for the 
#              decomposed MPCs 


# Required files: JuMP, Ipopt, LinearAlgebra, SparseArrays

# Required self written codes:
# colMatrix.jl, msUtils.jl, myParams.jl



## First function: The subproblem definitions:

# Write two versions, for two different types of initialization
# 1. For initialization with SS data
# 2. For initialization with previous iterations dynamic data (does this in the centralized mpc allready)

function subProbNlp(ps,par,xi,ui,uSS,xSS,zSS)
    col = [0, 0.155051,0.644949,1.0] #colpoint 1 is the zeroth colpoint in this formulation
    colMat,colCont,colObj = colMatrix.collocationMatrix(col) #colMat is square but first row will be ignored   


    αo =  20.    # [$/m³] Oil price per volume
    αgl = 1.2   # [$/m³] Gas lift cost per volume
    wgmax = 12.
    wgmax_track = wgmax - 1.8
    γwmax = 100


    # Simulation parameters
    Nx  = 3*par.Nw+2
    Nu  = par.Nw
    Nz  = 13*par.Nw+7
    Np  = par.Nw
    NCP = 3
    NFE = 24
    NPR = 2^Np + 1
    NRH = 1
    NSC = NPR^NRH
    dt  = 300. #mpcTime/NFE

    γr = 1.2e3


    optimizer = optimizer_with_attributes(Ipopt.Optimizer, 
                        "linear_solver" => "ma57",
                        "print_level" => 1)


    nlp = Model(optimizer)

    # Declare variables:
    #input variables
    @variable(nlp, u[nu=1:Nu,-1:NFE-1], start = uSS[nu]) #input variable
    #Differential state variables
    @variable(nlp, x[nx=1:Nx,1:NFE,0:NCP], start = xSS[nx])
    @variable(nlp, xend[nx=1:Nx], start = xSS[nx])
    #Algebraic state variables:
    @variable(nlp, z[nz=1:Nz,1:NFE,1:NCP], start = zSS[nz])

    #Add the parameter as a variable (for nonlinear programming)
    @variable(nlp, l[nu=1:Nu,1:NRH], start = uSS[nu])

    ## Add parameters to the model: The coordinator variables, parameters and initial conditions
    # Coordinator variables
    @NLparameter(nlp, lp[1:Nu,1:NRH] == 0.)

    # Uncertain parameters:
    @NLparameter(nlp, ΔGOR[nw=1:par.Nw,nrh=1:NRH] == ps[nrh][nw])

    # Initial condition:
    @NLparameter(nlp, x0p[i=1:Nx] == xi[i])
    @NLparameter(nlp, uprev[i=1:Nu] == ui[i])


    # Unpack the decission variables:
    # Inputs:
    wgl = [u[i,:] for i in 1:par.Nw]
    # Dynamic States: 
    mga = [x[i,:,:] for i in 1:par.Nw]
    mgt = [x[i+par.Nw,:,:] for i in 1:par.Nw]
    mot = [x[i+2*par.Nw,:,:] for i in 1:par.Nw]
    mgR = x[3*par.Nw+1,:,:]
    moR = x[3*par.Nw+2,:,:]
    # Algebraic States: 
    ρa  = [z[i,:,:] for i in 1:par.Nw]
    ρm  = [z[i+par.Nw,:,:] for i in 1:par.Nw]
    ρw  = [z[i+2*par.Nw,:,:] for i in 1:par.Nw]
    pa  = [z[i+3*par.Nw,:,:] for i in 1:par.Nw]
    pwh = [z[i+4*par.Nw,:,:] for i in 1:par.Nw]
    pwi = [z[i+5*par.Nw,:,:] for i in 1:par.Nw]
    pbh = [z[i+6*par.Nw,:,:] for i in 1:par.Nw]
    wiv = [z[i+7*par.Nw,:,:] for i in 1:par.Nw]
    wpc = [z[i+8*par.Nw,:,:] for i in 1:par.Nw]
    wpg = [z[i+9*par.Nw,:,:] for i in 1:par.Nw]
    wpo = [z[i+10*par.Nw,:,:] for i in 1:par.Nw]
    wro = [z[i+11*par.Nw,:,:] for i in 1:par.Nw]
    wrg = [z[i+12*par.Nw,:,:] for i in 1:par.Nw]
    ρRm = z[13*par.Nw+1,:,:]
    ρRH = z[13*par.Nw+2,:,:]
    pRh = z[13*par.Nw+3,:,:]
    pm  = z[13*par.Nw+4,:,:]
    wRH = z[13*par.Nw+5,:,:]
    wtg = z[13*par.Nw+6,:,:]
    wto = z[13*par.Nw+7,:,:]



    # Set up Boundary constraints
    # Inputs: 
    for i in 1:par.Nw
        for nfe in NRH:NFE-1
            set_lower_bound(wgl[i][nfe],0.)
            #set_upper_bound(wgl[i],10.)
        end
    end

    # Dynamic States:
    for i in 1:par.Nw
        set_lower_bound.(mga[i][:,0],0.)
        #set_upper_bound(mga[i],10.)

        set_lower_bound.(mgt[i][:,0],0.)
        #set_upper_bound(mgt[i],10.)

        set_lower_bound.(mot[i][:,0],0.)
        #set_upper_bound(mot[i],10.)
    end
    set_lower_bound.(mgR[:,0],0)
    #set_upper_bound(mgR,0)

    set_lower_bound.(moR[:,0],0)
    #set_upper_bound(moR,0)

    # Calculate Areas of the different piping:
    Ar = [par.Dr[nw]^2*pi/4 for nw in 1:par.Nw]
    Aw = [par.Dw[nw]^2*pi/4 for nw in 1:par.Nw]
    Aa = [par.Da[nw]^2*pi/4 - Aw[nw] for nw in 1:par.Nw] 
    AR = par.DR^2*pi/4



    # Set up the model equations:
    # Algebraic Equations:
    for nw in 1:par.Nw
        @NLconstraints(nlp, begin
            [nfe=1:NFE,ncp=1:NCP], ρa[nw][nfe,ncp] - par.Mw*pa[nw][nfe,ncp]*1e5/par.Ta[nw]/par.R == 0

            [nfe=1:NFE,ncp=1:NCP], ρm[nw][nfe,ncp] - (mgt[nw][nfe,ncp] + mot[nw][nfe,ncp])*pwh[nw][nfe,ncp]*par.Mw*par.ρo[nw]*1e5/(mot[nw][nfe,ncp]*pwh[nw][nfe,ncp]*par.Mw*1e5+par.ρo[nw]*par.R*par.Tw[nw]*mgt[nw][nfe,ncp]) == 0

            [nfe=1:NFE,ncp=1:NCP], ρw[nw][nfe,ncp] - (mgt[nw][nfe,ncp] + mot[nw][nfe,ncp]-par.ρo[nw]*par.Lr[nw]*Ar[nw])/(par.Lw[nw]*Aw[nw]) == 0

            [nfe=1:NFE,ncp=1:NCP], 1e-5*(pa[nw][nfe,ncp]*1e5 - (1 + par.g*par.Ha[nw]*par.Mw/par.R/par.Ta[nw])*(par.Ta[nw]*par.R + par.g*par.Ha[nw]*par.Mw)/(par.Mw*par.La[nw]*Aa[nw])*mga[nw][nfe,ncp]) == 0

            [nfe=1:NFE,ncp=1:NCP], 1e-5*(pwh[nw][nfe,ncp]*1e5 - (par.Tw[nw]*par.R/par.Mw)*mgt[nw][nfe,ncp]/(par.Lw[nw]*Aw[nw] + par.Lr[nw]*Ar[nw] - mot[nw][nfe,ncp]/par.ρo[nw]) + (mgt[nw][nfe,ncp] + mot[nw][nfe,ncp])/(2*par.Lw[nw]*Aw[nw])*par.g*par.Hw[nw]) == 0

            [nfe=1:NFE,ncp=1:NCP], 1e-5*(pwi[nw][nfe,ncp]*1e5 - pwh[nw][nfe,ncp]*1e5 - ρw[nw][nfe,ncp]*par.g*par.Hw[nw] - par.μ0*par.Lw[nw]*wpc[nw][nfe,ncp]/(pi*par.Dw[nw]^4*ρm[nw][nfe,ncp])) == 0

            [nfe=1:NFE,ncp=1:NCP], 1e-5*(pbh[nw][nfe,ncp]*1e5 - pwi[nw][nfe,ncp]*1e5 - par.ρo[nw]*par.g*par.Hr[nw] - par.μ0*par.Lr[nw]*wro[nw][nfe,ncp]/(pi*par.Dr[nw]^4*par.ρo[nw])) == 0

            [nfe=1:NFE,ncp=1:NCP], 1e-5*(wiv[nw][nfe,ncp]/par.Civ[nw])^2 - (ρa[nw][nfe,ncp]*(pa[nw][nfe,ncp] - pwi[nw][nfe,ncp])) == 0

            [nfe=1:NFE,ncp=1:NCP], 1e-5*(wpc[nw][nfe,ncp]/par.Cpc[nw])^2 - (ρm[nw][nfe,ncp]*(pwh[nw][nfe,ncp] - pm[nfe,ncp])) == 0

            [nfe=1:NFE,ncp=1:NCP], wpg[nw][nfe,ncp] - mgt[nw][nfe,ncp]/(mgt[nw][nfe,ncp]+mot[nw][nfe,ncp])*wpc[nw][nfe,ncp] == 0

            [nfe=1:NFE,ncp=1:NCP], wpo[nw][nfe,ncp] - wpc[nw][nfe,ncp] + wpg[nw][nfe,ncp] == 0

            [nfe=1:NFE,ncp=1:NCP], wro[nw][nfe,ncp] - par.PI[nw]*(par.pr[nw] - pbh[nw][nfe,ncp]) == 0

            [nfe=1:NRH,ncp=1:NCP], wrg[nw][nfe,ncp] - (par.GOR[nw]+ΔGOR[nw,nfe])*wro[nw][nfe,ncp] == 0

            [nfe=NRH+1:NFE,ncp=1:NCP], wrg[nw][nfe,ncp] - (par.GOR[nw]+ΔGOR[nw,NRH])*wro[nw][nfe,ncp] == 0
        end)
    end

    @NLconstraints(nlp, begin
        [nfe=1:NFE,ncp=1:NCP], ρRm[nfe,ncp] - (mgR[nfe,ncp]+moR[nfe,ncp])*pRh[nfe,ncp]*par.Mw*par.ρmo*1e5/(moR[nfe,ncp]*pRh[nfe,ncp]*par.Mw*1e5 + par.ρmo*par.R*par.TR*mgR[nfe,ncp]) == 0

        [nfe=1:NFE,ncp=1:NCP], ρRH[nfe,ncp] - (mgR[nfe,ncp]+moR[nfe,ncp])/(par.LR*AR) == 0

        [nfe=1:NFE,ncp=1:NCP], 1e-5*(pRh[nfe,ncp]*1e5 - par.TR*par.R*mgR[nfe,ncp]/(par.Mw*par.LR*AR) + ρRH[nfe,ncp]*par.g*par.HR/2) == 0

        [nfe=1:NFE,ncp=1:NCP], 1e-5*(pm[nfe,ncp]*1e5 - pRh[nfe,ncp]*1e5 - ρRH[nfe,ncp]*par.g*par.HR - par.μ0*par.LR*wRH[nfe,ncp]/(pi*par.DR^4*ρRm[nfe,ncp])) == 0
        
        [nfe=1:NFE,ncp=1:NCP], 1e-5*(wRH[nfe,ncp]/par.CRH)^2 - ρRm[nfe,ncp]*(pRh[nfe,ncp]-par.ps) == 0

        [nfe=1:NFE,ncp=1:NCP], wtg[nfe,ncp] - mgR[nfe,ncp]/(mgR[nfe,ncp]+moR[nfe,ncp])*wRH[nfe,ncp] == 0

        [nfe=1:NFE,ncp=1:NCP], wto[nfe,ncp] + wtg[nfe,ncp] - wRH[nfe,ncp] == 0
    end)



    for nw in 1:par.Nw
        @NLconstraints(nlp, begin
            [nfe=1:NFE,ncp=1:NCP], sum(colMat[ncp+1,h+1]*mga[nw][nfe,h] for h in 0:NCP) - dt*(wgl[nw][nfe-1] - wiv[nw][nfe,ncp]) == 0
            [nfe=1:NFE,ncp=1:NCP], sum(colMat[ncp+1,h+1]*mgt[nw][nfe,h] for h in 0:NCP) - dt*(wiv[nw][nfe,ncp] - wpg[nw][nfe,ncp] + wrg[nw][nfe,ncp]) == 0
            [nfe=1:NFE,ncp=1:NCP], sum(colMat[ncp+1,h+1]*mot[nw][nfe,h] for h in 0:NCP) - dt*(wro[nw][nfe,ncp] - wpo[nw][nfe,ncp]) == 0
        end)
    end


    @NLconstraints(nlp, begin
        [nfe=1:NFE,ncp=1:NCP], sum(colMat[ncp+1,h+1]*mgR[nfe,h] for h in 0:NCP) - dt*(sum(wpg[nw][nfe,ncp] for nw in 1:par.Nw) - wtg[nfe,ncp]) == 0
        [nfe=1:NFE,ncp=1:NCP], sum(colMat[ncp+1,h+1]*moR[nfe,h] for h in 0:NCP) - dt*(sum(wpo[nw][nfe,ncp] for nw in 1:par.Nw) - wto[nfe,ncp]) == 0
    end)

    # Multiple shooting/continuity - constraints
    @NLconstraint(nlp, [nx=1:Nx,nfe=2:NFE], x[nx,nfe,0] - sum(colCont[h+1]*x[nx,nfe-1,h] for h in 0:NCP) == 0)
    @NLconstraint(nlp, [nx=1:Nx], x[nx,1,0] - x0p[nx] == 0)
    @NLconstraint(nlp, [nu=1:Nu], u[nu,-1] - uprev[nu] == 0)
    @NLconstraint(nlp, [nx=1:Nx], xend[nx] - sum(colCont[h+1]*x[nx,NFE,h] for h in 0:NCP) == 0)


    # Add inequality constraints:
    # Operational ones:
    #@NLconstraint(nlp, [nfe=1:NFE],  wtg[nfe,NCP] - wgmax <= 0) #Gas production constraint
    #@NLconstraint(nlp, [nu=1:Nu,nfe=1:NFE], u[nu,nfe-1]-u[nu,nfe-2] - 0.8 <= 0) #MV changing inhibition
    #@NLconstraint(nlp, [nu=1:Nu,nfe=1:NFE], u[nu,nfe-2]-u[nu,nfe-1] - 0.8 <= 0) #MV changing inhibition


    # Physical ones: (makes sure the model equation does not brake appart)
    #@NLconstraint(nlp, [nw=1:par.Nw,nfe=1:NFE], pwi[nw][nfe,NCP] - pa[nw][nfe,NCP] <= 0)
    #@NLconstraint(nlp, [nw=1:par.Nw,nfe=1:NFE], pm[nfe,NCP] - pwh[nw][nfe,NCP] <= 0)
    #@NLconstraint(nlp, [nw=1:par.Nw,nfe=1:NFE], pbh[nw][nfe,NCP] - par.pr[nw] <= 0)
    #@NLconstraint(nlp, [nfe=1:NFE], par.ps - pRh[nfe,NCP] <= 0)
    #@NLconstraint(nlp, [nw=1:par.Nw,nfe=1:NFE], -ρw[nw][nfe,NCP] <= 0)

    ### Add parametric/coupling/NAC constraints
    @NLconstraint(nlp, [nu=1:Nu,nrh=1:NRH], u[nu,nrh-1] - l[nu,nrh] == 0) # NAC
    @NLconstraint(nlp, [nu=1:Nu,nrh=1:NRH], l[nu,nrh] - lp[nu,nrh] == 0)  # coupling 



    # Constrained Objective:
    #=@NLobjective(nlp,Min, (1/NSC)*(-αo*sum(sum(sum(wpo[nw][nfe,ncp]*colObj[ncp]*dt for ncp in 1:NCP) for nfe in 1:NFE) for nw in 1:par.Nw) 
                            + αgl*sum(sum(wgl[nw][nfe-1]*dt for nfe in 1:NFE) for nw in 1:par.Nw) 
                            + sum(sum(γr*(u[nu,nfe-1]-u[nu,nfe-2])^2 for nfe in 1:NFE) for nu in 1:Nu)))=#
    # Tracking Objective: 
    @NLobjective(nlp,Min, (1/NSC)*(-αo*sum(sum(sum(wpo[nw][nfe,ncp]*colObj[ncp]*dt for ncp in 1:NCP) for nfe in 1:NFE) for nw in 1:par.Nw) 
                            + αgl*sum(sum(wgl[nw][nfe-1]*dt for nfe in 1:NFE) for nw in 1:par.Nw) 
                            + sum(sum(γr*(u[nu,nfe-1]-u[nu,nfe-2])^2 for nfe in 1:NFE) for nu in 1:Nu)
                            + γwmax*sum(sum(colObj[ncp]*dt*(wtg[nfe,ncp] - wgmax_track)^2  for ncp in 1:NCP) for nfe in 1:NFE)
                            ))

    return nlp
end

function subProbNlp(ps,par,xi,ui,uG,xG,xendG,zG,lG)
    col = [0, 0.155051,0.644949,1.0] #colpoint 1 is the zeroth colpoint in this formulation
    colMat,colCont,colObj = colMatrix.collocationMatrix(col) #colMat is square but first row will be ignored   


    αo =  20.    # [$/m³] Oil price per volume
    αgl = 1.2   # [$/m³] Gas lift cost per volume
    wgmax = 12.
    wgmax_track = wgmax - 1.8
    γwmax = 100


    # Simulation parameters
    Nx  = 3*par.Nw+2
    Nu  = par.Nw
    Nz  = 13*par.Nw+7
    Np  = par.Nw
    NCP = 3
    NFE = 24
    NPR = 2^Np + 1
    NRH = 1
    NSC = NPR^NRH
    dt  = 300. #mpcTime/NFE

    γr = 1.2e3


    optimizer = optimizer = optimizer_with_attributes(Ipopt.Optimizer, 
                        "linear_solver" => "ma57",
                        "print_level" => 1)


    nlp = Model(optimizer)

    # Declare variables:
    #input variables
    @variable(nlp, u[nu=1:Nu,nfe=-1:NFE-1], start = uG[nu,nfe]) #input variable
    #Differential state variables
    @variable(nlp, x[nx=1:Nx,nfe=1:NFE,ncp=0:NCP], start = xG[nx,nfe,ncp])
    @variable(nlp, xend[nx=1:Nx], start = xendG[nx])
    #Algebraic state variables:
    @variable(nlp, z[nz=1:Nz,nfe=1:NFE,ncp=1:NCP], start = zG[nz,nfe,ncp])
    #@variable(nlp, 0 <= ϵ[1:NFE], start = 0)

    #Add the parameter as a variable (for nonlinear programming)
    @variable(nlp, l[nu=1:Nu,nrh=1:NRH], start = lG[nu,nrh])

    ## Add parameters to the model: The coordinator variables, parameters and initial conditions
    # Coordinator variables
    @NLparameter(nlp, lp[1:Nu,1:NRH] == 0.)

    # Uncertain parameters:
    @NLparameter(nlp, ΔGOR[nw=1:par.Nw,nrh=1:NRH] == ps[nrh][nw])

    # Initial condition:
    @NLparameter(nlp, x0p[i=1:Nx] == xi[i])
    @NLparameter(nlp, uprev[i=1:Nu] == ui[i])


    # Unpack the decission variables:
    # Inputs:
    wgl = [u[i,:] for i in 1:par.Nw]
    # Dynamic States: 
    mga = [x[i,:,:] for i in 1:par.Nw]
    mgt = [x[i+par.Nw,:,:] for i in 1:par.Nw]
    mot = [x[i+2*par.Nw,:,:] for i in 1:par.Nw]
    mgR = x[3*par.Nw+1,:,:]
    moR = x[3*par.Nw+2,:,:]
    # Algebraic States: 
    ρa  = [z[i,:,:] for i in 1:par.Nw]
    ρm  = [z[i+par.Nw,:,:] for i in 1:par.Nw]
    ρw  = [z[i+2*par.Nw,:,:] for i in 1:par.Nw]
    pa  = [z[i+3*par.Nw,:,:] for i in 1:par.Nw]
    pwh = [z[i+4*par.Nw,:,:] for i in 1:par.Nw]
    pwi = [z[i+5*par.Nw,:,:] for i in 1:par.Nw]
    pbh = [z[i+6*par.Nw,:,:] for i in 1:par.Nw]
    wiv = [z[i+7*par.Nw,:,:] for i in 1:par.Nw]
    wpc = [z[i+8*par.Nw,:,:] for i in 1:par.Nw]
    wpg = [z[i+9*par.Nw,:,:] for i in 1:par.Nw]
    wpo = [z[i+10*par.Nw,:,:] for i in 1:par.Nw]
    wro = [z[i+11*par.Nw,:,:] for i in 1:par.Nw]
    wrg = [z[i+12*par.Nw,:,:] for i in 1:par.Nw]
    ρRm = z[13*par.Nw+1,:,:]
    ρRH = z[13*par.Nw+2,:,:]
    pRh = z[13*par.Nw+3,:,:]
    pm  = z[13*par.Nw+4,:,:]
    wRH = z[13*par.Nw+5,:,:]
    wtg = z[13*par.Nw+6,:,:]
    wto = z[13*par.Nw+7,:,:]



    # Set up Boundary constraints
    # Inputs: 
    for i in 1:par.Nw
        for nfe in NRH:NFE-1
            set_lower_bound(wgl[i][nfe],0.)
            #set_upper_bound(wgl[i],10.)
        end
    end

    # Dynamic States:
    for i in 1:par.Nw
        set_lower_bound.(mga[i][:,0],0.)
        #set_upper_bound(mga[i],10.)

        set_lower_bound.(mgt[i][:,0],0.)
        #set_upper_bound(mgt[i],10.)

        set_lower_bound.(mot[i][:,0],0.)
        #set_upper_bound(mot[i],10.)
    end
    set_lower_bound.(mgR[:,0],0)
    #set_upper_bound(mgR,0)

    set_lower_bound.(moR[:,0],0)
    #set_upper_bound(moR,0)

    # Calculate Areas of the different piping:
    Ar = [par.Dr[nw]^2*pi/4 for nw in 1:par.Nw]
    Aw = [par.Dw[nw]^2*pi/4 for nw in 1:par.Nw]
    Aa = [par.Da[nw]^2*pi/4 - Aw[nw] for nw in 1:par.Nw] 
    AR = par.DR^2*pi/4



    # Set up the model equations:
    # Algebraic Equations:
    for nw in 1:par.Nw
        @NLconstraints(nlp, begin
            [nfe=1:NFE,ncp=1:NCP], ρa[nw][nfe,ncp] - par.Mw*pa[nw][nfe,ncp]*1e5/par.Ta[nw]/par.R == 0

            [nfe=1:NFE,ncp=1:NCP], ρm[nw][nfe,ncp] - (mgt[nw][nfe,ncp] + mot[nw][nfe,ncp])*pwh[nw][nfe,ncp]*par.Mw*par.ρo[nw]*1e5/(mot[nw][nfe,ncp]*pwh[nw][nfe,ncp]*par.Mw*1e5+par.ρo[nw]*par.R*par.Tw[nw]*mgt[nw][nfe,ncp]) == 0

            [nfe=1:NFE,ncp=1:NCP], ρw[nw][nfe,ncp] - (mgt[nw][nfe,ncp] + mot[nw][nfe,ncp]-par.ρo[nw]*par.Lr[nw]*Ar[nw])/(par.Lw[nw]*Aw[nw]) == 0

            [nfe=1:NFE,ncp=1:NCP], 1e-5*(pa[nw][nfe,ncp]*1e5 - (1 + par.g*par.Ha[nw]*par.Mw/par.R/par.Ta[nw])*(par.Ta[nw]*par.R + par.g*par.Ha[nw]*par.Mw)/(par.Mw*par.La[nw]*Aa[nw])*mga[nw][nfe,ncp]) == 0

            [nfe=1:NFE,ncp=1:NCP], 1e-5*(pwh[nw][nfe,ncp]*1e5 - (par.Tw[nw]*par.R/par.Mw)*mgt[nw][nfe,ncp]/(par.Lw[nw]*Aw[nw] + par.Lr[nw]*Ar[nw] - mot[nw][nfe,ncp]/par.ρo[nw]) + (mgt[nw][nfe,ncp] + mot[nw][nfe,ncp])/(2*par.Lw[nw]*Aw[nw])*par.g*par.Hw[nw]) == 0

            [nfe=1:NFE,ncp=1:NCP], 1e-5*(pwi[nw][nfe,ncp]*1e5 - pwh[nw][nfe,ncp]*1e5 - ρw[nw][nfe,ncp]*par.g*par.Hw[nw] - par.μ0*par.Lw[nw]*wpc[nw][nfe,ncp]/(pi*par.Dw[nw]^4*ρm[nw][nfe,ncp])) == 0

            [nfe=1:NFE,ncp=1:NCP], 1e-5*(pbh[nw][nfe,ncp]*1e5 - pwi[nw][nfe,ncp]*1e5 - par.ρo[nw]*par.g*par.Hr[nw] - par.μ0*par.Lr[nw]*wro[nw][nfe,ncp]/(pi*par.Dr[nw]^4*par.ρo[nw])) == 0

            [nfe=1:NFE,ncp=1:NCP], 1e-5*((wiv[nw][nfe,ncp]/par.Civ[nw])^2 - (ρa[nw][nfe,ncp]*(pa[nw][nfe,ncp] - pwi[nw][nfe,ncp])*1e5)) == 0

            [nfe=1:NFE,ncp=1:NCP], 1e-5*((wpc[nw][nfe,ncp]/par.Cpc[nw])^2 - (ρm[nw][nfe,ncp]*(pwh[nw][nfe,ncp] - pm[nfe,ncp])*1e5)) == 0

            [nfe=1:NFE,ncp=1:NCP], wpg[nw][nfe,ncp] - mgt[nw][nfe,ncp]/(mgt[nw][nfe,ncp]+mot[nw][nfe,ncp])*wpc[nw][nfe,ncp] == 0

            [nfe=1:NFE,ncp=1:NCP], wpo[nw][nfe,ncp] - wpc[nw][nfe,ncp] + wpg[nw][nfe,ncp] == 0

            [nfe=1:NFE,ncp=1:NCP], wro[nw][nfe,ncp] - par.PI[nw]*(par.pr[nw] - pbh[nw][nfe,ncp]) == 0

            [nfe=1:NRH,ncp=1:NCP], wrg[nw][nfe,ncp] - (par.GOR[nw]+ΔGOR[nw,nfe])*wro[nw][nfe,ncp] == 0
            [nfe=NRH+1:NFE,ncp=1:NCP], wrg[nw][nfe,ncp] - (par.GOR[nw]+ΔGOR[nw,NRH])*wro[nw][nfe,ncp] == 0
        end)
    end

    @NLconstraints(nlp, begin
        [nfe=1:NFE,ncp=1:NCP], ρRm[nfe,ncp] - (mgR[nfe,ncp]+moR[nfe,ncp])*pRh[nfe,ncp]*par.Mw*par.ρmo*1e5/(moR[nfe,ncp]*pRh[nfe,ncp]*par.Mw*1e5 + par.ρmo*par.R*par.TR*mgR[nfe,ncp]) == 0

        [nfe=1:NFE,ncp=1:NCP], ρRH[nfe,ncp] - (mgR[nfe,ncp]+moR[nfe,ncp])/(par.LR*AR) == 0

        [nfe=1:NFE,ncp=1:NCP], 1e-5*(pRh[nfe,ncp]*1e5 - par.TR*par.R*mgR[nfe,ncp]/(par.Mw*par.LR*AR) + ρRH[nfe,ncp]*par.g*par.HR/2) == 0

        [nfe=1:NFE,ncp=1:NCP], 1e-5*(pm[nfe,ncp]*1e5 - pRh[nfe,ncp]*1e5 - ρRH[nfe,ncp]*par.g*par.HR - par.μ0*par.LR*wRH[nfe,ncp]/(pi*par.DR^4*ρRm[nfe,ncp])) == 0
        
        [nfe=1:NFE,ncp=1:NCP], 1e-5*((wRH[nfe,ncp]/par.CRH)^2 - ρRm[nfe,ncp]*(pRh[nfe,ncp]-par.ps)*1e5) == 0

        [nfe=1:NFE,ncp=1:NCP], wtg[nfe,ncp] - mgR[nfe,ncp]/(mgR[nfe,ncp]+moR[nfe,ncp])*wRH[nfe,ncp] == 0

        [nfe=1:NFE,ncp=1:NCP], wto[nfe,ncp] + wtg[nfe,ncp] - wRH[nfe,ncp] == 0
    end)



    for nw in 1:par.Nw
        @NLconstraints(nlp, begin
            [nfe=1:NFE,ncp=1:NCP], sum(colMat[ncp+1,h+1]*mga[nw][nfe,h] for h in 0:NCP) - dt*(wgl[nw][nfe-1] - wiv[nw][nfe,ncp]) == 0
            [nfe=1:NFE,ncp=1:NCP], sum(colMat[ncp+1,h+1]*mgt[nw][nfe,h] for h in 0:NCP) - dt*(wiv[nw][nfe,ncp] - wpg[nw][nfe,ncp] + wrg[nw][nfe,ncp]) == 0
            [nfe=1:NFE,ncp=1:NCP], sum(colMat[ncp+1,h+1]*mot[nw][nfe,h] for h in 0:NCP) - dt*(wro[nw][nfe,ncp] - wpo[nw][nfe,ncp]) == 0
        end)
    end


    @NLconstraints(nlp, begin
        [nfe=1:NFE,ncp=1:NCP], sum(colMat[ncp+1,h+1]*mgR[nfe,h] for h in 0:NCP) - dt*(sum(wpg[nw][nfe,ncp] for nw in 1:par.Nw) - wtg[nfe,ncp]) == 0
        [nfe=1:NFE,ncp=1:NCP], sum(colMat[ncp+1,h+1]*moR[nfe,h] for h in 0:NCP) - dt*(sum(wpo[nw][nfe,ncp] for nw in 1:par.Nw) - wto[nfe,ncp]) == 0
    end)

    # Multiple shooting/continuity - constraints
    @NLconstraint(nlp, [nx=1:Nx,nfe=2:NFE], x[nx,nfe,0] - sum(colCont[h+1]*x[nx,nfe-1,h] for h in 0:NCP) == 0)
    @NLconstraint(nlp, [nx=1:Nx], x[nx,1,0] - x0p[nx] == 0)
    @NLconstraint(nlp, [nu=1:Nu], u[nu,-1] - uprev[nu] == 0)
    @NLconstraint(nlp, [nx=1:Nx], xend[nx] - sum(colCont[h+1]*x[nx,NFE,h] for h in 0:NCP) == 0)


    # Add inequality constraints:
    # Operational ones:
    #@NLconstraint(nlp, [nfe=1:NFE],  wtg[nfe,NCP] - wgmax - <= 0) #ϵ[nfe] <= 0) #Gas production constraint
    #@NLconstraint(nlp, [nu=1:Nu,nfe=1:NFE], u[nu,nfe-1]-u[nu,nfe-2] - 0.8 <= 0) #MV changing inhibition
    #@NLconstraint(nlp, [nu=1:Nu,nfe=1:NFE], u[nu,nfe-2]-u[nu,nfe-1] - 0.8 <= 0) #MV changing inhibition


    # Physical ones: (makes sure the model equation does not brake appart) #! Don't add, they cause LICQ to fail. 
    #@NLconstraint(nlp, [nw=1:par.Nw,nfe=1:NFE], pwi[nw][nfe,NCP] - pa[nw][nfe,NCP] <= 0)
    #@NLconstraint(nlp, [nw=1:par.Nw,nfe=1:NFE], pm[nfe,NCP] - pwh[nw][nfe,NCP] <= 0)
    #@NLconstraint(nlp, [nw=1:par.Nw,nfe=1:NFE], pbh[nw][nfe,NCP] - par.pr[nw] <= 0)
    #@NLconstraint(nlp, [nfe=1:NFE], par.ps - pRh[nfe,NCP] <= 0)
    #@NLconstraint(nlp, [nw=1:par.Nw,nfe=1:NFE], -ρw[nw][nfe,NCP] <= 0)

    ### Add parametric/coupling/NAC constraints
    @NLconstraint(nlp, [nu=1:Nu,nrh=1:NRH], u[nu,nrh-1] - l[nu,nrh] == 0) # NAC
    @NLconstraint(nlp, [nu=1:Nu,nrh=1:NRH], l[nu,nrh] - lp[nu,nrh] == 0)  # coupling 


    # Constrained Objective:
    #=@NLobjective(nlp,Min, (1/NSC)*(-αo*sum(sum(sum(wpo[nw][nfe,ncp]*colObj[ncp]*dt for ncp in 1:NCP) for nfe in 1:NFE) for nw in 1:par.Nw) 
                            + αgl*sum(sum(wgl[nw][nfe-1]*dt for nfe in 1:NFE) for nw in 1:par.Nw) 
                            + sum(sum(γr*(u[nu,nfe-1]-u[nu,nfe-2])^2 for nfe in 1:NFE) for nu in 1:Nu)
                            ))=#
    # Tracking Objective:
    @NLobjective(nlp,Min, (1/NSC)*(-αo*sum(sum(sum(wpo[nw][nfe,ncp]*colObj[ncp]*dt for ncp in 1:NCP) for nfe in 1:NFE) for nw in 1:par.Nw) 
                            + αgl*sum(sum(wgl[nw][nfe-1]*dt for nfe in 1:NFE) for nw in 1:par.Nw) 
                            + sum(sum(γr*(u[nu,nfe-1]-u[nu,nfe-2])^2 for nfe in 1:NFE) for nu in 1:Nu)
                            + γwmax*sum(sum(colObj[ncp]*dt*(wtg[nfe,ncp] - wgmax_track)^2  for ncp in 1:NCP) for nfe in 1:NFE)
                            ))

    return nlp
end

### Functions related to solving the nlps and updating them for the next iteration
function solveNlp(nlp)
    optimize!(nlp)

    # Staus on the solution
    exit_flag = termination_status(nlp)
    local_sol = exit_flag == JuMP.MOI.TerminationStatusCode(4) # LOCALLY_SOLVED
    sol_time = solve_time(nlp)

    return sol_time,local_sol,exit_flag
end

# Write a function which takes takes the nlp and gives better start guesses for next iteration
function updateGuessNlp(nlp)
    # Get the previous iterations solution
    uv = value.(nlp[:u])
    xv = value.(nlp[:x])
    xendv = value.(nlp[:xend])
    zv = value.(nlp[:z])
    lv = value.(nlp[:l])

    # Update the initial guesses
    set_start_value.(nlp[:u],uv)
    set_start_value.(nlp[:x],xv)
    set_start_value.(nlp[:xend],xendv)
    set_start_value.(nlp[:z],zv)
    set_start_value.(nlp[:l],lv)
end


#################################################################################################
# New help function: Feasible initial condition:
#* Don't think i need this one, [Due to slack variables]
#? Instead add a max step length in my cordinator step in addition to 


#################################################################################################
# New help functions:
# Do the NLP sensitivity ones:

function nlpSensitivity(nlp::Model,
                        p::Array{Float64,1},
                        ω::Array{Float64,1},
                        λ::Array{Float64,1},
                        νₗ::Array{Float64,1},
                        νᵤ::Array{Float64,1},
                        γ::Array{Float64,1})

    # Dimensions of the problem
    Nω = length(ω)
    Np = length(p)
    Nh = length(λ)

    #Divide the NLconstraints into equality and inequality constraint set 
    #c_index_set = [i for i in 1:Nh if nlp.nlp_data.nlconstr[i].lb == nlp.nlp_data.nlconstr[i].ub]#! old syntax
    c_index_set = [i for i in 1:Nh if nlp.nlp_model.constraints[nlp.nlp_model.constraints.keys[i]].set == JuMP.MOI.EqualTo{Float64}(0.0)]
    g_index_set = setdiff(1:Nh,c_index_set)

    Nc = length(c_index_set)
    Ng = Nh - Nc


    #Find the boundary constriant index sets
    gₗ_index_set = [i.value for i in JuMP.MOI.get(nlp,JuMP.MOI.ListOfConstraintIndices{JuMP.MOI.VariableIndex,JuMP.MOI.GreaterThan{Float64}}())] 
    gᵤ_index_set = [i.value for i in JuMP.MOI.get(nlp,JuMP.MOI.ListOfConstraintIndices{JuMP.MOI.VariableIndex,JuMP.MOI.LessThan{Float64}}())]

    Ngₗ = length(gₗ_index_set)
    Ngᵤ = length(gᵤ_index_set)

    # Find the parameter dummy variable selection matrix
    Sₚ = [spzeros(Np,Nω-Np) I(Np)]
    cₚ = ω[end-Np+1:end] - p # should be 0/or ϵ = precision of solution


    # Set up a NLPevaluator
    NLPeval = JuMP.NLPEvaluator(nlp)
    JuMP.MOI.initialize(NLPeval,[:Grad,:Jac,:Hess])

    # Evaluate NLconstraints 
    hold = Array{Float64,1}(undef,Nh+Np)
    MOI.eval_constraint(NLPeval,hold,ω)
    h = hold[1:Nh]

    c = h[c_index_set]
    g = h[g_index_set]

    #λc = λ[c_index_set]
    λg = λ[g_index_set]

    # Evaluate the boundary constraints
    gₗ = lower_bound.(JuMP.all_variables(nlp)[gₗ_index_set]) - ω[gₗ_index_set]
    gᵤ = ω[gᵤ_index_set] - upper_bound.(JuMP.all_variables(nlp)[gᵤ_index_set])

    # Evaluate the gradient of the objective function
    ∇F = Array{Float64,1}(undef,Nω)
    MOI.eval_objective_gradient(NLPeval,∇F,ω)

    # Evaluate the Jacobian of the NLconstraints
    J_struckt = MOI.jacobian_structure(NLPeval)
    J_h = Array{Float64,1}(undef,length(J_struckt))
    MOI.eval_constraint_jacobian(NLPeval,J_h,ω)

    ∇h = sparse(last.(J_struckt),first.(J_struckt),J_h)[:,1:Nh] #The reason for calling last before first is because we want the gradient and not the jacobian
    ∇c = ∇h[:,c_index_set]
    #∇g = ∇h[:,g_index_set]

    # Evaluate the Hessian 
    H_struckt = MOI.hessian_lagrangian_structure(NLPeval)
    H_a = Array{Float64,1}(undef,length(H_struckt))
    MOI.eval_hessian_lagrangian(NLPeval,H_a,ω,1.0,[λ;γ])

    push!(H_struckt,(Nω,Nω))
    push!(H_a,0)
    H = sparse(Symmetric(sparse(first.(H_struckt),last.(H_struckt),H_a),:L))



    # Calculate the gradient of the hessian_lagrangian #? Might not be necessary
    #∇L = ∇F + ∇h*λ + Sᵤ'*νᵤ - Sₗ'*νₗ + Sₚ'*γ
    ∇L = ∇F + ∇h*λ
    ∇L[gₗ_index_set] -= νₗ
    ∇L[gᵤ_index_set] += νᵤ
    ∇L[end-Np+1:end] += γ

    # Calculte the optimality residual
    opt_res = [norm(∇L,Inf),
        norm(c,Inf),
        norm(cₚ,Inf),
        norm(min.(-min.(g,0),max.(λg,0)),Inf),
        norm(min.(-gₗ,νₗ),Inf),
        norm(min.(-gᵤ,νᵤ),Inf)]

    η = maximum(opt_res)


    if η < 1e-4
        # Divide the NLconstraints into strongly active constraint set and other
        g₊_index_set = [i for i in g_index_set if λ[i] > η]
        g₀_index_set = setdiff(g_index_set,g₊_index_set)
        g₋_index_set = [i for i in g₀_index_set if h[i] > -η] #* Probably never going to see this one

        ∇g₊ = ∇h[:,g₊_index_set]

        Ng₊ = length(g₊_index_set)
        Ng₀ = Ng - Ng₊

        # Divide boundry constraints into active and inactive sets

        gₗ₊_index_set = [i for i in 1:Ngₗ if νₗ[i] > η]
        gₗ₀_index_set = setdiff(1:Ngₗ,gₗ₊_index_set)
        gₗ₋_index_set = [i for i in gₗ₀_index_set if gₗ[i] > -η] #* Probably never going to see this one

        ∇gₗ₊ = -I(Nω)[:,gₗ_index_set[gₗ₊_index_set]]

        Ngₗ₊ = length(gₗ₊_index_set)
        Ngₗ₀ = Ngₗ - Ngₗ₊

        gᵤ₊_index_set = [i for i in 1:Ngᵤ if νᵤ[i] > η]
        gᵤ₀_index_set = setdiff(1:Ngᵤ,gᵤ₊_index_set)
        gᵤ₋_index_set = [i for i in gᵤ₀_index_set if gᵤ[i] > -η] #* Probably never going to see this one

        ∇gᵤ₊ = I(Nω)[:,gᵤ_index_set[gᵤ₊_index_set]]

        Ngᵤ₊ = length(gᵤ₊_index_set)
        Ngᵤ₀ = Ngᵤ - Ngᵤ₊
    else
        # Divide the Inequalities into B-active and B-inactive constraints
        # Divide the NLconstraints into the 2 sets
        g₊_index_set = [i for i in g_index_set if λ[i] + h[i] > 0]
        g₋_index_set = setdiff(g_index_set,g₊_index_set)

        g₊  = h[g₊_index_set]
        ∇g₊ = ∇h[:,g₊_index_set]
        λg₊ = λ[g₊_index_set]

        Ng₊ = length(g₊_index_set)
        Ng₋ = Ng - Ng₊

        # Divide boundry constraints into the 2 sets 

        gₗ₊_index_set = [i for i in 1:Ngₗ if νₗ[i] + gₗ[i] > 0]
        gₗ₋_index_set = setdiff(1:Ngₗ,gₗ₊_index_set)

        gₗ₊ = gₗ[gₗ₊_index_set]
        ∇gₗ₊ = -I(Nω)[:,gₗ_index_set[gₗ₊_index_set]]
        νₗ₊ = νₗ[gₗ₊_index_set]

        Ngₗ₊ = length(gₗ₊_index_set)
        Ngₗ₋ = Ngₗ - Ngₗ₊

        gᵤ₊_index_set = [i for i in 1:Ngᵤ if νᵤ[i] + gᵤ[i] > 0]
        gᵤ₋_index_set = setdiff(1:Ngᵤ,gᵤ₊_index_set)

        gᵤ₊ = gᵤ[gᵤ₊_index_set]
        ∇gᵤ₊ = I(Nω)[:,gᵤ_index_set[gᵤ₊_index_set]]
        νᵤ₊ = νᵤ[gᵤ₊_index_set]

        Ngᵤ₊ = length(gᵤ₊_index_set)
        Ngᵤ₋ = Ngᵤ - Ngᵤ₊
    end


    #* NEW SECTION
    # Set up/build the KKT matrix for the sensitivity calculation:
    MXY = [∇c ∇g₊ ∇gₗ₊ ∇gᵤ₊ Sₚ']
    nx,ny = size(MXY)

    M = [H MXY;MXY' zeros(ny,ny)]
    #dM = Array(M)

    #LDLTM = LinearAlgebra.lu(M)
    #LDLMA57M = Ma57(M)
    #ma57_factorize(LDLMA57M)

    # Find the gradient of kkt conditions with respect to the parameters
    ∇ₚρ = Array([zeros(nx+ny-Np,Np);-I(Np)])

    # Solve the linear system of equations to obtain the nlp sensitivity with respect to the primal and dual all_variables
    #∇ₚz =  -(LDLTM \ ∇ₚρ)
    #∇ₚz = -ma57_solve(LDLMA57M,∇ₚρ)
    #∇ₚz = -inv(M)*∇ₚρ
    timeFac = @elapsed Mf = lu(M)
    timeBac = @elapsed ∇ₚz = -(Mf\∇ₚρ)
    #∇ₚz = -M\∇ₚρ


    #unpack the z vector
    ∇ₚω   = ∇ₚz[1:Nω,:]
    ∇ₚλc  = ∇ₚz[Nω+1:end,:][1:Nc,:]
    ∇ₚλg₊ = ∇ₚz[Nω+1:end,:][Nc+1:end,:][1:Ng₊,:]
    ∇ₚνₗ₊ = ∇ₚz[Nω+1:end,:][Nc+1:end,:][Ng₊+1:end,:][1:Ngₗ₊,:]
    ∇ₚνᵤ₊ = ∇ₚz[Nω+1:end,:][Nc+1:end,:][Ng₊+1:end,:][Ngₗ₊+1:end,:][1:Ngᵤ₊,:]
    ∇ₚγ   = ∇ₚz[Nω+1:end,:][Nc+1:end,:][Ng₊+1:end,:][Ngₗ₊+1:end,:][Ngᵤ₊+1:end,:]


    # Further modify the the results
    ∇ₚλ = zeros(Nh,Np)
    ∇ₚλ[c_index_set,:]  = ∇ₚλc
    ∇ₚλ[g₊_index_set,:] = ∇ₚλg₊

    ∇ₚνₗ = zeros(Ngₗ,Np)
    ∇ₚνₗ[gₗ₊_index_set,:] = ∇ₚνₗ₊

    ∇ₚνᵤ = zeros(Ngᵤ,Np)
    ∇ₚνᵤ[gᵤ₊_index_set,:] = ∇ₚνᵤ₊

    return ∇ₚω,∇ₚλ,∇ₚνₗ,∇ₚνᵤ,∇ₚγ,timeFac,timeBac
end

function nlpSensitivity(nlp::Model,p::Array{Float64,1})
    lb_MOI_index = JuMP.MOI.get(nlp,JuMP.MOI.ListOfConstraintIndices{JuMP.MOI.VariableIndex,JuMP.MOI.GreaterThan{Float64}}())#* New syntax
    ub_MOI_index = JuMP.MOI.get(nlp,JuMP.MOI.ListOfConstraintIndices{JuMP.MOI.VariableIndex,JuMP.MOI.LessThan{Float64}}())#* New syntax

    Np = length(p)
    ωS = JuMP.all_variables(nlp)
    #Nω = length(ωS)

    ω = [JuMP.MOI.get(nlp,JuMP.MOI.VariablePrimal(),i) for i in ωS]

    λAll = -JuMP.MOI.get(nlp,JuMP.MOI.NLPBlockDual())
    λ = λAll[1:end-Np]
    γ = λAll[end-Np+1:end]

    if length(lb_MOI_index)>0
        νₗ = JuMP.MOI.get(nlp.moi_backend,JuMP.MOI.ConstraintDual(),lb_MOI_index)
    else
        νₗ = Array{Float64,1}(undef,0)
    end

    if length(ub_MOI_index)>0
        νᵤ = -JuMP.MOI.get(nlp.moi_backend,JuMP.MOI.ConstraintDual(),ub_MOI_index)
    else
        νᵤ = Array{Float64,1}(undef,0)
    end

    ∇ω,∇λ,∇νₗ,∇νᵤ,∇γ,timeFac,timeBac = nlpSensitivity(nlp,p,ω,λ,νₗ,νᵤ,γ)

    return ω,∇ω,λ,∇λ,νₗ,∇νₗ,νᵤ,∇νᵤ,γ,∇γ,timeFac,timeBac
end

function nlpSensitivity(nlp::Model,Np::Int64) # Finds the sensitivity of the current point
    p0 = nlp.nlp_model.parameters[1:Np] 

    ω,∇ω,λ,∇λ,νₗ,∇νₗ,νᵤ,∇νᵤ,γ,∇γ,timeFac,timeBac = nlpSensitivity(nlp,p0)

    return ω,∇ω,λ,∇λ,νₗ,∇νₗ,νᵤ,∇νᵤ,γ,∇γ,timeFac,timeBac
end


##################################################
##################################################
# Algorithm for predictor and corrector steps. 
function nlpPredictor(nlp::Model,
                    p::Array{Float64,1},
                    dp::Array{Float64,1},
                    ω::Array{Float64,1},
                    λ::Array{Float64,1},
                    νₗ::Array{Float64,1},
                    νᵤ::Array{Float64,1},
                    γ::Array{Float64,1})

    # Dimensions of the problem
    Nω = length(ω)
    Np = length(p)
    Nh = length(λ)

    #Divide the NLconstraints into equality and inequality constraint set 
    c_index_set = [i for i in 1:Nh if nlp.nlp_model.constraints[nlp.nlp_model.constraints.keys[i]].set == JuMP.MOI.EqualTo{Float64}(0.0)]
    g_index_set = setdiff(1:Nh,c_index_set)

    Nc = length(c_index_set)
    Ng = Nh - Nc


    #Find the boundary constriant index sets
    gₗ_index_set = [i.value for i in JuMP.MOI.get(nlp,JuMP.MOI.ListOfConstraintIndices{JuMP.MOI.VariableIndex,JuMP.MOI.GreaterThan{Float64}}())] 
    gᵤ_index_set = [i.value for i in JuMP.MOI.get(nlp,JuMP.MOI.ListOfConstraintIndices{JuMP.MOI.VariableIndex,JuMP.MOI.LessThan{Float64}}())]

    Ngₗ = length(gₗ_index_set)
    Ngᵤ = length(gᵤ_index_set)

    # Find the parameter dummy variable selection matrix
    Sₚ = [spzeros(Np,Nω-Np) I(Np)]
    cₚ = ω[end-Np+1:end] - (p+dp) # should be 0/or ϵ = precision of solution


    # Set up a NLPevaluator
    NLPeval = JuMP.NLPEvaluator(nlp)
    JuMP.MOI.initialize(NLPeval,[:Grad,:Jac,:Hess])

    # Evaluate NLconstraints 
    hold = Array{Float64,1}(undef,Nh+Np)
    MOI.eval_constraint(NLPeval,hold,ω)
    h = hold[1:Nh]

    c = h[c_index_set]
    g = h[g_index_set]

    #λc = λ[c_index_set]
    λg = λ[g_index_set]

    # Evaluate the boundary constraints
    gₗ = lower_bound.(JuMP.all_variables(nlp)[gₗ_index_set]) - ω[gₗ_index_set]
    gᵤ = ω[gᵤ_index_set] - upper_bound.(JuMP.all_variables(nlp)[gᵤ_index_set])

    # Evaluate the gradient of the objective function
    ∇F = Array{Float64,1}(undef,Nω)
    MOI.eval_objective_gradient(NLPeval,∇F,ω)

    # Evaluate the Jacobian of the NLconstraints
    J_struckt = MOI.jacobian_structure(NLPeval)
    J_h = Array{Float64,1}(undef,length(J_struckt))
    MOI.eval_constraint_jacobian(NLPeval,J_h,ω)

    ∇h = sparse(last.(J_struckt),first.(J_struckt),J_h)[:,1:Nh] #The reason for calling last before first is because we want the gradient and not the jacobian
    ∇c = ∇h[:,c_index_set]
    #∇g = ∇h[:,g_index_set]

    # Evaluate the Hessian 
    H_struckt = MOI.hessian_lagrangian_structure(NLPeval)
    H_a = Array{Float64,1}(undef,length(H_struckt))
    MOI.eval_hessian_lagrangian(NLPeval,H_a,ω,1.0,[λ;γ])

    push!(H_struckt,(Nω,Nω))
    push!(H_a,0)
    H = sparse(Symmetric(sparse(first.(H_struckt),last.(H_struckt),H_a),:L))



    # Calculate the gradient of the hessian_lagrangian #? Might not be necessary
    #∇L = ∇F + ∇h*λ + Sᵤ'*νᵤ - Sₗ'*νₗ + Sₚ'*γ
    ∇L = ∇F + ∇h*λ
    ∇L[gₗ_index_set] -= νₗ
    ∇L[gᵤ_index_set] += νᵤ
    ∇L[end-Np+1:end] += γ

    # Calculte the optimality residual
    opt_res = [norm(∇L,Inf),
    norm(c,Inf),
    norm(cₚ,Inf),
    norm(min.(-min.(g,0),max.(λg,0)),Inf),
    norm(min.(-gₗ,νₗ),Inf),
    norm(min.(-gᵤ,νᵤ),Inf)]

    η = maximum(opt_res)


    # Divide the NLconstraints into strongly active constraint set and other
    g₊_index_set = [i for i in g_index_set if λ[i] > η]
    g₀_index_set = setdiff(g_index_set,g₊_index_set)
    #g₋_index_set = [i for i in g₀_index_set if h[i] > -η] #* Probably never going to see this one

    ∇g₊ = ∇h[:,g₊_index_set]

    Ng₊ = length(g₊_index_set)
    #Ng₀ = Ng - Ng₊

    # Divide boundry constraints into active and inactive sets

    gₗ₊_index_set = [i for i in 1:Ngₗ if νₗ[i] > η]
    gₗ₀_index_set = setdiff(1:Ngₗ,gₗ₊_index_set)
    #gₗ₋_index_set = [i for i in gₗ₀_index_set if gₗ[i] > -η] #* Probably never going to see this one

    ∇gₗ₊ = -I(Nω)[:,gₗ_index_set[gₗ₊_index_set]]

    Ngₗ₊ = length(gₗ₊_index_set)
    Ngₗ₀ = Ngₗ - Ngₗ₊

    gᵤ₊_index_set = [i for i in 1:Ngᵤ if νᵤ[i] > η]
    gᵤ₀_index_set = setdiff(1:Ngᵤ,gᵤ₊_index_set)
    #gᵤ₋_index_set = [i for i in gᵤ₀_index_set if gᵤ[i] > -η] #* Probably never going to see this one

    ∇gᵤ₊ = I(Nω)[:,gᵤ_index_set[gᵤ₊_index_set]]

    Ngᵤ₊ = length(gᵤ₊_index_set)
    #Ngᵤ₀ = Ngᵤ - Ngᵤ₊



    # Set up/build the KKT matrix for the sensitivity calculation:
    MXY = [∇c ∇g₊ ∇gₗ₊ ∇gᵤ₊ Sₚ']
    nx,ny = size(MXY)

    M = [H MXY;MXY' zeros(ny,ny)]
    #dM = Array(M)


    # Find the directional derivative 
    Δρ = [zeros(nx+ny-Np);-dp]

    # Solve the linear system of equations to obtain the nlp sensitivity with respect to the primal and dual all_variables
    timeFac = @elapsed Mf = lu(M)
    timeBac = @elapsed Δz = -(Mf\Δρ)


    #unpack the z vector
    Δω   = Δz[1:Nω]
    Δλc  = Δz[Nω+1:end][1:Nc]
    Δλg₊ = Δz[Nω+1:end][Nc+1:end][1:Ng₊]
    Δνₗ₊ = Δz[Nω+1:end][Nc+1:end][Ng₊+1:end][1:Ngₗ₊]
    Δνᵤ₊ = Δz[Nω+1:end][Nc+1:end][Ng₊+1:end][Ngₗ₊+1:end][1:Ngᵤ₊]
    Δγ   = Δz[Nω+1:end][Nc+1:end][Ng₊+1:end][Ngₗ₊+1:end][Ngᵤ₊+1:end]


    # Further modify the the results
    Δλ = zeros(Nh)
    Δλ[c_index_set]  = Δλc
    Δλ[g₊_index_set] = Δλg₊

    Δνₗ = zeros(Ngₗ)
    Δνₗ[gₗ₊_index_set] = Δνₗ₊

    Δνᵤ = zeros(Ngᵤ)
    Δνᵤ[gᵤ₊_index_set] = Δνᵤ₊

    return Δω,Δλ,Δνₗ,Δνᵤ,Δγ,timeFac,timeBac
end 


function nlpPredictor(nlp::Model,p::Array{Float64,1},dp::Array{Float64,1})
    lb_MOI_index = JuMP.MOI.get(nlp,JuMP.MOI.ListOfConstraintIndices{JuMP.MOI.VariableIndex,JuMP.MOI.GreaterThan{Float64}}())
    ub_MOI_index = JuMP.MOI.get(nlp,JuMP.MOI.ListOfConstraintIndices{JuMP.MOI.VariableIndex,JuMP.MOI.LessThan{Float64}}())

    Np = length(p)
    ωS = JuMP.all_variables(nlp)
    #Nω = length(ωS)

    ω = [JuMP.MOI.get(nlp,JuMP.MOI.VariablePrimal(),i) for i in ωS]

    λAll = -JuMP.MOI.get(nlp,JuMP.MOI.NLPBlockDual())
    λ = λAll[1:end-Np]
    γ = λAll[end-Np+1:end]

    if length(lb_MOI_index)>0
        νₗ = JuMP.MOI.get(nlp.moi_backend,JuMP.MOI.ConstraintDual(),lb_MOI_index)
    else
        νₗ = Array{Float64,1}(undef,0)
    end

    if length(ub_MOI_index)>0
        νᵤ = -JuMP.MOI.get(nlp.moi_backend,JuMP.MOI.ConstraintDual(),ub_MOI_index)
    else
        νᵤ = Array{Float64,1}(undef,0)
    end

    dω,dλ,dνₗ,dνᵤ,dγ,timeFac,timeBac = nlpPredictor(nlp,p,dp,ω,λ,νₗ,νᵤ,γ)

    return ω,λ,νₗ,νᵤ,γ,dω,dλ,dνₗ,dνᵤ,dγ,timeFac,timeBac
end


# Corrector steps
function nlpCorrector(nlp::Model,
                    p::Array{Float64,1},
                    dp::Array{Float64,1},
                    ω::Array{Float64,1},
                    λ::Array{Float64,1},
                    νₗ::Array{Float64,1},
                    νᵤ::Array{Float64,1},
                    γ::Array{Float64,1})

    # Dimensions of the problem
    Nω = length(ω)
    Np = length(p)
    Nh = length(λ)

    #Divide the NLconstraints into equality and inequality constraint set 
    c_index_set = [i for i in 1:Nh if nlp.nlp_model.constraints[nlp.nlp_model.constraints.keys[i]].set == JuMP.MOI.EqualTo{Float64}(0.0)]
    g_index_set = setdiff(1:Nh,c_index_set)

    Nc = length(c_index_set)
    Ng = Nh - Nc


    #Find the boundary constriant index sets
    gₗ_index_set = [i.value for i in JuMP.MOI.get(nlp,JuMP.MOI.ListOfConstraintIndices{JuMP.MOI.VariableIndex,JuMP.MOI.GreaterThan{Float64}}())] 
    gᵤ_index_set = [i.value for i in JuMP.MOI.get(nlp,JuMP.MOI.ListOfConstraintIndices{JuMP.MOI.VariableIndex,JuMP.MOI.LessThan{Float64}}())]

    Ngₗ = length(gₗ_index_set)
    Ngᵤ = length(gᵤ_index_set)

    # Find the parameter dummy variable selection matrix
    Sₚ = [spzeros(Np,Nω-Np) I(Np)]


    # Set up a NLPevaluator
    NLPeval = JuMP.NLPEvaluator(nlp)
    JuMP.MOI.initialize(NLPeval,[:Grad,:Jac,:Hess])

    # Evaluate NLconstraints 
    hold = Array{Float64,1}(undef,Nh+Np)
    MOI.eval_constraint(NLPeval,hold,ω)
    h = hold[1:Nh]

    c = h[c_index_set]
    g = h[g_index_set]

    λc = λ[c_index_set]
    λg = λ[g_index_set]

    # Evaluate the boundary constraints
    gₗ = lower_bound.(JuMP.all_variables(nlp)[gₗ_index_set]) - ω[gₗ_index_set]
    gᵤ = ω[gᵤ_index_set] - upper_bound.(JuMP.all_variables(nlp)[gᵤ_index_set])

    # Evaluate the gradient of the objective function
    ∇F = Array{Float64,1}(undef,Nω)
    MOI.eval_objective_gradient(NLPeval,∇F,ω)

    # Evaluate the Jacobian of the NLconstraints
    J_struckt = MOI.jacobian_structure(NLPeval)
    J_h = Array{Float64,1}(undef,length(J_struckt))
    MOI.eval_constraint_jacobian(NLPeval,J_h,ω)

    ∇h = sparse(last.(J_struckt),first.(J_struckt),J_h)[:,1:Nh] #The reason for calling last before first is because we want the gradient and not the jacobian
    ∇c = ∇h[:,c_index_set]
    #∇g = ∇h[:,g_index_set]

    # Evaluate the Hessian 
    H_struckt = MOI.hessian_lagrangian_structure(NLPeval)
    H_a = Array{Float64,1}(undef,length(H_struckt))
    MOI.eval_hessian_lagrangian(NLPeval,H_a,ω,1.0,[λ;γ])

    push!(H_struckt,(Nω,Nω))
    push!(H_a,0)
    H = sparse(Symmetric(sparse(first.(H_struckt),last.(H_struckt),H_a),:L))



    # Divide the Inequalities into B-active and B-inactive constraints
    # Divide the NLconstraints into the 2 sets
    g₊_index_set = [i for i in g_index_set if λ[i] + h[i] > 0]
    g₋_index_set = setdiff(g_index_set,g₊_index_set)

    g₊  = h[g₊_index_set]
    ∇g₊ = ∇h[:,g₊_index_set]
    λg₊ = λ[g₊_index_set]

    Ng₊ = length(g₊_index_set)
    Ng₋ = Ng - Ng₊

    # Divide boundry constraints into the 2 sets 

    gₗ₊_index_set = [i for i in 1:Ngₗ if νₗ[i] + gₗ[i] > 0]
    gₗ₋_index_set = setdiff(1:Ngₗ,gₗ₊_index_set)

    gₗ₊ = gₗ[gₗ₊_index_set]
    ∇gₗ₊ = -I(Nω)[:,gₗ_index_set[gₗ₊_index_set]]
    νₗ₊ = νₗ[gₗ₊_index_set]

    Ngₗ₊ = length(gₗ₊_index_set)
    Ngₗ₋ = Ngₗ - Ngₗ₊

    gᵤ₊_index_set = [i for i in 1:Ngᵤ if νᵤ[i] + gᵤ[i] > 0]
    gᵤ₋_index_set = setdiff(1:Ngᵤ,gᵤ₊_index_set)

    gᵤ₊ = gᵤ[gᵤ₊_index_set]
    ∇gᵤ₊ = I(Nω)[:,gᵤ_index_set[gᵤ₊_index_set]]
    νᵤ₊ = νᵤ[gᵤ₊_index_set]

    Ngᵤ₊ = length(gᵤ₊_index_set)
    Ngᵤ₋ = Ngᵤ - Ngᵤ₊



    # Set up/build the KKT matrix for the Newton step:
    MXY = [∇c ∇g₊ ∇gₗ₊ ∇gᵤ₊ Sₚ']
    nx,ny = size(MXY)

    M = [H MXY;MXY' zeros(ny,ny)]


    # Calculate the reduced residual
    ∇L = ∇F + ∇c*λc + ∇g₊*λg₊ 
    ∇L[gₗ_index_set[gₗ₊_index_set]] -= νₗ₊
    ∇L[gᵤ_index_set[gᵤ₊_index_set]] += νᵤ₊
    ∇L[end-Np+1:end] += γ

    ρ = [∇L;c;g₊;gₗ₊;gᵤ₊;ω[end-Np+1:end]-(p+dp)]

    # Solve the linear system of equations to obtain the Newton step to the primal and dual all_variables
    timeFac = @elapsed Mf = lu(M)
    timeBac = @elapsed Δz = -(Mf\ρ)


    #unpack the z vector
    Δω   = Δz[1:Nω]
    Δλc  = Δz[Nω+1:end][1:Nc]
    Δλg₊ = Δz[Nω+1:end][Nc+1:end][1:Ng₊]
    Δνₗ₊ = Δz[Nω+1:end][Nc+1:end][Ng₊+1:end][1:Ngₗ₊]
    Δνᵤ₊ = Δz[Nω+1:end][Nc+1:end][Ng₊+1:end][Ngₗ₊+1:end][1:Ngᵤ₊]
    Δγ   = Δz[Nω+1:end][Nc+1:end][Ng₊+1:end][Ngₗ₊+1:end][Ngᵤ₊+1:end]


    # Further modify the the results
    Δλ = zeros(Nh)
    Δλ[c_index_set]  = Δλc
    Δλ[g₊_index_set] = Δλg₊
    Δλ[g₋_index_set] = -λ[g₋_index_set]

    Δνₗ = zeros(Ngₗ)
    Δνₗ[gₗ₊_index_set] = Δνₗ₊
    Δνₗ[gₗ₋_index_set] = -νₗ[gₗ₋_index_set] 

    Δνᵤ = zeros(Ngᵤ)
    Δνᵤ[gᵤ₊_index_set] = Δνᵤ₊
    Δνᵤ[gᵤ₋_index_set] = -νᵤ[gᵤ₋_index_set]

    return Δω,Δλ,Δνₗ,Δνᵤ,Δγ,timeFac,timeBac
end

# Predictor and Corrector:
function nlpPredictorCorrector(nlp::Model,
                            p::Array{Float64,1},
                            dp::Array{Float64,1},
                            ω::Array{Float64,1},
                            λ::Array{Float64,1},
                            νₗ::Array{Float64,1},
                            νᵤ::Array{Float64,1},
                            γ::Array{Float64,1})

    # Dimensions of the problem
    Nω = length(ω)
    Np = length(p)
    Nh = length(λ)

    #Divide the NLconstraints into equality and inequality constraint set 
    c_index_set = [i for i in 1:Nh if nlp.nlp_model.constraints[nlp.nlp_model.constraints.keys[i]].set == JuMP.MOI.EqualTo{Float64}(0.0)]
    g_index_set = setdiff(1:Nh,c_index_set)

    Nc = length(c_index_set)
    Ng = Nh - Nc


    #Find the boundary constriant index sets
    gₗ_index_set = [i.value for i in JuMP.MOI.get(nlp,JuMP.MOI.ListOfConstraintIndices{JuMP.MOI.VariableIndex,JuMP.MOI.GreaterThan{Float64}}())] 
    gᵤ_index_set = [i.value for i in JuMP.MOI.get(nlp,JuMP.MOI.ListOfConstraintIndices{JuMP.MOI.VariableIndex,JuMP.MOI.LessThan{Float64}}())]

    Ngₗ = length(gₗ_index_set)
    Ngᵤ = length(gᵤ_index_set)


    # Set up a NLPevaluator
    NLPeval = JuMP.NLPEvaluator(nlp)
    JuMP.MOI.initialize(NLPeval,[:Grad,:Jac,:Hess])

    # Evaluate NLconstraints 
    hold = Array{Float64,1}(undef,Nh+Np)
    MOI.eval_constraint(NLPeval,hold,ω)
    h = hold[1:Nh]

    c = h[c_index_set]
    g = h[g_index_set]

    λc = λ[c_index_set]
    λg = λ[g_index_set]

    # Evaluate the boundary constraints
    gₗ = lower_bound.(JuMP.all_variables(nlp)[gₗ_index_set]) - ω[gₗ_index_set]
    gᵤ = ω[gᵤ_index_set] - upper_bound.(JuMP.all_variables(nlp)[gᵤ_index_set])

    # Evaluate the gradient of the objective function
    ∇F = Array{Float64,1}(undef,Nω)
    MOI.eval_objective_gradient(NLPeval,∇F,ω)

    # Evaluate the Jacobian of the NLconstraints
    J_struckt = MOI.jacobian_structure(NLPeval)
    J_h = Array{Float64,1}(undef,length(J_struckt))
    MOI.eval_constraint_jacobian(NLPeval,J_h,ω)

    ∇h = sparse(last.(J_struckt),first.(J_struckt),J_h)[:,1:Nh] #The reason for calling last before first is because we want the gradient and not the jacobian
    ∇c = ∇h[:,c_index_set]
    #∇g = ∇h[:,g_index_set]

    # Evaluate the Hessian 
    H_struckt = MOI.hessian_lagrangian_structure(NLPeval)
    H_a = Array{Float64,1}(undef,length(H_struckt))
    MOI.eval_hessian_lagrangian(NLPeval,H_a,ω,1.0,[λ;γ])

    push!(H_struckt,(Nω,Nω))
    push!(H_a,0)
    H = sparse(Symmetric(sparse(first.(H_struckt),last.(H_struckt),H_a),:L))



    # Divide the Inequalities into B-active and B-inactive constraints
    # Divide the NLconstraints into the 2 sets
    g₊_index_set = [i for i in g_index_set if λ[i] + h[i] > 0]
    g₋_index_set = setdiff(g_index_set,g₊_index_set)

    g₊  = h[g₊_index_set]
    ∇g₊ = ∇h[:,g₊_index_set]
    λg₊ = λ[g₊_index_set]

    Ng₊ = length(g₊_index_set)
    Ng₋ = Ng - Ng₊

    # Divide boundry constraints into the 2 sets 

    gₗ₊_index_set = [i for i in 1:Ngₗ if νₗ[i] + gₗ[i] > 0]
    gₗ₋_index_set = setdiff(1:Ngₗ,gₗ₊_index_set)

    gₗ₊ = gₗ[gₗ₊_index_set]
    ∇gₗ₊ = -I(Nω)[:,gₗ_index_set[gₗ₊_index_set]]
    νₗ₊ = νₗ[gₗ₊_index_set]

    Ngₗ₊ = length(gₗ₊_index_set)
    Ngₗ₋ = Ngₗ - Ngₗ₊

    gᵤ₊_index_set = [i for i in 1:Ngᵤ if νᵤ[i] + gᵤ[i] > 0]
    gᵤ₋_index_set = setdiff(1:Ngᵤ,gᵤ₊_index_set)

    gᵤ₊ = gᵤ[gᵤ₊_index_set]
    ∇gᵤ₊ = I(Nω)[:,gᵤ_index_set[gᵤ₊_index_set]]
    νᵤ₊ = νᵤ[gᵤ₊_index_set]

    Ngᵤ₊ = length(gᵤ₊_index_set)
    Ngᵤ₋ = Ngᵤ - Ngᵤ₊



    # Set up/build the KKT matrix for the Newton step:
    MXY = [∇c ∇g₊ ∇gₗ₊ ∇gᵤ₊ Sₚ']
    nx,ny = size(MXY)

    M = [H MXY;MXY' zeros(ny,ny)]


    # Calculate the reduced residual
    ∇L = ∇F + ∇c*λc + ∇g₊*λg₊ 
    ∇L[gₗ_index_set[gₗ₊_index_set]] -= νₗ₊
    ∇L[gᵤ_index_set[gᵤ₊_index_set]] += νᵤ₊
    ∇L[end-Np+1:end] += γ

    ρ = [∇L;g₊;gₗ₊;gᵤ₊;ω[end-Np+1:end]-(p)]

    # Calculate the sensitivity
    Δρ = [zeros(nx+ny-Np);-dp]


    # Factorize the KKT matrix and the do two backsolves to get the predictor and Corrector Steps
    timeFac  = @elapsed Mf = lu(M)
    timeBacC = @elapsed ΔzC = -(Mf\ρ)
    timeBacP = @elapsed ΔzP = -(Mf\Δρ)


    #unpack the z vector
    ΔωC   = ΔzC[1:Nω]
    ΔλcC  = ΔzC[Nω+1:end][1:Nc]
    Δλg₊C = ΔzC[Nω+1:end][Nc+1:end][1:Ng₊]
    Δνₗ₊C = ΔzC[Nω+1:end][Nc+1:end][Ng₊+1:end][1:Ngₗ₊]
    Δνᵤ₊C = ΔzC[Nω+1:end][Nc+1:end][Ng₊+1:end][Ngₗ₊+1:end][1:Ngᵤ₊]
    ΔγC   = ΔzC[Nω+1:end][Nc+1:end][Ng₊+1:end][Ngₗ₊+1:end][Ngᵤ₊+1:end]

    #unpack the z vector
    ΔωP   = ΔzP[1:Nω]
    ΔλcP  = ΔzP[Nω+1:end][1:Nc]
    Δλg₊P = ΔzP[Nω+1:end][Nc+1:end][1:Ng₊]
    Δνₗ₊P = ΔzP[Nω+1:end][Nc+1:end][Ng₊+1:end][1:Ngₗ₊]
    Δνᵤ₊P = ΔzP[Nω+1:end][Nc+1:end][Ng₊+1:end][Ngₗ₊+1:end][1:Ngᵤ₊]
    ΔγP   = ΔzP[Nω+1:end][Nc+1:end][Ng₊+1:end][Ngₗ₊+1:end][Ngᵤ₊+1:end]


    # Further modify the the results
    ΔλC = zeros(Nh)
    ΔλC[c_index_set]  = ΔλcC
    ΔλC[g₊_index_set] = Δλg₊C
    ΔλC[g₋_index_set] = -λ[g₋_index_set]

    ΔνₗC = zeros(Ngₗ)
    ΔνₗC[gₗ₊_index_set] = Δνₗ₊C
    ΔνₗC[gₗ₋_index_set] = -νₗ[gₗ₋_index_set] 

    ΔνᵤC = zeros(Ngᵤ)
    ΔνᵤC[gᵤ₊_index_set] = Δνᵤ₊C
    ΔνᵤC[gᵤ₋_index_set] = -νᵤ[gᵤ₋_index_set]

    #######
    ΔλP = zeros(Nh)
    ΔλP[c_index_set]  = ΔλcP
    ΔλP[g₊_index_set] = Δλg₊P

    ΔνₗP = zeros(Ngₗ)
    ΔνₗP[gₗ₊_index_set] = Δνₗ₊P

    ΔνᵤP = zeros(Ngᵤ)
    ΔνᵤP[gᵤ₊_index_set] = Δνᵤ₊P


    return ΔωC,ΔλC,ΔνₗC,ΔνᵤC,ΔγC,ΔωP,ΔλP,ΔνₗP,ΔνᵤP,ΔγP,timeFac,timeBacC,timeBacP
end 


function residualEquation(nlp::Model,
                        p::Array{Float64,1},
                        dp::Array{Float64,1},
                        ω::Array{Float64,1},
                        λ::Array{Float64,1},
                        νₗ::Array{Float64,1},
                        νᵤ::Array{Float64,1},
                        γ::Array{Float64,1})

    # Dimensions of the problem
    Nω = length(ω)
    Np = length(p)
    Nh = length(λ)

    #Divide the NLconstraints into equality and inequality constraint set 
    c_index_set = [i for i in 1:Nh if nlp.nlp_model.constraints[nlp.nlp_model.constraints.keys[i]].set == JuMP.MOI.EqualTo{Float64}(0.0)]
    g_index_set = setdiff(1:Nh,c_index_set)

    Nc = length(c_index_set)
    Ng = Nh - Nc


    #Find the boundary constriant index sets
    gₗ_index_set = [i.value for i in JuMP.MOI.get(nlp,JuMP.MOI.ListOfConstraintIndices{JuMP.MOI.VariableIndex,JuMP.MOI.GreaterThan{Float64}}())] 
    gᵤ_index_set = [i.value for i in JuMP.MOI.get(nlp,JuMP.MOI.ListOfConstraintIndices{JuMP.MOI.VariableIndex,JuMP.MOI.LessThan{Float64}}())]

    Ngₗ = length(gₗ_index_set)
    Ngᵤ = length(gᵤ_index_set)

    # Find the parameter dummy variable selection matrix
    Sₚ = [spzeros(Np,Nω-Np) I(Np)]
    cₚ = ω[end-Np+1:end] - (p+dp) # should be 0/or ϵ = precision of solution


    # Set up a NLPevaluator
    NLPeval = JuMP.NLPEvaluator(nlp)
    JuMP.MOI.initialize(NLPeval,[:Grad,:Jac,:Hess])

    # Evaluate NLconstraints 
    hold = Array{Float64,1}(undef,Nh+Np)
    MOI.eval_constraint(NLPeval,hold,ω)
    h = hold[1:Nh]

    c = h[c_index_set]
    g = h[g_index_set]

    #λc = λ[c_index_set]
    λg = λ[g_index_set]

    # Evaluate the boundary constraints
    gₗ = lower_bound.(JuMP.all_variables(nlp)[gₗ_index_set]) - ω[gₗ_index_set]
    gᵤ = ω[gᵤ_index_set] - upper_bound.(JuMP.all_variables(nlp)[gᵤ_index_set])

    # Evaluate the gradient of the objective function
    ∇F = Array{Float64,1}(undef,Nω)
    MOI.eval_objective_gradient(NLPeval,∇F,ω)

    # Evaluate the Jacobian of the NLconstraints
    J_struckt = MOI.jacobian_structure(NLPeval)
    J_h = Array{Float64,1}(undef,length(J_struckt))
    MOI.eval_constraint_jacobian(NLPeval,J_h,ω)

    ∇h = sparse(last.(J_struckt),first.(J_struckt),J_h)[:,1:Nh] #The reason for calling last before first is because we want the gradient and not the jacobian
    ∇c = ∇h[:,c_index_set]
    #∇g = ∇h[:,g_index_set]




    # Calculate the gradient of the hessian_lagrangian #? Might not be necessary
    #∇L = ∇F + ∇h*λ + Sᵤ'*νᵤ - Sₗ'*νₗ + Sₚ'*γ
    ∇L = ∇F + ∇h*λ
    ∇L[gₗ_index_set] -= νₗ
    ∇L[gᵤ_index_set] += νᵤ
    ∇L[end-Np+1:end] += γ

    # Calculte the optimality residual
    ρ = [∇L;c;min.(-g,λg);min.(-gₗ,νₗ);min.(-gᵤ,νᵤ);cₚ]

    return ρ
end


# Pathfollow algorithm 
function pathFollowAlg1(nlp::Model,p::Array{Float64,1},dp::Array{Float64,1},D=1e-6,Dinter=1e6)
    t0 = 0.
    dt = 1.
    α1 = 0.5
    iterCount = 0
    resList = Array{Float64,1}(undef,0)
    timeFacL = Array{Float64,1}(undef,0)
    timeBacL = Array{Float64,1}(undef,0)

    tList = [t0]

    #dpk = (t0+dt)*dp

    
    # Predictor Step: (only just after the nlpsolution)
    ωk,λk,νₗk,νᵤk,γk,dωk,dλk,dνₗk,dνᵤk,dγk,tFac,tBac = nlpPredictor(nlp,p,(t0+dt)*dp)
    iterCount += 1
    push!(timeFacL,tFac)
    push!(timeBacL,tBac)

    for i in 1:40
        ωu  = ωk +(t0+dt)*dωk
        λu  = λk +(t0+dt)*dλk
        νₗu = νₗk+(t0+dt)*dνₗk
        νᵤu = νᵤk+(t0+dt)*dνᵤk
        γu  = γk +(t0+dt)*dγk

        r = residualEquation(nlp,p,(t0+dt)*dp,ωu,λu,νₗu,νᵤu,γu)
        η = norm(r,Inf)
        if η <= Dinter
            push!(resList,η)
            t0 += dt
            dt = 1.0 - t0
            push!(tList,t0)

            ωk  = ωu
            λk  = λu
            νₗk = νₗu
            νᵤk = νᵤu
            γk  = γu
            break
        end
        if i == 40
            t0 += dt
            dt = 1.0 - t0
            push!(tList,t0)

            ωk  = copy(ωu)
            λk  = copy(λu)
            νₗk = copy(νₗu)
            νᵤk = copy(νᵤu)
            γk  = copy(γu)
            break
        end
        dt *= α1
    end
    #println("Iter ",iterCount,",t0 =",t0,",η =",resList[end])


    # Homotopy Phase
    function homoPhase(nlp,p,t0,dt,ωk,λk,νₗk,νᵤk,γk,iterCount,resList,timeFacL,timeBacL)
        while t0 < 1.0
            dωk,dλk,dνₗk,dνᵤk,dγk,tFac,tBac = nlpCorrector(nlp,p,(t0+dt)*dp,ωk,λk,νₗk,νᵤk,γk)
            iterCount += 1
            push!(timeFacL,tFac)
            push!(timeBacL,tBac)

            for i in 1:40
                ωu  = ωk +(t0+dt)*dωk
                λu  = λk +(t0+dt)*dλk
                νₗu = νₗk+(t0+dt)*dνₗk
                νᵤu = νᵤk+(t0+dt)*dνᵤk
                γu  = γk +(t0+dt)*dγk
        
                r = residualEquation(nlp,p,(t0+dt)*dp,ωu,λu,νₗu,νᵤu,γu)
                η = norm(r,Inf)
                if η <= Dinter
                    push!(resList,η)
                    t0 += dt
                    dt = 1.0 - t0
                    push!(tList,t0)
                
                    ωk  = ωu
                    λk  = λu
                    νₗk = νₗu
                    νᵤk = νᵤu
                    γk  = γu
                    break
                end
                if i == 40
                    t0 += dt
                    dt = 1.0 - t0
                    push!(tList,t0)
                
                    ωk  = ωu
                    λk  = λu
                    νₗk = νₗu
                    νᵤk = νᵤu
                    γk  = γu
                    break
                end
                dt *= α1
            end
            #println("Iter ",iterCount,",t0 =",t0,",η =",resList[end])
        end
        return ωk,λk,νₗk,νᵤk,γk,iterCount,resList,timeFacL,timeBacL
    end

    ωk,λk,νₗk,νᵤk,γk,iterCount,resList,timeFacL,timeBacL = homoPhase(nlp,p,t0,dt,ωk,λk,νₗk,νᵤk,γk,iterCount,resList,timeFacL,timeBacL)
    
    # Corrector Phase
    for i in 1:40
        r = residualEquation(nlp,p,dp,ωk,λk,νₗk,νᵤk,γk)
        η = norm(r,Inf)
        push!(resList,η)
        #println("Iter ",iterCount,",t0 =",t0+dt,",η =",resList[end])
        if η <= D
            break
        end
        dωk,dλk,dνₗk,dνᵤk,dγk,tFac,tBac = nlpCorrector(nlp,p,dp,ωk,λk,νₗk,νᵤk,γk)
        iterCount += 1
        push!(timeFacL,tFac)
        push!(timeBacL,tBac)

        ωk  += dωk
        λk  += dλk
        νₗk += dνₗk
        νᵤk += dνᵤk
        γk  += dγk
    end

    return ωk,λk,νₗk,νᵤk,γk,tList,iterCount,resList,timeFacL,timeBacL
end

function pathFollowAlg1(nlp::Model,
                        p::Array{Float64,1},
                        dp::Array{Float64,1},
                        ωk::Array{Float64,1},
                        λk::Array{Float64,1},
                        νₗk::Array{Float64,1},
                        νᵤk::Array{Float64,1},
                        γk::Array{Float64,1},
                        D=1e-6,
                        Dinter=1e6)
    t0 = 0.
    dt = 1.
    α1 = 0.5
    iterCount = 0
    resList = Array{Float64,1}(undef,0)
    timeFacL = Array{Float64,1}(undef,0)
    timeBacL = Array{Float64,1}(undef,0)

    tList = [t0]


    # Homo phase
    function homoPhase(nlp,p,t0,dt,ωk,λk,νₗk,νᵤk,γk,iterCount,resList,timeFacL,timeBacL)
        while t0 < 1.0
            dωk,dλk,dνₗk,dνᵤk,dγk,tFac,tBac = nlpCorrector(nlp,p,(t0+dt)*dp,ωk,λk,νₗk,νᵤk,γk)
            iterCount += 1
            push!(timeFacL,tFac)
            push!(timeBacL,tBac)

            for i in 1:40
                ωu  = ωk +(t0+dt)*dωk
                λu  = λk +(t0+dt)*dλk
                νₗu = νₗk+(t0+dt)*dνₗk
                νᵤu = νᵤk+(t0+dt)*dνᵤk
                γu  = γk +(t0+dt)*dγk
        
                r = residualEquation(nlp,p,(t0+dt)*dp,ωu,λu,νₗu,νᵤu,γu)
                η = norm(r,Inf)
                if η <= Dinter
                    push!(resList,η)
                    t0 += dt
                    dt = 1.0 - t0
                    push!(tList,t0)
                
                    ωk  = ωu
                    λk  = λu
                    νₗk = νₗu
                    νᵤk = νᵤu
                    γk  = γu
                    break
                end
                if i == 40
                    t0 += dt
                    dt = 1.0 - t0
                    push!(tList,t0)
                
                    ωk  = ωu
                    λk  = λu
                    νₗk = νₗu
                    νᵤk = νᵤu
                    γk  = γu
                    break
                end
                dt *= α1
            end
            #println("Iter ",iterCount,",t0 =",t0,",η =",resList[end])
        end
        return ωk,λk,νₗk,νᵤk,γk,iterCount,resList,timeFacL,timeBacL
    end

    ωk,λk,νₗk,νᵤk,γk,iterCount,resList,timeFacL,timeBacL = homoPhase(nlp,p,t0,dt,ωk,λk,νₗk,νᵤk,γk,iterCount,resList,timeFacL,timeBacL)

    # Corrector Phase
    for i in 1:40
        r = residualEquation(nlp,p,dp,ωk,λk,νₗk,νᵤk,γk)
        η = norm(r,Inf)
        push!(resList,η)
        #println("Iter ",iterCount,",t0 =",t0+dt,",η =",resList[end])
        if η <= D
            break
        end
        dωk,dλk,dνₗk,dνᵤk,dγk,tFac,tBac = nlpCorrector(nlp,p,dp,ωk,λk,νₗk,νᵤk,γk)
        iterCount += 1
        push!(timeFacL,tFac)
        push!(timeBacL,tBac)

        ωk  += dωk
        λk  += dλk
        νₗk += dνₗk
        νᵤk += dνᵤk
        γk  += dγk
    end

    return ωk,λk,νₗk,νᵤk,γk,tList,iterCount,resList,timeFacL,timeBacL
end

# Pathfollowing without the path
function pathFollowAlg2(nlp::Model,p::Array{Float64,1},dp::Array{Float64,1},D=1e-6)
    iterCount = 0
    resList = Array{Float64,1}(undef,0)
    timeFacL = Array{Float64,1}(undef,0)
    timeBacL = Array{Float64,1}(undef,0)
    timeTotL = Array{Float64,1}(undef,0)

    
    # Predictor Step: (only just after the nlpsolution)
    th = @elapsed ωk,λk,νₗk,νᵤk,γk,dωk,dλk,dνₗk,dνᵤk,dγk,tFac,tBac = nlpPredictor(nlp,p,dp)
    iterCount += 1
    push!(timeFacL,tFac)
    push!(timeBacL,tBac)
    push!(timeTotL,th)
    ωk  += dωk
    λk  += dλk
    νₗk += dνₗk
    νᵤk += dνᵤk
    γk  += dγk

    
    # Corrector Phase
    for i in 1:40
        r = residualEquation(nlp,p,dp,ωk,λk,νₗk,νᵤk,γk)
        η = norm(r,Inf)
        push!(resList,η)
        #println("Iter ",iterCount,",t0 =",t0+dt,",η =",resList[end])
        if η <= D
            break
        end
        th = @elapsed dωk,dλk,dνₗk,dνᵤk,dγk,tFac,tBac = nlpCorrector(nlp,p,dp,ωk,λk,νₗk,νᵤk,γk)
        iterCount += 1
        push!(timeFacL,tFac)
        push!(timeBacL,tBac)
        push!(timeTotL,th)

        ωk  += dωk
        λk  += dλk
        νₗk += dνₗk
        νᵤk += dνᵤk
        γk  += dγk
    end

    return ωk,λk,νₗk,νᵤk,γk,iterCount,resList,timeFacL,timeBacL,timeTotL
end

function pathFollowAlg2(nlp::Model,
                        p::Array{Float64,1},
                        dp::Array{Float64,1},
                        ωk::Array{Float64,1},
                        λk::Array{Float64,1},
                        νₗk::Array{Float64,1},
                        νᵤk::Array{Float64,1},
                        γk::Array{Float64,1},
                        D=1e-6)
    iterCount = 0
    resList = Array{Float64,1}(undef,0)
    timeFacL = Array{Float64,1}(undef,0)
    timeBacL = Array{Float64,1}(undef,0)
    timeTotL = Array{Float64,1}(undef,0)


    # Corrector Phase
    for i in 1:40
        th = @elapsed dωk,dλk,dνₗk,dνᵤk,dγk,tFac,tBac = nlpCorrector(nlp,p,dp,ωk,λk,νₗk,νᵤk,γk)
        iterCount += 1
        push!(timeFacL,tFac)
        push!(timeBacL,tBac)
        push!(timeTotL,th)

        ωk  += dωk
        λk  += dλk
        νₗk += dνₗk
        νᵤk += dνᵤk
        γk  += dγk

        r = residualEquation(nlp,p,dp,ωk,λk,νₗk,νᵤk,γk)
        η = norm(r,Inf)
        push!(resList,η)
        #println("Iter ",iterCount,",t0 =",t0+dt,",η =",resList[end])
        if η <= D
            break
        end
    end

    return ωk,λk,νₗk,νᵤk,γk,iterCount,resList,timeFacL,timeBacL,timeTotL
end



#######
# Do some debugging
#=
ΔGOR1Set = [-0.04,0,0.04]
ΔGOR2Set = [-0.01,0,0.01]

pBounds = [ΔGOR1Set,ΔGOR2Set]
pSc,pSet = scenarioGenerationBox(pBounds,1,1)

pSc[3]
Nu = 2
Nx = 8
Nz = 33
par = gasLiftPar()
include("/home/simen/PhD/utils/colMatrix.jl")

# Get initial condtions and guesses for mpc
ω0 = readdlm("ss_data/init_cond/noGasLift.txt")[:,1]
u0 = ω0[1:Nu]
x0 = ω0[Nu+1:end][1:Nx]
z0 = ω0[Nu+1:end][Nx+1:end]

u0 = round.(u0,digits=4)
x0 = round.(x0,digits=4)
z0 = round.(z0,digits=4)

ωSS = readdlm("ss_data/nominalOpt.txt")[:,1]
uSS = ωSS[1:Nu]
xSS = ωSS[Nu+1:end][1:Nx]
zSS = ωSS[Nu+1:end][Nx+1:end]


nlp = subProbNlp(pSc[3],par,x0,u0,uSS,xSS,zSS)



nlp.nlp_model.parameters[1:2] = uSS

optimize!(nlp)

p = uSS

dp = [0.5,1.0]*0.01

ωk,λk,νₗk,νᵤk,γk,tList,iterCount,resList = pathFollowAlg1(nlp,p,dp,0.1,1e4)
=#