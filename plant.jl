# Author: Simen Bjorvand
# Date: 18.10.22
# Description:
# This file containts a model and plant function for simulating the system:

# Required packages:
# DifferentialEquations, Sundials

#Required files:
# ~, 

# Required variables in scope:
# ~

function DAEmodel(out,dω,ω,p,t)
    # Unpacking parameters and inputs first:
    wgl,ΔGOR,par = p


    Nw = par.Nw
    Nx = 3Nw+2

    #Unpack outputs
    f = @view(out[1:Nx])
    g = @view(out[Nx+1:end])

    # Unpacking differentials:
    dx = dω[1:Nx]
    dmga = dx[1:Nw]
    dmgt = dx[Nw+1:2Nw]
    dmot = dx[2Nw+1:3Nw]
    dmgR,dmoR = dx[3Nw+1:end]


    # Unpacking states:
    x = ω[1:Nx]
    z = ω[Nx+1:end]

    # Differential states
    mga = x[1:Nw]
    mgt = x[Nw+1:2Nw]
    mot = x[2Nw+1:3Nw]
    mgR,moR = x[3Nw+1:end]

    # Algebraic states
    ρa  = [z[i] for i in 1:Nw]
    ρm  = [z[i+Nw] for i in 1:Nw]
    ρw  = [z[i+2*Nw] for i in 1:Nw]
    pa  = [z[i+3*Nw] for i in 1:Nw]
    pwh = [z[i+4*Nw] for i in 1:Nw]
    pwi = [z[i+5*Nw] for i in 1:Nw]
    pbh = [z[i+6*Nw] for i in 1:Nw]
    wiv = [z[i+7*Nw] for i in 1:Nw]
    wpc = [z[i+8*Nw] for i in 1:Nw]
    wpg = [z[i+9*Nw] for i in 1:Nw]
    wpo = [z[i+10*Nw] for i in 1:Nw]
    wro = [z[i+11*Nw] for i in 1:Nw]
    wrg = [z[i+12*Nw] for i in 1:Nw]
    ρRm = z[13*Nw+1]
    ρRH = z[13*Nw+2]
    pRh = z[13*Nw+3]
    pm  = z[13*Nw+4]
    wRH = z[13*Nw+5]
    wtg = z[13*Nw+6]
    wto = z[13*Nw+7]



    # Differential model equations:
    for nw in 1:Nw
        f[nw]     = wgl[nw] - wiv[nw] - dmga[nw]
        f[Nw+nw]  = wiv[nw] - wpg[nw] + wrg[nw] - dmgt[nw]
        f[2Nw+nw] = wro[nw] - wpo[nw] - dmot[nw]
    end
    f[3Nw+1] = sum(wpg[nw] for nw in 1:Nw) - wtg - dmgR
    f[3Nw+2] = sum(wpo[nw] for nw in 1:Nw) - wto - dmoR



    Ar = [par.Dr[nw]^2*pi/4 for nw in 1:par.Nw]
    Aw = [par.Dw[nw]^2*pi/4 for nw in 1:par.Nw]
    Aa = [par.Da[nw]^2*pi/4 - Aw[nw] for nw in 1:par.Nw] 
    AR = par.DR^2*pi/4
    # Algebraic model equations

    # Model equations
    for nw in 1:Nw
        g[nw]     = ρa[nw] - par.Mw*pa[nw]*1e5/par.Ta[nw]/par.R
        g[Nw+nw]  = ρm[nw] - (mgt[nw] + mot[nw])*pwh[nw]*par.Mw*par.ρo[nw]*1e5/(mot[nw]*pwh[nw]*par.Mw*1e5+par.ρo[nw]*par.R*par.Tw[nw]*mgt[nw]) 
        g[2Nw+nw] = ρw[nw] - (mgt[nw] + mot[nw]-par.ρo[nw]*par.Lr[nw]*Ar[nw])/(par.Lw[nw]*Aw[nw]) 
        g[3Nw+nw] = pa[nw]*1e5 - (1 + par.g*par.Ha[nw]*par.Mw/par.R/par.Ta[nw])*(par.Ta[nw]*par.R + par.g*par.Ha[nw]*par.Mw)/(par.Mw*par.La[nw]*Aa[nw])*mga[nw] 
        g[4Nw+nw] = pwh[nw]*1e5 - (par.Tw[nw]*par.R/par.Mw)*mgt[nw]/(par.Lw[nw]*Aw[nw] + par.Lr[nw]*Ar[nw] - mot[nw]/par.ρo[nw]) + (mgt[nw] + mot[nw])/(2*par.Lw[nw]*Aw[nw])*par.g*par.Hw[nw] 
        g[5Nw+nw] = pwi[nw]*1e5 - pwh[nw]*1e5 - ρw[nw]*par.g*par.Hw[nw] - par.μ0*par.Lw[nw]*wpc[nw]/(pi*par.Dw[nw]^4*ρm[nw])   
        g[6Nw+nw] = pbh[nw]*1e5 - pwi[nw]*1e5 - par.ρo[nw]*par.g*par.Hr[nw] - par.μ0*par.Lr[nw]*wro[nw]/(pi*par.Dr[nw]^4*par.ρo[nw]) 
        g[7Nw+nw] = wiv[nw]/par.Civ[nw] - sqrt(max(0,ρa[nw]*(pa[nw] - pwi[nw])*1e5)) 
        g[8Nw+nw] = wpc[nw]/par.Cpc[nw] - sqrt(max(0,ρm[nw]*(pwh[nw] - pm)*1e5)) 
        g[9Nw+nw] = wpg[nw] - mgt[nw]/(mgt[nw]+mot[nw])*wpc[nw] 
        g[10Nw+nw] = wpo[nw] - wpc[nw] + wpg[nw] 
        g[11Nw+nw] = wro[nw] - par.PI[nw]*(par.pr[nw] - pbh[nw]) 
        g[12Nw+nw] = wrg[nw] - (par.GOR[nw]+ΔGOR[nw])*wro[nw] 
    end

    g[13Nw+1] = ρRm - (mgR+moR)*pRh*par.Mw*par.ρmo*1e5/(moR*pRh*par.Mw*1e5 + par.ρmo*par.R*par.TR*mgR) 
    g[13Nw+2] = ρRH - (mgR+moR)/(par.LR*AR) 
    g[13Nw+3] = pRh*1e5 - par.TR*par.R*mgR/(par.Mw*par.LR*AR) + ρRH*par.g*par.HR/2 
    g[13Nw+4] = pm*1e5 - pRh*1e5 - ρRH*par.g*par.HR - par.μ0*par.LR*wRH/(pi*par.DR^4*ρRm) 
    g[13Nw+5] = (wRH/par.CRH)^2 - ρRm*(pRh-par.ps)*1e5 
    g[13Nw+6] = wtg - mgR/(mgR+moR)*wRH 
    g[13Nw+7] = wto + wtg - wRH 
end


function Plant(x0,z0,u,ΔGOR,par,dt)
    Nx = 3*par.Nw+2
    Nz = 13*par.Nw+7

    tspan = (0.0,dt)

    ω0 = [x0;z0]
    dω = zeros(Nx+Nz)
    diff_var = [[true for i in 1:Nx];[false for i in 1:Nz]]
    p = u,ΔGOR,par
    prob = DAEProblem(DAEmodel,dω,ω0,tspan,p,differential_vars = diff_var)
    sol = DifferentialEquations.solve(prob,IDA(),reltol=1e-8, abstol=1e-8)#IDA(),reltol=1e-8, abstol=1e-8)

    xk = sol.u[end][1:Nx]
    zk = sol.u[end][Nx+1:end]
    return xk,zk
end