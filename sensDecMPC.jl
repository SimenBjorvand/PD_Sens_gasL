# Date: 25.09 2023
# Author: Simen Bjorvand
# Description: A script where the decomposed mpcs will be implemented
#              Both based on Primal Decomposition:
#              1. Solved using only stepest decent 
#              2. Solved using second derivative information 
#              Use NLP sensitivities instead of resolving the nlps at each iteration 


# Required Libraries:
# JuMP, Ipopt, LinearAlgebra, SparseArrays

# Required self-written codes: 
# colMatrix.jl, msUtils.jl, myParams.jl, decMPCUtil.jl, msUtils.jl



#! Part 0: for debugging: get init conditions:
#=ωSS = readdlm("ss_data/nominalOpt.txt")[:,1]
ω0 = readdlm("ss_data/init_cond/noGasLift.txt")[:,1]



uSS = ωSS[1:Nu] 
xSS = ωSS[Nu+1:end][1:Nx] 
zSS = ωSS[Nu+1:end][Nx+1:end] 

u0 = ω0[1:Nu]
x0 = ω0[Nu+1:end][1:Nx] 

l0 = [uSS]

l0 = [[2.24,1.6]]
l0 = [[0.7,0.7]]
#l0 = [[1.6991246588750275,0.47053334951076875]]
=#



# First Function should be the basic SD decent algorithm
function sensDistMPCSD(x0,u0,l0,par,uSS,xSS,zSS,maxIter)
    #* PART 1: Build the subproblems
    Np = 2
    NPR = 2^Np + 1
    NRH = 1
    NSC = NPR^NRH


    ΔGOR1Set = [-0.04,0,0.04]
    ΔGOR2Set = [-0.01,0,0.01]

    pBounds = [ΔGOR1Set,ΔGOR2Set]
    pSc,pSet = scenarioGenerationBox(pBounds,NRH,NRH)

    #par = gasLiftPar()
    #Nx  = 3*par.Nw+2
    Nu  = par.Nw
    #Nz  = 13*par.Nw+7

    nlpS = Array{Model,1}(undef,NSC)
    for i in 1:NSC
        nlpS[i] = subProbNlp(pSc[i],par,x0,u0,uSS,xSS,zSS)
    end


    #~,~,~,~,~,~,~,~,γ5,∇γ5 = nlpSensitivity(nlpS[5],2)

    #* PART 2: Initialize the coordination algorithm
    # Step length
    α0 = 0.00003

    # Set up the primal decomposition algorithm: 
    l = Array{Array{Float64,1},1}(undef,sum(NPR^i for i in 0:(NRH-1)))
    for i in 1:length(l)
        l[i] = 0.7*ones(Nu)
    end

    #make an index array indicating which decompostion variable each variable belongs to?
    lIndexList = Array{Array{Int64,1}}(undef,NRH)
    lIndexList[1] = [1]
    if NRH > 1
        for i in 2:NRH
            lIndexList[i] = Array((1:NPR^(i-1)).+lIndexList[i-1][end])
        end
    end


    #Need to construct these selection matrices: 
    Sl = Array{Array{Float64,2},1}(undef,NPR^(NRH-1))
    for i in 1:length(Sl)
        Sl[i] = zeros(Nu*NRH,Nu*lIndexList[end][end])
    end

    for j in 1:NRH
        S1 = NPR^(j-1)
        for s1 in 1:S1
            S3 = NPR^(NRH-j)
            for s3 in 1:S3
                lI = lIndexList[j][s1]
                SlIndex = (s1-1)*S3+s3
                Sl[SlIndex][(j-1)*Nu+1:j*Nu,(lI-1)*Nu+1:lI*Nu] = I(Nu)
            end
        end   
    end

    # give initial values to the coordinator variables
    for j in 1:NRH
        for s1 in 1:NPR^(j-1)
            lI = lIndexList[j][s1]
            for nu = 1:Nu
                l[lI][nu] = l0[lI][nu]
            end
        end
    end

    #Update the coupled variables in the NLPs based on the NACs
    #Non-Anticapativity-Constraints:
    for j in 1:NRH
        for s1 in 1:NPR^(j-1)
            for s2 in 1:NPR^(NRH-j+1)
                for nu in 1:Nu
                    scenarioNumber = (s1-1)*(NPR^(NRH-j+1))+s2
                    lI = lIndexList[j][s1]
                    nlpS[scenarioNumber].nlp_model.parameters[(j-1)*Nu+nu] = l[lI][nu]
                end
            end
        end   
    end

    #Set up the l and dl variable related to each scenario
    lS = Array{Array{Float64,1},1}(undef,NSC)
    dlS = Array{Array{Float64,1},1}(undef,NSC)
    #Non-Anticapativity-Constraints:
    for i in 1:NSC
        lS[i] = Array{Float64,1}(undef,NRH*Nu)
        dlS[i] = zeros(NRH*Nu)
    end

    for j in 1:NRH
        for s1 in 1:NPR^(j-1)
            for s2 in 1:NPR^(NRH-j+1)
                scenarioNumber = (s1-1)*(NPR^(NRH-j+1))+s2
                lI = lIndexList[j][s1]
                lS[scenarioNumber][(j-1)*Nu+1:j*Nu] = l[lI]
                #dlS[scenarioNumber][(j-1)*Nu+1:j*Nu] = l[lI]
            end
        end   
    end 


    #* PART 3: Set up empty containers to store solution data
    # Set up variables which stores date related to the iterative algorithm 
    nlpSolveTimes = Array{Array{Float64,1},1}(undef,NSC)
    for i in 1:NSC
        nlpSolveTimes[i] = Array{Float64,1}(undef,0)
    end

    nlpIterList = Array{Array{Float64,1},1}(undef,NSC)
    for i in 1:NSC
        nlpIterList[i] = Array{Float64,1}(undef,0)
    end

    nlpFacTimeList = Array{Array{Float64,1},1}(undef,NSC)
    for i in 1:NSC
        nlpFacTimeList[i] = Array{Float64,1}(undef,0)
    end

    nlpBacTimeList = Array{Array{Float64,1},1}(undef,NSC)
    for i in 1:NSC
        nlpBacTimeList[i] = Array{Float64,1}(undef,0)
    end

    limSolTime = Array{Float64,1}(undef,0) #An array with the limiting solution times per iteration (due to parallelization)


    lList = [vcat(l...)]
    costList = Array{Float64,1}(undef,0)

    # Make sure the iterate with the best cost is feed to the system:
    ∇ϕList = Array{Array{Float64,1},1}(undef,0)



    # Set up error tolerance variables:
    Dtilde = 1e-2
    Dmin  = 1e-2
    D = Dmin

    # Store the Primal dual solutions
    ωS  = Array{Array{Float64,1},1}(undef,NSC)
    λS  = Array{Array{Float64,1},1}(undef,NSC)
    νₗS = Array{Array{Float64,1},1}(undef,NSC)
    νᵤS = Array{Array{Float64,1},1}(undef,NSC)
    γS  = Array{Array{Float64,1},1}(undef,NSC)


    #* PART 4: Iterations, solving the problem

    iterations = 0
    while true
        iterations +=1
        

        #* PART 4.1: Solve each of the nlps
        #=if iterations == 1
            tot_solve = 0.
            for i in 1:NSC
                print(i)
                hold_time,~ = solveNlp(nlpS[i])
                push!(nlpSolveTimes[i],hold_time)
                tot_solve += hold_time
            end

        elseif iterations == 2 # use the sensitivity
            for i in 1:NSC
                print(i)
                # Use the predictor step here
                ωs,λs,νₗs,νᵤs,γs,dωs,dλs,dνₗs,dνᵤs,dγs = nlpPredictor(nlpS[i],lS[i]-dlS[i],dlS[i])
                ωS[i]  = ωs  + dωs
                λS[i]  = λs  + dλs
                νₗS[i] = νₗs + dνₗs
                νᵤS[i] = νᵤs + dνᵤs
                γS[i]  = γs  + dγs
            end
        else 
            for i in 1:NSC
                print(i)
                dωs,dλs,dνₗs,dνᵤs,dγs = nlpCorrector(nlpS[i],lS[i]-dlS[i],dlS[i],ωS[i],λS[i],νₗS[i],νᵤS[i],γS[i])
                ωS[i]  += dωs
                λS[i]  += dλs
                νₗS[i] += dνₗs
                νᵤS[i] += dνᵤs
                γS[i]  += dγs
            end
        end=#
        if iterations == 1
            tot_solve = 0.
            for i in 1:NSC
                print(i)
                hold_time,~ = solveNlp(nlpS[i])
                push!(nlpSolveTimes[i],hold_time)
                tot_solve += hold_time
            end

        elseif iterations == 2 # use the sensitivity
            println("D = ",D)
            tot_solve = 0.
            tot_solve_alt = 0.
            for i in 1:NSC
                print(i)
                ωs,λs,νₗs,νᵤs,γs,iterCount,resList,tFac,tBac,tEx = pathFollowAlg2(nlpS[i],lS[i]-dlS[i],dlS[i],D)
                ωS[i]  = ωs
                λS[i]  = λs
                νₗS[i] = νₗs
                νᵤS[i] = νᵤs
                γS[i]  = γs
                push!(nlpSolveTimes[i],sum(tFac)+sum(tBac))
                tot_solve += sum(tFac)+sum(tBac)
                tot_solve_alt += sum(tEx)
                push!(nlpIterList[i],iterCount)
                push!(nlpFacTimeList[i],tFac...)
                push!(nlpBacTimeList[i],tBac...)
                #println(i,", ",iterCount,", ",resList[end],", Solve time: ",sum(tFac)+sum(tBac))
            end
        else
            println("D = ",D)
            tot_solve = 0.
            tot_solve_alt = 0.
            for i in 1:NSC
                print(i)
                ωs,λs,νₗs,νᵤs,γs,iterCount,resList,tFac,tBac,tEx = pathFollowAlg2(nlpS[i],lS[i]-dlS[i],dlS[i],ωS[i],λS[i],νₗS[i],νᵤS[i],γS[i],D)
                ωS[i]  = ωs
                λS[i]  = λs
                νₗS[i] = νₗs
                νᵤS[i] = νᵤs
                γS[i]  = γs
                push!(nlpSolveTimes[i],sum(tFac)+sum(tBac))
                tot_solve += sum(tFac)+sum(tBac)
                tot_solve_alt += sum(tEx)
                push!(nlpIterList[i],iterCount)
                push!(nlpFacTimeList[i],tFac...)
                push!(nlpBacTimeList[i],tBac...)
                #println(i,", ",iterCount,", ",resList[end],", Solve time: ",sum(tFac)+sum(tBac))
            end
        end

        push!(limSolTime,maximum([i[end] for i in nlpSolveTimes]))
        
        if iterations == 1
            println("   Iteration: ",iterations,", NLP solve Time: ",tot_solve,", ",limSolTime[end])
        else
            println("   Iteration: ",iterations,", NLP solve Time: ",tot_solve,", ",limSolTime[end],", total backsolves: ",sum(i[end] for i in nlpIterList))
        end


        ∇ϕ = zeros(length(l)*Nu)
        

        # Collect the dual variables of the complicating constraints:
        Λ = Array{Array{Float64,1},1}(undef,NSC)
        for i in 1:NSC
            if iterations == 1
                Λ[i] = -JuMP.MOI.get(nlpS[i],JuMP.MOI.NLPBlockDual())[end+1-NRH*Nu:end]
            else
                Λ[i] = γS[i]
            end                
        end



        #Updata the ∇ϕ by retriving the assosiated λs 
        for j in 1:NRH
            for s1 = 1:NPR^(j-1)
                for s2 = 1:NPR^(NRH-j+1)
                    scenarioNumber = (s1-1)*(NPR^(NRH-j+1))+s2
                    li = lIndexList[j][s1]
                    for nu in 1:Nu
                        ∇ϕ[(li-1)*Nu + nu] += -Λ[scenarioNumber][(j-1)*Nu+nu]
                        #Δl[li][nu] += Λ[scenarioNumber][(j-1)*par.Nu+nu]
                    end
                end
            end   
        end
        
        

        #Calculate Δl
        Δl = -∇ϕ            # SD

        #Find an appropriate step length
        αk = α0         # SD

        #Update the ls, such that the Input constraints are not broken
        for li in 1:length(l)
            l[li] += αk*Δl[(li-1)*Nu+1:li*Nu]
        end

        #Input constraint #* A better method should be adopted
        for i in 1:length(l)
            for nu in 1:Nu
                if l[i][nu] < 0.0
                    l[i][nu] = 0.0
                end
            end
        end


        # Update the nlps, with the new coupling constraint:
        for j in 1:NRH
            for s1 = 1:NPR^(j-1)
                for s2 = 1:NPR^(NRH-j+1)
                    scenarioNumber = (s1-1)*(NPR^(NRH-j+1))+s2
                    lI = lIndexList[j][s1]

                    dlS[scenarioNumber][(j-1)*Nu+1:j*Nu] = l[lI] - lS[scenarioNumber][(j-1)*Nu+1:j*Nu]
                    lS[scenarioNumber][(j-1)*Nu+1:j*Nu] = l[lI]
                end
            end   
        end

        # Update error tolerance:
        D = max(norm(Δl,Inf)*Dtilde,Dmin)

        push!(lList,vcat(l...))
        #push!(costList,costk)
        push!(∇ϕList,∇ϕ)



        #check inf norm of Δl and compare to some measure ϵ 
        Δla = lList[end]-lList[end-1]
        if LinearAlgebra.norm(vcat(Δla...),Inf)<1e-4
            break
        end
        if iterations >= maxIter
            break
        end
    end





    
    Δla = lList[end]-lList[end-1]
    nΔl = LinearAlgebra.norm(vcat(Δla...),Inf)
    nΔϕ = LinearAlgebra.norm(∇ϕList[end],Inf)

    uk = l[1]


    return uk, iterations,nlpSolveTimes,l,nΔl,nΔϕ,nlpIterList,nlpFacTimeList,nlpBacTimeList,limSolTime
end


# Second Function: Newton method
function sensDistMPCN(x0,u0,l0,par,uSS,xSS,zSS,maxIter)
    #* PART 1: Build the subproblems
    Np = 2
    NPR = 2^Np + 1
    NRH = 1
    NSC = NPR^NRH


    ΔGOR1Set = [-0.04,0,0.04]
    ΔGOR2Set = [-0.01,0,0.01]

    pBounds = [ΔGOR1Set,ΔGOR2Set]
    pSc,pSet = scenarioGenerationBox(pBounds,NRH,NRH)

    #par = gasLiftPar()
    #Nx  = 3*par.Nw+2
    Nu  = par.Nw
    #Nz  = 13*par.Nw+7

    nlpS = Array{Model,1}(undef,NSC)
    for i in 1:NSC
        nlpS[i] = subProbNlp(pSc[i],par,x0,u0,uSS,xSS,zSS)
    end


    #* PART 2: Initialize the coordination algorithm
    # Step length
    α0 = 1.0

    # Set up the primal decomposition algorithm: 
    l = Array{Array{Float64,1},1}(undef,sum(NPR^i for i in 0:(NRH-1)))
    for i in 1:length(l)
        l[i] = 0.7*ones(Nu)
    end

    #make an index array indicating which decompostion variable each variable belongs to?
    lIndexList = Array{Array{Int64,1}}(undef,NRH)
    lIndexList[1] = [1]
    if NRH > 1
        for i in 2:NRH
            lIndexList[i] = Array((1:NPR^(i-1)).+lIndexList[i-1][end])
        end
    end


    #Need to construct these selection matrices: 
    Sl = Array{Array{Float64,2},1}(undef,NPR^(NRH-1))
    for i in 1:length(Sl)
        Sl[i] = zeros(Nu*NRH,Nu*lIndexList[end][end])
    end

    for j in 1:NRH
        S1 = NPR^(j-1)
        for s1 in 1:S1
            S3 = NPR^(NRH-j)
            for s3 in 1:S3
                lI = lIndexList[j][s1]
                SlIndex = (s1-1)*S3+s3
                Sl[SlIndex][(j-1)*Nu+1:j*Nu,(lI-1)*Nu+1:lI*Nu] = I(Nu)
            end
        end   
    end

    # give initial values to the coordinator variables
    for j in 1:NRH
        for s1 in 1:NPR^(j-1)
            lI = lIndexList[j][s1]
            for nu = 1:Nu
                l[lI][nu] = l0[lI][nu]
            end
        end
    end

    #Update the coupled variables in the NLPs based on the NACs
    #Non-Anticapativity-Constraints:
    for j in 1:NRH
        for s1 in 1:NPR^(j-1)
            for s2 in 1:NPR^(NRH-j+1)
                for nu in 1:Nu
                    scenarioNumber = (s1-1)*(NPR^(NRH-j+1))+s2
                    lI = lIndexList[j][s1]
                    nlpS[scenarioNumber].nlp_model.parameters[(j-1)*Nu+nu] = l[lI][nu]
                end
            end
        end   
    end

    #Set up the l and dl variable related to each scenario
    lS = Array{Array{Float64,1},1}(undef,NSC)
    dlS = Array{Array{Float64,1},1}(undef,NSC)
    #Non-Anticapativity-Constraints:
    for i in 1:NSC
        lS[i] = Array{Float64,1}(undef,NRH*Nu)
        dlS[i] = zeros(NRH*Nu)
    end

    for j in 1:NRH
        for s1 in 1:NPR^(j-1)
            for s2 in 1:NPR^(NRH-j+1)
                scenarioNumber = (s1-1)*(NPR^(NRH-j+1))+s2
                lI = lIndexList[j][s1]
                lS[scenarioNumber][(j-1)*Nu+1:j*Nu] = l[lI]
                #dlS[scenarioNumber][(j-1)*Nu+1:j*Nu] = l[lI]
            end
        end   
    end 


    #* PART 3: Set up empty containers to store solution data
    # Set up variables which stores date related to the iterative algorithm 
    nlpSolveTimes = Array{Array{Float64,1},1}(undef,NSC)
    for i in 1:NSC
        nlpSolveTimes[i] = Array{Float64,1}(undef,0)
    end

    nlpIterList = Array{Array{Float64,1},1}(undef,NSC)
    for i in 1:NSC
        nlpIterList[i] = Array{Float64,1}(undef,0)
    end

    nlpFacTimeList = Array{Array{Float64,1},1}(undef,NSC)
    for i in 1:NSC
        nlpFacTimeList[i] = Array{Float64,1}(undef,0)
    end

    nlpBacTimeList = Array{Array{Float64,1},1}(undef,NSC)
    for i in 1:NSC
        nlpBacTimeList[i] = Array{Float64,1}(undef,0)
    end

    limSolTime = Array{Float64,1}(undef,0) #An array with the limiting solution times per iteration (due to parallelization)


    lList = [vcat(l...)]

    # Make sure the iterate with the best cost is feed to the system:
    #bestCost = Inf
    ∇ϕList = Array{Array{Float64,1},1}(undef,0)
    #bestl = copy(l)
    HList = Array{Array{Float64,2},1}(undef,0)

    # Alt stopping criteria: Assume all nlps could/are solved in parallel and stop when max solve time is supersceded 
    timeSpent = 0.0
    #timeSpent += st

    # Set up error tolerance variables:
    Dtilde = 1e-2
    Dmin  = 1e-2
    D = Dmin

    # Store the Primal dual solutions
    ωS  = Array{Array{Float64,1},1}(undef,NSC)
    λS  = Array{Array{Float64,1},1}(undef,NSC)
    νₗS = Array{Array{Float64,1},1}(undef,NSC)
    νᵤS = Array{Array{Float64,1},1}(undef,NSC)
    γS  = Array{Array{Float64,1},1}(undef,NSC)


    #* PART 4: Iterations, solving the problem

    iterations = 0
    while true
        iterations +=1
        

        #* PART 4.1: Solve each of the nlps
        if iterations == 1
            tot_solve = 0.
            for i in 1:NSC
                print(i)
                hold_time,~ = solveNlp(nlpS[i])
                push!(nlpSolveTimes[i],hold_time)
                #println(i,", Solve time: ",hold_time)
                tot_solve += hold_time
            end

        elseif iterations == 2 # use the sensitivity
            println("D = ",D)
            tot_solve = 0.
            tot_solve_alt = 0.
            for i in 1:NSC
                print(i)
                ωs,λs,νₗs,νᵤs,γs,iterCount,resList,tFac,tBac,tEx = pathFollowAlg2(nlpS[i],lS[i]-dlS[i],dlS[i],D)
                ωS[i]  = ωs
                λS[i]  = λs
                νₗS[i] = νₗs
                νᵤS[i] = νᵤs
                γS[i]  = γs
                push!(nlpSolveTimes[i],sum(tFac)+sum(tBac))
                tot_solve += sum(tFac)+sum(tBac)
                tot_solve_alt += sum(tEx)
                push!(nlpIterList[i],iterCount)
                push!(nlpFacTimeList[i],tFac...)
                push!(nlpBacTimeList[i],tBac...)
                #println(i,", ",iterCount,", ",resList[end],", Solve time: ",sum(tFac)+sum(tBac))
            end
        else
            println("D = ",D)
            tot_solve = 0.
            tot_solve_alt = 0.
            for i in 1:NSC
                print(i)
                ωs,λs,νₗs,νᵤs,γs,iterCount,resList,tFac,tBac,tEx = pathFollowAlg2(nlpS[i],lS[i]-dlS[i],dlS[i],ωS[i],λS[i],νₗS[i],νᵤS[i],γS[i],D)
                ωS[i]  = ωs
                λS[i]  = λs
                νₗS[i] = νₗs
                νᵤS[i] = νᵤs
                γS[i]  = γs
                push!(nlpSolveTimes[i],sum(tFac)+sum(tBac))
                tot_solve += sum(tFac)+sum(tBac)
                tot_solve_alt += sum(tEx)
                push!(nlpIterList[i],iterCount)
                push!(nlpFacTimeList[i],tFac...)
                push!(nlpBacTimeList[i],tBac...)
                #println(i,", ",iterCount,", ",resList[end],", Solve time: ",sum(tFac)+sum(tBac))
            end
        end

        push!(limSolTime,maximum([i[end] for i in nlpSolveTimes]))

        if iterations == 1
            println("   Iteration: ",iterations,", NLP solve Time: ",tot_solve)
        else
            println("   Iteration: ",iterations,", NLP solve Time: ",tot_solve,", ",limSolTime[end],", total backsolves: ",sum(i[end] for i in nlpIterList))
        end


        ∇ϕ = zeros(length(l)*Nu)
        

        # Collect the dual variables of the complicating constraints:
        Λ = Array{Array{Float64,1},1}(undef,NSC)
        for i in 1:NSC
            if iterations == 1
                Λ[i] = -JuMP.MOI.get(nlpS[i],JuMP.MOI.NLPBlockDual())[end+1-NRH*Nu:end]
            else
                Λ[i] = γS[i]
            end  
        end


        #Updata the ∇ϕ by retriving the assosiated λs 
        for j in 1:NRH
            for s1 = 1:NPR^(j-1)
                for s2 = 1:NPR^(NRH-j+1)
                    scenarioNumber = (s1-1)*(NPR^(NRH-j+1))+s2
                    li = lIndexList[j][s1]
                    for nu in 1:Nu
                        ∇ϕ[(li-1)*Nu + nu] += -Λ[scenarioNumber][(j-1)*Nu+nu]
                        #Δl[li][nu] += Λ[scenarioNumber][(j-1)*par.Nu+nu]
                    end
                end
            end   
        end

        
        # The Hessian calculation   
        ∇Λ = Array{Array{Float64,2},1}(undef,NSC)
        for i in 1:NSC
            if iterations == 1
                ~,~,~,~,~,~,~,~,~,∇γ,tFac,tBac = nlpSensitivity(nlpS[i],length(l)*Nu)
                ∇Λ[i] = ∇γ
            else
                ~,~,~,~,∇γ,tFac,tBac = nlpSensitivity(nlpS[i],lS[i],ωS[i],λS[i],νₗS[i],νᵤS[i],γS[i])
                ∇Λ[i] = ∇γ
            end
        end

        H = zeros(length(l)*Nu,length(l)*Nu)
        for i in 1:(NSC ÷ NPR)
            for j in 1:NPR
                sN = (i-1)*NPR+j
                H -= Sl[i]'*∇Λ[sN]*Sl[i]
            end
        end
        
        #println(eigen(H).values)
        println(H)

        #Calculate Δl
        Δl = -inv(H)*∇ϕ     # N


        #αk = stepLength(x0,l,vcat(Δl...),par,lIndexList,pSet)
        αk = 1.0       # N

        #Update the ls, such that the Input constraints are not broken
        for li in 1:length(l)
            l[li] += αk*Δl[(li-1)*Nu+1:li*Nu]
        end

        #Input constraint #* A better method should be adopted
        for i in 1:length(l)
            for nu in 1:Nu
                if l[i][nu] < 0.0
                    l[i][nu] = -0.0
                end
            end
        end


        # Update the nlps, with the new coupling constraint:
        for j in 1:NRH
            for s1 = 1:NPR^(j-1)
                for s2 = 1:NPR^(NRH-j+1)
                    scenarioNumber = (s1-1)*(NPR^(NRH-j+1))+s2
                    lI = lIndexList[j][s1]

                    dlS[scenarioNumber][(j-1)*Nu+1:j*Nu] = l[lI] - lS[scenarioNumber][(j-1)*Nu+1:j*Nu]
                    lS[scenarioNumber][(j-1)*Nu+1:j*Nu] = l[lI]
                end
            end   
        end

        # Update error tolerance:
        D = max(norm(∇ϕ,Inf)*Dtilde,Dmin)


        push!(lList,vcat(l...))
        push!(∇ϕList,∇ϕ)
        push!(HList,H)



        #check inf norm of Δl and compare to some measure ϵ 
        Δla = lList[end]-lList[end-1]
        if LinearAlgebra.norm(vcat(Δla...),Inf)<1e-4#1e-4
            break
        end
        if iterations>=maxIter
            break
        end
    end




    Δla = lList[end]-lList[end-1]
    nΔl = LinearAlgebra.norm(vcat(Δla...),Inf)
    nΔϕ = LinearAlgebra.norm(∇ϕList[end],Inf)

    uk = l[1]


    return uk, iterations,nlpSolveTimes,l,nΔl,nΔϕ,nlpIterList,nlpFacTimeList,nlpBacTimeList,limSolTime
end



#=
#Plot the solution
using PyPlot
fig,ax = subplots(2,figsize=[12,10])
ax[1].plot(1:length(lList),[i[1] for i in lList])
ax[2].plot(1:length(lList),[i[2] for i in lList])
fig



fig,ax = subplots(2,figsize=[12,8])
ax[1].semilogy(1:length(lList),[abs(i[1] - lList[end][1]) for i in lList])
ax[2].semilogy(1:length(lList),[abs(i[2] - lList[end][2]) for i in lList])
fig=#