# Date: 21.10 2022
# Modified: 18.09.2023
# Author: Simen Bjorvand
# Description: A script where the decomposed mpcs will be implemented
#              Both based on Primal Decomposition:
#              1. Trained using only stepest decent
#              2. Trained using second derivative information 


# Required Libraries:
# JuMP, Ipopt, LinearAlgebra, SparseArrays

# Required self-written codes: 
# colMatrix.jl, msUtils.jl, myParams.jl, decMPCUtil.jl, msUtils.jl


#=
#! Part 0: for debugging: get init conditions:
ωSS = readdlm("ss_data/nominalOpt.txt")[:,1]
ω0 = readdlm("ss_data/init_cond/noGasLift.txt")[:,1]



uSS = ωSS[1:Nu] 
xSS = ωSS[Nu+1:end][1:Nx] 
zSS = ωSS[Nu+1:end][Nx+1:end] 

u0 = ω0[1:Nu]
x0 = ω0[Nu+1:end][1:Nx] 

l0 = [uSS]

#l0 = [[2.24,1.6]]

#l0 = [[1.6991246588750275,0.47053334951076875]]
=#



# First Function should be the basic SD decent algorithm
function distMPCSD(x0,u0,l0,par,uSS,xSS,zSS,maxIter)
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

    #* PART 3: Set up empty containers to store solution data
    # Set up variables which stores date related to the iterative algorithm 
    nlpSolveTimes = Array{Array{Float64,1},1}(undef,NSC)
    for i in 1:NSC
        nlpSolveTimes[i] = Array{Float64,1}(undef,0)
    end

    lList = [vcat(l...)]
    costList = Array{Float64,1}(undef,0)

    # Make sure the iterate with the best cost is feed to the system:
    bestCost = Inf
    ∇ϕList = Array{Array{Float64,1},1}(undef,0)
    bestl = copy(l)

    # Alt stopping criteria: Assume all nlps could/are solved in parallel and stop when max solve time is supersceded 
    timeSpent = 0.0
    #timeSpent += st


    #* PART 4: Iterations, solving the problem

    iterations = 0
    while true
        iterations +=1
        

        #* PART 4.1: Solve each of the nlps
        tot_solve = 0.
        for i in 1:NSC
            print(i)
            hold_time,~ = solveNlp(nlpS[i])
            push!(nlpSolveTimes[i],hold_time)
            tot_solve += hold_time
        end
        costk = sum(objective_value.(nlpS))

        #holdTimek = [i[end] for i in nlpSolveTimes]
        #timeSpent += maximum(holdTimek)

        #if timeSpent>=max_solvetime
        #    iterations -= 1
        #    break
        #end
        println("   Iteration: ",iterations,", NLP solve Time: ",tot_solve," s, Cost = ",costk)

        # Best cost/length
        if costk < bestCost
            bestCost = costk
            bestl = copy(l)
        end

        ∇ϕ = zeros(length(l)*Nu)
        

        # Collect the dual variables of the complicating constraints:
        Λ = Array{Array{Float64,1},1}(undef,NSC)
        for i in 1:NSC
            Λ[i] = -JuMP.MOI.get(nlpS[i],JuMP.MOI.NLPBlockDual())[end+1-NRH*Nu:end]
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

        #=
        # The Hessian calculation   
        ∇Λ = Array{Array{Float64,2},1}(undef,NSC)
        for i in 1:NSC
            ~,~,~,~,~,~,~,~,~,∇γ = nlpSensitivity(nlpS[i],length(l)*Nu)
            ∇Λ[i] = ∇γ
        end

        H = zeros(length(l)*Nu,length(l)*Nu)
        for i in 1:(NSC ÷ NPR)
            for j in 1:NPR
                sN = (i-1)*NPR+j
                H -= Sl[i]'*∇Λ[sN]*Sl[i]
            end
        end
        =#
        
        
        

        #Calculate Δl
        Δl = -∇ϕ            # SD
        #Δl = -inv(H)*∇ϕ     # N

        #Find an appropriate step length #! need this
        #αk = stepLength(x0,l,vcat(Δl...),par,lIndexList,pSet)
        αk = α0         # SD
        #αk = 1.0       # N

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

        # Update the guesses of the nlps
        #Calculate Δl
        updateGuessNlp.(nlpS)



        # Update the nlps, with the new coupling constraint:
        for j in 1:NRH
            for s1 = 1:NPR^(j-1)
                for s2 = 1:NPR^(NRH-j+1)
                    scenarioNumber = (s1-1)*(NPR^(NRH-j+1))+s2
                    lI = lIndexList[j][s1]
                    for nu in 1:Nu
                        nlpS[scenarioNumber].nlp_model.parameters[(j-1)*Nu+nu] = l[lI][nu]
                    end
                end
            end   
        end

        push!(lList,vcat(l...))
        push!(costList,costk)
        push!(∇ϕList,∇ϕ)



        #check inf norm of Δl and compare to some measure ϵ 
        Δla = lList[end]-lList[end-1]
        if LinearAlgebra.norm(vcat(Δla...),Inf)<1e-4
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


    return uk, iterations,nlpSolveTimes,l,bestCost,nΔl,nΔϕ
end


# Second Function: Newton method
function distMPCN(x0,u0,l0,par,uSS,xSS,zSS,maxIter)
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

    #* PART 3: Set up empty containers to store solution data
    # Set up variables which stores date related to the iterative algorithm 
    nlpSolveTimes = Array{Array{Float64,1},1}(undef,NSC)
    for i in 1:NSC
        nlpSolveTimes[i] = Array{Float64,1}(undef,0)
    end

    lList = [vcat(l...)]
    costList = Array{Float64,1}(undef,0)

    # Make sure the iterate with the best cost is feed to the system:
    bestCost = Inf
    ∇ϕList = Array{Array{Float64,1},1}(undef,0)
    bestl = copy(l)
    HList = Array{Array{Float64,2},1}(undef,0)

    # Alt stopping criteria: Assume all nlps could/are solved in parallel and stop when max solve time is supersceded 
    timeSpent = 0.0
    #timeSpent += st


    #* PART 4: Iterations, solving the problem

    iterations = 0
    while true
        iterations +=1
        

        #* PART 4.1: Solve each of the nlps
        tot_solve = 0.
        for i in 1:NSC
            print(i)
            hold_time,~ = solveNlp(nlpS[i])
            push!(nlpSolveTimes[i],hold_time)
            tot_solve += hold_time
        end
        costk = sum(objective_value.(nlpS))

        #holdTimek = [i[end] for i in nlpSolveTimes]
        #timeSpent += maximum(holdTimek)

        #if timeSpent>=max_solvetime
        #    iterations -= 1
        #    break
        #end
        println("   Iteration: ",iterations,", NLP solve Time: ",tot_solve," s, Cost = ",costk)

        # Best cost/length
        if costk < bestCost
            bestCost = costk
            bestl = copy(l)
        end

        ∇ϕ = zeros(length(l)*Nu)
        

        # Collect the dual variables of the complicating constraints:
        Λ = Array{Array{Float64,1},1}(undef,NSC)
        for i in 1:NSC
            Λ[i] = -JuMP.MOI.get(nlpS[i],JuMP.MOI.NLPBlockDual())[end+1-NRH*Nu:end]
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
            ~,~,~,~,~,~,~,~,~,∇γ = nlpSensitivity(nlpS[i],length(l)*Nu)
            ∇Λ[i] = ∇γ
        end

        H = zeros(length(l)*Nu,length(l)*Nu)
        for i in 1:(NSC ÷ NPR)
            for j in 1:NPR
                sN = (i-1)*NPR+j
                H -= Sl[i]'*∇Λ[sN]*Sl[i]
            end
        end
        
        

        #Calculate Δl
        Δl = -inv(H)*∇ϕ     # N

        #Find an appropriate step length#! This is vital
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

        # Update the guesses of the nlps
        #Calculate Δl
        updateGuessNlp.(nlpS)



        # Update the nlps, with the new coupling constraint:
        for j in 1:NRH
            for s1 = 1:NPR^(j-1)
                for s2 = 1:NPR^(NRH-j+1)
                    scenarioNumber = (s1-1)*(NPR^(NRH-j+1))+s2
                    lI = lIndexList[j][s1]
                    for nu in 1:Nu
                        nlpS[scenarioNumber].nlp_model.parameters[(j-1)*Nu+nu] = l[lI][nu]
                    end
                end
            end   
        end

        push!(lList,vcat(l...))
        push!(costList,costk)
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


    return uk, iterations,nlpSolveTimes,l,bestCost,nΔl,nΔϕ
end