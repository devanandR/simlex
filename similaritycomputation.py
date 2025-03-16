# -*- coding: utf-8 -*-
"""
Created on Sun Dec 19 17:07:31 2021

@author: devanand
"""

import cplex
import sys
import time

#Brief: 
"""At each lp call in SimLex, we first take our previous basis,
 columns of nonbasic vars , non basic vars, current LPs cost vector and get fractional similarity scores
 for checking whether current LP is similar to previous LP or not
 
 """


def dotProduct(K,L):
    
   if len(K) != len(L):
      return 0

   return sum(i[0]*i[1] for i in zip(K, L))

def computedelta_BBinverse(model:cplex.Cplex, deltaB):
    
    print("Inverse operation start")
    #print(deltaB,"deltaB vector")
    #print(len(deltaB),"length of deta B")
    #print(sum(deltaB),"sum of deltaB")
    deltaBinverse = model.solution.advanced.btran(deltaB)
    #print(deltaBinverse,"deltaBinverse")
    #print(len(deltaBinverse)," size of deltaBinverse")
    #print()
    #print()
    #print("Inverse operation done")
    
    return deltaBinverse
    

def ComputeRCComponent(cost,prevBasisI, prevBasisStatus): # C_b^TB^-1
    
    # Cost vector associated to basic vars
    CB = []
    for varname,i in zip(prevBasisStatus.keys(),range(len(prevBasisStatus.keys()))):
        if prevBasisStatus[varname]==1 :
            CB.append(cost[i])
    
    #print(CB,"CB", len(CB)," its length")
    #print(prevBasisI,"prevBasisI", len(prevBasisI),"its length")
    
    # Compute Cb^T.B^-1
    CbBinv = {}
    
    for a_clmn in prevBasisI:
        CbBinv[a_clmn] = cost[a_clmn] - dotProduct(CB, prevBasisI[a_clmn])
    return CbBinv 
    
def computeBndSim(model, prevSolution):
    
    return 0.
def computeDualAj(prevColj,dual_val, prev_cost,costj):
    
    nonzindex = prevColj.ind
    nonzval = prevColj.val
                 
    c_bB_invAj =0.
    for indx in range(len(nonzindex)):
                
        rowId = nonzindex[indx]
        if prev_cost[rowId] == 0:
            ratio = costj[rowId]*1.0E+4 # Some high value
        else:
            ratio = costj[rowId]/prev_cost[rowId]
            
        c_bB_invAj = c_bB_invAj + dual_val[rowId]*nonzval[indx]*(ratio)            
                
    #c_bar = costj - c_bB_invAj 
    
    return c_bB_invAj

def computeObjSimDual(model, prev_cost, DUAL_VALUE, prevCol, prevBasisStatus):
    
    varNames = model.variables.get_names() 
    varCoef = model.objective.get_linear()
    cost ={}
    for i in range(len(varNames)):
        cost[varNames[i]] = varCoef[i]
    #c_bB_inv = ComputeRCComponent(cost,DUAL_VALUE, prevBasisStatus)
    #c_bB_inv = DUAL_VALUE
    
    violation =0.
    totalrc =0.
    # compute reduced cost of each nonbasic vars
    for varname,i in zip(prevBasisStatus.keys(),range(len(prevBasisStatus.keys()))):
        if prevBasisStatus[varname]==0: # it means it is non basic and at its lower bound
            c_bar = computeDualAj(prevCol[i],DUAL_VALUE,
                                  prev_cost,cost[varname])
            #c_bar = cost[varname] - ()),c_bB_inv,prev_cost[varname],cost[varname])
            
            if c_bar < 0.:
                violation += abs(c_bar)
            
            totalrc +=abs(c_bar)
                
        if prevBasisStatus[varname]==2: # it means it is non basic and at its upper bound
            
            c_bar = computeDualAj(prevCol[i],DUAL_VALUE,
                                  prev_cost,cost[varname])
            
            if c_bar > 0.:
                violation += c_bar
        
            totalrc +=c_bar
                
    return violation/totalrc          
    
def  computeBndSimDual(model, prevSolution):
    
    return 0.


def computeObjSim(model:cplex.Cplex, BINVERSE_A, prev_cost, prevCol, prevBasisStatus):
    
    varNames = model.variables.get_names() 
    varCoef = model.objective.get_linear()
    cost ={}
    for i in range(len(varNames)):
        cost[varNames[i]] = varCoef[i]
    current_rc = ComputeRCComponent(varCoef,BINVERSE_A, prevBasisStatus)
    
    violation =0.
    totalrc =0.
    
    for varname,i in zip(prevBasisStatus.keys(),range(len(prevBasisStatus.keys()))):
        if prevBasisStatus[varname]==0: # it means it is non basic and at its lower bound
            
            
            c_bar = current_rc[i]
            if c_bar < 0.:
                #violation += abs(c_bar)
                violation += abs(prev_cost[varname] - cost[varname])
            
            totalrc +=abs(c_bar)
                
        if prevBasisStatus[varname]==2: # it means it is non basic and at its upper bound
            
            c_bar = current_rc[i]
            
            if c_bar > 0.:
                #violation += abs(c_bar)
                violation += abs(prev_cost[varname] - cost[varname])
            
            totalrc +=abs(prev_cost[varname] - cost[varname])
     
    if totalrc ==0:
        return 0
    else:
        return violation/totalrc          
    
def settingCplexOption(model:cplex.Cplex,fracScoreObj, fracScoreBnd,logging,similarity_type,PRIMAL_FLAG):
    
    storeSense = False
        
    if similarity_type ==1:
        
        logging.info("similar type = 1")    
        thrshold_obj = 0.5
        threshold_bnd = 0.4
        
    if similarity_type ==2:
        
        logging.info("similar type = 2")
        thrshold_obj = 0.001
        threshold_bnd = 0.2
        
    logging.info("obj_frac_dev: " + str(fracScoreObj) + "and threshold obj is: " + str(thrshold_obj))
    if fracScoreObj <= thrshold_obj and fracScoreBnd<= threshold_bnd:
        
        ## Similar: It means keep the method primal and adv =1
        model.parameters.advance.set(1)
        model.parameters.preprocessing.presolve.set(0)
        
        if PRIMAL_FLAG:
            logging.info("current problem is similar to prev: adv =1, presolve =0, method = primal")
        
            model.parameters.lpmethod.set(model.parameters.lpmethod.values.primal)
        else:
            logging.info("current problem is similar to prev: adv =1, presolve =0, method = dual")
        
            model.parameters.lpmethod.set(model.parameters.lpmethod.values.dual)
        
    else:
        ## Dissimilar: It means start from scratch , keep concurrent and adv = 0
        
        model.parameters.advance.set(0)
        model.parameters.preprocessing.presolve.set(1)
        
        if PRIMAL_FLAG:
            logging.info("current problem is not similar to prev: adv =0, presolve =1, method = concurrent")
            model.parameters.lpmethod.set(model.parameters.lpmethod.values.primal)
        else:
            logging.info("current problem is not similar to prev: adv =0, presolve =1, method = concurrent")
        
            model.parameters.lpmethod.set(model.parameters.lpmethod.values.dual)
            
        model.parameters.lpmethod.set(model.parameters.lpmethod.values.concurrent)
        storeSense = True
    
    return storeSense
        
def computeSimRCBAsed(model:cplex.Cplex, BINVERSE_A, prev_cost, prevCol, prevBasisStatus,
                      prevSolution,logging, PRIMAL_FLAG):
    
    #global storeSense
    similarity_type = 2
    storeSense = False
    #prevBasisI = model.solution.advanced.binvcol(itr)
    fracScoreObj = computeObjSim(model, BINVERSE_A, prev_cost, prevCol, prevBasisStatus)
    fracScoreBnd = computeBndSim(model, prevSolution)

    storeSense = settingCplexOption(model, fracScoreObj, fracScoreBnd, 
                                    logging, similarity_type,PRIMAL_FLAG)
    return storeSense
def computeSimDualBAsed(model:cplex.Cplex, prev_cost, dual_val, prevCol, prevBasisStatus,
                      prevSolution,logging, PRIMAL_FLAG):
    
    #global storeSense
    similarity_type = 2
    storeSense = False
    #prevBasisI = model.solution.advanced.binvcol(itr)
    fracScoreObj = computeObjSimDual(model, prev_cost, dual_val, prevCol, prevBasisStatus)
    fracScoreBnd = computeBndSimDual(model, prevSolution)

    storeSense = settingCplexOption(model, fracScoreObj, fracScoreBnd, 
                                    logging, similarity_type,PRIMAL_FLAG)
    return storeSense

def computeSimSimple(model:cplex.Cplex, prev_cost,logging, PRIMAL_FLAG):
    
    # Present cost 
    similarity_type = 1
    present_cost = {}
    varNames = model.variables.get_names() 
    varCoef = model.objective.get_linear()

    for i in range(len(varNames)):
        present_cost[varNames[i]] = varCoef[i]
    
    deviation = 0.
    dev_val =0.
    total_val =0.
    total = 0.
    for name in prev_cost:
        
        if present_cost[name]==0 and prev_cost[name]== 0:
            continue
        if present_cost[name]*prev_cost[name] <= 0:
            
            dev_val  = dev_val + abs(present_cost[name] - prev_cost[name])
            deviation+=1
        
        total_val = total_val + max(present_cost[name],prev_cost[name])
        total += 1 
    
    if total_val==0:
        fractional_dev_obj = 0 
    #fractional_dev_obj = deviation/total
    else:
        fractional_dev_obj = dev_val/(2*total_val)
    ## for bound similarity always true
    fracScoreBnd = 0. ## 0 means always similar
    
    storeSense = settingCplexOption(model, fractional_dev_obj, fracScoreBnd, 
                                    logging, similarity_type,PRIMAL_FLAG)
    
    return storeSense

def computeSimScore(new_problem, logging,PREV_COST, prevRCForSim):
    
    startTime = time.time()
    solveStatus = new_problem.solve()
    TotalTime = time.time() - startTime
        
    if solveStatus == 0:
        print("Probem is infeasible")
        
    print(TotalTime,"Time Taken in Solving dummy Model")
    logging.info("Time Taken in Solving dummy Model" +str(TotalTime))    
    #obj_value = new_problem.solution.get_objective_value()
    
    varReducedCosts = new_problem.solution.get_reduced_costs()
    
    varNames = new_problem.variables.get_names() 
    
    varCoef = new_problem.objective.get_linear()
    
    prevobj_val = 0
    
    violation =0.
    total =0.
    similarity_score =0
    new_violation = 0
    v_count = 0
    
    for varname, var_new,var_id in zip(PREV_COST,varNames,range(len(varCoef))):
        if varname != var_new:
            print("index based computation is not possible")
            sys.exit(1)
        else:
            if prevRCForSim[varname] <0 and varReducedCosts[var_id] >=0:
                v_count+=1
                new_violation += abs(prevRCForSim[varname] - varReducedCosts[var_id])
                violation = violation + abs(PREV_COST[varname] - varCoef[var_id])
                #violation = violation + abs(PREV_COST[varname] - varCoef[var_id])
            
            if prevRCForSim[varname] >0 and varReducedCosts[var_id] <=0:
                v_count+=1
                violation = violation + abs(PREV_COST[varname] - varCoef[var_id])
                new_violation += abs(prevRCForSim[varname] - varReducedCosts[var_id])
                
            if prevRCForSim[varname] ==0 and varReducedCosts[var_id] !=0:
                violation = violation + abs(PREV_COST[varname] - varCoef[var_id])
                new_violation += abs(prevRCForSim[varname] + varReducedCosts[var_id])
                v_count+=1
                
                
            total = total + abs(PREV_COST[varname] - varCoef[var_id])
            #total_new_violation+= abs(abs(prevRCForSim[varname]) - abs(varReducedCosts[var_id]))
      
    if total==0:
        similarity_score = 1
    else:
        similarity_score = 1 - (violation/total)
    
    print("similarity score is ", similarity_score)
    logging.info(str(similarity_score))    
    
    return similarity_score, new_violation, v_count

    

def basicObjTric(model:cplex.Cplex, prevRCForSim, PREV_COST, problemFileWithBasis,\
                 problemFileWithoutBasis, logging,PRIMAL_FLAG):
    
    
    storeSense = False
    threshold_val = 0.9
    # Generate Random objective value (from first feasible solution)
    print("similarity score computation for ", problemFileWithBasis)
    logging.info("similarity score computation for "+ problemFileWithBasis)    
    
    #model.solution.basis.write("basis.bas")
    new_problem = cplex.Cplex(problemFileWithBasis)
    new_problem.parameters.simplex.limits.iterations.set(0) # without iteration
    new_problem.parameters.advance.set(1)
    #new_problem.parameters.preprocessing.presolve.set(0)
    
    new_problem.start.read_basis("basis.bas") # basis.bas was saved is the basis of last optimal basis
        
    if PRIMAL_FLAG:
        new_problem.parameters.lpmethod.set(new_problem.parameters.lpmethod.values.primal)
    else:
        new_problem.parameters.lpmethod.set(new_problem.parameters.lpmethod.values.dual)
     
    sim_using_basis, new_violation_basis, v_count1 = computeSimScore(new_problem,logging,PREV_COST, prevRCForSim)
    
    # Similarly for the case without basis
    # Generate Random objective value (from first feasible solution)
    print("similarity score computation for ", problemFileWithoutBasis)
    logging.info("similarity score computation for "+ problemFileWithoutBasis)    
    
    new_problem1 = cplex.Cplex(problemFileWithoutBasis)
    new_problem1.parameters.advance.set(0)
    new_problem1.parameters.simplex.limits.lowerobj.set(1.0E+75) # Setting to get feasible sol
    #new_problem1.parameters.simplex.limits.iterations.set(1) # without iteration
    #new_problem1.parameters.preprocessing.presolve.set(0)
    
    if PRIMAL_FLAG:
        new_problem1.parameters.lpmethod.set(new_problem1.parameters.lpmethod.values.primal)
    else:
        new_problem1.parameters.lpmethod.set(new_problem1.parameters.lpmethod.values.dual)
    
    sim_no_basis, new_violation_no_basis, v_count2 = computeSimScore(new_problem1,logging,PREV_COST, prevRCForSim)
    #new_problem1.solution.basis.write("temp.bas")
    print()
    print()
    print("simscore_using_basis: ",sim_using_basis)
    print("new_violation_basis: ",new_violation_basis)
    print()
    print("simscore_without_basis", sim_no_basis)
    print("new_violation_no_basis", new_violation_no_basis)
    print()
    print()
    if new_violation_no_basis > new_violation_basis:# or sim_using_basis >=  sim_no_basis:
    #if v_count2 >v_count1:# or sim_using_basis >=  sim_no_basis:
    #if sim_using_basis > sim_no_basis:
        
        print("Found a situation where B is similar: useful")
        logging.info("Found a situation where B is similar: useful")    
    
        #sys.exit(1)
        ## Similar: It means keep the method primal and adv =1
        
        if PRIMAL_FLAG:
            logging.info("current problem is similar to prev: adv =1, presolve =0, method = primal")
            
            model.parameters.lpmethod.set(model.parameters.lpmethod.values.primal)
        else:
            logging.info("current problem is similar to prev: adv =1, presolve =0, method = dual")
        
            model.parameters.lpmethod.set(model.parameters.lpmethod.values.dual)
        
        #model.start.read_basis("basis.bas") # basis.bas was saved is the basis of last optimal basis
        
        model.parameters.advance.set(1)
        
        return storeSense, model
    
    else:
        print("Found a situation where B is not similar: B is not useful")
        logging.info("Found a situation where B is not similar: B is not useful")    
    
        ## Dissimilar: It means start from scratch , keep concurrent and adv = 0
        #new_problem1 = cplex.Cplex(problemFileWithoutBasis)
        #new_problem.parameters.lpmethod.set(new_problem.parameters.lpmethod.values.concurrent)
        #new_problem1.start.read_basis("temp.bas") 
        
        #new_problem1.parameters.advance.set(1)
        #model.parameters.preprocessing.presolve.set(0)
        
        if PRIMAL_FLAG:
            logging.info("current problem is not similar to prev: adv =0, presolve =1, method = concurrent")
            model.parameters.lpmethod.set(model.parameters.lpmethod.values.primal)
        else:
            logging.info("current problem is not similar to prev: adv =0, presolve =1, method = concurrent")
        
            model.parameters.lpmethod.set(model.parameters.lpmethod.values.dual)
            
        #model.parameters.lpmethod.set(model.parameters.lpmethod.values.concurrent)
        model.parameters.advance.set(0)
        
        storeSense = True
    
        return storeSense, model
    
    
    