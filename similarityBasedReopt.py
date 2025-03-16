# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 10:24:47 2021

@author: 1021399
"""

import os
import cplex
import logging
import time
import sys
import numpy as np
import similaritycomputation as SimLexRC

### DynamicParamSelect.py



### Fix the variables based on the last solution 
def fix_variables(basedictVarValues, lastProblem : cplex.Cplex , lastDictProbVarIndex, lastDictProbLBs , lastDictProbUBs , lastDictVarReducedCosts , lastDictVarBasisStats ,baseVarNames , lastLBFixedStatus , lastUBFixedStatus, debugVariableFixing):
    
    
    upperBoundChangeList = []
    lowerBoundChangeList = []
    
    numChangedUB = 0
    numChangedLB = 0
    fixdvar = 0
    for varName in baseVarNames:        
            reducedCost = lastDictVarReducedCosts[varName] #Assume that the key is available in the reduced cost. Otherwise will use dict.get
            basisStatus = lastDictVarBasisStats[varName]

            lbFixStatus  =  False
            ubFixStatus  =  False

            try:
                lbFixStatus  =  lastLBFixedStatus[varName]
            except BaseException:
                None
            try:
                ubFixStatus  =  lastUBFixedStatus[varName]
            except BaseException:
                None
            
            if(ubFixStatus or lbFixStatus):
                continue
            #if(reducedCost > 0.000000 and basisStatus == 0 ):
            if(reducedCost > 10**-5 and basisStatus == 0 ):
            #if(reducedCost > 10**(-thrsld) and basisStatus == 0 ):
                #fix this variable to the lower bound
                varLB = lastDictProbLBs[varName]
                fixdvar+=1
                upperBoundChangeList.append((lastProblem.variables.get_indices(varName), varLB))
                #lastProblem.variables.set_upper_bounds([(varName , varLB)])
                lastLBFixedStatus[varName] = True
                lastDictProbUBs[varName] =  basedictVarValues[varName]
                #lastDictProbUBs[varName] =  varLB
                numChangedUB = numChangedUB + 1
                if(debugVariableFixing):
                    logging.info(" Variable "  + varName + " , reducedCost : " + str(reducedCost) + ", LB : " + str(lastDictProbLBs[varName]) + " , UB : " +str(lastDictProbUBs[varName])  + " , fixedToLB : " + str(lastDictProbLBs[varName]))
                    
            #if(reducedCost < -0.000000  and basisStatus == 2):
            if(reducedCost < -10**-5  and basisStatus == 2):
            #if(reducedCost < -10**(-thrsld) and basisStatus == 2):
                #fix teh variable to the upper bound
                varUB = lastDictProbUBs[varName]
                lowerBoundChangeList.append((lastProblem.variables.get_indices(varName), varUB))
                #lastProblem.variables.set_lower_bounds([(varName , varUB)])
                fixdvar+=1
                lastUBFixedStatus[varName] = True
                lastDictProbLBs[varName] =  basedictVarValues[varName]
                #lastDictProbLBs[varName] =  varUB
                numChangedLB = numChangedLB + 1
                if(debugVariableFixing):
                    logging.info(" Variable "  + varName + " , reducedCost : " + str(reducedCost) + ", LB : " + str(lastDictProbLBs[varName]) + " , UB : " +str(lastDictProbUBs[varName])  + " , fixedToUB : " + str(lastDictProbUBs[varName]))
            
               
    if(numChangedUB > 0) :
        lastProblem.variables.set_upper_bounds(upperBoundChangeList)
    if(numChangedLB > 0):
        lastProblem.variables.set_lower_bounds(lowerBoundChangeList)
    #print(str(time.ctime())+fixdvar,"fixedVar")
    return (lastDictProbLBs , lastDictProbUBs , lastLBFixedStatus , lastUBFixedStatus)

def solve_orig_problems_with_Var_Fixes(files_to_solve , debug=False, fileTypeCplex = False, fix_var= True):
    dictProblemWiseInformation = {}
    numProbs = len(files_to_solve)
    lastProblem = None
    ## First Problem
    baseNumVars = 0
    globalUBFixedStatus = None
    globalUBFixedStatus = None
    if fileTypeCplex:
        baseProblem  = files_to_solve[0]
    else:
        baseProblem  = cplex.Cplex(files_to_solve[0])
    
    
    #baseProblem = changeModelToMin(baseProblem)
    # if baseProblem.objective.sense[baseProblem.objective.get_sense()] =='maximize':
    #     v = baseProblem.objective.get_linear()
    #     baseProblem.objective.set_linear([(i,-v[i]) for i in range(len(v))])
    #     baseProblem.objective.set_sense(baseProblem.objective.sense.minimize)
    
    lastProblem = baseProblem    
    if(debug):
        if fileTypeCplex:
            baseProblem.write("matrix.orig."+files_to_solve[0].get_problem_name()+"Obcrunch.LP")    
        else:
            baseProblem.write("matrix.orig."+files_to_solve[0]+".LP")
    
    iterationN = 0
    objValue ,baseVarNames, baseDictVarValues  , baseDictVarIndex , baseDictVarLBs , baseDictVarUBs , baseDictObjCoefs ,solutionTime , numIteration, baseDictVarReducedCosts , baseDictBasisStatus = solve_and_collect_stats(baseProblem,iterationN)
    dictVarValues = baseDictVarValues
    dictVarLBs  = baseDictVarLBs
    dictVarUBs = baseDictVarUBs
    dictVarReducedCosts  = baseDictVarReducedCosts
    dictVarBasisStatus =  baseDictBasisStatus
    # intialize it:Dev
    dictLBFixedStatus = dict.fromkeys(baseVarNames)
    dictUBFixedStatus = dict.fromkeys(baseVarNames)
    
    if fix_var:
        dictVarLBs , dictVarUBs , dictLBFixedStatus , dictUBFixedStatus  = fix_variables(dictVarValues,lastProblem ,baseDictVarIndex, dictVarLBs , dictVarUBs , dictVarReducedCosts , dictVarBasisStatus ,baseVarNames , dictLBFixedStatus , dictUBFixedStatus , debug and DEBUG_VAR_FIX)    
    
    dictProblemWiseInformation[0] = {}
    
    dictProblemWiseInformation[0]["objValue"] = objValue 
    dictProblemWiseInformation[0]["solutionTime"] = solutionTime 
    dictProblemWiseInformation[0]["iteration"] = numIteration
    
    #dictProblemWiseInformation[0]["dictVarValues"] = baseDictVarValues 
    
    for i in range(1,numProbs):
      iterationN+=1
      fileName = files_to_solve[i]
      print("\nfile to solve", fileName)
      print()
      objValue ,dictVarValues  , dictObjCoefs  , solutionTime  , numIteration, dictVarReducedCosts , dictVarBasisStatus = solve_next_problem(lastProblem ,fileTypeCplex, baseVarNames ,   baseDictVarIndex, dictLBFixedStatus , dictUBFixedStatus,fileName , i ,iterationN, debug , matrixFileNamePrefix="matrix.orig.")

      dictProblemWiseInformation[i] = {}
      dictProblemWiseInformation[i]["objValue"] = objValue 
      dictProblemWiseInformation[i]["solutionTime"] = solutionTime 
      dictProblemWiseInformation[i]["iteration"] = numIteration
      #dictProblemWiseInformation[i]["dictVarValues"] = dictVarValues 
      
      if fix_var:
          dictProbLBs , dictProbUBs , dictLBFixedStatus , dictUBFixedStatus = fix_variables(dictVarValues, lastProblem , baseDictVarIndex,dictVarLBs , dictVarUBs , dictVarReducedCosts , dictVarBasisStatus ,baseVarNames , dictLBFixedStatus , dictUBFixedStatus , debug and DEBUG_VAR_FIX)
      
    return dictProblemWiseInformation


def solve_and_collect_stats(problem: cplex.Cplex,iterationN ):
        varNames , dictVarIndex , dictVarLBs, dictVarUBs , dictObjCoefs = collect_problem_stats(problem)
        solutionReturnCode , solutionTime, numIteration = solve_problem(problem,iterationN)
        objValue , varValues  , varReducedCosts , basisStats= collect_solution_stats(problem, varNames)
        logging.info("--------------- Solution for " + str(problem))
        logging.info("Solution Status "  + str(solutionReturnCode) + " ObjValue " + str(objValue))
        logging.info("---------------")
        return (objValue ,varNames ,varValues , dictVarIndex , dictVarLBs , dictVarUBs , dictObjCoefs ,solutionTime , numIteration, varReducedCosts , basisStats)


    
def computeNONZandFixedVars(model):
    
    U = model.variables.get_upper_bounds()
    L = model.variables.get_lower_bounds()
    FixedVars = [i for i in range(len(U)) if U[i] == L[i]]

    #### Abt nonzero Pair
    ObjCoefInModel = model.objective.get_linear()
    VarNamesInModel = model.variables.get_names()
    dictCoefPerVars = dict.fromkeys(VarNamesInModel)
    
    numVars = len(VarNamesInModel)
    for i in range(numVars):
        varName = VarNamesInModel[i]
        dictCoefPerVars[varName] = ObjCoefInModel[i]
        
    return len(FixedVars),dictCoefPerVars

    

def nonZeroDev(nonZeroDicCurrent, nonZeroDicPrev):
    
    nonZUnique = 0
    
    if nonZeroDicPrev=={}:
        return 0
    for vName in nonZeroDicCurrent:
        if abs(nonZeroDicCurrent[vName]) > 0. :
            if not abs(nonZeroDicPrev[vName]) > 0.:
                nonZUnique+=1            
    
    ## nonZUnique = zero coef in prev run now becomes nonzero in current run.
    ## counterNonPrev = nonzero coef in prev run now becomes zero in current run
    counterNonPrev = len([1 for vName in nonZeroDicPrev if abs(nonZeroDicPrev[vName])>0. 
                          and (not abs(nonZeroDicCurrent[vName]) > 0)])    
    logging.info("--------------- current unique nonzero and previous nonzero in ObjFunction" + 
                 str(nonZUnique) +" and " + str(counterNonPrev))
    
    coefDiff = (counterNonPrev + nonZUnique)/len(nonZeroDicCurrent.keys())
    return coefDiff # dissimilar
    #return 100*(min(nonZUnique,counterNonPrev))/(nonZUnique + counterNonPrev)
        
    
#return solution status and the time taken
def solve_problem(problem : cplex.Cplex,iterationN):
    
    logging.info("Solving Problem ...")
    global PRIMAL_FLAG
    global basisAndrcShouldStoreForSim
    global PREV_VAR_NAMES
    global PREV_FIXED_VARS
    global PREV_COST
    global storeSense
    global CONCURRENT_FLAG
    global PREVBASISINV
    global DUAL_VALUE
    global COLUMNS
    global prevBasisStatusForSim
    global prevRCForSim
    global SOLVED_TIME, TIMELIM
    global BINVERSE_A
    global BASIC_VAR
    global prevObjSol
    global WITH_BASIS
    
    logging.info("----------------------------------------------")
    #problem.parameters.simplex.tolerances.markowitz.set(0.001)
    #problem.parameters.simplex.tolerances.optimality.set(1E-6)
    #problem.parameters.simplex.tolerances.feasibility.set(1E-6)
    #if iterationN == 1 or iterationN == 2 or iterationN == 3 or iterationN == 4 or iterationN == 6 or iterationN == 7 or iterationN == 8 or iterationN == 9:
    if iterationN == 0:
        storeSense = True
        logging.info("solving first problem")
        problem.parameters.advance.set(0)
        #problem.parameters.preprocessing.presolve.set(0)
        
        #problem.parameters.preprocessing.dual.set(1)
        current_fixed_vars,currentNonZ = computeNONZandFixedVars(problem)
        PREV_FIXED_VARS = current_fixed_vars
        PREV_COST = currentNonZ
        
        if not CONCURRENT_FLAG:
        
            if PRIMAL_FLAG:
            
                problem.parameters.lpmethod.set(problem.parameters.lpmethod.values.primal)
                logging.info("solving method: primal ")
            else:
                
                problem.parameters.lpmethod.set(problem.parameters.lpmethod.values.dual)
                logging.info("solving method: dual ")
        else:
            #if DYNAMIC_FLAG:
            #    logging.info("exiting due to flag mismatch: concurrent can not be enabled with dynamic ")
            #    print("exiting due to flag mismatch: concurrent can not be enabled with dynamic ")
            #    sys.exit(1)
                
            problem.parameters.lpmethod.set(problem.parameters.lpmethod.values.concurrent)
            logging.info("solving method: concurrent ")
            
        ### fixed variables in current problem. 
        #problem.set_results_stream(None)
        #problem.set_log_stream(None)
        
        startTime = time.time()
        problem.solve()
        endTime = time.time()
        problem.solution.basis.write("basis.bas")
        
        numIteration = problem.solution.progress.get_num_iterations()

    else:
        if DYNAMIC_FLAG:  ### For Dynamic flag enabled, concurrent flag has to disabled 
            
            
                #logging.info("exiting due to flag mismatch: concurrent can not be enabled with dynamic ")
                #print("exiting due to flag mismatch: concurrent can not be enabled with dynamic ")
                #sys.exit(1)
        
            if similarityType == 1:        
                ### Set cplex parameters as per current_fixed_vars
                storeSense = SimLexRC.computeSimSimple(problem, PREV_COST,logging, PRIMAL_FLAG)
                #updateCplexParamS1(current_fixed_vars,currentNonZ, problem)
            
            elif similarityType==2: # RC based similarity
                storeSense = SimLexRC.computeSimRCBAsed(problem, BINVERSE_A, PREV_COST, COLUMNS,
                                           prevBasisStatusForSim,[],logging, PRIMAL_FLAG)
            elif similarityType == 3:
                storeSense = SimLexRC.computeSimDualBAsed(problem,PREV_COST, DUAL_VALUE, COLUMNS,
                                           prevBasisStatusForSim,[],logging, PRIMAL_FLAG)
            elif similarityType == 4:
                problem.write("WithOptBasis.MPS.gz")
                problem.write("WithoutOptBasis.MPS.gz")
                storeSense, problem = SimLexRC.basicObjTric(problem, prevRCForSim,PREV_COST,\
                                                          "WithOptBasis.MPS.gz","WithoutOptBasis.MPS.gz",logging, PRIMAL_FLAG)
            if CONCURRENT_FLAG:
                problem.parameters.lpmethod.set(problem.parameters.lpmethod.values.concurrent)
                logging.info("solving method: concurrent ")
            
        else:
            if not WITH_BASIS:
                problem.parameters.advance.set(0)
            else:
                problem.parameters.advance.set(1)
            #problem.parameters.preprocessing.presolve.set(1)
            #problem.parameters.lpmethod.set(problem.parameters.lpmethod.values.primal)
                    
            if PRIMAL_FLAG:
                problem.parameters.lpmethod.set(problem.parameters.lpmethod.values.primal)
                logging.info("solving method: primal ")
            else:
                problem.parameters.lpmethod.set(problem.parameters.lpmethod.values.dual)
                logging.info("solving method: dual ")
            
            if CONCURRENT_FLAG:
                problem.parameters.lpmethod.set(problem.parameters.lpmethod.values.concurrent)
                logging.info("solving method: concurrent ")
            
        #if not CONCURRENT_FLAG:
        #    
        #    if PRIMAL_FLAG:
        #        problem.parameters.lpmethod.set(problem.parameters.lpmethod.values.primal)
        #        logging.info("solving method: primal ")
        #    else:
        #        problem.parameters.lpmethod.set(problem.parameters.lpmethod.values.dual)
        #        logging.info("solving method: dual ")
        
        #problem.set_results_stream(None)
        #problem.set_log_stream(None)
        #else:
        #    problem.parameters.preprocessing.presolve.set(1)
            
        startTime = time.time()
        problem.solve()
        endTime = time.time()
        numIteration = problem.solution.progress.get_num_iterations()
        problem.solution.basis.write("basis.bas")
        
        if not PRIMAL_FLAG:
            PRIMAL_FLAG = True
    
    print("time taken in solving problem number ", iterationN, " is", endTime - startTime)
    print()
    print("....model solved....")
    
    ### Record parameters for Similarity computation 
    SOLVED_TIME += (endTime - startTime)
    
    if SOLVED_TIME >= TIMELIM:
        logging.info("-------------------------------------------")
        logging.info("time out while solving problem with iteretion no. " + str(iterationN))
        print("time out while solving " + str(iterationN))    
        sys.exit(1)        
    
    if storeSense:    
        
        if  similarityType == 4:
            prevObjSol = problem.solution.get_values()
             #### Abt nonzero Pair
            ObjCoefInModel = problem.objective.get_linear()
        
            PREV_VAR_NAMES = problem.variables.get_names() 
            numVars = len(PREV_VAR_NAMES)
        
            for j in range(numVars):
                varName = PREV_VAR_NAMES[j]
                PREV_COST[varName] = ObjCoefInModel[j]
        else:
            prevObjSol = problem.solution.get_values()
             #### Abt nonzero Pair
            ObjCoefInModel = problem.objective.get_linear()
            PREV_VAR_NAMES = problem.variables.get_names() 
            numVars = len(PREV_VAR_NAMES)
        
            # get dual info
            #DUAL_VALUE = problem.solution.get_dual_values()        
            COLUMNS = fetchColumns(problem)                
            #PREVBASISINV = problem.solution.advanced.binvcol()                
            PrevBinvPerCol = {}
            varSolutionBasis = problem.solution.basis.get_basis()[0]
            for j in range(numVars):
                varName = PREV_VAR_NAMES[j]
                if varSolutionBasis[j]==1:
                    BASIC_VAR.append(j)
                PREV_COST[varName] = ObjCoefInModel[j]
                if varSolutionBasis[j] ==0 or varSolutionBasis[j]==2: #for A_j
                    BINVERSE_A[j], PrevBinvPerCol = getBinverseAj(problem, COLUMNS[j], 
                                                                  PrevBinvPerCol)
            PrevBinvPerCol = {}
            
            #rc_comp = time.time()
            #CbTransBinv  = getcBTransposeAinverse(problem, ObjCoefInModel,BASIC_VAR)
            #print(time.time() - rc_comp, "is the time in CbTransBinv computation")
            #Binv = problem.solution.advanced.binvrow()
            #problem.solution.advanced.
            #nonbasic =0
            
        #sys.exit(1)
        basisAndrcShouldStoreForSim = True
        storeSense = False
    
    #print("Cost C and Cb",PREV_COST)
    print()
    #print("index of B",BASIC_VAR)
    #problem.solution.quality_metric.max_dual_residual;
    logging.info("Done Solving Problem , Solution Time  " + str(endTime -startTime) + " : objValue  "+ str(problem.solution.get_objective_value()))
    
    return problem.solution.get_status(), (endTime - startTime), numIteration

def getBinverseAj(problem:cplex.Cplex, prevColj,PrevBinvPerCol):
    
    #total_var = problem.variables.get_num()
    total_const = problem.linear_constraints.get_num()
        
    nonzindex =  prevColj.ind
    nonzval = prevColj.val
            
    B_invAj = [0]*total_const
    for indx in range(len(nonzindex)):
        rowId = nonzindex[indx]
        
        B_invAjtemp =[]
        
        if rowId in PrevBinvPerCol:
            B_invAjtemp = PrevBinvPerCol[rowId]
        else:
            B_invAjtemp = problem.solution.advanced.binvcol(rowId)
            PrevBinvPerCol[rowId] = B_invAjtemp 
        B_invAjtemp = [x*nonzval[indx] for x in B_invAjtemp]
        
        B_invAj = [x + y for x, y in zip(B_invAjtemp,B_invAj) ]

    return B_invAj, PrevBinvPerCol



def getcBTransposeAinverse(problem:cplex.Cplex, PrevCost,BASIC_VAR):
    
    #total_var = problem.variables.get_num()
    total_const = problem.linear_constraints.get_num()
        
            
    CbB_inv = [0]*total_const
    
    for indx ,i in zip(BASIC_VAR, range(len(BASIC_VAR))):
        CB_invtemp = []
        temp = PrevCost[indx]
        
        print("Hi")
        if temp !=0:
            CB_invtemp  = [temp*x for x in problem.solution.advanced.binvrow(i)]
        
        CbB_inv = [x + y for x, y in zip(CB_invtemp,CbB_inv)]

    return CbB_inv



def fetchColumns(problem:cplex.Cplex):
    
    COLUMNS = {}
    columns_temp = problem.variables.get_cols()
    return columns_temp
    row_indx = [pair.ind for pair in columns_temp]
    col_coef = [pair.val for pair in columns_temp]
    
    for var,i in zip(PREV_VAR_NAMES, range(len(PREV_VAR_NAMES))):
        
        COLUMNS[var] = []
        temp_col_coef = col_coef[i]
        #temp_row_indx  = row_indx[i]
        
        for r in range(problem.linear_constraints.get_num()):
            if r in row_indx[i]:
                temp_val = temp_col_coef[row_indx[i].index(r)]
                COLUMNS[var].append(temp_val)
            else:
                COLUMNS[var].append(0)
    return COLUMNS
    
def collect_solution_stats(problem : cplex.Cplex , varNames ):
    
    global basisAndrcShouldStoreForSim
    global prevBasisStatusForSim
    global prevRCForSim
    
    if(varNames == None):
        varNames = problem.variables.get_names()
    objValue = problem.solution.get_objective_value()
#    logging.info("TSDEBUG  collect_solution_stats" )
#    logging.info("objCoef " + str(objValue))
#    logging.info("varNames length " + str(len(varNames)))

    varValues = problem.solution.get_values()
    varReducedCosts = problem.solution.get_reduced_costs()
    varSolutionBasis = problem.solution.basis.get_basis()[0]
    
    
    dictVarValues = dict.fromkeys(varNames)
    dictVarReducedCosts = dict.fromkeys(varNames)
    dictVarBasisStatus = dict.fromkeys(varNames)
    
    numVars = len(varNames)
    for i in range(numVars):
        varName = varNames[i]
        dictVarValues[varName] = varValues[i]
        dictVarReducedCosts[varName] = varReducedCosts[i]
        dictVarBasisStatus[varName] = varSolutionBasis[i]
        if basisAndrcShouldStoreForSim and similarityType != 1:
            prevRCForSim[varName] = varReducedCosts[i]
            prevBasisStatusForSim[varName] = varSolutionBasis[i]
    if prevBasisStatusForSim:
        basisAndrcShouldStoreForSim = False
    
    
    return objValue , dictVarValues , dictVarReducedCosts , dictVarBasisStatus
    
 
# return varNames , dictVarIndex , dictVarLBs, dictVarUBs , dictObjCoefs
def collect_problem_stats(problem : cplex.Cplex , varNames = None, onlyObjCoefs = False):
    logging.info("collect_problem_stats" )
    problemVars = problem.variables
    numVars = problemVars.get_num()
    
    #if(varNames == None):
    #Commented out the above line, this performance trick did not work as the indices would be different for different problems
    varNames = problemVars.get_names()
    logging.info("\t collected varNames" )
    
    
    varIndices = list( i for i in range(0,numVars))

    if(not onlyObjCoefs):
        #varIndices = problemVars.get_indices()
        #contlogging.info("\t collected varIndices" )
    
        varLBs = problemVars.get_lower_bounds()
        logging.info("\t collected varLBs" )
    
        varUBs = problemVars.get_upper_bounds()
        logging.info("\t collected varUBs" )

    
    varObjCoefs = problem.objective.get_linear()
    logging.info("\t collected varObjCoefs" )
  
    dictObjCoefs = dict.fromkeys(varNames)
    
    dictVarIndex = None
    dictVarLBs = None
    dictVarUBs = None

    if(not onlyObjCoefs):
        dictVarIndex = dict.fromkeys(varNames)
        dictVarLBs = dict.fromkeys(varNames)
        dictVarUBs = dict.fromkeys(varNames)

    for i in range(numVars):
        name = varNames[i]
        objCoef = varObjCoefs[i]
        dictObjCoefs[name] = objCoef
        if(not onlyObjCoefs):
            lb = varLBs[i]
            ub = varUBs[i]
            index = varIndices[i]
            dictVarIndex[name] = index
            dictVarLBs[name] = lb
            dictVarUBs[name] = ub
    logging.info("\t collected UB,LB,Index" )
    return varNames , dictVarIndex , dictVarLBs, dictVarUBs , dictObjCoefs

def changeModelToMin(model):
    
    Change = False
    if model.objective.sense[model.objective.get_sense()] =='maximize':
        v = model.objective.get_linear()
        model.objective.set_linear([(i,-v[i]) for i in range(len(v))])
        model.objective.set_sense(model.objective.sense.minimize)
        Change = True
    
    return model, Change

###SOLVE THE NEXT PROBLEM WITH THIS AS THE BASE PROBLEM
#return lastProblem Stats , lastSolution Stats
def solve_next_problem(lastProblem : cplex.Cplex , problemFileCplex, baseVarNames , baseVarIndices,dictLBFixedStatus , dictUBFixedStatus,  nextFileName , problemIndex  ,iterationN, debug= False , matrixFileNamePrefix = "matrix.orig."):

    
    if problemFileCplex:
        newProblem = nextFileName
        newProbName = nextFileName.get_problem_name()
        logging.info("#SOLVE_NEXT_PROBLEM for " + newProbName)
    else:
        newProblem = cplex.Cplex(nextFileName)
        logging.info("#SOLVE_NEXT_PROBLEM for " + nextFileName)  
    
    
    #newProblem = changeModelToMin(newProblem)
    
    problemVars = newProblem.variables
    aggregatedVarNames  = problemVars.get_names()  
    varNames , dictVarIndex , dictVarLBs, dictVarUBs , dictObjCoefs  = collect_problem_stats(newProblem , varNames = aggregatedVarNames, onlyObjCoefs=True)
    if problemFileCplex:
        logging.info("\tDone calling collect_problem_stats " + newProbName)        
        logging.info("\tDone calling collect_problem_stats " + newProbName)        
    else:
        logging.info("\tDone calling collect_problem_stats " + nextFileName)        
        logging.info("\tDone calling collect_problem_stats " + nextFileName)        

    #newProblem.end() ## Dev Commented
    if problemFileCplex:
        logging.info("\tEnded the newProblem for " + newProbName)        
    else:
        logging.info("\tEnded the newProblem for " + nextFileName)        
    varCounter = 0
    newObjCoefs = list()
    zeroList = list()        

    for varName in baseVarIndices:
        
        newObjCoef = 0.
        zeroList.append( (varCounter , 0.0))
        
        #if varName in dictLBFixedStatus or varName in dictUBFixedStatus:
        #    if dictLBFixedStatus[varName]==True or dictUBFixedStatus[varName]==True:
        #        varCounter = varCounter + 1
        #        continue
        
        try:
            newObjCoef = dictObjCoefs[varName]
        except BaseException:
            None
        if(newObjCoef != 0 and newObjCoef!=None):
            #newObjCoefs.append((varName , newObjCoef)) 
            #varIndex = baseVarIndices[varName]
            ## make the new objcoef index based rather than the name based, it was taking a lot of time
            #newObjCoefs.append((varName , newObjCoef))
            newObjCoefs.append((lastProblem.variables.get_indices(varName), newObjCoef))
        varCounter = varCounter + 1
    if problemFileCplex:
        logging.info("\tCollected changes to be made for next " + newProbName)        
        logging.info("\tGoing to set zeroList " + newProbName)        
        lastProblem.objective.set_linear(zeroList)
        logging.info("\tGoing to set new ObjectiveFunction " + newProbName)
    else:
        logging.info("\tCollected changes to be made for next " + nextFileName)        
        logging.info("\tGoing to set zeroList " + nextFileName)        
        lastProblem.objective.set_linear(zeroList)
        logging.info("\tGoing to set new ObjectiveFunction " + nextFileName)
    
    if(len(newObjCoefs)==0):
	    return (0,{},{}, 0,{},{})
    lastProblem.objective.set_linear(newObjCoefs) ###VILLAIN LINE
    """
    for newCoefEntry in newObjCoefs:
        lastProblem.objective.set_linear(newCoefEntry[0], newCoefEntry[1])
    """
    if problemFileCplex:
        logging.info("\tDone setting new ObjectiveFunction " + newProbName)
        logging.info("#\tSolving Problem for " + newProbName)
    else:
        logging.info("\tDone setting new ObjectiveFunction " + nextFileName)
        logging.info("#\tSolving Problem for " + nextFileName)
    
    if(debug):
        if problemFileCplex:
            lastProblem.write(matrixFileNamePrefix + newProbName +str(problemIndex)+ "Obcrunch.LP")
        else:
            lastProblem.write(matrixFileNamePrefix + nextFileName + "Obcrunch.LP")
    solutionReturnCode ,  solveTime, numIteration = solve_problem(lastProblem,iterationN)
    
    if problemFileCplex:
        logging.info("#\tSolved Problem " + newProbName)

        logging.info("#\tCollecting Solution Stats for  " + newProbName)
    else:
        logging.info("#\tSolved Problem " + nextFileName)

        logging.info("#\tCollecting Solution Stats for  " + nextFileName)
    
    objValue , dictVarValues  ,  dictVarReducedCosts , dictVarBasisStatus = collect_solution_stats(lastProblem, baseVarNames)        
    
    if problemFileCplex:
        logging.info("#\tCollected Solution Stats for  " + newProbName)
        logging.info("#\t#DONE SOLVE_NEXT_PROBLEM  " + newProbName)
    else:
        logging.info("#\tCollected Solution Stats for  " + nextFileName)
        logging.info("#\t#DONE SOLVE_NEXT_PROBLEM  " + nextFileName)

    return (objValue ,dictVarValues  , dictObjCoefs  , solveTime  , numIteration, dictVarReducedCosts , dictVarBasisStatus)



if __name__ == "__main__": 
    
    ## Set current working directory:
    CWD = r"D:\JDA\SimilarityBasedLex\MOLPInstances\MOLPS\3"
    os.chdir(CWD)    
    
    ### Type of Similarity based Lexicographic method
    similarityType = 4
    
    PRIMAL_FLAG = False
    WITH_BASIS = False
    DYNAMIC_FLAG = False # Simlex : True   
    fix_var = True
    CONCURRENT_FLAG = True
    
    
    PREV_FIXED_VARS = 0
    PREV_COST ={}
    DEBUG_VAR_FIX = True
    PREV_VAR_NAMES =[]
    SOLVED_TIME = 0.
    TIMELIM = 7200.
    basisAndrcShouldStoreForSim = False
    storeSense  = False
    prevRCForSim = {}
    prevBasisStatusForSim = {}
    PREVBASISINV = []
    prevObjSol =0.
    COLUMNS = {}
    DUAL_VALUE = []
    BINVERSE_A = {}
    BASIC_VAR = []
    if similarityType ==1:
        CONSOLE_LOG_FILE_NAME = "dynamicSettings_sim_type1.log"
    elif similarityType ==2:
        CONSOLE_LOG_FILE_NAME = "dynamicSettings_sim_type2.log"
    elif similarityType ==3:
        CONSOLE_LOG_FILE_NAME = "dynamicSettings_sim_type3.log"
    elif similarityType ==4:
        CONSOLE_LOG_FILE_NAME = "dynamicSettings_sim_type4.log"
    elif DYNAMIC_FLAG==False:
        CONSOLE_LOG_FILE_NAME = "nonDynamic"
    FORMAT = "%(asctime)s %(message)s"
    logging.basicConfig(format=FORMAT , filename=CONSOLE_LOG_FILE_NAME , filemode='w' , level=logging.INFO)

    logging.info("STARTED THE Dynamic Parameter select and solve")
    logging.info("SETUP INFORMATION")
    logging.info("------------------------------")
    #logging.info(str(FILES_TO_SOLVE))
    logging.info("------------------------------")


    FILES_TO_SOLVE=  ["LP0.LP.LP", "LP1.LP.LP", 
                      "LP2.LP.LP"]
    #, "LP3.LP.LP",
     #                 "LP4.LP.LP", "LP5.LP.LP",
      #                "LP6.LP.LP", "LP7.LP.LP",
       #               "LP8.LP.LP", "LP9.LP.LP",
        #              "LP10.LP.LP"]
                      # "molp_10_900_60_assignment2.LP",
                      # "molp_10_900_60_assignment3.LP",
                      # "molp_10_900_60_assignment4.LP",
                      # "molp_10_900_60_assignment5.LP",
                      # "molp_10_900_60_assignment6.LP",
                      # "molp_10_900_60_assignment7.LP",
                      # "molp_10_900_60_assignment8.LP",
                      # "molp_10_900_60_assignment9.LP",
                      # "molp_10_900_60_assignment10.LP"]
    
    for file in FILES_TO_SOLVE:
        model, change_flag = changeModelToMin(cplex.Cplex(file))
        if change_flag:
            model.write(file)
    
    debug = False
    fileTypeCplex = False
    
    InfoDict = solve_orig_problems_with_Var_Fixes(FILES_TO_SOLVE,debug,fileTypeCplex, fix_var)
    
    # Print the header
    print(f"{'Index':<10} {'objValue':<20} {'solutionTime':<20} {'iteration':<20}")
    #print("-" * 50)
    # Print the data in column format
    for index, values in InfoDict.items():
        print(f"{index:<10} {values['objValue']:<20} {values['solutionTime']:<20} {values['iteration']:<20}")
    
    #print(InfoDict)
    print()
    total_time = sum([InfoDict[i]["solutionTime"] for i in range(len(InfoDict))]) 
    print("total time taken ",total_time)
    print()
    total_iteration = sum([InfoDict[i]["iteration"] for i in range(len(InfoDict))]) 
    print("total_iteration ",total_iteration)
    
    
    
