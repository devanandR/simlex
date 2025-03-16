# -*- coding: utf-8 -*-
"""
Created on Sat Sep  4 18:32:25 2021

@author: 1021399
"""


### 
import os
import logging
import sys
from datetime import date
import re
import cplex
import time
from cplex.exceptions import CplexSolverError





def multiobjex1(filenameDict, paramfiles=None):
    
    if debugLex:
        listOfModels = [cplex.Cplex(filenameDict[i]) for i in filenameDict]
    else:
        listOfModels = [filenameDict[i] for i in filenameDict]
    """Solve a multi-objective model."""
    with listOfModels as c:
        if c.multiobj.get_num() == 1:
            print('Model is not multi-objective')
            return
        try:
            if paramfiles:
                # We need to read-in all the parameter sets
                paramsets = [c.create_parameter_set() for i in paramfiles]
                _ = [p.read(f) for p, f in zip(paramsets, paramfiles)]
                c.solve(paramsets=paramsets)
            else:
                c.solve()
        except CplexSolverError:
            print("Exception raised during solve")
            return

        # solution.get_status() returns an integer code
        status = c.solution.get_status()
        print(c.solution.status[status])
        if status == c.solution.status.multiobj_unbounded:
            print("Model is unbounded")
            return
        if status == c.solution.status.multiobj_infeasible:
            print("Model is infeasible")
            return
        if status == c.solution.status.multiobj_inforunbd:
            print("Model is infeasible or unbounded")
            return
        if status == c.solution.status.multiobj_stopped:
            print("Optimization was not finished")
            return


        print("Solution status = ", status, ":", end=' ')
        # the following line prints the status as a string
        print(c.solution.status[status])

        # Now print the values of the various objective for the
        # solution
        print('Objective values...')
        for i in range(c.multiobj.get_num()):
            print("Objective {} value = {}".
                  format(i, c.solution.multiobj.get_objective_value(i)))

        print()

        # Now print the objective values by priorities
        priorities = sorted(set([c.multiobj.get_priority(i)
                                 for i in range(c.multiobj.get_num())]),
                            reverse=True)
        print('Objective values by priorities...')
        objval_by_priority = c.solution.multiobj.get_objval_by_priority
        for p in priorities:
            print("Objective priority {} value = {}".
                  format(p, objval_by_priority(p)))



### Key point : HOw to change variable type  c.variables.set_types(0, c.variables.type.binary)

    

def add_next_objective(prob:cplex.Cplex, newObjLin,indx,fileName,ProblemSize):
    
    prob.multiobj.set_num(indx+1)
    
    ### Set priority
    # Adjust priorities, we want to optimize the first objective
    # first so we give it a higher priority.
    
    # Set the objective name of the second objective (0-indexed).
    prob.multiobj.set_name(indx, fileName)
    
    prob.multiobj.set_linear(indx, newObjLin) ##
    prob.multiobj.set_priority(indx, ProblemSize - indx)
    prob.multiobj.set_abstol(indx, 1e-6)   ### Set abs tolerance 
    prob.multiobj.set_reltol(indx, .000001)
    
    prob.multiobj.set_weight(indx, 1)
    


def solve_orig_problems_with_Var_Fixes(files_to_solve):
    
    dictProblemWiseInformation = {}
    numProbs = len(files_to_solve)
    ## First Problem
    baseProblem  = cplex.Cplex(files_to_solve[0])
    
    ### Important function where we do lot sizing based model modification
    dictProblemWiseInformation[0] = baseProblem
    
    VarNamesInModel = baseProblem.variables.get_names()
    dictCoefPerVars = dict.fromkeys(VarNamesInModel)
    
    baseProblem.multiobj.set_name(0, files_to_solve[0])
    baseProblem.multiobj.set_priority(0, numProbs)
    baseProblem.multiobj.set_abstol(0, 1e-6)   ### Set abs tolerance 
    baseProblem.multiobj.set_reltol(0, .000001)
    
    
    numVars = len(VarNamesInModel)
    for i in range(numVars):
        varName = VarNamesInModel[i]
        dictCoefPerVars[varName] = i
    
    for i in range(1,numProbs):
        fileName = files_to_solve[i]
        tempModel = cplex.Cplex(fileName)
        
        
        #updateLotOpVarInObj(tempModel, lotSizingOperations,OpLotAry)
    
        
        objLin = tempModel.objective.get_linear()
        
        ### variables names
        newObjNames = tempModel.variables.get_names()
        
        
        newObjLin = [(dictCoefPerVars[nm],val) for nm,val in zip(newObjNames,objLin)]
        
        #Add one more objective, multiobj
        add_next_objective(baseProblem, newObjLin,i,files_to_solve[i],len(files_to_solve))
        ##
    if debugLex: 
        baseProblem.write("MultiObj.LP")
      
    return baseProblem


def solveAndlogMultiObj(problem:cplex.Cplex):
    
    
    
    if DEFAULT_RUN:
        logging.info("current problem is solved with default; cplex.solve() no change in solving parameters")
    
    else:
        if WITH_BASIS:
            problem.parameters.advance.set(1)
            problem.parameters.preprocessing.presolve.set(1)
                
            if  CONCURRENT:
                logging.info("current problem is solved with: adv =1, presolve =1, method = concurrent")
                problem.parameters.lpmethod.set(problem.parameters.lpmethod.values.concurrent)
            elif PRIMAL_FLAG:
                if CONCURRENT:
                    logging.info("exiting due to mismatch in param: CONCURRENT and PRIMOPT can not be one together")
                    sys.exit(1)
                logging.info("current problem is solved with: adv =1, presolve =1, method = primal")
                problem.parameters.lpmethod.set(problem.parameters.lpmethod.values.primal)
            else:
                logging.info("current problem is solved with: adv =1, presolve =1, method = dual")
        
                problem.parameters.lpmethod.set(problem.parameters.lpmethod.values.dual)
            
        else:
            problem.parameters.advance.set(0)
            problem.parameters.preprocessing.presolve.set(1)
            
            if  CONCURRENT:
                logging.info("current problem is solved with: adv =0, presolve =1, method = concurrent")
                problem.parameters.lpmethod.set(problem.parameters.lpmethod.values.concurrent)
                
            elif PRIMAL_FLAG:
                if CONCURRENT:
                    logging.info("exiting due to mismatch in param: CONCURRENT and PRIMOPT can not be one together")
                    sys.exit(1)
                        
                logging.info("current problem is solved with: adv =0, presolve =1, method = primal")
                problem.parameters.lpmethod.set(problem.parameters.lpmethod.values.primal)
            else:
                logging.info("current problem is solved with: adv =0, presolve =1, method = dual")
        
                problem.parameters.lpmethod.set(problem.parameters.lpmethod.values.dual)
                
    ### Now solve it 
    
    #Set a time limit
    problem.parameters.timelimit.set(TIMELIM)
    
    #problem.parameters.lpmethod.set(problem.parameters.lpmethod.values.primal)
        
    
    startTime = time.time()
    problem.solve()
    endTime = time.time()
    cplextotal = endTime - startTime

    #logging
    
    #cplex logging 
    #problem.set_results_stream(cplexlog)
    #problem.set_warning_stream(cplexlog)
    #problem.set_error_stream(cplexlog)
    #problem.set_log_stream(cplexlog)

    logging.info("Raw Objective Values : ")
    logging.info("Time Taken : " + str(cplextotal))
    for i in range(0, problem.solution.multiobj.get_num_solves()) :
        logging.info("Objective : " + str(problem.solution.multiobj.get_objective_value(i)))
        
    ### Get the rawValue
    logging.info("And Regular Objective Values : ")
    x = problem.solution.get_values() ## Solution of the problem
    for i in range(0, problem.solution.multiobj.get_num_solves()) :
        objCoef = problem.multiobj.get_linear(i)
        regularValue = sum([x[i] for i in range(len(objCoef)) if objCoef[i]!=0])
        logging.info(" Regular Objective : " + str(regularValue))




if __name__ == "__main__": 
    
    ## Set current working directory:
    CWD = r"D:\JDA\SimilarityBasedLex\MOLPInstances\MOLPS\1"
    os.chdir(CWD)    
    
    debugLex =True
    TIMELIM = 1500.0
    
    PRIMAL_FLAG = True
    CONCURRENT = False
    WITH_BASIS = False
    DEFAULT_RUN = True
    
    
    CONSOLE_LOG_FILE_NAME = "lexicographCplex.log"
    FORMAT = "%(asctime)s %(message)s"
    logging.basicConfig(format=FORMAT , filename=CONSOLE_LOG_FILE_NAME , filemode='w' , level=logging.INFO)

    logging.info("STARTED Cplex Standard Lexicographic Method for Hierarchical MOLP")
    logging.info("SETUP INFORMATION")
    logging.info("------------------------------")
    #logging.info(str(FILES_TO_SOLVE))
    

    FILES_TO_SOLVE=  ["LP0.LP.LP", "LP1.LP.LP", 
                      "LP2.LP.LP", "LP3.LP.LP",
                      "LP4.LP.LP", "LP5.LP.LP",
                      "LP6.LP.LP", "LP7.LP.LP",
                      "LP8.LP.LP", "LP9.LP.LP",
                      ]
    
    
    ### Parsing lpopt_inputfile and option file
    #logging.info("Parsing operation: Using input file and opMap and find lot op")
    logging.info("------------------------------")
    
    time1 = time.time()
    
    finalModel = solve_orig_problems_with_Var_Fixes(FILES_TO_SOLVE)
    
    solveAndlogMultiObj(finalModel)
    
    print("Total Time it takes to finish end to end: ", str(time.time() - time1))
    ###
    #multiobjex1(InfoDict)
    
