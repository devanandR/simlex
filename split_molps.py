# -*- coding: utf-8 -*-
"""
Created on Sun Mar 16 22:02:00 2025

@author: devii
"""
import cplex
import os


def model_LPs(model:cplex.Cplex, linear_obj_coef):
    
    no_obj = model.multiobj.get_num()
    
    vname_set = model.variables.get_names()
    const_set = model.linear_constraints.get_names()
    vname_up_set = model.variables.get_upper_bounds()
    vname_lb_set = model.variables.get_lower_bounds()
    var_num = model.variables.get_num()
    const_num = model.linear_constraints.get_num()
    const_rhs = model.linear_constraints.get_rhs()
    const_sense = model.linear_constraints.get_senses()
    c_set = model.variables.get_cols()
    
    # Create an empty LP
    lp_model  = cplex.Cplex()
    
    # add variables
    lp_model.linear_constraints.add(names = [const_set[i] for i in range(const_num)], senses= const_sense, 
                                    rhs = const_rhs)
    
    lp_model.variables.add(obj=[0 for i in range(var_num)],names = vname_set, lb = vname_lb_set, ub = vname_up_set 
                           ,columns = c_set # We need to add constraint 
                           )
    #lp_model.variables.add(columns= c_set)
    lp_model.write("dummy_lp.LP")
    
    #print(x_set,"x_set")
    #print(r_set, "r_set")
    #sys.exit(1)
    
    
    LP_Directory_List = ["LP"+str(i)+".LP" for i in range(no_obj)]
    for i in range(no_obj):
        
        obj_coef_list = linear_obj_coef[i]
        n = len(obj_coef_list)
        lp_model.objective.set_linear([(j,obj_coef_list[j]) for j in range(n)])
        lp_model.write(LP_Directory_List[i]+".LP", filetype= "LP")
    
    return LP_Directory_List

def collectSeriesLPs(file):
    
    
    empty_model = cplex.Cplex(file)
    linear_obj_coef = {}
    
    no_obj = empty_model.multiobj.get_num()
    
    for i in range(no_obj):
        
        linear_obj_coef[i] = empty_model.multiobj.get_linear(i)
    
    return model_LPs(empty_model, linear_obj_coef)


File_Location = r"D:\JDA\SimilarityBasedLex\MOLPInstances\MOLPS"
solving_method = 'cplex'
# Update directory
os.chdir(File_Location)


collectSeriesLPs("molp_27_28_218_entropy.LP")