"""

Second version of MOABC Algorithm based on paper by Reza Akbari et al.

Created by Rafael Hunt-Stokes @06/08/18

"""
import random 
import numpy as np  





store_address = None                         #directory in which output data will be stored
completed_iteration = 0                      #number of completed iterations
completed_percentage = 0.0                   #fraction of optimisation completed
pareto_front = ()                            #current pareto-front with the format (((param1,param2,...),(obj1,obj2,...),(err1,err2,...), (std1, std2, ...)),...)



def nothing_function(data):
    pass

class Optimiser(object):

    def __init__(self, settings_dict, interactor, store_location, param_min, param_max, progress_handler = None):
   
        self.interactor = interactor                                                                  #interactor with dls_optimizer util.py 
        self.store_location = store_location                                                          #location of output file 
        self.pop_size = settings_dict['pop_size']                                                     #number of employed + onlooker bees (2N, N employed, N onlooker) 
        self.max_iter = settings_dict['max_iter']                                                     #number of iterations                                                   
        self.param_min = param_min                                                                    #minimum bounds of parameters
        self.param_max = param_max                                                                    #maximum bounds of parameters
        self.max_trial = settings_dict['max_trials']                                                  #maximum number of iterations food source quality not improved before abandonment 
        self.param_count = len(interactor.param_var_groups)                                           #number of parameters 
        self.obj_count = len(interactor.measurement_vars)                                             #number of objectives
        self.employed_weight = settings_dict['w1']                                                    #importance of randomly selected other food source in calculation of updated positions 
        self.onlooker_weight = settings_dict['w2']                                                    #importance of information given by employed bees to onlookers when onlookers attempt food source improvements  



        if progress_handler == None:
            progress_handler = nothing_function

        self.progress_handler = progress_handler                                                      #window showing progress plots 

        self.pause = False                                                                            #ability to pause and cancel optimisations 
        self.cancel = False
            
        self.add_current_to_individuals = settings_dict['add_current_to_individuals']                 #ability to set current machine status to initial point 
        if self.add_current_to_individuals == True:
            self.initParams = interactor.get_ap()
        else:
            self.initParams = []

    def optimise(self):
        
        self.population = [] 
        for i in range(self.pop_size):

            self.population.append(Bee(self.param_count, self.param_min, self.param_max))

           


class Bee(object): 
    
    def __init__(self, param_count, param_min, param_max):

        self.fitness = None
        self.location = []
        self.initial_locale(param_count, param_min, param_max)
        self.abandonment_count = 0 
        
    
    def initial_locale(param_count, param_min, param_max):
        
                 
        for i in range(param_count):
            x = param_min[i] + random.uniform(0,1) * (param_max[i] - param_min[i]) 
            self.location.append(x)     
               

    








   
 

      
          
        
                 

            
            
                            
                
 

































#---------------------------------------settings window-------------------------------------------------------#


class import_algo_frame(Tkinter.Frame):
    """
    This class produces the MOABC options window that appears after defining parameters, objectives etc.
    """

    def __init__(self, parent):

        Tkinter.Frame.__init__(self, parent)

        self.parent = parent

        self.initUi()

    def initUi(self):
        """
        This initialises the GUI for this window, including text entries for each algorithm option
        """

        self.add_current_to_individuals = Tkinter.BooleanVar(self)
        self.add_current_to_individuals.set(True)

        Tkinter.Label(self, text="Population size:").grid(row=0, column=0, sticky=Tkinter.E)
        self.i0 = Tkinter.Entry(self)
        self.i0.grid(row=0, column=1, sticky=Tkinter.E+Tkinter.W)

        Tkinter.Label(self, text="Max. iterations:").grid(row=1, column=0, sticky=Tkinter.E)
        self.i1 = Tkinter.Entry(self)
        self.i1.grid(row=1, column=1, sticky=Tkinter.E+Tkinter.W)

        Tkinter.Label(self, text="Abandonment count:").grid(row=2, column=0, sticky=Tkinter.E)
        self.i2 = Tkinter.Entry(self)
        self.i2.grid(row=2, column=1, sticky=Tkinter.E+Tkinter.W)

        Tkinter.Label(self, text="Employee update coefficient:").grid(row=3, column=0, sticky=Tkinter.E)
        self.i3 = Tkinter.Entry(self)
        self.i3.grid(row=3, column=1, sticky=Tkinter.E+Tkinter.W)

        Tkinter.Label(self, text="Onlooker update coefficient:").grid(row=4, column=0, sticky=Tkinter.E)
        self.i4 = Tkinter.Entry(self)
        self.i4.grid(row=4, column=1, sticky=Tkinter.E+Tkinter.W)

        self.c0 = Tkinter.Checkbutton(self, text="Use current machine state", variable=self.add_current_to_individuals)
        self.c0.grid(row=5, column=1)

        #Tkinter.Label(self, text="Recommended:\nSwarm Size: 50\nMax. Iterations: 5\nParticle Inertia: 0.5\nSocial Parameter: 1.5\nCognitive Parameter: 2.0", justify=Tkinter.LEFT).grid(row=6, column=0, columnspan=2, sticky=Tkinter.W)

        #self.i0.insert(0, "50")     #defaults are added in case user does not want to decide
        #self.i1.insert(0, "5")
        #self.i2.insert(0, "0.5")
        #self.i3.insert(0, "1.5")
        #self.i4.insert(0, "2.0")



    def get_dict(self):
        """
        Upon clicking OK, errors are lifted if any settings don't work
        """
        setup = {}

        try:
            setup['pop_size'] = int(self.i0.get())
        except:
            raise ValueError("The value for \"pop Size\": \"{0}\", could not be converted to an int".format(self.i0.get()))
        try:
            setup['max_iter'] = int(self.i1.get())
        except:
            raise ValueError("The value for \"Max. Iterations\": \"{0}\", could not be converted to an int".format(self.i1.get()))
        try:
            setup['max_trials'] = float(self.i2.get())
        except:
            raise ValueError("The value for \"Abandonment count\": \"{0}\", could not be converted to a float".format(self.i2.get()))
        try:
            setup['w1'] = float(self.i3.get())
        except:
            raise ValueError("The value for \"Employee update coefficient\": \"{0}\", could not be converted to a float".format(self.i3.get()))
        try:
            setup['w2'] = float(self.i4.get())
        except:
            raise ValueError("The value for \"Onlooker update coefficient\": \"{0}\", could not be converted to a float".format(self.i4.get()))

        if self.add_current_to_individuals.get() == 0:
            setup['add_current_to_individuals'] = False
        elif self.add_current_to_individuals.get() == 1:
            setup['add_current_to_individuals'] = True

        return setup        












