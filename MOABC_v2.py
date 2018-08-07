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
        self.pop_size = settings_dict['pop_size']                                                     #number of employed bees (N employed, N onlooker) 
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





    def measure_nectar(self, bee_locations):
        """Function that performs objective measurement of Bees"""
            
        self.interactor.set_ap(bee_locations)                                     #configure machine for measurement
        all_data = self.interactor.get_ar()                                       #perform measurement
        all_results = [i.mean for i in all_data]                                  #retrieve mean from measurement
                        

            
            
        return all_results
 


    def dominance_check(self,a,b):
    """
        Function determines which of two points is the dominant in objective space.

        Args:
            a: list of objective values [obj1,obj2,...]
            b: list of objective values [obj1,obj2,...]

        Returns:
            Function will return the point that dominates the other. If neither dominates, return is False.
    """
        if all(a_i > b_i for (a_i,b_i) in zip(a,b)):         #does a dominate b?
            return a
        if all(a_i < b_i for (a_i,b_i) in zip(a,b)):         #does b dominate b?
            return b
        if all(a_i == b_i for (a_i,b_i) in zip(a,b)):        #are the points the same?
            return b
        else:
            return False
            
    def num_dominated(self, index):
        """
        Function to iterate over dominance check function for each point and return number of points indexed position dominates.
        Used in onlooker methods to calculate roulette wheel selection probabilities.
        """
        
        numb_dominated = 0 
        curr_objectives = list(np.loadtxt("./external_archive")) 
        
        point_selected = curr_objectives[index]
        
        for i in range(len(self.population)): 
            test_against = curr_objectives[i]
            
            if test_against == point_selected: 
                break 
            else: 
                dominator = dominance_check(point_selected, test_against) 

            if dominator == point_selected:
                numb_dominated += 1 
            
        return numb_dominated 
        
    def calc_fitness(self, num_bees): 
        """
        Calculates fitness values of each point in archive and returns list. Fitness is number of sources dominated 
        by source i / number of sources.
        """
        tot_fit = [] 
        for i in range(len(self.population)):
            number_dominated = num_dominated(i) 
            fit = number_dominated / num_bees 
            tot_fit.append(fit) 


        return tot_fit 


    def calc_probability(self, index): 
        """
        calculates pk fitness values using roulette wheel selection method. Used within onlooker bee methods
        """
        all_fit = calc_fitness(len(self.population)) 
          
        pk = all_fit[index] / sum(all_fit)

        return pk  


    def send_scouts(self, index): 
    """
        Called by send_employees() & send_onlookers if abandonment count exceeds max_trials.
        Initialises new random solution within search space.  

    """
       
        
        #randomly produce new food source
        self.population[index] = Bee(self.param_count, self.param_min, self.param_max)

        
        






    def send_employees(self, curr_nectar):
            
        current_pos = [Bee.location for Bee in self.population] 
             
        for i in range(len(current_positions)):
            rand_param = randint(0, self.param_count -1) 
        
                
            neighbour_index = randint(0, len(current_pos) - 1)
            while current_pos[neighbour_index] == current_pos[i]:
                    
                neighbour_index = randint(0, len(current_pos) - 1)

            #takes randomly selected parameter from Bee and randomly selected neighbour and alters it  
            current_pos[i][rand_param] = current_pos[i][rand_param] + self.employed_weight * random.uniform(0,1) * (current_pos[i][rand_param] - current_pos[neighbour_index][rand_param]) 

            #next, check whether nectar amount has been improved. If not, increase abandonment count
            new_nectar = measure_nectar(current_pos[i]) 
                
            old_nectar = curr_nectar[i] 

            x = dominance_check(new_nectar, old_nectar)
               
            if x == new_nectar:
                
                curr_nectar[i] = new_nectar 
                
                update_archive(curr_nectar)       
  
            else: 
                #increase abandoment count; if exceeds max_trials send scouts  
                self.population.abandonment_count[i] += 1
                
                if self.population.abandonment_count[i] > self.max_trial:
                    send_scouts(i)  
            
            #finally, update bee's position in population 
            self.population.location[i] = current_pos[i]

    def roulette_selection(self, probs):
        """
        Called by onlookers function to select neighbour in local search eqn.
        """

        max_prob = sum(probs) 
        pick = random.uniform(0, max_prob) 
        current = 0 
        for i in range(len(probs)): 
            current += probs[i]
            if current > pick:
                return i  
        
        
    def send_onlookers(self):
        """
            Function essentially does the same as employed bees, but prioritises local search about 
            food source positions with a higher fitness value. 
            Fitness of food source is calculated using:

            pk = fit(xk) / SUM{fit(xm)}       
                
            where fit(xm) = dom(m) / pop_size & dom(m) = function returning number of sources dominated 
            by source m. 
        """
        
        current_pos = [Bee.location for Bee in self.population]
        
        probs_points = [] 
        for i in range(len(self.population)):
            
            pk = calc_probability(i) 
            probs_points.append(pk) 

        #now we have all the probabilities in a list, we can pass this to roulette selection function
        for i in range(len(current_pos)):
            neighbour_index = roulette_selection(probs_points) 
            rand_param = randint(0, self.param_count - 1) 
            while neighbour_index == i:
                rand_param = randint(0, self.param_count - 1) 

            current_pos[i][rand_param] = current_pos[i][rand_param] + self.onlooker_weight * random.uniform(0,1) * (current_pos[i][rand_param] - current_pos[neighbour_index][rand_param])

            new_nectar = measure_nectar(current_pos[i]) 
                
            old_nectar = curr_nectar[i] 

            x = dominance_check(new_nectar, old_nectar)
               
            if x == new_nectar:
                
                curr_nectar[i] = new_nectar 
                
                update_archive(curr_nectar)       
  
            else: 
                #increase abandoment count; if exceeds max_trials send scouts  
                self.population.abandonment_count[i] += 1
                
                if self.population.abandonment_count[i] > self.max_trial:
                    send_scouts(i)  
            
            #finally, update bee's position in population 
            self.population.location[i] = current_pos[i]


    def update_archive(self, initial_config = False, init_nectar):   
        
        if initial_config = True:
            np.savetxt("./external_archive", init_nectar)    
        

        else: 
            #will add the comparison functions later#
            np.savetxt("./external_archive", init_nectar)
            pass 


        #next, add function to check crowding/diversity in archive and remove points according to epsilon method 





    def optimise(self):
        
        self.population = []                                                                          #initialise the population of employed bees and save Bee objects in list 
        for i in range(self.pop_size):

            self.population.append(Bee(self.param_count, self.param_min, self.param_max))

        if self.add_current_to_individuals:                                                           #can have current machine status as one of the bee's initial location
            current_ap = self.interactor.get_ap()
            self.population[0].choose_position(current_ap)
        

        init_nectar = measure_nectar(self.population.location)                                        #list of initial nectar amounts - used for comparison with local improvements
        initial_config = True                                                                         #initially populates archive with initial positions 
        update_archive(initial_config, init_nectar) 
        for i in range(1, self.max_iters):                                                            #begin main loop 
            #these are the functions which need to be called to run main body of algorithm. Will write each one individually and gradually add them when ready     
            send_employees(init_nectar) 
            send_onlookers()
            #send scouts removed from here - only send them if abandonment count requires them to be sent#
            #update archive removed pending "investingataion"      
            
        




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
               

    def choose_position(self, x0):
        """
        Function that allows a specific particle in swarm to have a specific location (used for 'use current' option)

        Args:
            coords of new position in parameter space

        Returns:
            None, but updates positions of particle
        """
        self.location = list(x0)





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

        Tkinter.Label(self, text="Number of employed bees:").grid(row=0, column=0, sticky=Tkinter.E)
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












