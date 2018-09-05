"""

Second version of MOABC Algorithm based on paper by Reza Akbari et al.

Created by Rafael Hunt-Stokes @06/08/18

"""
from __future__ import division
import random
import numpy as np
import Tkinter
import ttk
import os
from scipy import spatial
from dlsoo import plot
import copy
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import datetime 

store_address = None                         #directory in which output data will be stored
completed_iteration = 0                      #number of completed iterations
completed_percentage = 0.0                   #fraction of optimisation completed
pareto_front = ()                            #current pareto-front with the format (((param1,param2,...),(obj1,obj2,...),(err1,err2,...), (std1, std2, ...)),...)



def nothing_function(data):
    pass

class Optimiser(object):

    def __init__(self, settings_dict, interactor, store_location, a_min_var, a_max_var, progress_handler = None):

        self.interactor = interactor                                                                  #interactor with dls_optimizer util.py
        self.store_location = store_location                                                          #location of output file
        self.pop_size = settings_dict['pop_size']                                                     #number of employed bees (N employed, N onlooker)
        self.max_iter = settings_dict['max_iter']                                                     #number of iterations
        self.param_min = a_min_var                                                                    #minimum bounds of parameters
        self.param_max = a_max_var                                                                    #maximum bounds of parameters
        self.max_trial = settings_dict['max_trials']                                                  #maximum number of iterations food source quality not improved before abandonment
        self.param_count = len(interactor.param_var_groups)                                           #number of parameters
        self.obj_count = len(interactor.measurement_vars)                                             #number of objectives
        self.employed_weight = settings_dict['w1']                                                    #importance of randomly selected other food source in calculation of updated positions
        self.onlooker_weight = settings_dict['w2']                                                    #importance of information given by employed bees to onlookers when onlookers attempt food source improvements
        self.scouting_activated_number = 0
        self.archive_size = settings_dict['archive_size']

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



    def save_details_file(self):
        """
        Function writes a file containing details of algorithm run
        """
        file_return = ""

        file_return += "dlsoo_moabc.py algorithm\n"
        file_return += "=================\n\n"
        file_return += "Iterations: {0}\n".format(self.max_iter)
        file_return += "Pop size: {0}\n\n".format(self.pop_size)
        file_return += "Parameter count: {0}\n".format(self.param_count)
        file_return += "Objectives count: {0}\n\n".format(self.obj_count)
        file_return += "Minimum bounds: {0}\n".format(self.param_min)
        file_return += "Maximum bounds: {0}\n\n".format(self.param_max)
        file_return += "Employee Weight: {0}\n".format(self.employed_weight)
        file_return += "Onlooker weight: {0}\n".format(self.onlooker_weight)


        return file_return

    def logbook(self, population, nectar_list, objectives, result, deleted_index, front_pop):
        """
        Function used to record details of algorithm processes to aid with bug fixing.
        Args: 
                Completed Iteration 
                Population Passed
                Nectar List Used (in create_pareto func.) 
                Objectives Considered in pareto_remover
                Result of pareto_remover check 
                Indices to delete 
                Pop. passed to front

        Outputs: File in top directory containing details of each iteration run.
        """
       
        f = file("{0}/algo_run_data.{1}".format(self.store_location, completed_iteration), "w")
        f.write("\n\n***ITERATION: {}***\n\n".format(completed_iteration))
        f.write("Population to create_Pareto:\n\n{0}\n\nNectar list used:\n\n {1}\n\n***Results of Dominance Checks***".format(population, nectar_list)) 
        
        for i in range(len(objectives)): 
            f.write("Check between: {0}\n\n".format(objectives[i]))
            if result[i] == False:
                result[i] = 0 
            f.write("Objective removed: {}\n\n".format(result[i]))
              
        f.write("Indices to delete: {0}\n\nPop. passed to front: {1}".format(deleted_index, front_pop)) 
        f.close

    def measure_nectar(self, bee_location):
        """Function that performs objective measurement of Bees"""



        self.interactor.set_ap(bee_location)                                     #configure machine for measurement
        data = self.interactor.get_ar()                                          #perform measurement
        results = [i.mean for i in data]                                         #retrieve mean from measurement
        std = [i.dev for i in data]                                              #retrieve standard deviation
        err = [i.err for i in data]                                              #retrieve standard error of mean



        return results, std, err

    def dump_fronts(self, fronts, iteration):
        """
        Saves front information to file after each iteration. Used by plotting functions.
        """
        f = file("{0}/FRONTS/fronts.{1}".format(self.store_location, iteration), "w")
        f.write("fronts = ((\n")
        for i, data in enumerate(fronts):
            f.write("    ({0}, {1}, {2}, {3}), \n".format(data[0], tuple(data[1]), data[2], data[3]))
        f.write("),)\n")
        f.close

    def pareto_remover(self,a,b):
        """
        Function determines which of two points is the dominant in objective space.
        Used in create_pareto() function.  
        Args:
            a: list of objective values [obj1,obj2,...]
            b: list of objective values [obj1,obj2,...]

        Returns:
            Function will return the point that dominates the other. If neither dominates, return is False.
        """
        if all(a_i > b_i for (a_i,b_i) in zip(a,b)):         #does b dominate a?
            return a
        if all(a_i < b_i for (a_i,b_i) in zip(a,b)):         #does a dominate b?
            return b
        if all(a_i == b_i for (a_i,b_i) in zip(a,b)):        #are the points the same?
            return b
        else:
            return False                                    #points are non-dominant and should both remain in front

    def dominance_check(self,a,b):
        """
        Function determines which of two points is the dominant in objective space.
        Used in send_onlookers() selection methods. 
        Args:
            a: list of objective values [obj1,obj2,...]
            b: list of objective values [obj1,obj2,...]

        Returns:
            Function will return the point that dominates the other. If neither dominates, return is False.
        """
        if all(a_i > b_i for (a_i,b_i) in zip(a,b)):         #does b dominate a?
            return b
        if all(a_i < b_i for (a_i,b_i) in zip(a,b)):         #does a dominate b?
            return a

        else:
            return False                                     # a & b non-dominant

    def num_dominated(self, front_bee, curr_pop):
        """
        Function to iterate over dominance check function for each point and return number of points indexed position dominates.
        Used in onlooker methods to calculate roulette wheel selection probabilities.
        """

        numb_dominated = 0

        point_selected = front_bee[1]
        

        for i in range(len(curr_pop)):
            test_against = curr_pop[i].nectar

            if test_against == point_selected:

                continue
            else:
                dominator = self.dominance_check(point_selected, test_against)

            if dominator == point_selected:
                numb_dominated += 1

        return numb_dominated

    def create_pareto(self, pareto_data):
        """
        Updates global Pareto front with latest iterations information. Called by optimise() in main loop.
        Args: pareto_data - population data + previous front 
        Output: none, but updates pareto front     
        """
        global pareto_front
        
        #deepcopy pareto data - subsequent operations do not alter original data until the end (when it's needed)
        copy_pareto_data = copy.deepcopy(pareto_data)
        
        #remove possible duplicate points 
        no_duples = [] 
        for i in copy_pareto_data:
            if i not in no_duples:
                no_duples.append(i) 

        #for use in logbook() function 
        log_copy = copy.deepcopy(pareto_data)
        
        #extract objectives only from the pareto data
        nectar_list = list([Bee[1] for Bee in no_duples])
        
        #for use in logbook function 
        log_nectar = copy.deepcopy(nectar_list) 
        log_objectives = []
        log_results = []

        #create list to store index of dominated solutions for removal 
        delete_index = []
        
        #main loop of create_pareto() 
        for i in range(len(no_duples)):                                                          #cycle through swarm and compare objectives
            for j in range(len(no_duples)):

                if i==j:                                                                         #no need to compare solution with itself
                    continue
                log_objectives.append([log_nectar[i], log_nectar[j]])  
                particle_to_remove = self.pareto_remover(nectar_list[i], nectar_list[j])         #determine which solution is dominated
                log_results.append(particle_to_remove) 
                if particle_to_remove == False:                                                  #if neither are dominant, leave both in front
                    continue
                else:

                    delete_index.append(nectar_list.index(particle_to_remove))                   #store index of solution if it is dominated 

        delete_index = sorted(set(delete_index))
        log_index = delete_index[:] 

        print "indices to be deleted: {}".format(delete_index)

        #cycle through pareto data in no_duples, delete the dominated solutions
        #subtract 1 from index to delete after each loop to account for shrinking pareto data     
        for i in delete_index:
                                                                        
            del no_duples[i]
            for j in range(len(delete_index)):
                delete_index[j] = delete_index[j] - 1 
                

        #redefine global pareto front for latest iteration         
        pareto_front = list(no_duples)

        #call logbook() function to record data of function run in .txt file "algo_run_details.{iteration}"
        self.logbook(log_copy, log_nectar, log_objectives, log_results, log_index, no_duples)   



    def calc_fitness(self, num_front_bee, front, curr_pop):
        """
        Calculates fitness values of each point in archive and returns list. Fitness is number of sources dominated
        by source i / number of sources.
        """
        tot_fit = []
        number_dominated = []
        for i in range(num_front_bee):
            number_dominated.append(self.num_dominated(front[i], curr_pop))

            fit = number_dominated[i] / num_front_bee

            tot_fit.append(fit)


        print "Fitness of each bee: {0}\nNumber bees: {1}\nNumber dominated: {2}".format(tot_fit, num_front_bee, number_dominated)
        return tot_fit


    def calc_probability(self, index, fitness_calc):
        """
        calculates pk fitness values using roulette wheel selection method. Used within onlooker bee methods
        """

        if sum(fitness_calc) == 0:
            pk = 0
        else:
            pk = fitness_calc[index] / sum(fitness_calc)


        return pk


    def send_scouts(self, index, population):
        """
        Called by send_employees() & send_onlookers if abandonment count exceeds max_trials.
        Initialises new random solution within search space.

        """
        self.scouting_activated_number += 1
        
        #randomly produce new food source
        population[index] = Bee(self.param_count, self.param_min, self.param_max)
        #measure new bee location's nectar and update other attributes 
        new_nect, std, err = self.measure_nectar(population[index].location)
        population[index].nectar = new_nect
        population[index].std = std
        population[index].err = err
        print "new BEE: {}, {}, {}, {}".format(population[index].location, population[index].nectar, population[index].std, population[index].err)
        
        return population[index]



    def send_employees(self, population):
        """
        This function performs takes current population and performs local searches with random neighbour location, selected from current front, 
        and performs local searches. A single random parameter is changed and nectar measured. If new nectar is better than old, update position.
        """
        
        global pareto_front 
        global completed_iteration
        #deepcopy to avoid changing population/front before the end of function      
        copy_front = copy.deepcopy(pareto_front) 
        copy_population = copy.deepcopy(population)
		
        #arbritary check of front crowding - this calls the archiver code  
        
        while len(copy_front) > self.archive_size:
            copy_front = self.nearest_neighbours_calc(copy_front)
 
        current_pos = []

        #extract the current parameters of each bee in population and assign to list 
        for i in range(len(copy_population)):
            current_pos.append(copy_population[i].location)

	
        #select random parameter to change 	
        for i in range(len(copy_population)):
            rand_param = random.randint(0, self.param_count-1)

            #select a neighbouring bee location randomly from current front 
            neighbour_index = random.randint(0, len(copy_front)-1)
            
            #if bee and its neighbour are the same, select a new neighbour 
            while copy_population[i] == copy_front[neighbour_index]:
                neighbour_index = random.randint(0, len(copy_front)-1)
  
            #change the randomly selected parameter according to employee equation (see report)
            current_pos[i][rand_param] = current_pos[i][rand_param] + self.employed_weight * random.uniform(0,1) * (current_pos[i][rand_param] - copy_front[neighbour_index][0][rand_param]) 
            
            #ensure bee's do not exceed parameter bounds 
            for j in range(len(population[i].bounds)):
                if current_pos[i][rand_param] < population[i].bounds[0][j]:
                    current_pos[i][rand_param] = population[i].bounds[0][j] 

                if current_pos[i][rand_param] > population[i].bounds[1][j]:
                    current_pos[i][rand_param] = population[i].bounds[1][j]   
	   	
            #next, check whether nectar amount has been improved. If not, increase abandonment count
            new_nectar, std, err = self.measure_nectar(current_pos[i])
	 	     
            #assign previous nectar value to variable for comparison  
            old_nectar = population[i].nectar
	 	
		    #see whether old or new location dominates. x variable is dominating bee  
            x = self.dominance_check(new_nectar, old_nectar)


            #increase abandoment count; if exceeds max_trials send scouts
            if x == old_nectar:
                #increase abandoment count; if exceeds max_trials send scouts
                population[i].abandonment_count += 1

                if population[i].abandonment_count > self.max_trial:
                    population[i] = self.send_scouts(i, population)
                    pass
            else:
                #update bee's position, nectar, std and err in population only if new position is better
                population[i].std = std
                population[i].err = err
                population[i].nectar = new_nectar
                population[i].location = current_pos[i]
        	
    def roulette_selection(self, probs):
        """
        Called by onlookers function to select neighbour in local search eqn.
        """

        max_prob = sum(probs)
        pick = random.uniform(0, max_prob)
        current = 0
        for i in range(len(probs)):
            current += probs[i]
            if current >= pick:



                return i

    def normalised_front(self, front):
        """
        For a given front in objective space, this function will normalise the front to a unit square

        Args:
            front: a list of objectives for all solutions in front

        Returns:
            front_norm: a list of normalised (0.0-->1.0) objective coords
        """

        front_x = [i[0] for i in front]
        front_y = [i[1] for i in front]

        max_x = max(front_x)
        max_y = max(front_y)

        min_x = min(front_x)
        min_y = min(front_y)

        x_norm = [(i-min_x)/(max_x-min_x) for i in front_x]
        y_norm = [(i-min_y)/(max_y-min_y) for i in front_y]

        front_norm = zip(x_norm,y_norm)
        return front_norm


    def nearest_neighbours_calc(self, front):
        """
        Calculates the crowding of solutions in front for onlookers selection function.
        """

	    #don't want to normalise the actual front
	    normalised = copy.deepcopy(front)
        #extract objectives from front 
        front_obj = [i[1] for i in normalised] 
       
        
        #normalise the front objectives to a unit square in obj. space (??) 
        norm_obj = self.normalised_front(front_obj) 
                 
        #pass objective points to KDTree          
        kd_tree = spatial.KDTree(norm_obj) 

        #calculate the number of nearest neighbours of each point, create list of them 
        density = [len(kd_tree.query_ball_point(x = i, r = 0.05)) for i in norm_obj]  
        
        #next, find most crowded, delete that point, return "trimmed" front to onlookers  
        max_crowd_idx = density.index(max(density))
        del front[max_crowd_idx]
        print "trimmed front: {}".format(front_obj) 

        return front  

    def send_onlookers(self, population):
        """
            Function essentially does the same as employed bees, but neighbours are selected from current front 
            via roulette wheel selection method, based on fitness of each solution. 
            Fitness of food source is calculated using:

            pk = fit(xk) / SUM{fit(xm)}

            where fit(xm) = dom(m) / pop_size & dom(m) = function returning number of sources dominated
            by source m.
        """
        global pareto_front 
        
        copy_population = copy.deepcopy(population) 
        copy_front = copy.deepcopy(pareto_front)

        #arbritary check of front crowding - maybe mess with the size of this 
      
     	while len(copy_front) > self.archive_size: 
            copy_front = self.nearest_neighbours_calc(copy_front) 
                            

        current_pos = []

        #assign parameters of each bee in front to a list 
        for i in range(len(copy_population)):
            current_pos.append(copy_population[i].location
        
        #next, call roulette seelction functions to calculate fitness of each solution in front 
        probs_points = []
        fitness_calc = self.calc_fitness(len(copy_front), copy_front, copy_population)
        for i in range(len(copy_front)):

            pk = self.calc_probability(i, fitness_calc)
            probs_points.append(pk)
        print "props points: {}".format(probs_points)
        
        #now we have all the probabilities in a list, we can pass this to roulette selection function
        for i in range(len(copy_population)):
            neighbour_index = self.roulette_selection(probs_points)
            rand_param = random.randint(0, self.param_count - 1)
            curr_nectar = population[i].nectar
            
            #incase selected neighbour is the point itself, redo selection until it isn't    
              
            while copy_front[neighbour_index] == copy_population[i]:
                neighbour_index = self.roulette_selection(probs_points) 
                    
            #change single parameter according to onlooker equation (see report)
            current_pos[i][rand_param] = current_pos[i][rand_param] + self.onlooker_weight * random.uniform(0,1) * (current_pos[i][rand_param] - copy_front[neighbour_index][0][rand_param])
            
            #ensure bee's do not exceed parameter bounds 
            for j in range(len(population[i].bounds)):
                if current_pos[i][rand_param] < population[i].bounds[0][j]:
                    current_pos[i][rand_param] = population[i].bounds[0][j] 

                if current_pos[i][rand_param] > population[i].bounds[1][j]:
                    current_pos[i][rand_param] = population[i].bounds[1][j]
                    
            #measure new nectar amount         
            new_nectar, std, err = self.measure_nectar(current_pos[i])
	    
            #assign old nectar value to variable for comparison 
            old_nectar = curr_nectar
            
            #x is the dominating solution 
            x = self.dominance_check(new_nectar, old_nectar)

            
            #If not improved, increase abandoment count; if exceeds max_trials send scouts
            if x == old_nectar:
                population[i].abandonment_count += 1

                if population[i].abandonment_count > self.max_trial:
                    population[i] = self.send_scouts(i, population)

            else:
                #update bee's position, nectar, std and err in population, only if new position is better

                population[i].nectar = new_nectar
                population[i].location = current_pos[i]
                population[i].std = std
                population[i].err = err
        
    
    def optimise(self):
        """
        Main function in algorithm. This calls all other functions and runs the main loop. 
        """
        
        global store_address
        store_address = self.store_location
        global pareto_front
        global completed_percentage
        global completed_iteration
        population = []                                                                         


        percentage_interval = 1/self.max_iter
        
        #initialise the population of employed bees and save Bee objects in list
        for i in range(self.pop_size):

            population.append(Bee(self.param_count, self.param_min, self.param_max))



        #can have current machine status as one of the bee's initial location
        if self.add_current_to_individuals:                                                           
            current_ap = self.interactor.get_ap()
            population[0].choose_position(current_ap)


        #measure initial nectar amounts and assign each bee attributes 
        for i in range(len(population)):
            
            init_nectar, std, err = self.measure_nectar(population[i].location)                  
            #populate bee parameters for creation of initial pareto front
            population[i].nectar = init_nectar                                                   
            population[i].std = std
            population[i].err = err
        
        #extract information from population to create first front 
        pareto_data = [[j.location,j.nectar,j.err, j.std] for j in population]

        #create first front 
        self.create_pareto(pareto_data)
        
        front_to_dump = tuple(list(pareto_front))
        #dump front in fronts.0 file 
        self.dump_fronts(front_to_dump, 0)
        completed_iteration = 1                                                                              

        #begin main loop
        for i in range(1, self.max_iter):                                                           
            #call employee function 
            print "EMPLOYEE PHASE"
            self.send_employees(population)
            #call onlooker function 
            print "ONLOOKER PHASE"
            self.send_onlookers(population)
            
            #create pareto front after employees and onlooker phases 
            print "CREATING PARETO FRONT"

            #various code to govern pause, cancel and update progress bar 
            completed_percentage += percentage_interval
            print completed_percentage*100,'%'
            self.progress_handler(completed_percentage, completed_iteration)
            while self.pause:                                                         
                self.progress_handler(completed_percentage, completed_iteration)

            if self.cancel:
                break

            #extract information from population and previous front for processing 
            pareto_data = [[j.location, j.nectar, j.err, j.std] for j in population] + pareto_front
            
            #pass data to create_pareto to make new front 
            self.create_pareto(pareto_data)
            
            #dump front in front file 
            front_to_dump = tuple(list(pareto_front))
            self.dump_fronts(front_to_dump, completed_iteration)
            completed_iteration += 1

        #after optimisation is completed, print out the abandonment count of bees     
        abandonment_counters = []
        for i in range(len(population)):
            abandonment_counters.append(population[i].abandonment_count)
        print "OPTIMISATION COMPLETED\nNumber Scouts sent: {0}\nAbandonment_counts: {1}".format(self.scouting_activated_number, abandonment_counters)



class Bee(object):
    """
    Bee objects are individuals in population. They are initialised with attributes nectar, location, bounds, fitness and errors. 
    """
    
    def __init__(self, param_count, param_min, param_max):

        self.fitness = None
        self.location = []
        self.initial_locale(param_count, param_min, param_max)
        self.abandonment_count = 0
        self.nectar = []
        self.std = []
        self.err = []
        self.bounds = [param_min, param_max]
    
    def initial_locale(self, param_count, param_min, param_max):
        """
        Called by init. Will create initial parameter values according to initialisation equation (see report).
        """

        for i in range(param_count):
            x = param_min[i] + random.uniform(0,1) * (param_max[i] - param_min[i])
            self.location.append(x)


    def choose_position(self, x0):
        """
        Function that allows a specific bee in population to have a specific location (used for 'use current' option)

        Args:
            coords of new position in parameter space

        Returns:
            None, but updates positions of bee
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

        Tkinter.Label(self, text="Archive Size:").grid(row=5, column=0, sticky=Tkinter.E)
        self.i5 = Tkinter.Entry(self)
        self.i5.grid(row=5, column =1, sticky = Tkinter.E+Tkinter.W) 

        self.c0 = Tkinter.Checkbutton(self, text="Use current machine state", variable=self.add_current_to_individuals)
        self.c0.grid(row=6, column=1)

        Tkinter.Label(self, text="Recommended:\nNumber employed: 50\nMax. Iterations: 5\nEmployee Coefficient: 0.7\nOnlooker coefficient: 0.8\nAbandonment Count: 5\nArchive Size: pop / 3", justify=Tkinter.LEFT).grid(row=6, column=0, columnspan=2, sticky=Tkinter.W)

        self.i0.insert(0, "10")     #defaults are added in case user does not want to decide
        self.i1.insert(0, "3")
        self.i2.insert(0, "5")
        self.i3.insert(0, "0.7")
        self.i4.insert(0, "0.8")
     


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
        try:
            setup['archive_size'] = int(self.i5.get())
        except:
            raise ValueError("The value for \"archive_size\": \"{0}\", could not be converted to an integer".format(self.i5.get())) 
        if self.add_current_to_individuals.get() == 0:
            setup['add_current_to_individuals'] = False
        elif self.add_current_to_individuals.get() == 1:
            setup['add_current_to_individuals'] = True

        return setup


class import_algo_prog_plot(Tkinter.Frame):
    """
    This class sets up the MOPSO progress plot that shows during the optimisation, including a percentage bar and a plot of the current fronts.
    """

    def __init__(self, parent, axis_labels, signConverter):

        Tkinter.Frame.__init__(self, parent)

        self.parent = parent
        self.signConverter = signConverter
        self.axis_labels = axis_labels

        self.initUi()
        
    def initUi(self):
        """
        setup plots
        """
        
        self.fig = Figure(figsize=(5, 5), dpi=100)
        self.a = self.fig.add_subplot(111)

        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=Tkinter.BOTTOM, fill=Tkinter.BOTH, expand=True)

    def update(self):
        """
        after each iteration, plot the new front
        """
        global store_address
        global completed_iteration
        
        self.a.clear()
        file_names = []
        print "Store address: {0}\nCompleted iteration: {1}".format(store_address, completed_iteration)
        for i in range(completed_iteration):
            file_names.append("{0}/FRONTS/fronts.{1}".format(store_address, i))
        
        plot.plot_pareto_fronts(file_names, self.a, self.axis_labels, self.signConverter)

        self.canvas.show()

#--------------------------------------------------------------- CLASS FOR FINAL RESULTS WINDOW --------------------------------------------------------#

class import_algo_final_plot(Tkinter.Frame):
    """
    This class sets up the final results window that shows after the optimisation is complete. The actual plot of the final fronts is
    imported from the final_plot class.
    """

    def __init__(self, parent, pick_handler, axis_labels, signConverter, initial_config=None, post_analysis_store_address = None):

        global store_address
        Tkinter.Frame.__init__(self, parent)

        self.parent = parent
        self.signConverter = signConverter

        self.pick_handler = pick_handler
        self.axis_labels = axis_labels

        if initial_config is not None:
            self.initial_measurements = initial_config

        if post_analysis_store_address is not None:             #this class is also used for post_analysis
            store_address = post_analysis_store_address


    def initUi(self, initial_config_plot=True):
        """
        Setup window GUI
        """
        global store_address

        self.parent.title("MOABC results")

        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=0)

        self.rowconfigure(0, weight=1)

        self.view_mode = Tkinter.StringVar()
        self.view_mode.set('No focus')

        if initial_config_plot is True:
            self.plot_frame = final_plot(self, self.axis_labels, self.signConverter, initial_config=self.initial_measurements)
        else:
            self.plot_frame = final_plot(self, self.axis_labels, self.signConverter)
        self.plot_frame.grid(row=0, column=0, pady=20, padx=20, rowspan=1, sticky=Tkinter.N+Tkinter.S+Tkinter.E+Tkinter.W)

        Tkinter.Label(self, text="View mode:").grid(row=0, column=1)

        self.cbx_view_mode = ttk.Combobox(self, textvariable=self.view_mode, values=('No focus', 'Best focus'))
        self.cbx_view_mode.bind("<<ComboboxSelected>>", lambda x: self.plot_frame.initUi())
        self.cbx_view_mode.grid(row=0, column=2)

        self.grid(sticky=Tkinter.N+Tkinter.S+Tkinter.E+Tkinter.W)
        self.parent.columnconfigure(0, weight=1)
        self.parent.rowconfigure(0, weight=1)

    def on_pick(self, event):
        """
        This function gathers information from the saved files to allow the user to see the machine/algorithm parameters/results upon clicking
        on a solution on the Pareo front.
        """
        global store_address
        completed_iteration = len(os.listdir('{0}/FRONTS'.format(store_address)))

        my_artist = event.artist
        x_data = my_artist.get_xdata()
        y_data = my_artist.get_ydata()
        ind = event.ind
        point = tuple(zip(self.signConverter[0]*x_data[ind], self.signConverter[1]*y_data[ind]))

        print "Point selected, point: {0}".format(point)

        ''' By this point we have the ars, but not the aps. We get these next. '''

        file_names = []
        for i in range(completed_iteration):
            file_names.append("{0}/FRONTS/fronts.{1}".format(store_address, i))


        fs = []

        for file_name in file_names:
            execfile(file_name)

            fs.append(locals()['fronts'][0])

        aggregate_front_data = []
        for i in fs:
            for j in i:
                aggregate_front_data.append(j)
        aggregate_front_results = [i[1] for i in aggregate_front_data]
        point_number = aggregate_front_results.index(point[0])
        point_a_params = tuple(aggregate_front_data[point_number][0])

        print "ap: {0}".format(point_a_params)

        ''' By this point he have the aps, but not the mps. We don't find these in the algorithm. '''

        print "pick_handler passed this: ", point[0], point_a_params
        self.pick_handler(point[0], point_a_params)

#--------------------------------------------------------- CLASS FOR FINAL PLOT IN RESULTS WINDOW -------------------------------------------------#

class final_plot(Tkinter.Frame):
    """
    This class is called upon in the import_algo_final_plot class. It retrives the dumped fronts and then used dls_optimiser_plot to plot the
    Pareto fronts.
    """

    def __init__(self, parent, axis_labels, signConverter, initial_config=None):

        Tkinter.Frame.__init__(self, parent)

        self.parent = parent
        self.signConverter = signConverter
        self.axis_labels = axis_labels

        if initial_config is not None:
            self.initial_measurements = initial_config
            self.initUi(initial_config_plot=True)
        else:
            self.initUi()

    def initUi(self, initial_config_plot=False):
        global store_address
        completed_iteration = len(os.listdir('{0}/FRONTS'.format(store_address)))

        for widget in self.winfo_children():
            widget.destroy()

        fig = Figure(figsize=(5, 5), dpi=100)
        a = fig.add_subplot(111)
        fig.subplots_adjust(left=0.15)

        file_names = []
        for i in range(completed_iteration):                                            #gather fronts
            file_names.append("{0}/FRONTS/fronts.{1}".format(store_address, i))

        print 'file names', file_names

        if initial_config_plot is True:

            plot.plot_pareto_fronts_interactive(file_names,
                                                a,
                                                self.axis_labels,
                                                None,
                                                None,
                                                self.parent.view_mode.get(),
                                                self.signConverter,
                                                initial_measurements=self.initial_measurements)
        else:

            plot.plot_pareto_fronts_interactive(file_names,
                                                a,
                                                self.axis_labels,
                                                None,
                                                None,
                                                self.parent.view_mode.get(),
                                                self.signConverter)

        canvas = FigureCanvasTkAgg(fig, self)
        canvas.mpl_connect('pick_event', self.parent.on_pick)
        canvas.show()
        canvas.get_tk_widget().pack(side=Tkinter.BOTTOM, fill=Tkinter.BOTH, expand=True)

        toolbar = NavigationToolbar2TkAgg(canvas, self)
        toolbar.update()
        canvas._tkcanvas.pack(side=Tkinter.TOP, fill=Tkinter.BOTH, expand=True)









