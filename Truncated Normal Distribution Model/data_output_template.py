"""A simulation of coevolution of prey and predator
( birth、death、competition and mutation) using Weini's Gillespie Algorithm* """
# * See https://github.com/WeiniHuangBW/DynamicalTradeoffandCoevolution

"""Truncated Normal Distribution Model"""

# This version is for output the data I need for my research.
# 1000 runs and time period equals to 2000 per run.
# From line 115 can see the parameters I set.

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import truncnorm

class Prey():
    """ 
    A class to represent a type of prey with follwing attribute:
    ❶ Number of this species of prey: Nx
    ❷ Growth-defense trade off parameter: gi ∈ (0,1) 
      (truncated normal distribution model)
    """
    def __init__(self, Nx, gi):
        """initialize attribute to descibe a kind of prey"""
        self.gi = gi
        self.Nx = Nx
        self.population_size = []
        self.time = []

    def prey_mutation(self):
        """mutation of prey to generate a new species (gi following an truncated normal distribution)."""
        new_prey = Prey(Nx=1, gi=get_truncated_normal_random_number(self.gi, gi_sd, gi_low, gi_up, number=1))
        current_prey_type.append(new_prey)

class Predator():
    """
    A class to represent a type of predator with following attribute:
    ❶ Number of this species of predator: Ny
    ❷ Reproduction efficacy(defined as the ratio of the reproduction rate to 
    the predation rate of predators): kl ∈ (0, kmax], kmax=0.3
     (truncated normal distribution model)
    """
    def __init__(self, Ny, kl):
        """initialize attribute to descibe a kind of predator"""
        self.kl = kl
        self.Ny = Ny
        self.population_size = []
        self.time = []

    def predator_mutation(self):
        """mutation of predator to generate a new species (kl following an truncated normal distribution)."""
        new_predator = Predator(Ny=1, kl=get_truncated_normal_random_number(self.kl, kl_sd, kl_low, kl_up, number=1))
        current_predator_type.append(new_predator)

def get_truncated_normal_random_number(mean, sd, low, up, number=1):
    """to get truncated normal random number"""
    truncated_normal_generator = truncnorm(
        (low - mean)/sd, (up -mean)/sd, loc=mean, scale=sd)
    truncated_normal_random_number = truncated_normal_generator.rvs(number)
    return truncated_normal_random_number

def predation_rate(gi, kl):
    """Calculate the predation rate."""
    f = p * (gi ** (m*kl/kmax))
    return f

def reaction_time(reaction_rate):
    """Calculate the reaction time with random number."""
    reaction_time = - (1/reaction_rate) * np.log(np.random.random_sample())
    return reaction_time

def predation_without_reproduction(current_prey_type, current_predator_type): 
    """Calculate the sum of all predator reaction coefficient for each prey and store it into a list. """
    ls_0 = []   
    for i in current_prey_type:
        a = sum([i.Nx * l.Ny * (1-l.kl) * predation_rate(i.gi, l.kl) for l in current_predator_type])
        ls_0.append(a)

    return ls_0

def prey_extinction_check(current_prey_type, s=0):
    """If one kind of prey species extinct, remove it from current_type to extinct_type."""
    while s < len(current_prey_type):
        if current_prey_type[s].Nx == 0:
            del current_prey_type[s]

        else:
            s += 1

def predator_extinction_check(current_predator_type, s=0):
    """If one kind of predator species extinct, remove it from current_type to extinct_type."""
    while s < len(current_predator_type):
        if current_predator_type[s].Ny == 0:
            del current_predator_type[s]

        else:
            s += 1
 
def calculate_simpson_index(species_list, simpson_index_list):
    """Calculate the Simpson's Diversity Index: 1 - ∑pi**2,  pi is the proportion of characters belonging to the ith type of letter in the string of interest. And then append the index into simpson_index_list."""
    
    species_list_array = np.array(species_list) # transform it into array
    ratio = species_list_array / species_list_array.sum()
    simpson_index = 1 - sum(ratio**2)

    simpson_index_list.append(simpson_index)

def calculate_shannon_index(species_list_, shannon_index_list):
    """Calculate the Shannon index: -∑pilnpi,  pi is the proportion of characters belonging to the ith type of letter in the string of interest. And then append the index into shannon_index_list."""
    species_list_array_ = np.array(species_list_)
    ratio_ = species_list_array_ / species_list_array_.sum()

    shannon_index = - sum(ratio_ * np.log(ratio_))

    shannon_index_list.append(shannon_index)


# 1.Input Fixed Parameters
## Prey
bx = 1.0  # float; baseline growth rate of prey
dx = 0.1  # float; intrinsic death rate of prey
rc = 0.00005  # float; resource competition coefficient
ux = 0.0001  # mutation rate of prey per division

## Predation
p = 0.005  # float; scaling coefficient of the predation rate
dy = 0.5  # float; intrinsic death rate of predator
m = 0.5  # int; determines the initial shape of the growth-defense trade-off curve for the prey species
uy = 0.001  # mutation rate of predation per division
kmax = 0.3  # the maximum reproduction efficacy

## Time
t = 0.0  # float; start time
T = 2000.0 # float; maximum elapsed time

## For the mutation of gi and kl: Truncated Normal Distribution Parameters ☆☆☆
gi_sd = 0.01  # float; standard deviation of gi
gi_up = 1.0  # float; maximum value of gi
gi_low = 0.0  # float; minimum value of gi

kl_sd = 0.003  # float; standard deviation of gi
kl_up = kmax  # float; maximum value of kl
kl_low = 0.0  # float; minimum value of kl

# 2.List that collect all the data
outputlist = []
## NOTE:
## 0 stands for only predator go extinct.
## 1 stands for both prey and predator go extinct.
## 2 stands for coexsit.

check_list = [] # Using for check the maximum prey number to ensure it will not go extict.

co_exist_prey = []
co_exist_predator = []

prey_simpson_index_list = []
predator_simpson_index_list = []

prey_shannon_index_list =[]
predator_shannon_index_list = []


# 3.Outer Loop: the Gillespie Algorithem
loop_number = 0
while loop_number < 1000:
    loop_number = loop_number + 1

    # For each run, set the initial conditions
    ancestor_prey = Prey(Nx=1000, gi=1.0)
	ancestor_predator = Predator(Ny=100, kl=0.3)

    current_prey_type = [ancestor_prey]  # record current prey 
    current_predator_type = [ancestor_predator]  # record current predator

    # Inner loop(main): the Gillespie Algorithem
    while t <= T:
        
        # (1) Check the species extinction phenomenon.
        prey_extinction_check(current_prey_type)
        predator_extinction_check(current_predator_type)

        ## break the inner loop when all prey go extinct and make record
        if len(current_prey_type) == 0:
            outputlist.append(1)
            break
        
        ## break the inner loop when predator go extinct and make record
        if len(current_predator_type) == 0: 
            outputlist.append(0)  
            max_prey_number = max([i.Nx for i in current_prey_type])
            check_list.append(max_prey_number) 
            break

        # (2) Main Program
        total_prey_number = sum([i.Nx for i in current_prey_type])
    
        ## ❶ Define the reaction rates.
        ### PREY REACTION
        prey_birth_without_mutation = [bx * (1-ux) * i.gi * i.Nx for i in current_prey_type]
        prey_birth_with_mutation = [bx * ux * i.gi * i.Nx for i in current_prey_type]
        prey_competition_death = [rc * i.Nx * total_prey_number for i in current_prey_type]
        prey_intrinsic_death = [dx * i.Nx for i in current_prey_type]
        
        ### PREDATOR REACTION
        predation_without_birth = predation_without_reproduction(current_prey_type, current_predator_type)
        predation_with_birth_without_mutation = [i.Nx * l.Ny * l.kl * predation_rate(i.gi, l.kl) * (1-uy) for i in current_prey_type for l in current_predator_type] ## ☆ needed to be transformed into 2d_array
        predation_with_birth_with_mutation = [i.Nx * l.Ny * l.kl * predation_rate(i.gi, l.kl) * uy for i in current_prey_type for l in current_predator_type]  ## ☆ needed to be transformed into 2d_array
        predator_intrinsic_death = [dy * l.Ny for l in current_predator_type]

        ## ❷Calculate the reation time with a random number: (1/ai) * log(1/ri).
        ### PREY REACTION
        time_prey_birth_without_mutation = list(map(reaction_time, prey_birth_without_mutation))
        time_prey_birth_with_mutation = list(map(reaction_time, prey_birth_with_mutation))
        time_prey_competition_death = list(map(reaction_time, prey_competition_death))
        time_prey_intrinsic_death = list(map(reaction_time, prey_intrinsic_death))
        
        ### PREDATOR REACTION
        time_predation_without_birth = list(map(reaction_time, predation_without_birth))
        time_predation_with_birth_without_mutation = list(map(reaction_time, predation_with_birth_without_mutation))
        time_predation_with_birth_with_mutation = list(map(reaction_time, predation_with_birth_with_mutation ))
        time_predator_intrinsic_death = list(map(reaction_time, predator_intrinsic_death))

        ## ❸ Pick up the shortest time for the reaction.
        prey_array = np.array([time_prey_birth_without_mutation, time_prey_birth_with_mutation, time_prey_competition_death, time_prey_intrinsic_death, time_predation_without_birth])  ### here, I generate a n-darray for the convenience of later analysis and choice ;)
        prey_array_min = np.min(prey_array)

        prey_predator_without_mutation_array = np.array(time_predation_with_birth_without_mutation).reshape(len(current_prey_type), len(current_predator_type))
        prey_predator_without_mutation_array_min = np.min(prey_predator_without_mutation_array)

        prey_predator_with_mutation_array = np.array(time_predation_with_birth_with_mutation).reshape(len(current_prey_type), len(current_predator_type))
        prey_predator_with_mutation_array_min = np.min(prey_predator_with_mutation_array)

        time_predator_intrinsic_death_min = min(time_predator_intrinsic_death)

        ### ★tau for reaction and renew the time 
        tau = min(prey_array_min, prey_predator_without_mutation_array_min, prey_predator_with_mutation_array_min, time_predator_intrinsic_death_min)

        t = t + tau # ⚠️NOTE HERE! Things will go wrong when write 't += tau', still not knowing why...⚠️
        
        ## ❹ Select which reaction takes place and update the number.
        if tau == prey_array_min:
            prey_index_reation = np.where(prey_array == prey_array_min) 
            # prey_index_reation of the min value: a tuple——((array([number]), array([number])))
            # prey_index_reation[0] represent the reaction to happen; prey_index_reation[1] represent the species in prey_type list.
            if prey_index_reation[0][0] == 0: ## prey reproduction without mutation
                current_prey_type[prey_index_reation[1][0]].Nx += 1
            elif prey_index_reation[0][0] == 1: ## prey reproduction with mutation
                current_prey_type[prey_index_reation[1][0]].prey_mutation()
            else: ## prey die
                current_prey_type[prey_index_reation[1][0]].Nx -= 1
                
        elif tau == prey_predator_without_mutation_array_min:  ## prey be eaten and predator birth without mutation
            prey_predator_without_mutation_index = np.where(prey_predator_without_mutation_array == prey_predator_without_mutation_array_min) 

            prey_index_death_eaten = prey_predator_without_mutation_index[0][0]
            predator_index_birth = prey_predator_without_mutation_index[1][0]

            current_prey_type[prey_index_death_eaten].Nx -= 1
            current_predator_type[predator_index_birth].Ny += 1  

        elif tau == prey_predator_with_mutation_array_min:  ## prey be eaten and predator birth with mutation
            prey_predator_with_mutation_index = np.where(prey_predator_with_mutation_array == prey_predator_with_mutation_array_min)

            prey_index_death_eaten_ = prey_predator_with_mutation_index[0][0]
            predator_index_birth_mutant = prey_predator_with_mutation_index[1][0]

            current_prey_type[prey_index_death_eaten_].Nx -= 1
            current_predator_type[predator_index_birth_mutant].predator_mutation()
            
        else:  ## predator die naturally
            predator_index_death = time_predator_intrinsic_death.index(time_predator_intrinsic_death_min)
            current_predator_type[predator_index_death].Ny -= 1

    # when prey and predator can co-exist:
    if len(current_prey_type) != 0 and len(current_predator_type) != 0:  ## this judgment using for thwart extiction phenomenon 
        co_exist_prey.append(len(current_prey_type))
        co_exist_predator.append(len(current_predator_type))
        outputlist.append(2)

        # calculate the index we need:
        calculate_simpson_index([i.Nx for i in current_prey_type], prey_simpson_index_list)
        calculate_simpson_index([l.Ny for l in current_predator_type], predator_simpson_index_list)

        calculate_shannon_index([i.Nx for i in current_prey_type], prey_shannon_index_list)
        calculate_shannon_index([l.Ny for l in current_predator_type], predator_shannon_index_list)
    

# 4. Data processing
## Coexistence:
### array
co_exist_prey_array = np.array(co_exist_prey)
co_exist_predator_array = np.array(co_exist_predator)
prey_simpson_index_list_array = np.array(prey_simpson_index_list)
predator_simpson_index_list_array = np.array(predator_simpson_index_list)
prey_shannon_index_list_array = np.array(prey_shannon_index_list)
predator_shannon_index_list_array = np.array(predator_shannon_index_list)

mean_prey_types = co_exist_prey_array.mean() 
mean_predator_types = co_exist_predator_array.mean()

prey_simpson_index = prey_simpson_index_list_array.mean()
predator_simpson_index = predator_simpson_index_list_array.mean()
prey_shannon_index = prey_shannon_index_list_array.mean()
predator_shannon_index = predator_shannon_index_list_array.mean()

coexist_number = outputlist.count(2)

# Extinction:
predator_extinction_number = outputlist.count(0)
both_extinctin_number = outputlist.count(1)

# 5. Output data
with open('1_1.txt', 'w+') as fp1:
    fp1.write('m=0.5, sd=1%')
    fp1.write('\n' + 'mean_prey_types= ' + str(mean_prey_types))
    fp1.write('\n' + 'mean_predator_types= ' + str(mean_predator_types))
    fp1.write('\n' + 'prey_shannon_index=  ' + str(prey_shannon_index))
    fp1.write('\n' + 'predator_shannon_index= ' + str(predator_shannon_index))
    fp1.write('\n' + 'prey_simpson_index= ' + str(prey_simpson_index))
    fp1.write('\n' + 'predator_simpson_index= ' + str(predator_simpson_index))
    fp1.write('\n' + 'time that only predator extinct: ' + str(predator_extinction_number))
    fp1.write('\n' + 'time that both extinct: ' + str(both_extinctin_number))
    fp1.write('\n' + 'time that coexist:' + str(coexist_number))

with open('1_1_check.txt', 'w+') as fp:
    fp.write('m=0.5, sd=1%')
    fp.write('\n' + str(check_list))