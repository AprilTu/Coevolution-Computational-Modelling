"""A simulation of coevolution of prey and predator
( birth、death、competition and mutation) using Weini's Gillespie Algorithm* """
# * See https://github.com/WeiniHuangBW/DynamicalTradeoffandCoevolution

"""Truncated Normal Distribution Model"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import truncnorm

class Prey():
    """ 
    A class to represent a type of prey with follwing attribute:
    ❶ Number of this species of prey: Nx
    ❷ Growth-defense trade off parameter: gi ∈ (0,1) 
      (following an truncated normal distribution Model)
    """
    def __init__(self, Nx, gi):
        """initialize attribute to descibe a kind of prey"""
        self.gi = gi
        self.Nx = Nx
        self.population_size = []
        self.time = []

    def record_data(self, t):
        """to record the number of the prey in each reaction"""
        self.population_size.append(self.Nx)
        self.time.append(t)

    def prey_mutation(self):
        """mutation of prey to generate a new species."""
        new_prey = Prey(Nx=1, gi=get_truncated_normal_random_number(self.gi, gi_sd, gi_low, gi_up, number=1))
        current_prey_type.append(new_prey)

class Predator():
    """
    A class to represent a type of predator with following attribute:
    ❶ Number of this species of predator: Ny
    ❷ Reproduction efficacy(defined as the ratio of the reproduction rate to 
    the predation rate of predators): kl ∈ (0, kmax], kmax=0.3
     (following an truncated normal distribution Model)
    """
    def __init__(self, Ny, kl):
        """initialize attribute to descibe a kind of predator"""
        self.kl = kl
        self.Ny = Ny
        self.population_size = []
        self.time = []

    def record_data(self, t):
        """to record the number of predator in each reaction"""
        self.population_size.append(self.Ny)
        self.time.append(t) 

    def predator_mutation(self):
        """mutation of predator to generate a new species."""
        new_predator = Predator(Ny=1, kl=get_truncated_normal_random_number(self.kl, kl_sd, kl_low, kl_up, number=1))
        current_predator_type.append(new_predator)

def get_truncated_normal_random_number(mean, sd, low, up, number=1):
    """to get truncated normal random number"""
    truncated_normal_generator = truncnorm(
        (low - mean)/sd, (up -mean)/sd, loc=mean, scale=sd)
    truncated_normal_random_number = truncated_normal_generator.rvs(number)
    return truncated_normal_random_number

def prey_extinction_check(current_prey_type, extinct_prey_type, s=0):
    """if one kind of prey species extinct, remove it from current_type to extinct_type."""
    while s < len(current_prey_type):
        if current_prey_type[s].Nx == 0:
            extinct_prey_type.append(current_prey_type[s])
            del current_prey_type[s]

        else:
            s += 1

def predator_extinction_check(current_predator_type, extinct_predator_type, s=0):
    """if one kind of predator species extinct, remove it from current_type to extinct_type."""
    while s < len(current_predator_type):
        if current_predator_type[s].Ny == 0:
            extinct_predator_type.append(current_predator_type[s])
            del current_predator_type[s]

        else:
            s += 1

def predation_rate(gi, kl):
    f = p * (gi ** (m*kl/kmax))
    return f

def predation_without_reproduction(current_prey_type, current_predator_type): 
    ls_0 = []   
    for i in current_prey_type:
        a = sum([i.Nx * l.Ny * (1-l.kl) * predation_rate(i.gi, l.kl) for l in current_predator_type])
        ls_0.append(a)

    return ls_0

def reaction_time(reaction_rate):
    """Calculate the reaction time with random number"""
    reaction_time = - (1/reaction_rate) * np.log(np.random.random_sample())
    return reaction_time
   

# 1.Input Fixed Parameters and set them to numpy words
## Prey
bx = 1.0; bx = np.float64(bx)  # float; baseline growth rate of prey
dx = 0.1; dx = np.float64(dx)  # float; intrinsic death rate of prey
rc = 0.00005; rc = np.float64(rc)  # float; resource competition coefficient
ux = 0.0001; ux = np.float64(ux)  # mutation rate of prey per division

## Predation
p = 0.005; p = np.float64(p)  # float; scaling coefficient of the predation rate
dy = 0.5; dy = np.float64(dy)  # float; intrinsic death rate of predator
m = 2; m = np.int64(m)  # int; determines the initial shape of the growth-defense trade-off curve for the prey species
uy = 0.001; uy = np.float64(uy)  # mutation rate of predation per division
kmax = 0.3; kmax = np.float64(kmax)  # the maximum reproduction efficacy

## Time
t = 0.0; t = np.float64(t)  # float; start time
T = 1000.0; T = np.float64(T) # float; maximum elapsed time


## Normal Distribution Parameters ☆☆☆
gi_sd = 1; gi_sd = np.float64(gi_sd)
gi_up = 1.0; gi_up = np.float64(gi_up)
gi_low = 0.0; gi_low = np.float64(gi_low)

kl_sd = 0.3; kl_sd = np.float64(kl_sd)
kl_up = kmax; kl_up = np.float64(kl_up)
kl_low = 0.0; kl_low = np.float64(kl_low)

# 2.Set the initial conditions
ancestor_prey = Prey(Nx=1000, gi=1.0)
ancestor_predator = Predator(Ny=100, kl=0.3)

# 3.List that collect all the data
time_data = []  # record time

current_prey_type = [ancestor_prey]  # record current prey 
current_predator_type = [ancestor_predator]  # record current predator
extinct_prey_type = []  # record the extinct prey
extinct_predator_type = []  # record the extinct predator

# 4.Main Loop: the Gillespie Algorithem
while t <= T:
    # ★First check extinction phenomenon
    prey_extinction_check(current_prey_type, extinct_prey_type)
    predator_extinction_check(current_predator_type, extinct_predator_type)

    total_prey_number = sum([i.Nx for i in current_prey_type])
	
	# ❶Define the reaction rates.
    ## PREY REACTION
    prey_birth_without_mutation = [bx * (1-ux) * i.gi * i.Nx for i in current_prey_type]
    prey_birth_with_mutation = [bx * ux * i.gi * i.Nx for i in current_prey_type]
    prey_competition_death = [rc * i.Nx * total_prey_number for i in current_prey_type]
    prey_intrinsic_death = [dx * i.Nx for i in current_prey_type]
    
    ## PREDATOR REACTION
    predation_without_birth = predation_without_reproduction(current_prey_type, current_predator_type)
    predation_with_birth_without_mutation = [i.Nx * l.Ny * l.kl * predation_rate(i.gi, l.kl) * (1-uy) for i in current_prey_type for l in current_predator_type] ## ☆ needed to be transformed into 2d_array
    predation_with_birth_with_mutation = [i.Nx * l.Ny * l.kl * predation_rate(i.gi, l.kl) * uy for i in current_prey_type for l in current_predator_type]  ## ☆ needed to be transformed into 2d_array
    predator_intrinsic_death = [dy * l.Ny for l in current_predator_type]

    ## ❷calculate the reation time with a random number: (1/ai) * log(1/ri).
    time_prey_birth_without_mutation = list(map(reaction_time, prey_birth_without_mutation))
    time_prey_birth_with_mutation = list(map(reaction_time, prey_birth_with_mutation))
    time_prey_competition_death = list(map(reaction_time, prey_competition_death))
    time_prey_intrinsic_death = list(map(reaction_time, prey_intrinsic_death))

    time_predation_without_birth = list(map(reaction_time, predation_without_birth))
    time_predation_with_birth_without_mutation = list(map(reaction_time, predation_with_birth_without_mutation))
    time_predation_with_birth_with_mutation = list(map(reaction_time, predation_with_birth_with_mutation ))
    time_predator_intrinsic_death = list(map(reaction_time, predator_intrinsic_death))

    ## ❸Pick up the shortest time for the reaction and Select which reaction takes place and update the number.
    prey_array = np.array([time_prey_birth_without_mutation, time_prey_birth_with_mutation, time_prey_competition_death, time_prey_intrinsic_death, time_predation_without_birth])  ### here, I generate a n-darray for the convenience of later analysis ;)
    prey_array_min = np.min(prey_array)

    prey_predator_without_mutation_array = np.array(time_predation_with_birth_without_mutation).reshape(len(current_prey_type), len(current_predator_type))
    prey_predator_without_mutation_array_min = np.min(prey_predator_without_mutation_array)

    prey_predator_with_mutation_array = np.array(time_predation_with_birth_with_mutation).reshape(len(current_prey_type), len(current_predator_type))
    prey_predator_with_mutation_array_min = np.min(prey_predator_with_mutation_array)

    time_predator_intrinsic_death_min = min(time_predator_intrinsic_death)

    ### ★tau for reaction and renew the time 
    tau = min(prey_array_min, prey_predator_without_mutation_array_min, prey_predator_with_mutation_array_min, time_predator_intrinsic_death_min)

    t = t + tau # NOTE HERE! Things will go wrong when write 't += tau', still not knowing why...
    
    ### pick up the reaction:
    if tau == prey_array_min:
        prey_index_reation = np.where(prey_array == prey_array_min) 
        ## prey_index_reation of the min value: a tuple——((array([number]), array([number])))
        ## prey_index_reation[0] represent the reaction to happen; prey_index_reation[1] represent the species in prey_type list.
        if prey_index_reation[0][0] == 0: ## prey reproduction without mutation
            current_prey_type[prey_index_reation[1][0]].Nx += 1
        elif prey_index_reation[0][0] == 1: ## prey reproduction with mutation
            current_prey_type[prey_index_reation[1][0]].prey_mutation()
        else: ## prey die
            current_prey_type[prey_index_reation[1][0]].Nx -= 1

    elif tau == prey_predator_without_mutation_array_min:
        prey_predator_without_mutation_index = np.where(prey_predator_without_mutation_array == prey_predator_without_mutation_array_min)

        prey_index_death_eaten = prey_predator_without_mutation_index[0][0]
        predator_index_birth = prey_predator_without_mutation_index[1][0]

        current_prey_type[prey_index_death_eaten].Nx -= 1
        current_predator_type[predator_index_birth].Ny += 1        

    elif tau == prey_predator_with_mutation_array_min:
        prey_predator_with_mutation_index = np.where(prey_predator_with_mutation_array == prey_predator_with_mutation_array_min)

        prey_index_death_eaten_ = prey_predator_with_mutation_index[0][0]
        predator_index_birth_mutant = prey_predator_with_mutation_index[1][0]

        current_prey_type[prey_index_death_eaten_].Nx -= 1
        current_predator_type[predator_index_birth_mutant].predator_mutation()

    else:
        predator_index_death = time_predator_intrinsic_death.index(time_predator_intrinsic_death_min)
        current_predator_type[predator_index_death].Ny -= 1

    ## ❹Record the data.
    time_data.append(t)

    for i in current_prey_type:
        i.record_data(t)

    for l in current_predator_type:
        l.record_data(t)

# 5.Visualization of the population size.
all_prey_type = current_prey_type + extinct_prey_type
all_predator_type = current_predator_type + extinct_predator_type

# 画布1
plt.subplot(211)
plt.title("Predator Prey Coevolution model",fontsize=20)
plt.xlabel("Time", fontsize=14)
plt.ylabel("Number of prey individuals",fontsize=14)

for i in all_prey_type:
    plt.plot(i.time, i.population_size, linewidth=1)

plt.legend(['g = ' + str(i.gi) for i in all_prey_type], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

# 画布2
plt.subplot(212)
plt.xlabel("Time", fontsize=14)
plt.ylabel("Number of predator individuals",fontsize=14)

for l in all_predator_type:
    plt.plot(l.time, l.population_size, linewidth=1)

plt.legend(['k = ' + str(l.kl) for l in all_predator_type], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

plt.show()      
