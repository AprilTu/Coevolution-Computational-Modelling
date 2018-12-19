"""A simulation of coevolution of prey and predator
( birth、death、competition and mutation) using Weini's Gillespie Algorithm """

import numpy as np

import matplotlib.pyplot as plt

class Prey():
    """ 
    A class to represent a type of prey with follwing attribute:
    ❶ Number of this species of prey: Nx
    ❷ Growth-defense trade off parameter: gi ∈ (0,1) 
      (following an uniform distribution)
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

class Predator():
    """
    A class to represent a type of predator with following attribute:
    ❶ Number of this species of predator: Ny
    ❷ Reproduction efficacy(defined as the ratio of the reproduction rate to 
    the predation rate of predators): kl ∈ (0, kmax], kmax=0.3
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

def prey_mutation():
    """Mutation of prey to generate a new species."""
    new_prey = Prey(Nx=1, gi=np.random.random_sample())
    current_prey_type.append(new_prey) # 这个会不会覆盖掉前面的new prey呀？不会呀因为new prey只在这个方法里有效，不是全局变量

def predator_mutation():
    """Mutation of predator to generate a new species"""
    new_predator = Predator(Ny=1, kl=np.random.uniform(0,0.3))
    current_predator_type.append(new_predator)

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
        ls_1 = []

        for l in current_predator_type:
            element = i.Nx * l.Ny * (1-l.kl) * predation_rate(i.gi, l.kl)
            ls_1.append(element)
        
        ls_0.append(sum(ls_1))

    return ls_0

def predation_with_reproduction_with_mutation(current_prey_type, current_predator_type):
    ls_2 = []
    for i in current_prey_type:
        ls_3 = []
        for l in current_predator_type:
            element_1 = i.Nx * l.Ny * l.kl * predation_rate(i.gi, l.kl) *uy
            ls_3.append(element_1)

        ls_2.append(sum(ls_3))

    return ls_2

def reaction_time(reaction_rate):
    """Calculate the reaction time with random number"""
    reaction_time = - (1/reaction_rate) * np.log(np.random.random_sample())
    return reaction_time

# Fixed Parameters
## Prey
bx = 1  # float; baseline growth rate of prey
dx = 0.1  # float; intrinsic death rate of prey
rc = 0.00005  # float; resource competition coefficient
ux = 0.0001  # mutation rate of prey per division

## Predation
p = 0.005  # float; scaling coefficient of the predation rate
dy = 0.5  # float; intrinsic death rate of predator
m = 3  # int; determines the initial shape of the growth-defense trade-off curve for the prey species
uy = 0.001  # mutation rate of predation per division
kmax = 0.6

## Time
t = 0.0  # float; start time
T = 500.0  # # float; maximum elapsed time

# Set the initial conditions
ancestor_prey = Prey(Nx=1000, gi=1)
ancestor_predator = Predator(Ny=100, kl=kmax)

# List that collect all the data
time_data = []  # record time

current_prey_type = [ancestor_prey]  # record current prey 
current_predator_type = [ancestor_predator]  # record current predator
extinct_prey_type = []  # record the extinct prey
extinct_predator_type = []  # record the extinct predator

total_prey_population_data = []
total_predator_population_data = []

# Main Loop: the Gillespie Algorithem
while t <= T:
    # ★First check extinction phenomenon
    prey_extinction_check(current_prey_type, extinct_prey_type)
    predator_extinction_check(current_predator_type,extinct_predator_type)

    total_prey_number = sum([i.Nx for i in current_prey_type])
    total_predator_number = sum([i.Ny for i in current_predator_type])
    
    # ❶Define the reaction rates.
    # 列表解析法生成列表，顺序对应着 current_type 的顺序，便于后续检索： 很优秀！
    prey_birth_without_mutation = [(1-ux) * i.gi * i.Nx for i in current_prey_type]
    prey_birth_with_mutation = sum([ux * i.gi * i.Nx for i in current_prey_type])
    prey_competition_death = [rc * i.Nx * total_prey_number for i in current_prey_type]
    prey_intrinsic_death = [dx * i.Nx for i in current_prey_type]
    

    # 有predator的反应
    # A list 
    predation_without_birth = predation_without_reproduction(current_prey_type, current_predator_type)

    #### 处理一下最困难的二维数据：既需要知道i in Prey 也需要知道 l in Predator
    ### 掌握列表解析，好像不是很难呀——列表解析从左到右扫描
    ### 后面需要转化为数组
    predation_with_birth_without_mutation = [i.Nx * l.Ny * l.kl * predation_rate(i.gi, l.kl) * (1-uy) for i in current_prey_type for l in current_predator_type]

    # 可能会出问题 重写 >> ok
    predation_with_birth_with_mutation = predation_with_reproduction_with_mutation(current_prey_type, current_predator_type)

    predator_intrinsic_death = [dy * l.Ny for l in current_predator_type]

    ## ❷calculate the reation time with a random number: (1/ai) * log(1/ri).
    time_prey_birth_without_mutation = list(map(reaction_time, prey_birth_without_mutation))
    time_prey_birth_with_mutation = reaction_time(prey_birth_with_mutation)
    time_prey_competition_death = list(map(reaction_time, prey_competition_death))
    time_prey_intrinsic_death = list(map(reaction_time, prey_intrinsic_death))

    time_predation_without_birth = list(map(reaction_time, predation_without_birth))
    time_predation_with_birth_without_mutation = list(map(reaction_time, predation_with_birth_without_mutation))
    time_predation_with_birth_with_mutation = list(map(reaction_time, predation_with_birth_with_mutation ))
    time_predator_intrinsic_death = list(map(reaction_time, predator_intrinsic_death))

    ## ❸Pick up the shortest time for the reaction and Select which reaction takes place and update the number.
    prey_ndarray = np.array([time_prey_birth_without_mutation, time_prey_competition_death, time_prey_intrinsic_death, time_predation_without_birth, time_predation_with_birth_with_mutation])  ### generate a n-darray for latter analysis
    prey_ndarray_min = np.min(prey_ndarray)

    prey_predator_array = np.array(time_predation_with_birth_without_mutation).reshape(len(current_prey_type),len(current_predator_type))
    prey_predator_array_min = np.min(prey_predator_array)

    time_predator_intrinsic_death_min = min(time_predator_intrinsic_death)

    ### tau for reaction and renew the time 
    tau = min([prey_ndarray_min, prey_predator_array_min, time_prey_birth_with_mutation, time_predator_intrinsic_death_min])
    
    t += tau

    ### pick up the reaction:
    if tau == time_prey_birth_with_mutation:
        prey_mutation()

    elif tau == time_predator_intrinsic_death_min:
        predator_index_death = time_predator_intrinsic_death.index(time_predator_intrinsic_death_min)
        current_predator_type[predator_index_death].Ny -= 1

    elif tau == prey_ndarray_min:
        prey_index_reation = np.where(prey_ndarray == np.min(prey_ndarray)) 
        ## prey_index_reation of the min value: a tuple——((array([number]), array([number])))
        ## prey_index_reation[0] represent the reaction to happen; prey_index_reation[1] represent the species in prey_type list.
        if prey_index_reation[0][0] == 0: ## reproduction without mutation:
            current_prey_type[prey_index_reation[1][0]].Nx += 1

        elif prey_index_reation[0][0] == 1 or prey_index_reation[0][0] == 2 or prey_index_reation[0][0] == 3:
            current_prey_type[prey_index_reation[1][0]].Nx -= 1

        else:
            current_prey_type[prey_index_reation[1][0]].Nx -= 1
            predator_mutation()

    else:
        prey_predator_index = np.where(prey_predator_array == np.min(prey_predator_array))

        prey_index_death_eaten = prey_predator_index[0][0]
        predator_index_birth = prey_predator_index[1][0]

        current_prey_type[prey_index_death_eaten].Nx -= 1
        current_predator_type[predator_index_birth].Ny += 1

    ## Record the data.
    time_data.append(t)

    for i in current_prey_type:
        i.record_data(t)

    for l in current_predator_type:
        l.record_data(t)

# Visualization of the population size.
all_prey_type = current_prey_type + extinct_prey_type
all_predator_type = current_predator_type + extinct_predator_type

# 画布1
plt.subplot(211)

for i in all_prey_type:
    if max(i.population_size) >= 5:
        plt.plot(i.time, i.population_size, linewidth=1, label='g='+str(i.gi))
plt.legend(bbox_to_anchor=(1.05, 1), loc=2)

# 画布2
plt.subplot(212)
for l in all_predator_type:
    if max(l.population_size) >= 4:
        plt.plot(l.time, l.population_size, linewidth=1, label='k='+str(l.kl))
plt.legend(bbox_to_anchor=(1.05, 1), loc=2)

plt.show()  
