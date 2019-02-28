#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy.stats import truncnorm
import matplotlib.pyplot as plt

"""A simulation of co-evolution of prey and predator(birth、death、competition
and de novo mutation) using Gillespie Algorithm.

* This version is for drawing a picture of stochastic co-evolutionary dynamics
of the prey and predator species.

* The new g and k value drawn randomly from an Truncated Normal Distribution.

* Note: This computational model is based on Weini's work in 2017.
  See at: https://github.com/WeiniHuangBW/DynamicalTradeoffandCoevolution
  For more detailed information about this model, please visit:
  https://www.nature.com/articles/s41467-017-01957-8#Sec8
"""


class Prey:
    """
    A class to represent a type of prey species.

    Attributes:
        num: Number of this prey species.
        g_value: Growth-defense trade off parameter, g ∈ (0,1).
        num_list: A list that store the number of this species after a
            reaction happen.
        time_list: A list that store the time.
    """

    def __init__(self, num, g_value):
        """Initialize attributes to describe a kind of prey."""
        self.g_value = float('%.3f' % g_value)
        self.num = num
        self.num_list = []
        self.time_list = []

    def prey_reproduce(self):
        """Prey reproduces without mutation."""
        self.num += 1

    def prey_die(self):
        """Prey dies because of competition, predation or dies naturally."""
        self.num -= 1

    def prey_mutation(self, sd_value):
        """Mutation of parental prey to generate a new prey species.

        When a prey reproduces, a mutation may occur with a small probability
        μx. The mutant is characterized by a new g value(type specific rate)
        drawn randomly from a truncated normal distribution between 0 and 1.
        The mean of this distribution is the g value of the parental prey.
        The lower bound is 0 and the upper bound is 1.

        Args:
            sd_value: The standard deviation of the truncated normal
                distribution that the mutant prey's g value is drawn
                from.
        """
        new_prey = Prey(num=1,
                        g_value=get_truncated_normal_random_number(
                            mean=self.g_value, sd=sd_value, low=0.0, up=1.0))
        current_prey.append(new_prey)

    def record_data(self, time):
        """Record the number of prey and time after each reaction happen."""
        self.time_list.append(time)
        self.num_list.append(self.num)


class Predator:
    """
    A class to represent a type of predator species.

    Attributes:
        num: Number of this predator species.
        k_value: k value is the ratio of predator growth to predation
            and represents the reproduction efficacy of a predator type.
        num_list: A list that store the number of this species after a
            reaction happen.
        time_list: A list that store the time.
    """

    def __init__(self, num, k_value):
        """Initialize attributes to describe a kind of predator."""
        self.k_value = float('%.3f' % k_value)
        self.num = num
        self.num_list = []
        self.time_list = []

    def predator_reproduce(self):
        """Predator reproduces without mutation."""
        self.num += 1

    def predator_die(self):
        """Predator dies naturally."""
        self.num -= 1

    def predator_mutation(self, sd_value):
        """Mutation of parental predator to generate a new prey species.

        With a probability μy, a predator produces a mutant with a new
        k value drawn from a truncated normal distribution between 0 and
        kmax, the upper limit of the reproduction efficacy. The mean of
        this distribution is the k value of the parental predator. The
        lower bound is 0 and the upper bound is 0.3.

        Args:
            sd_value: The standard deviation of the truncated normal
                distribution.
        """
        new_predator = Predator(num=1,
                                k_value=get_truncated_normal_random_number(
                                    mean=self.k_value, sd=sd_value, low=0.0,
                                    up=1.0))
        current_predator.append(new_predator)

    def record_data(self, time):
        """Record the number of prey and time after each reaction happen."""
        self.time_list.append(time)
        self.num_list.append(self.num)


def get_truncated_normal_random_number(mean, sd, low, up):
    """Get a random number that follow a truncated normal distribution.

    The truncated normal distribution is the probability distribution derived
    from that of a normally distributed random variable by bounding the random
    variable from either below or above (or both).

    Args:
        mean: The mean or expectation of the distribution.
        sd: The standard deviation of the distribution.
        low: The lower bound of the distribution.
        up: The upper bound of the distribution.

    Returns:
        A random number that follow a truncated normal distribution.
    """
    truncated_normal_generator = truncnorm(
        (low - mean) / sd,
        (up - mean) / sd,
        loc=mean,
        scale=sd)
    truncated_normal_random_number = truncated_normal_generator.rvs()
    return truncated_normal_random_number


def extinction_check(current_type, extinct_type, s=0):
    """Check extinction phenomenon.

    If one kind of prey or predator species go extinct, remove it from
    current_type to extinct_type list.

    Args:
        current_type: A list that stores the all the current type of prey
            or predator.
        extinct_type: A list that stores the extinct type of prey or predator.
        s: s=0. It used for counting the loop.
    """
    while s < len(current_type):
        if current_type[s].num == 0:
            extinct_type.append(current_type[s])
            del current_type[s]
        else:
            s += 1


def reaction_time(reaction_rate_array):
    """Calculate the reaction time with random number.

    In Gillespie Algorithm, we need to calculate the (reaction_rate)*ln(1/r).
    r is a random number from uniform distribution between 0 and 1.

    Args:
        reaction_rate_array: An array that store the reaction rate.

    Returns:
        An array of reaction time that are calculated from reaction rate.
    """
    array_shape = reaction_rate_array.shape
    if len(array_shape) == 1:
        calculate_reaction_time = - (1 / reaction_rate_array) * \
                                  np.log(np.random.rand(array_shape[0]))
    else:
        calculate_reaction_time = - (1 / reaction_rate_array) * \
                                  np.log(np.random.rand(array_shape[0],
                                                        array_shape[1]))

    return calculate_reaction_time


# 1. GLOBAL PARAMETERS
bx = 1.0  # baseline growth rate of prey
dx = 0.1  # intrinsic death rate of prey
dy = 0.5  # intrinsic death rate of predator
rc = 0.00005  # resource competition coefficient
ux = 0.0001  # mutation rate of prey per division
uy = 0.001  # mutation rate of predation per division
p = 0.005  # scaling coefficient of the predation rate
m = 3  # determines the initial shape of the growth-defense trade-off

# 2. Set the initial conditions.
t = 0  # Star time
T = 500  # Time period
ancestor_prey = Prey(num=1000, g_value=1.0)
ancestor_predator = Predator(num=100, k_value=0.3)

current_prey = [ancestor_prey]  # Store the current prey types.
current_predator = [ancestor_predator]  # Store the current predator types.

extinct_prey = []  # Store the extinct prey type.
extinct_predator = []  # Store the extinct predator type.

# 3. MAIN LOOP: The Gillespie Algorithm to simulate the co-evolution system.
while t <= T:

    extinction_check(current_prey, extinct_prey)
    extinction_check(current_predator, extinct_predator)

    if len(current_prey) == 0 or len(current_predator) == 0:
        print('Extinction Phenomenon')
        break

    prey_g_array = np.array([prey.g_value for prey in current_prey])
    prey_num_array = np.array([prey.num for prey in current_prey])
    total_prey_num = prey_num_array.sum()

    prey_birth_nonmutant = bx * (1 - ux) * prey_num_array * prey_g_array
    prey_birth_mutant = bx * ux * prey_g_array * prey_num_array
    prey_competition_death = rc * total_prey_num * prey_num_array
    prey_intrinsic_death = dx * prey_num_array

    predator_k_array = np.array(
        [predator.k_value for predator in current_predator])
    predator_num_array = np.array(
        [predator.num for predator in current_predator])
    predation_rate = p * prey_g_array[:, np.newaxis] ** \
        (m * predator_k_array / 0.3)

    predation_no_birth = predator_num_array * (1 - predator_k_array) * \
        prey_num_array[:, np.newaxis] * predation_rate
    predation_birth_nonmutant = predator_num_array * predator_k_array * \
        (1 - uy) * prey_num_array[:, np.newaxis] * predation_rate
    predation_birth_mutant = predator_num_array * predator_k_array * \
        uy * prey_num_array[:, np.newaxis] * predation_rate
    predator_intrinsic_death = dy * predator_num_array

    # ❷Calculate the reaction time with a random number
    t_prey_birth_nonmutant = reaction_time(prey_birth_nonmutant)
    t_prey_birth_mutant = reaction_time(prey_birth_mutant)
    t_prey_competition_death = reaction_time(prey_competition_death)
    t_prey_intrinsic_death = reaction_time(prey_intrinsic_death)

    t_predation_no_birth = reaction_time(predation_no_birth)
    t_predation_birth_nonmutant = reaction_time(predation_birth_nonmutant)
    t_predation_birth_mutant = reaction_time(predation_birth_mutant)
    t_predator_intrinsic_death = reaction_time(predator_intrinsic_death)

    # ❸Pick up the shortest time for the reaction.
    min_t_prey_birth_nonmutant = t_prey_birth_nonmutant.min()
    min_t_prey_birth_mutant = t_prey_birth_mutant.min()
    min_t_prey_competition_death = t_prey_competition_death.min()
    min_t_prey_intrinsic_death = t_prey_intrinsic_death.min()

    min_t_predation_no_birth = t_predation_no_birth.min()
    min_t_predation_birth_nonmutant = t_predation_birth_nonmutant.min()
    min_t_predation_birth_mutant = t_predation_birth_mutant.min()
    min_t_predator_intrinsic_death = t_predator_intrinsic_death.min()

    tau = min(
        min_t_prey_birth_nonmutant, min_t_prey_birth_mutant,
        min_t_prey_competition_death, min_t_prey_intrinsic_death,
        min_t_predation_no_birth, min_t_predation_birth_nonmutant,
        min_t_predation_birth_mutant, min_t_predator_intrinsic_death)

    t = t + tau  # Renew the time.

    # ❹Select which reaction takes place and update the number.
    if tau == min_t_prey_birth_nonmutant:
        index_one = t_prey_birth_nonmutant.argmin()
        current_prey[index_one].prey_reproduce()

    elif tau == min_t_prey_birth_mutant:
        index_two = t_prey_birth_mutant.argmin()
        current_prey[index_two].prey_mutation(sd_value=1.0)

    elif tau == min_t_prey_competition_death:
        index_three = t_prey_competition_death.argmin()
        current_prey[index_three].prey_die()

    elif tau == min_t_prey_intrinsic_death:
        index_four = t_prey_intrinsic_death.argmin()
        current_prey[index_four].prey_die()

    elif tau == min_t_predation_no_birth:
        index_five = t_predation_no_birth.min(1).argmin()
        current_prey[index_five].prey_die()

    elif tau == min_t_predation_birth_nonmutant:
        index_six = t_predation_birth_nonmutant.min(1).argmin()
        index_seven = t_predation_birth_nonmutant.min(0).argmin()

        current_prey[index_six].prey_die()
        current_predator[index_seven].predator_reproduce()

    elif tau == min_t_predation_birth_mutant:
        index_eight = t_predation_birth_mutant.min(1).argmin()
        index_nine = t_predation_birth_mutant.min(0).argmin()

        current_prey[index_eight].prey_die()
        current_predator[index_nine].predator_mutation(sd_value=0.3)

    elif tau == min_t_predator_intrinsic_death:
        index_ten = t_predator_intrinsic_death.argmin()
        current_predator[index_ten].predator_die()

    # Record the time and number data after each reaction.
    for prey in current_prey:
        prey.record_data(t)
    for predator in current_predator:
        predator.record_data(t)

# 4. Visualization of the population dynamics.
if len(current_prey) != 0 and len(current_predator) != 0:
    all_prey = current_prey + extinct_prey
    all_predator = current_predator + extinct_predator

    plt.subplot(211)
    for prey in all_prey:
        if max(prey.num_list) >= 5:
            plt.plot(prey.time_list, prey.num_list,
                     linewidth=1, label='g='+str(prey.g_value))
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2)

    plt.subplot(212)
    for predator in all_predator:
        if max(predator.num_list) >= 3:
            plt.plot(predator.time_list, predator.num_list,
                     linewidth=1, label='k='+str(predator.k_value))
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2)

    plt.show()
