#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy.stats import truncnorm
import csv

"""A simulation of co-evolution of prey and predator(birth、death、competition
and de novo mutation) using Gillespie Algorithm. 

This computational model is based on Weini's work in 2017.
See at: https://github.com/WeiniHuangBW/DynamicalTradeoffandCoevolution
For more detailed information about this model, please visit:
https://www.nature.com/articles/s41467-017-01957-8#Sec8
"""


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


def new_prey_g(original_g, sd_value):
    """Get a new g value for a new prey.

    When a prey reproduces, a mutation may occur with a small probability μx.
    The mutant is characterized by a new g value(type specific rate) drawn
    randomly from a truncated normal distribution between 0 and 1. The mean
    of this distribution is the g value of the parental prey. The lower
    bound is 0 and the upper bound is 1.

    Args:
        original_g: The g value of the parental prey.
        sd_value: The standard deviation of the truncated normal distribution.

    Returns:
        A new g value. g ∈ (0,1).
    """
    new_g = get_truncated_normal_random_number(
        mean=original_g, sd=sd_value, low=0.0, up=1.0)
    return new_g


def new_predator_k(original_k, sd_value):
    """Get a new k value for a new predator.

    When a predator reproduces, a mutation may occur with a small probability
    μy. The mutant is characterized by a new k value drawn randomly from a
    truncated normal distribution between 0 and kmax. The mean of this
    distribution is the k value of the parental predator. The lower
    bound is 0 and the upper bound is 0.3.

    Args:
        original_k: The k value of the parental predator.
        sd_value: The standard deviation of the truncated normal distribution.

    Returns:
        A new k value. k ∈ (0,kmax).
    """
    new_k = get_truncated_normal_random_number(
        mean=original_k, sd=sd_value, low=0.0, up=0.3)
    return new_k


def extinction_check(current_type_number, current_type_attribute, s=0):
    """Delete extinct species.

    If one kind of prey or predator species go extinct, remove it from
    current_type. That means delete its number and g value in list
    [current_type_number] and [current_type_attribute].

    Args:
        current_type_number: A list that stores the population size of the
            species.
        current_type_attribute: A list that stores the attribute of the
        species. (g value for prey and k value for predator.)
        s: s=0. It used for counting the loop.
    """
    while s < len(current_type_number):
        if current_type_number[s] == 0:
            del current_type_number[s]
            del current_type_attribute[s]
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


def shannon_index(species_num_array):
    """Calculate the Shannon index: -∑pilnpi.

    The Shannon index is an information statistic index, which means it
    assumes all species are represented in a sample and that they are
    randomly sampled. p is the proportion (n/N) of individuals of one
    particular species found (n) divided by the total number of
    individuals found (N). The value of this index ranges between
    0 and 1, the greater the value, the greater the sample diversity.

    Args:
        species_num_array: An array that store the number of different kind
            of species.

    Returns:
        Shannon index of diversity of this population.
    """
    ratio = species_num_array / species_num_array.sum()
    shannon_index_diversity = - sum(ratio * np.log(ratio))
    return shannon_index_diversity


def simpson_index(species_num_array):
    """Calculate the Simpson's Diversity Index: 1 - ∑pi**2

    The Simpson index is a dominance index because it gives more weight to
    common or dominant species. In this case, a few rare species with only
    a few representatives will not affect the diversity.  p is the proportion
    (n/N) of individuals of one particular species found (n) divided by the
    total number of individuals found (N). The value of this index ranges
    between 0 and 1, the greater the value, the greater the sample
    diversity.

    Args:
        species_num_array: An array that store the number of different kind
            of species.

    Returns:
        Simpson's diversity index of this population.
    """
    ratio_ = species_num_array / species_num_array.sum()
    simpson_index_diversity = 1 - sum(ratio_**2)
    return simpson_index_diversity


# 1. GLOBAL PARAMETERS
bx = 1.0  # baseline growth rate of prey
dx = 0.1  # intrinsic death rate of prey
dy = 0.5  # intrinsic death rate of predator
rc = 0.00005  # resource competition coefficient
ux = 0.0001  # mutation rate of prey per division
uy = 0.001  # mutation rate of predation per division
p = 0.005  # scaling coefficient of the predation rate
m = 0.5  # determines the initial shape of the growth-defense trade-off

# 2. MAIN LOOP
frequency_count = []
# For each cycle, record 0, 1 or 2 in this list.
# NOTE:
# 0 stands for both prey and predator go extinct.
# 1 stands for only predator go extinct.
# 2 stands for coexist.

headers_1 = ['Prey Types', 'Predator Type',
             'Prey Shannon Index', 'Predator Shannon Index',
             'Prey Simpson Index', 'Predator Simpson Index']
with open('data1.csv', 'a') as file_object_1:
    file_csv_1 = csv.writer(file_object_1)
    file_csv_1.writerow(headers_1)

cycle_number_count = 0
cycle_number = 50

while cycle_number_count < cycle_number:
    cycle_number_count += 1

    # For each run, set the initial conditions
    ancestor_prey_g = 1.0
    ancestor_prey_num = 1000
    ancestor_predator_k = 0.3
    ancestor_predator_num = 100
    t = 0.0  # start time
    T = 2000  # time period

    # Store current prey's and  predator's attributes.
    # record g value of each prey species
    current_prey_g = [ancestor_prey_g]
    # record the population size of each prey species
    current_prey_nx = [ancestor_prey_num]
    # record k value of each predator species.
    current_predator_k = [ancestor_predator_k]
    # record the population size of each predator species
    current_predator_ny = [ancestor_predator_num]

    # For one main simulation cycle,
    # we use Gillespie Algorithm to simulate prey-predator dynamic.
    while t <= T:
        # Check the species extinction phenomenon.
        extinction_check(current_prey_nx, current_prey_g)
        extinction_check(current_predator_ny, current_predator_k)

        # Break the inner loop when all prey go extinct and make record.
        if len(current_prey_nx) == 0:
            frequency_count.append(0)
            break

        # Break the inner loop when predator go extinct and make record.
        if len(current_predator_ny) == 0:
            frequency_count.append(1)
            break

        # Gillespie Algorithm
        # ❶Define the reaction rates.
        prey_g_array = np.array(current_prey_g)
        prey_nx_array = np.array(current_prey_nx)
        total_prey_number = prey_nx_array.sum()

        prey_birth_nonmutant = bx * (1-ux) * prey_g_array * prey_nx_array
        prey_birth_mutant = bx * ux * prey_g_array * prey_nx_array
        prey_competition_death = rc * total_prey_number * prey_nx_array
        prey_intrinsic_death = dx * prey_nx_array

        predator_k_array = np.array(current_predator_k)
        predator_ny_array = np.array(current_predator_ny)
        predation_rate = p * prey_g_array[:, np.newaxis] ** \
            (m*predator_k_array/0.3)

        predation_no_birth = predator_ny_array * (1-predator_k_array) * \
            prey_nx_array[:, np.newaxis] * predation_rate
        predation_birth_nonmutant = predator_ny_array * predator_k_array * \
            (1-uy) * prey_nx_array[:, np.newaxis] * predation_rate
        predation_birth_mutant = predator_ny_array * predator_k_array * \
            uy * prey_nx_array[:, np.newaxis] * predation_rate
        predator_intrinsic_death = dy * predator_ny_array

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

        t = t + tau  # renew the time

        # ❹Select which reaction takes place and update the number.
        if tau == min_t_prey_birth_nonmutant:
            index_one = t_prey_birth_nonmutant.argmin()
            current_prey_nx[index_one] += 1

        elif tau == min_t_prey_birth_mutant:
            index_two = t_prey_birth_mutant.argmin()
            mutant_g_value = new_prey_g(current_prey_g[index_two], sd_value=1)

            current_prey_nx.append(1)
            current_prey_g.append(mutant_g_value)

        elif tau == min_t_prey_competition_death:
            index_three = t_prey_competition_death.argmin()
            current_prey_nx[index_three] -= 1

        elif tau == min_t_prey_intrinsic_death:
            index_four = t_prey_intrinsic_death.argmin()
            current_prey_nx[index_four] -= 1

        elif tau == min_t_predation_no_birth:
            index_five = t_predation_no_birth.min(1).argmin()
            current_prey_nx[index_five] -= 1

        elif tau == min_t_predation_birth_nonmutant:
            index_six = t_predation_birth_nonmutant.min(1).argmin()
            index_seven = t_predation_birth_nonmutant.min(0).argmin()

            current_prey_nx[index_six] -= 1
            current_predator_ny[index_seven] += 1

        elif tau == min_t_predation_birth_mutant:
            index_eight = t_predation_birth_mutant.min(1).argmin()
            index_nine = t_predation_birth_mutant.min(0).argmin()
            mutant_k_value = new_predator_k(current_predator_k[index_nine],
                                            sd_value=0.3)

            current_prey_nx[index_eight] -= 1
            current_predator_ny.append(1)
            current_predator_k.append(mutant_k_value)

        else:
            index_ten = t_predator_intrinsic_death.argmin()
            current_predator_ny[index_ten] -= 1

    # This judgment using for thwart extinction phenomenon.
    if len(current_prey_nx) != 0 and len(current_predator_ny) != 0:
        # record coexistence phenomenon
        frequency_count.append(2)

        # calculate the diversity index.
        # Richness
        coexist_prey_types = len(current_prey_nx)
        coexist_predator_types = len(current_predator_ny)
        #  Evenness
        prey_shannon = shannon_index(np.array(current_prey_nx))
        predator_shannon = shannon_index(np.array(current_predator_ny))

        prey_simpson = simpson_index(np.array(current_prey_nx))
        predator_simpson = simpson_index((np.array(current_predator_ny)))

        row_1 = [coexist_prey_types, coexist_predator_types,
                 prey_shannon, predator_shannon,
                 prey_simpson, predator_simpson]

        with open('data1.csv', 'a') as file_object_1:
            file_csv_1 = csv.writer(file_object_1)
            file_csv_1.writerow(row_1)

# 3. DATA PROCESSING
both_extinct_frequency = frequency_count.count(0)
predator_extinct_frequency = frequency_count.count(1)
coexistence_frequency = frequency_count.count(2)

# 4. DATE OUTPUT
headers_2 = ['Both Extinction', 'Predator Extinction', 'Coexistence']
row_2 = [both_extinct_frequency,
         predator_extinct_frequency,
         coexistence_frequency]

with open('data2.csv', 'a') as file_object_2:
    file_csv_2 = csv.writer(file_object_2)
    file_csv_2.writerow(headers_2)
    file_csv_2.writerow(row_2)
