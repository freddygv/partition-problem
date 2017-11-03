import sys
import random
import time

list_size = 100
trials = 50
iterations = 25000

file_name = sys.argv[1]

def randomSolution(in_list):
    seq = [0 for i in range(len(in_list))]

    # Initializing a random solution S by randomly assigning -1 or 1 to each index
    for i in range(len(in_list)):
        seq[i] = random.choice([-1,1])
        
    return seq

def randomPrepartitionSolution(in_list):
    size = len(in_list)
    
    P = [random.randint(0,size-1) for i in range(size)]
    
    return P


def reconcilePrepartition(p, in_list):
    new_list = [0 for i in range(len(in_list))]
    
    for j in range(len(in_list)):
        new_list[p[j]] = new_list[p[j]] + in_list[j]
    
    return new_list


def calculateResidue(sequence, in_list):
    residue = 0
    
    # Calculating the residue term by adding the product of the input and output sequence elements
    for i in range(len(in_list)):
        residue += sequence[i] * in_list[i]
        
    return abs(residue)


def KK(orig_list):
    in_list = [i for i in orig_list] # Deep copy of the array
    residue = 0
    
    while True:
        # Going to keep track of the highest and second highest values
        curr_max = 0
        next_max = 0

        curr_max_idx = -1
        next_max_idx = -1

        # Looping through to find the two highest values
        for i in range(len(in_list)):
            if in_list[i] > curr_max:
                next_max     = curr_max
                next_max_idx = curr_max_idx

                curr_max     = in_list[i]
                curr_max_idx = i

            elif in_list[i] > next_max:
                next_max     = in_list[i]
                next_max_idx = i
                
        if next_max == 0:
            residue = curr_max
            break
            
        # Residue is their difference
        residue = curr_max - next_max
        
        # Replacing the max with the difference, and the other value with 0
        in_list[next_max_idx] = 0
        in_list[curr_max_idx] = residue
        
    return residue


def repeatedRandom(original_seq, original_res, in_list, iterations):

    for j in range(iterations):
        # Randomly generating a new solution and updating the residue
        temp_seq = randomSolution(in_list)
        new_residue = calculateResidue(temp_seq, in_list)

        # If the residue was reduced, replace original solution
        if new_residue < original_res:
            original_seq = temp_seq
            original_res = new_residue
        else:
            continue
            
    return original_res


def repeatedRandomPP(original_p, original_res, in_list, iterations):
    """
    Prepartitioned Repeated Random
    """

    for j in range(iterations):
        # Randomly generating a completely new sequence P and updating the residue
        new_p = randomPrepartitionSolution(in_list)
        
        new_list    = reconcilePrepartition(new_p, in_list)
        new_residue = KK(new_list)

        # If the residue was reduced, replace original solution
        if new_residue < original_res:
            original_p   = new_p
            original_res = new_residue
        else:
            continue
            
    return original_res


def hillClimbing(original_seq, original_res, in_list, iterations):
    
    for j in range(iterations):
        temp_seq = original_seq

        # Picking two random indeces for the sequence
        k = random.randint(0,len(in_list)-1)
        l = k
        
        # The indices cannot be equal
        while l == k:
            l = random.randint(0,len(in_list)-1)
        
        # Converting a single value in the sequence
        temp_seq[k] = temp_seq[k] * -1 
        
        # Flip a coin to see if the second index should be moved
        if 0.5 > random.uniform(0, 1):
            temp_seq[l] = temp_seq[l] * -1 

        new_residue = calculateResidue(temp_seq, in_list)

        # If the residue was reduced, replace original solution
        if new_residue < original_res:
            original_seq = temp_seq
            original_res = new_residue
        else:
            continue
            
    return original_res


def hillClimbingPP(original_p, original_res, in_list, iterations):
    """
    Prepartitioned Hill Climbing
    """
    
    for j in range(iterations):
        new_p = original_p
        
        # Picking two random indeces for the sequence
        k = random.randint(0,len(in_list)-1)
        l = k
        
        while original_p[l] == k:
            l = random.randint(0,len(in_list)-1)
        
        # Converting a single value in the sequence and generating new residue
        new_p[l] = k
        new_list    = reconcilePrepartition(new_p, in_list)
        new_residue = KK(new_list)

        # If the residue was reduced, replace original solution
        if new_residue < original_res:
            original_p = new_p
            original_res = new_residue
        else:
            continue
            
    return original_res


def coolingSchedule(iteration):
    # Calculates the decreasing temperature for the simulatedAnnealing Algorithm
    return 10**10*((0.8)**math.floor(iteration/300))


def pMoveDown(new_res, old_res, temp):
    # Calculates the probability that we should move to a position that is worse off
    return math.exp(-abs(new_res-old_res)/temp)


def simulatedAnnealing(original_seq, original_res, in_list, iterations):
    
    for j in range(iterations):
        temperature = coolingSchedule(j)

        temp_seq = original_seq
        temp_res = original_res
        
        # Picking two random indeces for the sequence
        k = random.randint(0,len(in_list)-1)
        l = k
        
        # The indices cannot be equal
        while l == k:
            l = random.randint(0,len(in_list)-1)
            
        # Flip a coin to see if the second index should be moved
        if 0.5 > random.uniform(0, 1):
            temp_seq[l] = temp_seq[l] * -1 

        # Converting a single value in the sequence
        temp_seq[k] = temp_seq[k] * -1 
        new_residue = calculateResidue(temp_seq, in_list)

        # If the residue was reduced, replace original solution
        if new_residue < original_res:
            original_seq = temp_seq
            original_res = new_residue
        else:
            p = pMoveDown(new_residue, original_res, temperature)
            
            if p > random.uniform(0, 1):
                original_seq = temp_seq
                original_res = new_residue
            else:
                continue
        
    return original_res

def simulatedAnnealingPP(original_p, original_res, in_list, iterations):
    """
    Prepartitioned Simulated Annealing
    """
    
    optimal_p = original_p
    optimal_res = original_res
    
    for j in range(iterations):
        temperature = coolingSchedule(j)
        
        new_p = original_p
        temp_res = original_res
        
        # Picking a random index for the sequence
        k = random.randint(0,len(in_list)-1)
        l = k
        
        while original_p[l] == k:
            l = random.randint(0,len(in_list)-1)
        
        # Converting a single value in the sequence and generating new residue
        new_p[l]    = k
        new_list    = reconcilePrepartition(new_p, in_list)
        new_residue = KK(new_list)

        # If the residue was reduced, replace original solution
        if new_residue < original_res:
            optimal_p  = new_p
            optimal_res  = new_residue
            
            original_p = new_p
            original_res = new_residue
        else:
            p = pMoveDown(new_residue, original_res, temperature)
            
            if p > random.uniform(0, 1):
                original_p = new_p
                original_res = new_residue
            else:
                continue
        
    return optimal_res


def read_file(input_file):
    array = [int(item.rstrip()) for item in open(input_file)]

    return array


def main():
    input_list = read_file(file_name)

    print KK(input_list)

if __name__ == "__main__":
    main()
