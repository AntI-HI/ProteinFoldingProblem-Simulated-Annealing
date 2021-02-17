from pyrosetta import *
from pyrosetta.toolbox import *
import statistics
from copy import copy
from itertools import combinations
from operator import attrgetter, mod
import random
import math
import copy

init()

pose = pose_from_sequence('EAEDLQVGQVELGGGPGAGSLQPLALEGSLQ')
scorefxn = get_fa_scorefxn()

# pose.dump_pdb("sequence.pdb")

# config1 = pyrosetta.rosetta.core.pose.deep_copy(pose)
# relax = pyrosetta.rosetta.protocols.relax.FastRelax()
# relax.set_scorefxn(scorefxn)
# relax.apply(config1)

# print(scorefxn(config1))

# config1.dump_pdb("output_pdb_0.pdb")


# returns the lowest-energy pose fold you can find for p by running
# Simulated Annealing. 
def optimizeFold(pose):

    bestPose = pyrosetta.rosetta.core.pose.deep_copy(pose)     # Initial comformation with very high energy. Whenever we found a better conformation we update the bestPose.

    T_init = 10 * pose.total_residue() # Initial Temperature.
    r = 0.001   # Decay Factor for Temperature.
    T = T_init
    stopping_condition = pose.total_residue() * 1000000
    limit = 0.1
    
    currentPose = bestPose
    minEnergy = scorefxn(bestPose) # Initial Energy.
    j = 0 # Iteration number
    k = 0
    for i in range(stopping_condition):
        
        breakloop = False

        while breakloop is False:
            seq = random.randint(1, currentPose.total_residue())
            rotate = random.uniform(0, limit)
            #print(rotate)

            currentEnergy = scorefxn(currentPose)
            rand = random.randint(0,50)
            if rand < 45:   
                currentPose.set_phi(seq, currentPose.phi(seq) + rotate)
            elif rand >= 45 and rand < 50:
                currentPose.set_phi(seq, currentPose.phi(seq) + rotate)
                currentPose.set_psi(seq, currentPose.psi(seq) + rotate)
            else:
                currentPose.set_psi(seq, currentPose.psi(seq) + rotate)
                currentPose.set_omega(seq, currentPose.omega(seq) + rotate)
            
            foldedEnergy = scorefxn(currentPose)
            if foldedEnergy < minEnergy or random.random() < math.exp((currentEnergy/10 - foldedEnergy/10)/ T):
                if foldedEnergy < minEnergy:
                    minEnergy = foldedEnergy
                    bestPose = pyrosetta.rosetta.core.pose.deep_copy(currentPose)
                    k = 0
                    print(minEnergy)
                else:
                    k += 1

                breakloop = True
                T = T_init * math.exp(-j * r)
                if T < 1: # Reheat
                    j = j/10
                # print("{}, {}".format(foldedEnergy, i))
                j += 1
                    
            else:
                if rand < 45:   
                    currentPose.set_phi(seq, currentPose.phi(seq) - rotate)
                elif rand >= 45 and rand < 50:
                    currentPose.set_phi(seq, currentPose.phi(seq) - rotate)
                    currentPose.set_psi(seq, currentPose.psi(seq) - rotate)
                else:
                    currentPose.set_psi(seq, currentPose.psi(seq) - rotate)
                    currentPose.set_omega(seq, currentPose.omega(seq) - rotate)
                k += 1

                if k == 20:
                    limit = limit * 0.5    
                    k = 0 
                    breakloop = True  
                    #print(i)
                if limit < 0.35 :
                    print(limit)
                    limit = 1                

    return bestPose


config2 = pyrosetta.rosetta.core.pose.deep_copy(pose)
optimized_pose = optimizeFold(config2)
foldscore2 = scorefxn(optimized_pose)
print(foldscore2)

optimized_pose.dump_pdb("output_pdb.pdb")
