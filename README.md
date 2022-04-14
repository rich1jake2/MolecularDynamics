# MolecularDynamics
Showing the videos for molecular dynamics simulations in multiple dimensions

The first thing we've done is developed the code. Since we've developed this project before in other repositiories we've made a few changes to this code. First, we've chosen to use classes for the particles instead of just arrays and positions. This makes it a bit easier to track which object we update in each time step. We also use the velocity verlet method as suggested in the Computational Physics Textbook by Landau. 

##  - Primary For Loop - 

The most important part of the molecular dynamics simulations is the force accumulation and the subsequen position and velocity updates. Here is the class we used for the particles. 

'''
class particle:
    # intializing the particle attributes 
    def __init__(self, m0 = 1.0, pos0 = [0,0,0], v0 = [0,0,0], fDegs = 3 ):
        self.mass = m0 
        self.dim = fDegs
        self.pos = np.array(pos0)
        self.vel = np.array(v0)
        self.traj = [ pos0 ]
        self.energy = .5 * m0 * nr(v0)**2
''' 

Note our natrual degrees of freedom will be 3, so if we want to go in two dimensions then we have to set our z vectors to 0. Thus, we have many commented intializing conditions that we can present. 



Next, the bulk of the work is being done by the following triple nested for loop. Outside for loop is for each time step, the nested for loop in this is used to calculate the accumulated forces on each particle. 

'''


'''
