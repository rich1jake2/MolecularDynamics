# MolecularDynamics
Showing the videos for molecular dynamics simulations in multiple dimensions

The first thing we've done is developed the code. Since we've developed this project before in other repositiories we've made a few changes to this code. First, we've chosen to use classes for the particles instead of just arrays and positions. This makes it a bit easier to track which object we update in each time step. We also use the velocity verlet method as suggested in the Computational Physics Textbook by Landau. 

##  The Code  

The most important part of the molecular dynamics simulations is the force accumulation and the subsequen position and velocity updates. Here is the class we used for the particles. 



    class particle:
        # intializing the particle attributes 
        def __init__(self, m0 = 1.0, pos0 = [0,0,0], v0 = [0,0,0], fDegs = 3 ):
            self.mass = m0 
            self.dim = fDegs
            self.pos = np.array(pos0)
            self.vel = np.array(v0)
            self.traj = [ pos0 ]
            self.energy = .5 * m0 * nr(v0)**2
        


Note our natrual degrees of freedom will be 3, so if we want to go in two dimensions then we have to set our z vectors to 0. Thus, we have many commented intializing conditions that we can present and then initialize a list of particle class instances - here is a simple example of two particles which lie along the same line in 2 dimensions (this is what the d2 boolean variable is for when we get to graphing - so we only graph in 2D)

    randomPos = [ [bound*.25,1.,0.], [.75*bound,1.,0.]]
    randomVel = [ [3.,0.,0.], [-3.,0.,0.] ]
    nParticles = 2
    d2 =True
    
    particles = [particle(pos0 = randomPos[i], v0 = randomVel[i]) for i in range(nParticles)]

Then we do some initializing of the time step h, necessary lists to keep track of the trajectories, previous forces (for verlet integration), and the energies

    pETavg = 0.0 
    kETavg = 0.0
    eTotAvg = 0.0
    velocities = []
    h = .01
    particleTrajectories = [ [[*k]]  for k in randomPos]
    fPrev = np.zeros((nParticles,3))
    eList = []


To calculate the force we introduce the following function - which is the derivative of the Lennard-Jones potential with respect to r - the euclidean distance between 2 particles 
    
    def lennardJonesForce(r, sig = 1.0, epi = 1.0 ):
        if r > 5:
            return 0  
        sigma = sig
        epsi = epi
        sigR = sigma/r
    return 48*epsi*( sigR**13 - .5*sigR**7 )

Because we will be handling Periodic Boundary Conditions, the following functions will be called to make the code more readable. The first one is a replication of Landau's sign function to ensure contiuity. The 'ghostDist' function will call the positions of the two particles and ajdust the distances which should be used for the force calculations with Periodic Boundary conditions. The periodic boundary conditions function is called to update the particles position if it crosses a boundary. 


    def sign(a, b):
        if b>= 0.0:
            return abs(bound)
        else: 
            return -1*abs(bound)

    def ghostDist(r1,r2,bound0 = 0, boundf = bound):
        # updated ghosts 
        diff = r1 - r2
        for el in range(len(diff)): 
            if abs(diff[el])  > boundf/2:
                diff[el] = diff[el] - sign(boundf, diff[el])
        return diff, nr(diff)



    # function to update if the boundary conditions are satisfied 
    def periodicBoundary(r,bound0 = 0.0 ,boundf = bound):
        rf = r
        for ll in range(len(r)):
            if r[ll] > boundf: 
                rf[ll] -= boundf
            if r[ll] < bound0:
                rf[ll] += boundf
            else: 
                rf[ll] = rf[ll]
        return np.array(rf)
  
Next, the bulk of the work is being done by the following triple nested for loop. Outside for loop is for each time step, the nested for loop in this is used to calculate the accumulated forces on each particle. For a large number of timesteps I added a computer timeer to keep track of every two thousandths time step 

       for i in range(timeSteps):
        Fjk = np.zeros((nParticles, 3))
        energyTotal = 0.0
        pE = 0.0
        kE = 0.0
        velocities.append([])
        if i%2000 == 0:
            end = time.time()
            print(i, start - end )
        for j in range(nParticles - 1 ):
            # For ghost particle interactions w/ periodic boundaries
            rj = particles[j].pos 
            for k in range(j+1, nParticles):
                # Temp force 
                Ftemp = 0.0  
                rk = particles[k].pos
                sep, r = ghostDist(rj,rk)
                Ftemp = lennardJonesForce( r ) 
                Fjk[j] = Fjk[j] +  (Ftemp)*(sep/r)
                Fjk[k] = Fjk[k] -  (Ftemp)*(sep/r)

                pE += 4.*((1/(r**6))) * ( 1/((r**6)) - 1.)/nParticles

        pETavg += pE/(h*(i+1)) 
        energyTotal += pE
        
        # Updating Velocities 
        for n in range(nParticles):
            particles[n].vel = particles[n].vel + (Fjk[n] + fPrev[n])*(h/2)*(1/particles[n].mass) 
            kE +=(nr(particles[n].vel)**2 * (.5) * (particles[n].mass))/nParticles
            velocities[i].append( nr(particles[n].vel) )
        # Updating the Positions
        for m in range(nParticles):

            # updating positions with velocity verlet method
            particles[m].pos = particles[m].pos + particles[m].vel*h + .5*Fjk[m]*h**2
            particles[m].pos =  periodicBoundary([*(particles[m].pos)],boundf = bound)  


            # adding the positions to the trajectory
            particleTrajectories[m].append( particles[m].pos )

        

        energyTotal += kE
        kETavg += kE/(h*(i+1))
        eTotAvg += energyTotal/(h*(i+1) )
        eList.append([kE,pE,energyTotal])
        fPrev = Fjk

We recognize the compuational ineffeciences of the for loops - there ought to be a better way to implment this with python classes AND python's ability to handle lists - however, a lack of time and experience has made this task difficult for myself. After this for loop the primary goal is to create an animation - the following lines are made to make the plotting as simple as possible 

    eList = np.array(eList).T
    # getting the scatters for each particle 
    pltTrajs = np.array([ particleTrajectories[i]  for i in range(nParticles)])
    
This next portion is meant to create the histogram plot and line plots for the velocity distribution. Note: we won't actually use matplotlib's histrogram plot - we will use their bar plot becasue the histrogram density attribute changes the scales that we want to deal with and I wanted the velocity distrubtion to look more like a probability distribution where 1 is the maximum value 

    velocities.sort()
    binz = 20
    velXPlt , sep = np.linspace( 0, 20, binz, retstep = True )
    #print(len(velXPlt))
    #print(velXPlt)
    binNumsPlt = []
    separationList = []
    for j in range(len(velocities)) :
        binNumsPlt.append( [*np.zeros(binz)] )
        for i,k in enumerate(velXPlt):
            for ii, el2 in enumerate(velocities[j]):
               if el2 <= k and el2 > velXPlt[i-1]:
                   binNumsPlt[j][i] += 1/nParticles








