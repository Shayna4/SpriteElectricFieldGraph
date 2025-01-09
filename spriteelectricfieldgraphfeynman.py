import math
import matplotlib.pyplot as plt
import numpy as np
import scipy
#check retarded and real positions
'''
this program produced for graphs to help better understand if and how the velocity
of charged particles should contrbute to the calculations predicting 
upper-atmosphere lightning
the first graph compares the retarded position of the charge to the real postion
graph number two calculates the electric field of the charge without the velocity
contributing to the results through the method of images and the boundaries between
60 and 90 km
the third graph attempts to use the same methods of caluculations used in graph 
two to calculate the velocity's contribution to the electric field
Graph four calculates the velocity's contribution to the electric field through
the Feynman method

i need to add more to this also what should i do about the first four charts
'''
class Feynman:
    
    def __init__(self,dz,Q, L, time):
        self.dz=dz
        self.c = scipy.constants.c
        self.Q = Q
        self.easy = -self.Q/(4*math.pi*scipy.constants.epsilon_0)
        self.observeKm = list(range(60, 91,1))
        self.L = L
        self.time = time
        self.n = list(range(-100,100,1))
        self.zMeas = np.array(self.observeKm*1000)
    
    '''
    calculate electrostatic field of charge with image charges
    '''
    def electricFieldOne(self,Q,observe,dist):
        E=0
        for n in self.n:
            E+= (Q*(observe-(2*n*self.L*1e3+dist)))/((observe-(2*n*self.L*1e3+\
                                                                dist))**2)**(3/2)\
                -(Q*(observe-(2*n*self.L*1e3-dist)))/((observe-(2*n*self.L*1e3\
                                                                -dist))**2)**(3/2)
        
        return E
    '''
    #calculate the total electrostatic field using method of images on velocity 
    component
    '''
    def electricFieldTwo(self,observe,dist,a,bList,cList,aReal,time):
        
        E= [0 for k in self.observeKm]
        for n in self.n:
           
            zImageChargeMinus = [2 * n *(90000) - d for d in dist]
            zImageChargePlus = [2 * n * (90000) + d for d in dist]
          

            rRetardedMinus = [(o*1e3) - zMinus for zMinus,o in \
                                 zip(zImageChargeMinus,observe)]
            rRetardedPlus = [(o*1e3) - zPlus for zPlus,o in \
                                zip(zImageChargePlus,observe)]
            
        
            E = [(l-1 / (rRetardedP**2) * (b**2 - 4.0 * a * c)**(-1.5) )for l,rRetardedP,b,c in zip(E,rRetardedPlus,bList,cList)]    

        return E    
   
    '''
    this function creates 4 graphs to help visually understand how the velocity 
    component of the electrostatic field should be included 
    '''
    def secondGraph(self,other,charge,posit):
        fif,ax= plt.subplots(2,2)
        time = .0002
        aReal = -2*1e10
        retardedP=[]
        realP=[]
        t=[]
        while time < .0011:
            realPosition = self.dz*1e3+.5*aReal*time**2
            realPos = [ self.dz*1e3+.5*aReal*time**2 for zMeas in self.observeKm]
            timeDelay = (80*1e3-realPosition)/self.c
            #should be able to get rid of at least one of the three above lines
            retardedTime = time+timeDelay
            a = .5*aReal/self.c**2
            
            bQuad = [1-aReal/self.c**2*(zMeas*1e3) for zMeas in self.observeKm] 
           
            cQuad = [.5*aReal/self.c**2*(zMeas*1e3)**2-(self.dz*1e3+.5*aReal*time**2)\
                     for zMeas in self.observeKm]
           
           
           
            retardedPos1 = [(.5*(-b+math.sqrt(b**2-4*a*c)))/a for b,c in zip(\
                bQuad,cQuad)]
            retardedPos = [zMeas*1e3 - r for zMeas, r in zip(self.observeKm, retardedPos1)]
                       
          
            retardedP+= [retardedPos1[-2]]#?
            realP+=[realPosition]
            t+=[time]
            
            #set up/plot the first graph:
            #realPosition vs time and retarded position bs time
            plt.subplot(2,2,1)
            plt.scatter(t,realP,marker ='*')
            plt.scatter(t,retardedP,marker='o')
            
            #set up/ plot the second graph:
            #plot the electrostatic field of all three charges and their 
            #image charges
            '''
            this graph plots the electric field of the accelerating charge in the cloud
            using the method of images. Images charges for the accelerating charge and 
            for the other two charges are included. It calculates the titol electrostatic
            based on the positions of the charges at each ms in time completely neglecting
            velocity
            '''
            ESecond =  [(self.electricFieldOne(self.Q,zMeas*1e3,r) + other.electricFieldOne\
                (other.Q,zMeas*1e3,other.dz*1e3)+charge.electricFieldOne(charge.Q,zMeas*1e3,\
                charge.dz*1e3))/(4*math.pi*scipy.constants.epsilon_0) \
                        for zMeas,r in zip(self.observeKm,retardedPos1)]  
            plt.subplot(2,2,2)
            plt.plot(ESecond,self.observeKm, linewidth=1)
            
            #set up/plot the third graph:
            #the velocity field (Feynman's second term) on the single accelerating
            #charge alone with no boundary conditions/ image charges
            #double check this
            realVelocity = aReal*time
            EThird = [(self.easy*2*(-aReal*time)*(b**2-4*a*c)**\
                        (-1.5))/(self.c*(retardedPosition)**2)\
                       for b,c,retardedPosition in\
                       zip(bQuad, cQuad,retardedPos)]
            EThird = np.array(EThird)
            observeKm = np.array(self.observeKm)
            plt.subplot(2,2,3)
            plt.plot(EThird,observeKm,linewidth=1)
            
            '''
            set up/ graph the fourth graph:
            the velocity electric field (Feynman's second term) on a single 
            accelerating charge but now with the velocity field of the image
            charges added in
            ''' 
            EFourth = [(k*(-self.easy) *(-aReal * time  )*\
            (2.0 / scipy.constants.c))/2 for k in self.electricFieldTwo\
                    (self.observeKm,retardedPos1,a,bQuad,cQuad,aReal,time)]
               
            plt.subplot(2,2,4)
            plt.plot(EFourth,self.observeKm,linewidth=1)
            
            time+=.0001
            
        time_labels = [f'$t$ = {i * 0.0001:.4f} s' for i in range(2,11)] 
        plt.subplot(2,2,1)
        plt.title('position of one accelerating charge with constant velocity'\
                  ,fontsize=8)
        plt.xlabel('time in s',fontsize=5)
        plt.ylabel('position in m',fontsize = 5)
        plt.subplot(2,2,2).legend(time_labels, loc='upper left', fontsize=5)
        plt.title('delayed coulomb field, 1 accelerating charge 2 static charges\
        with bounds', fontsize= 8)
        plt.ylabel('km',fontsize=5)
        plt.xlabel('E[Vm^-1]',fontsize = 5)
        plt.subplot(2,2,3).legend(time_labels, loc='upper left', fontsize=5)
        plt.title('delayed velocity elecric field, 1 accelerating charge,\
        no bounds', fontsize = 8)
        plt.ylabel('km',fontsize=5)
        plt.xlabel('E[v*m^-1]',fontsize=5)        
        plt.subplot(2,2,4).legend(time_labels, loc='upper left', fontsize=5)
        plt.title('delayed velocity electric field, one accelerating charge,\
        with boundaries', fontsize= 8)
        plt.ylabel('km',fontsize=5)
        plt.xlabel('E[v*m^-1]',fontsize=5)
        plt.savefig('spriteelectricfieldgraph.png')
        plt.show()
        
            
                   
            
            
#i need to make code look nicer/ orginize better and add comments            
r = Feynman(10,200,90,1)
t = Feynman(5,-200, 90,1)
j = Feynman(.5*(10-5)+10,-3/8*200,90,1)
positions= \
[10000,10000,9900,9600,9100,8400,7500,6400,5099.999,3599.99,1899.9999,-3.63797*1e-12]
#r.electricField(t,j,positions)
r.secondGraph(t,j,positions)
plt.show()
        
