#!/usr/bin/env python
# coding: utf-8


import utils2
import utils1
import matplotlib.pyplot as plt
from feynman import Diagram
import time
from math import pi

utils1

utils1 

class Process:
    """
    This class reads the input interaction and creates the diagrams.
    """
        
    def __init__(self, process):
        """
        Initialise the class with the scattering process. Here "process" is a string.
        """
        
        # getting particle data informations
        self.pdb = utils2.ParticleDatabase()
        self.stringcp = process
        
        
        # list of standardmodel particle ids
        standardmodel = [1,2,3,4,5,6, 11,12,13,14,15,16,18, 22,  21, 23, 24, 25]
        # extend with the antiparticles
        a = [-item for item in standardmodel]
        standardmodel.extend(a)
        
        # list of prt id who are not valid inputs. in this case 
        invalid = [21, 23, 24, 25]
        
        # list of flavors
        iflavors = {1:0, 2:0, 3:0, 4:0, 5:0, 6:0, 11:1, 12:1, 13:1, 14:1, 15:1, 16:1, 18:1, 22:2}
        oflavors = {1:0, 2:0, 3:0, 4:0, 5:0, 6:0, 11:1, 12:1, 13:1, 14:1, 15:1, 16:1, 18:1, 22:2}
        
        self.valid = False
        self.prts, process = [None]*4, process.split()

        # check if the input string has 5 elements
        if len(process) != 5: print("Wrong input Syntax"); return
        for key, prt in zip(process[0:2] + process[3:5], self.prts):
            # check if the particel is in the particle data base
            if not key in self.pdb: print("Particle is not in the ParticleDatabase"); return
            else:
                # Check if particle is a standardmodel prt. 
                if not self.pdb[key].pid in standardmodel:
                    print("invalid input." + key + " is not a standardmodel particle");return
                # Check if W, Z gluons, higgs.
                if self.pdb[key].pid in invalid:
                    print("invalid particle. Gluons, W Bosons, Z Bosons and higgs are excluded"); return 
            
        iflavors[abs(self.pdb[self.stringcp.split()[0]].pid)] += -1 if self.pdb[self.stringcp.split()[0]].pid<0 else +1
        iflavors[abs(self.pdb[self.stringcp.split()[1]].pid)] += -1 if self.pdb[self.stringcp.split()[1]].pid<0 else +1
        for key in iflavors:
            if iflavors[key] != oflavors[key]: print("flavor not conserved")
 
            prt = key

        self.valid = True                
        
        self.input_order_prt1 = self.pdb[self.stringcp.split()[0]]
        self.input_order_prt2 = self.pdb[self.stringcp.split()[1]]
        self.input_order_prt3 = self.pdb[self.stringcp.split()[3]]
        self.input_order_prt4 = self.pdb[self.stringcp.split()[4]]
        
        # ordering the particle. In1 will be the particle and In2 the antiparticle
        #incoming particles
        #if list_incoming[0] >= list_incoming[1]:
        if self.input_order_prt1.pid >= self.input_order_prt2.pid:
            self.In1 = self.input_order_prt1
            self.In2 = self.input_order_prt2
        else:
            self.In1 = self.input_order_prt2
            self.In2 = self.input_order_prt1

        #outgoing particles
        if self.input_order_prt3.pid >= self.input_order_prt4.pid:
            self.Out1 = self.input_order_prt3
            self.Out2 = self.input_order_prt4
        else:
            self.Out1 = self.input_order_prt4
            self.Out2 = self.input_order_prt3


    ########################
    # Checks if the order changed
    ########################
    def order_Incoming(self):
        """
        Checks if the order of the input of the incoming particles have changed.
        """
        if self.input_order_prt1.pid >= self.input_order_prt2.pid:
            return True
        else:
            return False
        
    def order_Outgoing(self):
        """
        Checks if the order of the input of the outgoing particles have changed.
        """
        if self.input_order_prt3.pid >= self.input_order_prt4.pid:
            return True
        else:
            return False
            
    ################
    # Check if the propagator is photon or eleectron
    ################
    def prop_is_photon(self):
        """
        Checks if the propagator is a photon.
        """
        if self.In1.pid != 22 and self.In2.pid != 22 and self.Out1.pid != 22 and self.Out2.pid != 22 :
            return True
        else:
            return False
    
    def prop_is_fermion(self):
        """
        Checks if the propagator is a fermion.
        """
        if self.In1.pid ==22 or self.Out1.pid == 22:
            return True
        else:
            return False
        
    #################
    # checks if the particle is a quark
    #################
    def prt_is_quark(self):
        """
        Method that ckecks if the particle is a quark
        """
        quark_id = [1,2,3,4,5,6]
        a = [-item for item in quark_id ]
        quark_id .extend(a)
        if self.Out1.pid in quark_id:
            return True
        else:
            False
                
    ##############################
    # s -chanel diagram
    ##############################
    def s_channel(self):
        """
        Returns the Feynman diagram for the s channel.
        """
        In1 = self.In1.name
        In2 = self.In2.name
        Out1 = self.Out1.name
        Out2 = self.Out2.name
        
        fig = plt.figure(figsize=(10.,10.))
        ax = fig.add_axes([0,0,1,1], frameon=False)

        diagram = Diagram(ax)

        ############################
        # Defining the vertices
        ############################

        # In1_vtx and In2_vtx are the vertex of at the left end of the in coming lines
        # we need to define them, so that we can draw the incoming lines
        # Used marker = '' to not show the vertex
        in1_vtx = diagram.vertex(xy=(.1,.75), marker='')
        in2_vtx= diagram.vertex(xy=(.1,.25), marker='')
        # v1 is the vertex on the left of the propagator and v2 is the vertex on the right of the propagator
        v1 = diagram.vertex(xy=(.35,.5))
        v2 = diagram.vertex(xy=(.65,.5))

        # By defining the vertix we can define where the lines starts and ends. 
        # for example to draw the first incoming line: line goes from vertix in1 to v1, style defines the form of the line
                               # second incoming line: line goes from vertix in2 to v1, style defines the form of the line



        # out1_vtx is the vertex at the right end of the "higgs line"
        out1_vtx = diagram.vertex(xy=(.9,.75), marker = '')
        # out2_vtx is the vertex at the right end  of the outgoing W/Z Boson line
        out2_vtx = diagram.vertex(xy=(.9,.25),marker='')

        ###########################
        # Defining the lines
        ###########################
        
        # incoming lines
        # in1_ln is the line which goes from vertex in1_vtx to vertex v1 (incoming particle one)
        if self.In1.pid != 22:
            in1_ln = diagram.line(in1_vtx, v1)
        else:
            in1_ln = diagram.line(in1_vtx, v1, style='wiggly')
        # in2_ln is the line which goes from vertex v1 to vertex in2_vtx (incoming particle two)
        if self.In2.pid != 22:
            in2_ln = diagram.line(v1, in2_vtx)
        else:
            in2_ln = diagram.line(v1, in2_vtx, style='wiggly')

        # "v1", "v2" in the bracket tells us that the line goes from vertex one to vertex two
        # "style" is the form of the line
        if self.prop_is_photon():
            prop = diagram.line(v1, v2, style='wiggly')
        else:
            prop = diagram.line(v1, v2) # prop is a fermion
         
        # out going lines
        # out1_ln is the line between vertex  v2 and vertex out1_vtx (outgoing Fermion 1) 
        if self.Out1.pid != 22:
            out1_ln = diagram.line(v2, out1_vtx)
        else:
            out1_ln = diagram.line(v2, out1_vtx, style='wiggly')
        # out2_ln is the  line between vertex  v2 and vertex out2_vtx (outgoing Fermion 2) 
        if self.Out2.pid != 22:
            out2_ln = diagram.line(out2_vtx, v2)
        else:
            out2_ln = diagram.line(out2_vtx, v2, style='wiggly')

        #########################
        # Naming the lines
        #########################

        #incoming lines
        in1_ln.text("%r" %In1,fontsize=30)
        in2_ln.text("%r" %In2,fontsize=30)
        #wz1.text("$Z/W^\pm$",fontsize=30)
        #wz2.text("$Z/W^\pm$",fontsize=30)
        
        # propagator
        if self.prop_is_photon():
            diagram.text(0.5,0.55,"$\gamma$",fontsize=30)
        else:
            diagram.text(0.5,0.55,"$f$",fontsize=30)
            
        # outgoing lines
        diagram.text(0.69,0.35,"%r" %Out2,fontsize=30)
        out1_ln.text("%r" %Out1,fontsize=30)

        diagram.plot()
        plt.show()
        return diagram
    
    ##################################
    # t-channel diagram
    ##################################    
    def t_channel(self):
        """
        Returns the Feynman diagram for the t-channel.
        """

        In1 = self.In1.name
        In2 = self.In2.name
        Out1 = self.Out1.name
        Out2 = self.Out2.name

        fig = plt.figure(figsize=(10.,10.))
        ax = fig.add_axes([0,0,1,1], frameon=False)

        diagram = Diagram(ax)
        # In1_vtx and In2_vtx are the vertex of at the left end of the in coming lines
        in1_vtx = diagram.vertex(xy=(.1,.65), marker='')
        in2_vtx= diagram.vertex(xy=(.1,.15), marker='')
        # v1 is the vertex on the left of the propagator and v2 is the vertex on the right of the propagator
        v1 = diagram.vertex(xy=(.35,.5))
        v2 = diagram.vertex(xy=(.35,.3))
        # out1_vtx is the vertex at the right end of the "higgs line"
        out1_vtx = diagram.vertex(xy=(.6,.65), marker = '')
        # out2_vtx is the vertex at the right end  of the outgoing W/Z Boson line
        out2_vtx = diagram.vertex(xy=(.6,.15),marker='')
        
        ###########################
        # Defining the lines
        ###########################
        # Incoming lines
        # in1_ln is the line which goes from vertex in1_vtx to vertex v1 (incoming particle one)
        if self.In1.pid != 22:
            in1_ln = diagram.line(in1_vtx, v1)
        else:
            in1_ln = diagram.line(in1_vtx, v1, style='wiggly')
        # in2_ln is the line which goes from vertex v1 to vertex in2_vtx (incoming particle two)
        if self.In2.pid != 22:
            in2_ln = diagram.line(v2, in2_vtx)
        else:
            in2_ln = diagram.line(v2, in2_vtx, style='wiggly')

        # "v1", "v2" in the bracket tells us that the line goes from vertex one to vertex two
        # "style" is the form of the line
        if self.prop_is_photon():
            prop = diagram.line(v1, v2, style='wiggly')
        else:
            prop = diagram.line(v1, v2)
        
        # Outgoing particles
        # out1_ln is the line between vertex  v2 and vertex out1_vtx (outgoing Fermion 1) 
        if self.Out1.pid != 22:
            out1_ln = diagram.line(v1, out1_vtx)
        else:
            out1_ln = diagram.line(v1, out1_vtx, style='wiggly')
        # out2_ln is the  line between vertex  v2 and vertex out2_vtx (outgoing Fermion 2) 
        if self.Out2.pid != 22:
            out2_ln = diagram.line(out2_vtx, v2)
        else:
            out2_ln = diagram.line(out2_vtx, v2, style='wiggly')
        
        #########################
        # Naming the lines
        #########################
        #incoming lines
        in1_ln.text("%r" %In1,fontsize=30)
        in2_ln.text("%r" %In2,fontsize=30)
        # propagator
        if self.prop_is_photon():
            prop.text("$\gamma$",fontsize=30)
        else:
            prop.text("f",fontsize=30)
        
        # outgoing lines
        out2_ln.text("%r" %Out2,fontsize=30)
        #diagram.text(0.69,0.35,"out2",fontsize=30)
        out1_ln.text("%r" %Out1,fontsize=30)

        diagram.plot()
        plt.show()
        return diagram

    ##############################
    # u-chanel diagram
    ##############################    
    def u_channel(self):
        """
        Returns the Feynman diagram for the u channel.
        """
        In1 = self.In1.name
        In2 = self.In2.name
        Out1 = self.Out1.name
        Out2 = self.Out2.name

        fig = plt.figure(figsize=(10.,10.))
        ax = fig.add_axes([0,0,1,1], frameon=False)

        diagram = Diagram(ax)
        # In1_vtx and In2_vtx are the vertex of at the left end of the in coming lines
        in1_vtx = diagram.vertex(xy=(.1,.65), marker='')
        in2_vtx= diagram.vertex(xy=(.1,.15), marker='')
        # v1 is the vertex on the left of the propagator and v2 is the vertex on the right of the propagator
        v1 = diagram.vertex(xy=(.35,.5))
        v2 = diagram.vertex(xy=(.35,.3))
        # out1_vtx is the vertex at the right end of the "higgs line"
        out1_vtx = diagram.vertex(xy=(.7,.65), marker = '')
        # out2_vtx is the vertex at the right end  of the outgoing W/Z Boson line
        out2_vtx = diagram.vertex(xy=(.7,.15),marker='')

        ###########################
        # Defining the lines
        ###########################
        # Incoming lines
        # in1_ln is the line which goes from vertex in1_vtx to vertex v1 (incoming particle one)
        if self.In1.pid !=22:
            in1_ln = diagram.line(in1_vtx, v1)
        else:
            in1_ln = diagram.line(in1_vtx, v1, style='wiggly')
        # in2_ln is the line which goes from vertex v1 to vertex in2_vtx (incoming particle two)
        if self.In2.pid !=22:
            in2_ln = diagram.line(v2, in2_vtx)
        else:
            in2_ln = diagram.line(v2, in2_vtx, style='wiggly')

        # "v1", "v2" in the bracket tells us that the line goes from vertex one to vertex two
        # "style" is the form of the line
        if self.prop_is_photon():
            prop = diagram.line(v1, v2, style='wiggly')
        else:
            prop = diagram.line(v1, v2)
        
        # outgoing Particles
        # out1_ln is the line between vertex  v2 and vertex out1_vtx (outgoing Fermion 1) 
        if self.Out1.pid !=22:
            out1_ln = diagram.line(v2, out1_vtx)
        else:
            out1_ln = diagram.line(v2, out1_vtx, style='wiggly')
        # out2_ln is the  line between vertex  v2 and vertex out2_vtx (outgoing Fermion 2) 
        if self.Out2.pid !=22:
            out2_ln = diagram.line(out2_vtx, v1)
        else:
            out2_ln = diagram.line(out2_vtx, v1, style='wiggly')
        
        #########################
        # Naming the lines
        #########################
        #incoming lines
        in1_ln.text("%r" %In1,fontsize=30)
        in2_ln.text("%r" %In2,fontsize=30)
        
        # propagator
        if self.prop_is_photon():
            prop.text("$\gamma$",fontsize=30)
        else:
            prop.text("f",fontsize=30)

        # outgoing lines
        # out2
        #out2_ln.text(r"$\bar{\mathrm{out2}}$",fontsize=30)
        diagram.text(0.69,0.25,"%r" %Out2,fontsize=30)
        # out1
        #out1_ln.text("out1",fontsize=30)
        diagram.text(0.69,0.55,"%r" %Out1,fontsize=30)

        diagram.plot()
        plt.show()
        return diagram
        
    ##########################
    # Checking if the diagrams exist.
    #########################  
    def t_channel_exist(self):
        """
        Checks if the t channel diagram exists.
        """
        if self.In1.pid == self.Out1.pid and self.In2.pid == self.Out2.pid:
            return True
        # photon case
        #compton scattering
        elif self.In1.pid ==22 and self.In2.pid !=22:
            return True
        # pair production
        elif self.In1.pid ==22 and self.In2.pid ==22:
            return True
        # pair annihilation
        elif self.Out1.pid ==22 and self.Out1.pid ==22:
            return True
        else:
            return False
        
    def u_channel_exist(self):
        """
        Checks if the u channel diagram exist.
        """
        if self.t_channel_exist() and self.Out1.pid == self.Out2.pid:
            return True
        else:
            return False

    def s_channel_exist(self):
        """
        Checks if the t channel diagram exists.
        """

        if self.In1.pid == -self.In2.pid and self.Out1.pid == -self.Out2.pid:
            return True
        #compton scattering
        elif self.In1.pid ==22 and self.In2.pid !=22:
            return True
        else:
            return False
        
    def signs_su(self):
        """
        Checks if we have to consider a relative minus sign for the matrix elements of the u and s channel.
        """
        
        In1 = self.In1
        In2 = self.In2
        Out1 =self.Out1
        Out2 = self.Out2
        
        if self.s_channel_exist() and self.u_channel_exist():
            self.In1 = Out2
            self.Out2 = In1
            if self.u_channel_exist():
                return True
#        elif self.s_channel_exist() and self.u_channel_exist():
            self.In1 = In2
            self.In2 = In1
            if self.u_channel_exist():
                return True
#        elif self.s_channel_exist() and self.u_channel_exist():
            self.Out1 = Out2
            self.Out2 = Out1
            if self.u_channel_exist():
                return True
            return False
        else:
            return False
                
    def signs_st(self):
        """
        Checks if we have to consider a relative minus sign for the matrix elements of the s and t channel.
        """
        
        In1 = self.In1
        In2 = self.In2
        Out1 = self.Out1
        Out2 = self.Out2
        
        if self.s_channel_exist() and self.t_channel_exist():
            self.In1 = Out2
            self.Out2 = In1
            if self.s_channel_exist():
                return True
            self.In2 = In1
            self.Out2 = Out2
            if self.t_channel_exist():
                return True
            self.Out1 = Out2
            self.Out2 = Out1
            if self.t_channel_exist():
                return True
            return False
        else:
            return False
        
    def signs_tu(self):
        """
        Checks if we have to consider a relative minus sign for the matrix elements of the u and t channel.
        """
        
        In1 = self.In1
        In2 = self.In2
        Out1 =self.Out1
        Out2 = self.Out2
        
        if self.u_channel_exist() and self.t_channel_exist():
            self.In1 = Out2
            self.Out2 = In1
            if self.t_channel_exist():
                return True
            self.In1 = In2
            self.In2 = In1
            if self.t_channel_exist():
                return True
            self.Out1 = Out2
            self.Out2 = Out1
            if self.t_channel_exist():
                return True
            return False
        else:
            return False
                
    def diagram(self):
        """
        Returns all possible feynman diagrams.
        """
        
        # first check if u, s and t chanel diagrams exists         
        list_diagrams = []
        if self.t_channel_exist():
            list_diagrams.append(self.t_channel())
        if self.u_channel_exist():
            list_diagrams.append(self.u_channel())
        if self.s_channel_exist():
            list_diagrams.append(self.s_channel())
        
        return list_diagrams #, end_time - start_time
            


class Scattering():
    """
    This class defines the differential cross section function, needed to calculate
    the integrated cross section for a given process
    """
    def __init__(self, Process, pp1, pp2, h, frame = "lab"):
        """
        Takes as inpute "Process", which is a process initialised by the Process class, the four momentum
        pp1 and pp2 of the two incoming particles and the helicity configuration h
        """
        self.frame = frame
        self.process = Process
        self.In1 = Process.In2
        self.In2 = Process.In1
        self.Out1 = Process.Out2
        self.Out2 = Process.Out1
        
        
        # changing order of p1 and p2 in the case we changed the order of the incoming particles internally
        if Process.order_Incoming():
            self.pp1 = pp2
            self.pp2 = pp1
        else:
            self.pp1 = pp1
            self.pp2 = pp2
            
        # changing order of the helicity in the case we changed the order of the particles internally
        self.h = [None,None,None,None]
        if Process.order_Incoming():
            self.h[0] = h[1]
            self.h[1] = h[0]
        else:
            self.h[0] = h[0]
            self.h[1] = h[1]
        
        if Process.order_Outgoing():
            self.h[2] = h[3]
            self.h[3] = h[2]
        else:
            self.h[2] = h[2]
            self.h[3] = h[3]
            
        # exchange the non intitialized helicities by [-1,1]
        for x, i in enumerate(self.h):
            if i ==None:
                self.h[x] = [-1,1]
            else: #self.h[x] = [self.h[x]]
                try: self.h[x][0] 
                except: self.h[x] = [self.h[x]]     

        #gamma matrices
        self.dmu = utils2.DiracMatrices()
        self.dml = ~self.dmu        
        
        #################################
        # preparation for the scattering corss section
        #################################

        # creationg the boost matrix bm 
        self.bm = utils2.BoostMatrix(self.pp1 + self.pp2)

        # creating the boost back matrix bbm
        self.bbm = utils2.BoostMatrix(~self.pp1 + ~self.pp2)
        
        # setting outgoing momenta to zero
        self.pp3 = utils2.FourVector(0,0,0,0)
        self.pp4 = utils2.FourVector(0,0,0,0)
        
        # calculating the CM momenta
        from math import sqrt

        self.pp1_cm = self.bm*self.pp1
        self.pp2_cm = self.bm*self.pp2
        self.q = sqrt(self.pp1_cm[0]**2 - 1/2*(self.Out1.mass**2 + self.Out2.mass**2))
    
                
        # Creating the incoming particles
        if self.frame == "lab":
            # bringing particles in the lab frame where prt2 is at rest (fixed target experiment)
         #   bm0 = utils2.BoostMatrix(self.pp2)
         #   pp10= bm0*self.pp1
         #   pp20= bm0*self.pp2
            self.p1 = utils2.Particle(self.In1, self.pp1, 1)
            self.p2 = utils2.Particle(self.In2, self.pp2, 1)  
         #   self.p1 = utils2.Particle(self.In1, pp10, 1)
         #   self.p2 = utils2.Particle(self.In2, pp20, 1)     
        else:
            self.p1 = utils2.Particle(self.In1, self.pp1_cm, 1)
            self.p2 = utils2.Particle(self.In2, self.pp2_cm, 1)
        
        ###################################
        # Calculating factors
        ###################################

        ## Calculate the cross-section prefactor  
        # lab frame
        # (hbar/(8 pi))^2
        # hbar =6.5821 × e-25 GeV s
        # c = 2.99792458 × e8 m s-1
        if self.frame == "lab":
            self.xspre = (6.5821e-25/(8*pi))**2
        # CM frame
        # ((hbar c)/(8 pi))^2 in units m^2 GeV^2.
        else:
            self.xspre = (1.97326979e-16/(8*pi))**2
        
        # Calculate the matrix-element prefactor (-4 pi alpha).
        self.mepre = -4*pi/137.0
        
    def me(self):
        """
        Calculates the matrix element for each possible diagram. 
        Returns the sum of the matrix elements.
        """
  
        c = 2.99792458e8 # speed of light
        # define the identity matrix
        I = utils2.Matrix([1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1])
        ##########################################
        # The case when the propagator is a photon
        ##########################################
        if self.process.prop_is_photon():
           # print("start me() method")
            
            ##############
            # t channel
            #############
            if self.process.t_channel_exist():
            #    print("t channel")
                p0 = self.p1.p - self.p3.p
                Mt = self.mepre/p0**2*sum([
                    (self.p3.wbar()*self.dmu[mu]*self.p1.w())*
                    (self.p4.wbar()*self.dml[mu]*self.p2.w())
                    for mu in range(0, 4)])

            else:
                Mt =0

            ##############
            # u channel
            ##############
            if self.process.u_channel_exist():
                p0 = self.p1.p - self.p4.p
                Mu = self.mepre/p0**2*sum([
                    (self.p4.wbar()*self.dmu[mu]*self.p1.w())*
                    (self.p3.wbar()*self.dml[mu]*self.p2.w())
                    for mu in range(0, 4)])

            else:
                Mu = 0

            ##############
            # s cannel
            ##############
            if self.process.s_channel_exist():
                p0 = self.p1.p + self.p2.p
                Ms = self.mepre/p0**2*sum([
                    (self.p3.wbar()*self.dmu[mu]*self.p4.w())*
                    (self.p2.wbar()*self.dml[mu]*self.p1.w())
                    for mu in range(0, 4)])
            else:
                Ms = 0
            #checks if the matrix elements have to be added or subtracted
            if self.process.signs_su():
                return Ms - Mu
            else:
                return Ms +Mu
            if self.process.signs_st():
                return Ms - Mt
            else:
                return Ms + Mt
            if self.process.signs_tu():
                return Mt - Mu  
            else:
                return Mt + Mu
            

        ###################################################################
        # in case the propagator is a fermion and one incoming prt is a photon
        ###################################################################
        elif self.process.prop_is_fermion() and self.In2.pid != 22 and self.Out2.pid != 22:
            
            #############
            # t channel
            ############
            if self.process.t_channel_exist():
                p0 = self.p1.p - self.p3.p
                m = self.In2.mass
                p1_slash =  (self.dmu[0]*self.p1.p[0] - self.dmu[1]*self.p1.p[1] - self.dmu[2]*self.p1.p[2] - self.dmu[3]*self.p1.p[3])
                p3_slash =  (self.dmu[0]*self.p3.p[0] - self.dmu[1]*self.p3.p[1] - self.dmu[2]*self.p3.p[2] - self.dmu[3]*self.p3.p[3])
                
                ep2_slash = (self.dmu[0]*self.p2.ep[0] - self.dmu[1]*self.p2.ep[1]- self.dmu[2]*self.p2.ep[2] - self.dmu[3]*self.p2.ep[3])
                ep3_slash_con = (self.dmu[0]*self.p2.ep[0].conjugate() - self.dmu[1]*self.p2.ep[1].conjugate()- self.dmu[2]*self.p2.ep[2].conjugate() - self.dmu[3]*self.p2.ep[3].conjugate())
                Mt = self.mepre/(p0**2 - m**2)*sum([
                    (self.p4.wbar()*self.dmu[mu]*ep2_slash)* # where ep is the polarization vector
                    ( p1_slash - p3_slash + m*I)*
                    (ep3_slash_con*self.dml[mu]*self.p1.w())
                    for mu in range(0, 4)])
            else:
                Mt = 0

            ###########
            # u channel
            ###########
            # there is no u chaennel diagram
            
            
            ###########
            # s channel
            ###########
            if self.process.s_channel_exist():
                p0 = self.p1.p + self.p2.p
                m = self.prt2.data.mass
                p1_slash =  (self.dmu[0]*self.p1.p[0] - self.dmu[1]*self.p1.p[1] - self.dmu[2]*self.p1.p[2] - self.dmu[3]*self.p1.p[3])
                p2_slash =  (self.dmu[0]*self.p2.p[0] - self.dmu[1]*self.p2.p[1] - self.dmu[2]*self.p2.p[2] - self.dmu[3]*self.p2.p[3])
                
                ep2_slash = (self.dmu[0]*self.p2.ep[0] - self.dmu[1]*self.p2.ep[1]- self.dmu[2]*self.p2.ep[2] - self.dmu[3]*self.p2.ep[3])
                ep3_slash_con = (self.dmu[0]*self.p2.ep[0].conjugate() - self.dmu[1]*self.p2.ep[1].conjugate()- self.dmu[2]*self.p2.ep[2].conjugate() - self.dmu[3]*self.p2.ep[3].conjugate())
 
                Ms = self.mepre/(p0**2 - m**2)*sum([
                    (self.p4.wbar()*self.dmu[mu]*(~self.p3.ep()))*
                    ( p1_slash + p2_slash + m*I)*
                    (self.p2.ep()*self.dml[mu]*self.p1.w())
                    for mu in range(0, 4)]) 
            else:
                Ms =0
            # checks the relative sign between the matrix elements
            if self.process.signs_st():
                return Ms - Mt
            else:
                return Ms + Mt



        #######################################################################################
        # in case the propagator is a fermion and the two outgoing prts are photons
        #######################################################################################
        elif self.process.prop_is_fermion() and self.Out1.pid == 22 and self.Out2.pid == 22:
            
            #############
            # t channel
            ############
            if self.process.t_channel_exist():
                p0 = self.p1.p - self.p3.p
                m = self.p2.data.mass
                p1_slash =  (self.dmu[0]*self.p1.p[0] - self.dmu[1]*self.p1.p[1] - self.dmu[2]*self.p1.p[2] - self.dmu[3]*self.p1.p[3])
                p3_slash =  (self.dmu[0]*self.p3.p[0] - self.dmu[1]*self.p3.p[1] - self.dmu[2]*self.p3.p[2] - self.dmu[3]*self.p3.p[3])
                
                ep3_slash = (self.dmu[0]*self.p3.ep()[0] - self.dmu[1]*self.p3.ep()[1]- self.dmu[2]*self.p3.ep()[2] - self.dmu[3]*self.p3.ep()[3])
                ep4_slash = (self.dmu[0]*self.p4.ep()[0] - self.dmu[1]*self.p4.ep()[1]- self.dmu[2]*self.p4.ep()[2] - self.dmu[3]*self.p4.ep()[3])
 
                Mt = self.mepre/(p0**2 - m**2*c**2)*sum([
                    (self.p2.wbar()*self.dmu[mu]*ep4_slash)*
                    ( p1_slash - p3_slash + m*I*c)*
                    (ep3_slash*self.dml[mu]*self.p1.w())
                    for mu in range(0, 4)])

            ###########
            # u channel
            ###########
            if self.process.u_channel_exist():
                p0 = self.p1.p - self.p4.p
                m = self.p1.data.mass
                p1_slash =  (self.dmu[0]*self.p1.p[0] - self.dmu[1]*self.p1.p[1] - self.dmu[2]*self.p1.p[2] - self.dmu[3]*self.p1.p[3])
                p4_slash =  (self.dmu[0]*self.p4.p[0] - self.dmu[1]*self.p4.p[1] - self.dmu[2]*self.p4.p[2] - self.dmu[3]*self.p4.p[3])

                ep3_slash = (self.dmu[0]*self.p3.ep()[0] - self.dmu[1]*self.p3.ep()[1]- self.dmu[2]*self.p3.ep()[2] - self.dmu[3]*self.p3.ep()[3])
                ep4_slash = (self.dmu[0]*self.p4.ep()[0] - self.dmu[1]*self.p4.ep()[1]- self.dmu[2]*self.p4.ep()[2] - self.dmu[3]*self.p4.ep()[3])
  
                Mu = self.mepre/(p0**2 - m**2*c**2)*sum([
                    (self.p2.wbar()*self.dmu[mu]*ep3_slash)*
                    ( p1_slash - p4_slash + m*I*c)*
                    (ep4_slash*self.dml[mu]*self.p1.w())
                    for mu in range(0, 4)])
                
           
            
            ###########
            # s channel
            ###########
            # no s channel exists
              
            #amplitudes have to be added here
            return Mt + Mu 
        
        #############################################
        # In case the propagator is a fermion and two incoming particles are photon
        #########################################
        elif self.process.prop_is_fermion() and self.In1.pid == 22 and self.In2.pid == 22:

            #############
            # t channel
            ############
            if self.process.t_channel_exist():
                p0 = self.p1.p - self.p3.p
                m = self.p2.data.mass
                p1_slash =  (self.dmu[0]*self.p1.p[0] - self.dmu[1]*self.p1.p[1] - self.dmu[2]*self.p1.p[2] - self.dmu[3]*self.p1.p[3])
                p3_slash =  (self.dmu[0]*self.p3.p[0] - self.dmu[1]*self.p3.p[1] - self.dmu[2]*self.p3.p[2] - self.dmu[3]*self.p3.p[3])

                Mt = self.mepre/(p0**2 - m)*sum([
                    (self.p4.wbar()*self.dmu[mu]*self.p2.w())*
                    ( p1_slash - p3_slash + m*I)*
                    (self.p1.ep()*self.dml[mu]*self.p3.w())
                    for mu in range(0, 4)])
            
        ###########
        # u channel
        ###########
        #no u channek 

        
        ###########
        # s channel
        ###########
        # no s channel exists

            return Mt 
        
        
        
    
    def dxs(self, phi, theta):
        """
        Returns the differential scattering cross section. Takes as input phi and theta.
        """
        

        # for the factor S, check if the outgoing particles are identical
        if self.Out1.pid == self.Out2.pid:
            S= 1/2
        else:
            S=1
                    
 
        
        from math import sqrt, cos, sin
        ct = cos(theta)
        st = sin(theta)

        # Calculating p3 and p4 in the CM frame
        self.pp3[0] = self.pp1_cm[0]
        self.pp3[1] = self.q*st*cos(phi)
        self.pp3[2] = self.q*st*sin(phi)
        self.pp3[3] = self.q*ct
        self.pp4 = ~self.pp3
        
        # boosting p3 and p4 in the labframe
        if self.frame == "lab":
#            bm = utils2.BoostMatrix(self.pp1 +self.pp2)
            self.pp3 = self.bbm*self.pp3
            self.pp4 = self.bbm*self.pp4   
        # creating the outgoing particles
        self.p3 = utils2.Particle(self.Out1, self.pp3,1)
        self.p4 = utils2.Particle(self.Out2, self.pp4,1)
        
        list_M = []
        for i in self.h[0]:
            for j in self.h[1]:
                for k in self.h[2]:
                    for l in self.h[3]:
                        self.p1.h = i
                        self.p2.h = j
                        self.p3.h = k
                        self.p4.h = l
                        M = self.me()
       #                 print(M)
                        try: Me2 = M.real**2 + M.imag**2
                        except: Me2 = M**2
                        list_M.append(Me2)
       # if len(self.h[0])!=None and len(self.h[1])!=None: # case 1
        if self.h[0]!=None and self.h[1]!=None: # case 1
            me2 = sum(list_M)
        else: # case 2
            me2 = sum(list_M)/4

        # colour scale factors
        if self.process.prt_is_quark():
            Nc =3
        else: Nc =1
            
        # Quark charge
        uct = [2,4,6]
        dsb = [1,3,5]
        if self.p3.data.pid in uct:
            Q = 2/3
        elif self.p3.data.pid in dsb:
            Q = 2/3
        else:
            Q=1
        
        p1_abs = sqrt(sum([p1i**2 for p1i in self.p1.p[1:]]))
        p3_abs = sqrt(sum([p3i**2 for p3i in self.p3.p[1:]]))
        # cross section in the lab
        if self.frame == "lab":
            c = 2.99792458e8 
            dxs = (Q**2*Nc*self.xspre*me2*S*p3_abs**2 *st/ 
                    (self.p2.data.mass *p1_abs*(p3_abs *(self.p1.p[0] + self.p2.data.mass*c**2)
                     - p1_abs* self.p3.p[0]* ct) ))
        #CM cross section
        else:     
            dxs= Q**2*Nc*self.xspre*me2*S*p3_abs/p1_abs*st/(self.p1.p[0] + self.p2.p[0])**2
        
        
        return dxs #, end_time - start_time
        
        

