
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""MHOQ implementation with switching frequency minimisation

@author: Bikash Adhikari
@date: 14.11.2024
@license: BSD 3-Clause
"""

import numpy as np
from scipy import linalg , signal
import sys
import random
import gurobipy as gp
from gurobipy import GRB
import tqdm



class MHOQ_RLIM:
    def __init__(self, Nb, Qstep, QMODEL,  A, B, C, D):
        """
        Constructor for the Model Predictive Controller.
        :param Nb: Number of bits 
        :param Qstep: Quantizer step size / Least significant bit (LSB) 
        :param N_PRED: Prediction horizon | int 
        :param Xcs: Reference/Test signal 
        :param QL: Quantization levels 
        :param A, B, C, D: Matrices; state space representation of the reconstruction filter
        """
        self.Nb = Nb
        self.Qstep = abs(Qstep)
        self.QMODEL = QMODEL
        self.A = A
        self.B = B
        self.C = C
        self.D = D


    def state_prediction(self, st, con):
        """
        Predict the state for the given initial condition and control
        """
        x_iplus1 = self.A @ st + self.B * con
        return x_iplus1


    def get_codes(self, N_PRED, Xcs, YQns, MLns, Step):

        match self.QMODEL:
            case 1:
                QL = YQns.squeeze()
            case 2:
                QL = MLns.squeeze()

        # Storage container for code
        C_Store = []

        # Loop length
        len_MPC = Xcs.size - N_PRED

        # State dimension
        x_dim =  int(self.A.shape[0]) 

        # Initial state
        init_state = np.zeros(x_dim).reshape(-1,1)

        u_kminus1_ind = (np.floor(Xcs[0]/self.Qstep +1/2)).astype(int) 

        # MPC loop
        for j in tqdm.tqdm(range(len_MPC)):
        # for j in range(1, len_MPC):
            m = gp.Model("MPC- INL")
            u = m.addMVar((2**self.Nb, N_PRED), vtype=GRB.BINARY, name= "u") # control variable
            x = m.addMVar((x_dim*(N_PRED+1),1), vtype= GRB.CONTINUOUS, lb = -GRB.INFINITY, ub = GRB.INFINITY, name = "x")  # State varible 


            # Add objective function
            Obj = 0

            # Set initial constraint
            m.addConstr(x[0:x_dim,:] == init_state)
            for i in range(N_PRED):
                k = x_dim * i
                st = x[k:k+x_dim]
                bin_con =  QL.reshape(1,-1) @ u[:,i].reshape(-1,1) 
                con = bin_con - Xcs[j+i]

                # Switching set 
                # Objective update
                e_t = self.C @ x[k:k+x_dim] + self.D * con
                Obj = Obj + e_t * e_t 

                # Constraints update
                f_value = self.A @ st + self.B * con
                st_next = x[k+x_dim:k+2*x_dim]
                m.addConstr(st_next == f_value)

                #% Limit search space
                ub = (u_kminus1_ind +  Step)
                if ub >= 2**self.Nb-1:
                    ub = 2**self.Nb-1
                lb = (u_kminus1_ind -  Step)
                if lb <= 0:
                    lb = 0        
                U1 = u[lb:ub+1,i].reshape(-1,1) 
                # Add constraint that sum of variables should be 1
                consU1 = gp.quicksum(U1[:,i]) == 1
                m.addConstr(consU1)

                # extract other variables that are not in the given bound
                mask = np.ones(u.size, dtype=bool)
                mask[np.arange(lb,ub+1,1)] = False
                U2 = u[mask]
                
                # Binary varialble constraint
                consU2 = U2== np.zeros([U2.shape[0],1])
                m.addConstr(consU2)

            m.update
            # Set Gurobi objective
            m.setObjective(Obj, GRB.MINIMIZE)

            # 0 - Supress log output, 1- Print log outputs
            m.Params.LogToConsole = 0

            # Gurobi setting for precision  
            m.Params.IntFeasTol = 1e-9
            m.Params.IntegralityFocus = 1

            # Optimization 
            m.optimize()

            # Extract variable values 
            allvars = m.getVars()
            values = m.getAttr("X",allvars)
            values = np.array(values)

            # Variable dimension
            nr, nc = u.shape
            u_val = values[0:nr*nc]
            u_val = (np.reshape(u_val, (2**self.Nb, N_PRED))).astype(int)

            # Extract Code
            C_MPC = []
            for i in range(N_PRED):
                c1 = np.nonzero(u_val[:,i])[0][0]
                c1 = int(c1)
                C_MPC.append(c1)
            C_MPC = np.array(C_MPC)
            C_Store.append(C_MPC[0])

            U_opt = QL[C_MPC[0]] 
            # State prediction 
            con = U_opt - Xcs[j]
            x0_new = self.state_prediction(init_state, con)
            # State update for subsequent prediction horizon 
            init_state = x0_new
            u_kminus1_ind = C_MPC[0]
        C_MHOQ  = np.array(C_Store).reshape(1,-1)
        return C_MHOQ