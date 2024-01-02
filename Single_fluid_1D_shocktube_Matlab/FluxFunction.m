if limitedvariables == 1
    [w]  = Cons2Prim(q,gamma);
    [Wf] = limitersPrim(N,w,limiter);
elseif limitedvariables == 2
    [Qf] = limitersCons(N,q,limiter);
    [Wf]  = Cons2Prim(Qf,gamma);
elseif limitedvariables == 3
    [ri]  = Cons2Riem(q,gamma);
    [RIf] = limitersRiem(N,ri,limiter);
    [Wf]  = Riem2Prim(RIf,gamma);
end
if flux==1
    [F,Lambda] = ExactFlux(N,Wf,gamma,Lambda);
elseif flux==2
    [F,Lambda] = LinearOsherFlux(N,Wf,gamma,Lambda);
elseif flux==3
    [F,Lambda] = OsherFlux(N,Wf,gamma,Lambda);
elseif flux==4
    [F,Lambda] = HLLFlux(N,Wf,gamma,Lambda);
elseif flux==5
    [F,Lambda] = HLLCFlux(N,Wf,gamma,Lambda);
elseif flux==6
    [F,Lambda] = RoeFlux(N,Wf,gamma,Lambda);
elseif flux==7
    [F,Lambda] = PVRSFlux(N,Wf,gamma,Lambda);
elseif flux==8
    [F,Lambda] = OsherFlux_generalEOS(N,Wf,gamma,Lambda);
end
[F] = BoundaryConditions(F,Wf,gamma,N,BC_L,BC_R);
