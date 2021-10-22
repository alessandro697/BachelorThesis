#!/usr/bin/env python3

import sys
import os
import re
import os.path
import math
from rkf45 import r8_rkf45       # necessary: it provides RK integration
import numpy as np
# to install this lib, run as root: apt-get install python-numpy
from matplotlib.pyplot import *  # for plots: comment out if plots not needed
# to install this lib, run as root: apt-get install python-matplotlib

"""
"""	

def derivs(t, y): # here is all the physics: the left side of the equation
    neq=len(y)
    nhalf=neq//2
    deriv=[0.]*neq  # just to initialize the new array
    for i in range(nhalf):  # the second half of y is velocities
        deriv[i]=y[i+nhalf] # this enables the mapping of Newton to 1st order

    global derpotenz
    
    derpotenz=p*(6*pow(y[1]-y[0],5)-6*pow(y[1]-y[0],3)*(R3*R3+R2*R2)+6*(y[1]-y[0])*R3*R3*R2*R2)  #force pot.interazione
           
    deriv[2]=force-force_amplitude*math.sin(twopi*y[0])+derpotenz-gamma*deriv[0]

    deriv[3]=force-force_amplitude*math.sin(twopi*y[1])-derpotenz-gamma2*deriv[1]

    return deriv
   
def vcmcompute(tlist,xcmlist):
    nstep=len(tlist)
    vcmaveGood=(xcmlist[nstep-1]-xcmlist[nstep//2])/(tlist[nstep-1]-tlist[nstep//2])

    return vcmaveGood

def wholecalculation(npart,timefin,time_step,\
                     gamma,gamma2,pot_amplitude,forc,forc_at_end,R,R2,R3):
    global force, f_at_end
    force=forc
    f_at_end=forc_at_end
    if gamma == gamma2:
        label='__n_'+str(npart)+'__gam_'+str(gamma)+'__U0_'+str(pot_amplitude)+'__F_'+str(force)+'__R_'+str(R)+'__R2_'+str(R2)+'__R3_'+str(R3)+'__c_'+str(p)
    else:
        label='__n_'+str(npart)+'__gam1_'+str(gamma)+'__gam2_'+str(gamma2)+'__U0_'+str(pot_amplitude)+'__F_'+str(force)+'__R_'+str(R)+'__R2_'+str(R2)+'__R3_'+str(R3)+'__c_'+str(p)
    
    
    print("# starting a mol2 calculation for")
    print("#	"+label)
    filen="output_"+label+".dat"
    #   end of setting up stuff
    
#   spaced particles as initial condition:
    x0=(R*np.asarray(range(npart))).tolist() 
    v0=[0.]*npart	# all 0 speeds at start
    y=x0+v0
    if debug==1:
        print("initial condition:",x0,v0,y)
    neq=len(y)
    yp=derivs(0.,y)

    if debug==1:
        print("#starting point:",y)
        print("#starting deriv:",yp)

    relerr=1.e-6
    abserr=1.e-10
    flag=1

    tlist=[]
    xcmlist=[]
    vcmlist=[]
    fil = open(filen, 'w')
    o=""
    print(0, o.join((" "+str(i)) for i in y), file=fil)
    nstep=int(round(1.*timefin/time_step))
    for it in range(nstep):
        ti=it*time_step
        tf=(it+1)*time_step
        tlist.append(tf)
        y, yp, t, flag = r8_rkf45( derivs, neq, y, yp, ti, tf, relerr, abserr, flag )
        if flag!=2:
            print("Warning! flag =",flag,".... trying to keep on going")
            flag=2
        o=""
        print(tf, o.join((" "+str(i)) for i in y), file=fil)
        xcm=0
        vcm=0
        for i in range(npart):
            xcm+=y[i]
        xcm/=npart
        for i in range(npart,2*npart):
            vcm+=y[i]
        vcm/=npart
        print("integrated from",ti,"to",tf,"xcm=",xcm,"vcm=",vcm)
        xcmlist.append(xcm)
        vcmlist.append(vcm)

    endvcm=(abs(vcmlist[nstep-1])+abs(vcmlist[nstep-2])+abs(vcmlist[nstep-3]))/3.
    print("#end of the calculation, final |vel|=",endvcm)
    print("#energia stato eccitato e barriera di potenziale:",energy,barrier)
    fil.close()
    return endvcm,tlist,xcmlist,vcmlist

def bisection(f,otherarguments,a,b,tol):
    fa=f(a,otherarguments)
    fb=f(b,otherarguments)
    if(fb*fa>0):
        print("bisection: same sign at",a,":",fa,"and at",b,":",fb)
        exit()
    c = (a+b)/2.0
    fc=f(c,otherarguments)
    while (b-a)/2.0 > tol:
        if fc == 0:
            return c
        elif fa*fc < 0:
            b = c
            fb=fc
        else:
            a = c
            fa=fc
        c = (a+b)/2.0
        fc=f(c,otherarguments)
    return c

def zerofunc(x,otherarguments):
    [npart,timefin,time_step,gamma,gamma2,pot_amplitude,f_at_end,R,R2,R3]=otherarguments
    force=x
    endvcm,tlist,xcmlist,vcmlist = \
        wholecalculation(npart,timefin,time_step,\
                         gamma,gamma2,pot_amplitude,force,f_at_end,R,R2,R3)
    return endvcm-1.e-6

def mol(args):

    commandname=args.prog
    global debug
    debug = args.debug
    global gamma, gamma2,force_amplitude,R,R2,R3,p
    global pi,twopi,energy,barrier
    pi=math.pi
    twopi=2*pi

    
    force=args.force
    fin=args.fin
    ffin=args.ffin
    deltaf=args.deltaf
    f_at_end=args.f_at_end
    finddepforce=args.finddepforce
    gamma=args.gamma
    gamma2=args.gamma2
    if gamma2 == -1:
        gamma2 = gamma
        
    npart=args.npart
    doplot=args.doplot
    timefin=args.timefin
    time_step=args.time_step
    pot_amplitude=args.pot_amplitude
    R=args.R
    R2=args.R2
    R3=args.R3
    p=args.p
    
    
    force_amplitude=pot_amplitude*pi # derivative of Vext(x)=-U0/2 cos(2 pi x)
    energy=(p/2)*(3*pow(R2,2)*pow(R3,4)-pow(R3,6))
    barrier=(p/2.)*(pow(R3,6)-pow(R2,6)+3*pow(R3,2)*pow(R2,4)-3*pow(R3,4)*pow(R2,2))
    
      
    if finddepforce:
        tol=1.e-4  
        force_min=0
        force_max=force
        force_crit=bisection(zerofunc,[npart,timefin,time_step,\
                        gamma,gamma2,pot_amplitude,f_at_end,R,R2,R3],\
                        force_min,force_max,tol)
        print("# critical force =",force_crit)
        force=force_crit
    if deltaf==0: # single calculation:
        endvcm,tlist,xcmlist,vcmlist = wholecalculation(npart,\
          timefin,time_step,gamma,gamma2,pot_amplitude,force,f_at_end,R,R2,R3)
        print("#F= ",force," <vcm>=",vcmcompute(tlist,xcmlist))
    else:
        for force in np.linspace(fin,ffin,int(round((ffin-fin)/deltaf+1))):
            endvcm,tlist,xcmlist,vcmlist = wholecalculation(npart,\
              timefin,time_step,gamma,gamma2,pot_amplitude,force,f_at_end,R,R2,R3)
            print("F= ",force," <vcm>=",vcmcompute(tlist,xcmlist),"\n")

    sys.stdout.flush()
    if doplot:
# This graphics part below here is entirely optional.
# It may be worth commenting it out if you are unwilling to install matplotlib
#    figure( 1 )
        subplot( 2, 1, 1 )
        plot( tlist, xcmlist, 'b-o')
        xlabel( '$t$' )
        ylabel( '$x_{cm}$' )
        title( 'mol2 $n='+str(npart)+'$, $R = '+str(R)+'$, $R2 = '+str(R2)+'$, $R3 = '+str(R3)+'$, $c = '+str(p)+'$, $\gamma='+str(gamma)+'$, $F='+str(force)+'$')
#    legend( ( '$x_{cm}$' ), loc='upper left' )
        legend( ( "xcm" ), loc='upper left' )

#    figure( 2 )
        subplot( 2, 1, 2 )
        plot( tlist, vcmlist, 'b-o')
        xlabel( '$t$' )
        ylabel( '$v_{cm}$' )
        legend( ( '$v_{cm}$' ), loc='lower right' )

        show()




# the following function is only exectuted when the code is run as a script, and its purposes is to parse
# the command line and to generate a meaningful parsed args list to the actual function doing the job:
if __name__ == "__main__":
    import argparse

    desc="""simulate the molecular bistable model using accurate RKF integration

						
"""

    epil="""\t\t\t ultima modifica 09/06/2020"""

    ##  Argument Parser definition
    parser = argparse.ArgumentParser( formatter_class=argparse.ArgumentDefaultsHelpFormatter  #RawTextHelpFormatter
                                      , description=desc, epilog=epil)

    parser.add_argument( '-d', '--debug', action='store_true',
                         dest='debug', default=0,
                         help='activate debug mode -- WARNING! output is affected/spoiled!' )    
    parser.add_argument( '-f',
                         dest='force', type=float, default=0.0, action='store',
                         help='external driving force')
    parser.add_argument( '--Fin',
                         dest='fin', type=float, default=0.0, action='store',
                         help='initial driving force')
    parser.add_argument( '--Ffin',
                         dest='ffin', type=float, default=0.0, action='store',
                         help='final driving force')
    parser.add_argument( '--Fde',
                         dest='deltaf', type=float, default=0.0, action='store',
                         help='driving force increment')
    parser.add_argument( '--f-at-end-force',
                         dest='f_at_end', type=float, default=0.0, action='store',
                         help='external force applied only to the leftmost particle')
    parser.add_argument( '--find-depinning-force',
                         dest='finddepforce', action='store_true',
                         help='execute a bisection calculation of the depinning force')
    parser.add_argument( '-g',
                         dest='gamma', type=float, default=10., action='store',
                         help='viscous damping coefficient')
    parser.add_argument( '-g2',
                         dest='gamma2', type=float, default=-1, action='store',
                         help='viscous damping coefficient for the second particle')
    parser.add_argument( '-n',
                         dest='npart', type=int, default=2, action='store',
                         help='the number of particles')
    parser.add_argument( '--np',
                         dest='doplot', action='store_false',
                         help='generate no plot at the end')    
    parser.add_argument( '-t',
                         dest='timefin', type=float, default=3000.0, action='store',
                         help='the total simulation time')
    parser.add_argument( '--dt',
                         dest='time_step', type=float, default=0.1, action='store',
                         help='the time interval, used only for writing')
    parser.add_argument( '-U',
                         dest='pot_amplitude', type=float, default=1.0, action='store',
                         help='peak-peak potential amplitude U0 in -(U0/2)*cos(2pi x)')
    parser.add_argument( '--R',
                         dest='R', type=float, default=1., action='store',
                         help='spaziatura iniziale molecola, R=R3')
    parser.add_argument( '--R2',
                         dest='R2', type=float, default=1./(math.sqrt(3.)), action='store',
                         help='massimo')
    parser.add_argument( '--R3',
                         dest='R3', type=float, default=1., action='store',
                         help='secondo minimo')
    parser.add_argument( '-p',
                         dest='p', type=float, default=1., action='store',
                         help='prefattore')
    
    ## End arg parser definition
#    args = argparse.parser()
    args=parser.parse_args(sys.argv[1:])
    d = vars(args)	# adding prog in args, for unknown reasons it's not there...
    d['prog']=parser.prog
    
#    print('outside, in parser:',parser.prog)
#    print('outside, in parser:',args.prog)
#    print('outside, in parser:',args)
#   here the actual function doing the job is called:

    stored=mol(args)
# and the results are printed out:
#    for i in stored:
#        print(' '.join(map(str,i)))
