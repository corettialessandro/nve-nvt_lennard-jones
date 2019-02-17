def NH_mdlj(nstep=1000, rho=0.84, m=6, kt_targ=.684, tausq=.1, kt_init=.684, dt=0.005, freq=1, mode=0):
    """ Okinawa 3/04/2017 - LJMD  3D-version
    Molecular Dynamics program for N interacting LJ-particles
    in a 3 dimensional space using reduced LJ units (mass,sigma,epsilon):
    Dynamics: canonical - single nose-hoover algorithm;
    INPUT:
      initial configuration from file 'conf_in' (binary):
      - 1st row:  x_1 y_1 z_1 (initial coordinates of 1st particle)
      -  . . .
      - Nth row   x_N y_N z_N (initial coordinates of Nth particle)
      - N+1th row:  px_1 py_1 pz_1 (initial momenta of 1st particle)
      -  . . .
      - 2*Nth row   px_N py_N pz_N (initial momenta of Nth particle)
      if mode is 0 initial positions are sampled from FCC lattice and previous input/output files are removed
      if mode is not 0 initial positions are read from configuration file
      if mode is even initial velocities are extracted from Maxwell Distribution and previous output files are removed
    OUTPUT:
    - enep          total potential energy  E_p
    - enek          total kinetic energy E_k
    - vir_p         virial coefficient
    - enh           contribution to energy from nose-hoover thermostat
    - enet          total energy Ep+Ek (not constant for canonical dynamics)
    - enht          total energy plus nose-hoover contribution (constant of motion for N-H algorithm)
    - vcm(x,y,z)    center of mass momentum (constant of motion)
    - inst_T        instantaneous kinetic temperature
    - Delta_T       difference between inst and target temperatures
    PLOT:
    - RDF           Radial distribution function as a function of interparticle distance
    - Temperature   Instantaneous, mean and target temperatures vs. timestep
    """
    from numpy import random, sqrt, sum
    from MDLJ_lib import MDLJ_lib
#
    a=(4/rho)**(1./3.)
    Lref=a*m
    N=4*m**3
    md=MDLJ_lib(N, Lref)
    g=3*N-3;
    gdr_out='gdrmd.out'
    print( "# initial (kinetic) temperature kt = %8.4f" % kt_init )
    print( "# integration time step dt = %8.4f" % dt )
    print( "# number of time steps for this run ns = %d" % nstep )
    print( "# reference side of the cubic MD-box L(0) = %8.4f " % Lref )
    print( "# number of particles  N = %d" % N )
    print( "# density rho = %8.4f"  % rho )
    if mode :
        md.read_input(N, conf_in='conf_in.b')
        print("# initial configuration read from config file ")
    else :
        # initial positions mode=0 from FCC lattice
        md.clean_run()
        md.fcc(m, verbose=False)
        # initial velocities: realization of random gaussian process with zero mean
        # so that initial momentum is exactly zero at the initial time
        # this results in a stationary center of mass for the N particle system
    if mode%2 == 0:
        # if mode is even (re)samples initial velocities from gaussian
        md.clean_run(keepinput=True)
        md.t=0
        pstd=sqrt(md.mass*kt_init)
        md.px = random.normal(0., pstd, N)
        md.py = random.normal(0., pstd, N)
        md.pz = random.normal(0., pstd, N)
        vcmx  = sum(md.px)
        vcmy  = sum(md.py)
        vcmz  = sum(md.pz)
        md.px   -= vcmx/N
        md.py   -= vcmy/N
        md.pz   -= vcmz/N
        # reduced coordinates !
        md.px *= md.L
        md.py *= md.L
        md.pz *= md.L
        print("# initial momenta sampled from maxwellian at temperature %8.4f " % kt_init)
    print( "# starting with box side at initial time L(t=0) = %8.4f " % md.L )
    # ofstream eout("outmd.txt");
    tt, ekt, ept, vir, enht, ektsq, eptsq, etsq, enhtsq, enhttsq = md.nose_hoover(N, nstep, dt, kt_targ, tausq, freq)
    DEiE = sqrt(etsq/nstep - ((ekt + ept)/nstep)**2)/abs((ekt+ept)/nstep)
    DENHiENH = sqrt(enhttsq/nstep - ((ekt + ept + enht)/nstep)**2)/abs((ekt+ept+enht)/nstep)
    pres = (2*ekt + vir)/(3.*md.L**3)
    print( "# END OF THE RUN: Dumping average quantities")
    print( "# average potential energy  ep   = %10.5g " % (ept/nstep) )
    print( "# average kinetic   energy  ek   = %10.5g " % (ekt/nstep) )
    print( "# average N-H       energy  enh  = %10.5g " % (enht/nstep) )
    print( "# average total     energy  et   = %10.5g " % ((ekt+ept)/nstep) )
    print( "# average total N-H energy  enht = %10.5g " % ((ekt+ept+enht)/nstep) )
    print( "# average pressure          pres = %10.5g" %  (pres/nstep) )
    print( "# average temperature       kT   = %10.5g" %(2.*ekt/(nstep*g)) )
    print( "# average total     energy relative fluctuations: DE/E       = %.4g" % DEiE)
    print( "# average total N-H energy relative fluctuations: DE_NH/E_NH = %.4g" % DENHiENH)
    # gdr final printout - use actual value of density rho
    md.write_gdr( N, tt, rho, gdr_out )
#   end of md - visualization of G(r)
    from matplotlib.pyplot import plot, show
    from numpy import loadtxt,ones,size
    r,gdr = loadtxt( gdr_out, unpack=True, skiprows=1 )
    plot(r,gdr,'b.',r,ones(size(r)),'k-')
    show()
    #   Mean Temperature vs. timestep
    t,itemp,avgtemp = loadtxt('temperature.out', unpack=True)
    plot(t,itemp,'r.',t,avgtemp,'b-',t,kt_targ*ones(size(t)),'g-')
    show()
