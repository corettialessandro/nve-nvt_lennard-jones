class MDLJ_lib :

  def __init__(self, N, Lr):
    from numpy import zeros
    #-start-------------------------
    self.t = 0.
    self.mass= 1.
    self.L   = Lr
    self.NHchi = 0.0
    self.NHxi = 0.0
    self.itemp = 0.0
    self.PL  = 0.
    self.Q   = 100.*N
    #eps     = 1.
    #sig     = 1.
    self.rx      = zeros( N )
    self.ry      = zeros( N )
    self.rz      = zeros( N )
    self.px      = zeros( N )
    self.py      = zeros( N )
    self.pz      = zeros( N )
    self.fax     = zeros( N )
    self.fay     = zeros( N )
    self.faz     = zeros( N )
    self.fdx     = zeros( N )
    self.fdy     = zeros( N )
    self.fdz     = zeros( N )
    self.kg      = 512
    self.gcount  = zeros( self.kg )
    self.ekin    = 0.0
    self.ene     = 0.0
    self.etot    = 0.0
    self.pext    = 1.
    # L can fluctuate rmax for gdr limited to initial L_ref/2.5
    rmax = Lr/2.5
    self.r2max = rmax * rmax
    self.ldel = rmax/self.kg

  def clean_run(self, keepinput=False):
      from glob import glob
      from os import remove

      if (keepinput):
          extlist = ['*.out', '*.xyz']
      else:
          extlist = ['*.out', '*.xyz', '*.b']
      for ext in extlist:
          filelist = glob(ext)
          for filename in filelist:
              remove(filename)
      print("# cleaning previous run")

  def read_input(self, N, conf_in='conf_in.b'):
      import pickle
      with open(conf_in, 'rb') as ftrj:
         (Nr, t) = pickle.load(ftrj)
         self.t = t
         if N!=Nr :
             print(' ??? reading %d particle from step %d configuration expected %d' % (Nr,t,N) )
         ( self.rx, self.ry, self.rz, self.L ) = pickle.load( ftrj)
         ( self.px, self.py, self.pz, self.PL) = pickle.load( ftrj)

  def write_input(self, N, t, conf_out='conf_in.b'):
      import pickle
      with open(conf_out, 'wb') as ftrj:
          pickle.dump( (N, t) , ftrj, pickle.HIGHEST_PROTOCOL)
          pickle.dump( ( self.rx, self.ry, self.rz, self.L ), ftrj, pickle.HIGHEST_PROTOCOL)
          pickle.dump( ( self.px, self.py, self.pz, self.PL), ftrj, pickle.HIGHEST_PROTOCOL)

  def fcc(self, m, verbose=False):
      from numpy import  random, rint
      print( "# number of lattice cells m^3 = %d" % (m**3) )
      a = self.L/m
      print( "# lattice parameter a = %f" % a )
      natom = 4*m**3
      print( "# number of particles N = %d " % natom )
      print( "# sides of md-box L = [ %.2f %.2f %.2f ]" % (a*m, a*m, a*m) )
      j  = 0
      xi = 0.
      yi = 0.
      zi = 0.
      delta=0.025
      rrx = random.normal(0., delta, natom)
      rry = random.normal(0., delta, natom)
      rrz = random.normal(0., delta, natom)
      for nx in range(m) :
          for ny in range(m) :
             for nz in range(m) :
                 self.rx[j] = xi + a*nx + rrx[j]
                 self.ry[j] = yi + a*ny + rry[j]
                 self.rz[j] = zi + a*nz + rrz[j]
                 if verbose: print( "  %d   %8.3f   %8.3f   %8.3f " % (j, self.rx[j], self.ry[j], self.rz[j]) )
                 # reduced box coordinates in [-0.5:0.5]^3
                 self.rx[j]/= self.L
                 self.rx[j]-= rint(self.rx[j])
                 self.ry[j]/= self.L
                 self.ry[j]-= rint(self.ry[j])
                 self.rz[j]/= self.L
                 self.rz[j]-= rint(self.rz[j])
                 j +=1
                 self.rx[j] = xi + a*nx + rrx[j] + 0.5*a
                 self.ry[j] = yi + a*ny + rry[j] + 0.5*a
                 self.rz[j] = zi + a*nz + rrz[j]
                 if verbose: print( "  %d   %8.3f   %8.3f   %8.3f " % (j, self.rx[j], self.ry[j], self.rz[j]) )
                 # reduced box coordinates in [-0.5:0.5]^3
                 self.rx[j]/= self.L
                 self.rx[j]-= rint(self.rx[j])
                 self.ry[j]/= self.L
                 self.ry[j]-= rint(self.ry[j])
                 self.rz[j]/= self.L
                 self.rz[j]-= rint(self.rz[j])
                 j +=1
                 self.rx[j] = xi + a*nx + rrx[j] + 0.5*a
                 self.ry[j] = yi + a*ny + rry[j]
                 self.rz[j] = zi + a*nz + rrz[j] + 0.5*a
                 if verbose: print( "  %d   %8.3f   %8.3f   %8.3f " % (j, self.rx[j], self.ry[j], self.rz[j]) )
                 # reduced box coordinates in [-0.5:0.5]^3
                 self.rx[j]/= self.L
                 self.rx[j]-= rint(self.rx[j])
                 self.ry[j]/= self.L
                 self.ry[j]-= rint(self.ry[j])
                 self.rz[j]/= self.L
                 self.rz[j]-= rint(self.rz[j])
                 j +=1
                 self.rx[j] = xi + a*nx + rrx[j]
                 self.ry[j] = yi + a*ny + rry[j] + 0.5*a
                 self.rz[j] = zi + a*nz + rrz[j] + 0.5*a
                 if verbose: print( "  %d   %8.3f   %8.3f    %8.3f " % (j, self.rx[j], self.ry[j], self.rz[j]) )
                 # reduced box coordinates in [-0.5:0.5]^3
                 self.rx[j]/= self.L
                 self.rx[j]-= rint(self.rx[j])
                 self.ry[j]/= self.L
                 self.ry[j]-= rint(self.ry[j])
                 self.rz[j]/= self.L
                 self.rz[j]-= rint(self.rz[j])
                 j +=1

      print( "# initial atom disposition built on fcc lattice")

  def vel_verlet(self, N, nstep, dt, freq=1):
    from numpy import sum, rint
    from calcener import calcener
    #import pickle
    tt   = 0
    ept  = 0.
    ekt  = 0.
    eptsq = 0.
    ektsq = 0.
    etsq = 0.
    vir = 0.
    itempavg = 0.
    ftempout = 'temperature.out'
    dth=0.5*dt
    dtm=dt/(self.mass*self.L**2)
    # initial call to compute starting energies, forces and virial
    epa,epd,vira,vird,self.fax,self.fdx,self.fay,self.fdy,self.faz,self.fdz = calcener(self.rx,self.ry,self.rz,N,self.L)
    self.fax/= self.L**12
    self.fdx/= self.L**6
    self.fay/= self.L**12
    self.fdy/= self.L**6
    self.faz/= self.L**12
    self.fdz/= self.L**6
    enep = epa/self.L**12+epd/self.L**6
    virp = vira/self.L**12 + vird/self.L**6
    enek = 0.5*sum( self.px*self.px + self.py*self.py + self.pz*self.pz )/(self.mass*self.L**2)
    self.itemp = self.calc_temp( N )
    itempavg = self.write_avg(self.itemp, itempavg, tt+1, ftempout)
    vcmx = sum(self.px)
    vcmy = sum(self.py)
    vcmz = sum(self.pz)
    data = 'run.out'
    out_data = open(data, 'a')
    out_data.write("#   'time'    'enep'    'enek'  'vir_p'   'enet'    'vcmx'   'vcmy'   'vcmz'    'temp'\n")
    out_data.write(" %8.3f %9.4g %9.4g %9.4g %10.7f %7.2e %7.2e %7.2e %9.4g\n" % (self.t, enep/N, enek/N, virp, enep+enek, vcmx, vcmy, vcmz, self.itemp))
    print( "   'time'    'enep'    'enek'  'vir_p'   'enet'    'vcmx'   'vcmy'   'vcmz'   'temp'")
    print (" %8.3f %9.4g %9.4g %9.4g %10.7f %7.2e %7.2e %7.2e %9.4g" % (self.t, enep/N, enek/N, virp, enep+enek, vcmx, vcmy, vcmz, self.itemp))
    for pas in range(nstep) :
        vcmx = 0.
        vcmy = 0.
        vcmz = 0.
        self.t += dt
        # advance one step
        # momenta first
        self.px += (self.fax+self.fdx)*dth
        self.py += (self.fay+self.fdy)*dth
        self.pz += (self.faz+self.fdz)*dth
        # positions second
        self.rx += dtm*self.px
        self.ry += dtm*self.py
        self.rz += dtm*self.pz
        self.rx -= rint(self.rx)
        self.ry -= rint(self.ry)
        self.rz -= rint(self.rz)
        # compute forces
        epa,epd,vira,vird,self.fax,self.fdx,self.fay,self.fdy,self.faz,self.fdz = calcener(self.rx,self.ry,self.rz,N,self.L)
        enep = epa/self.L**12 + epd/self.L**6
        virp = vira/self.L**12 + vird/self.L**6
        self.fax/= self.L**12
        self.fdx/= self.L**6
        self.fay/= self.L**12
        self.fdy/= self.L**6
        self.faz/= self.L**12
        self.fdz/= self.L**6
        # momenta thrid
        self.px += (self.fax+self.fdx)*dth
        self.py += (self.fay+self.fdy)*dth
        self.pz += (self.faz+self.fdz)*dth
        vcmx = sum(self.px)
        vcmy = sum(self.py)
        vcmz = sum(self.pz)
        enek = 0.5*sum( self.px*self.px + self.py*self.py + self.pz*self.pz )/(self.mass*self.L**2)
        # computing gdr and single step printout ...
        ekt += enek
        ept += enep
        vir += virp
        ektsq += enek*enek
        eptsq += enep*enep
        etsq += (enek+enep)*(enek+enep)
        if (pas+1)%freq==0 :
           # compute g(R)
           self.calc_gdr( N )
           # compute and write out running averages (for example temperature)
           self.itemp = self.calc_temp(N)
           itempavg = self.write_avg(self.itemp, itempavg, tt+1, ftempout)
           # save configuration for VMD in xyz format
           self.write_xyz( N )
           # save Phase-Space configuration for correlation functions
           self.write_PStraj( N )
           out_data.write(" %8.3f %9.4g %9.4g %9.4g %10.7f %7.2e %7.2e %7.2e %9.4g\n" % (self.t, enep/N, enek/N, virp, enep+enek, vcmx, vcmy, vcmz, self.itemp))
           print (" %8.3f %9.4g %9.4g %9.4g %10.7f %7.2e %7.2e %7.2e %9.4g" % (self.t, enep/N, enek/N, virp, enep+enek, vcmx, vcmy,vcmz, self.itemp) )
           tt += 1
        # end of md run
        # final configuration
    self.write_input(N, self.t, conf_out='conf_in.b')
    return (tt, ekt, ept, vir, ektsq, eptsq, etsq)

  def nose_hoover(self, N, nstep, dt, kt, tausq=0.1, freq=1):
    from numpy import sum, rint, exp
    from calcener import calcener
    #import pickle
    tt = 0
    enht=0.
    ept=0.
    ekt=0.
    vir=0.
    ektsq=0.
    eptsq=0.
    etsq=0.
    enhtsq=0.
    enhttsq=0.
    itempavg=0.
    G = 3.*N-3.
    ftempout = 'temperature.out'
    dth=0.5*dt
    dthh = .5*dth
    dtm=dt/(self.mass*self.L**2)
    # initial call to compute starting energies, forces and virial
    epa,epd,vira,vird,self.fax,self.fdx,self.fay,self.fdy,self.faz,self.fdz = calcener(self.rx,self.ry,self.rz,N,self.L)
    self.fax/= self.L**12
    self.fdx/= self.L**6
    self.fay/= self.L**12
    self.fdy/= self.L**6
    self.faz/= self.L**12
    self.fdz/= self.L**6
    enh = G*kt*(.5*tausq*self.NHchi**2 + self.NHxi)
    enep = epa/self.L**12+epd/self.L**6
    virp = vira/self.L**12 + vird/self.L**6
    enek = 0.5*sum( self.px*self.px + self.py*self.py + self.pz*self.pz )/(self.mass*self.L**2)
    self.itemp = self.calc_temp( N )
    fnh = (self.itemp/kt - 1.)/tausq
    itempavg = self.write_avg(self.itemp, itempavg, tt+1, ftempout)
    vcmx = sum(self.px)
    vcmy = sum(self.py)
    vcmz = sum(self.pz)
    data = 'NH_output_data.out'
    out_data = open(data, 'w+')
    out_data.write("#   'time'    'enep'    'enek'  'vir_p'   'enh'   'enet'        'enht'        'vcmx'   'vcmy'   'vcmz'   'inst_T'   'Delta_T'\n")
    out_data.write(" %8.3f %9.4g %9.4g %9.4g %9.4g %10.7f %10.7f %7.2e %7.2e %7.2e %9.4g %9.4g\n" % (self.t, enep/N, enek/N, virp, enh, enep+enek, enep+enek+enh, vcmx, vcmy, vcmz, self.itemp, self.itemp-kt))
    print( "   'time'    'enep'    'enek'  'vir_p'   'enh'   'enet'        'enht'        'vcmx'   'vcmy'   'vcmz'   'inst_T'   'Delta_T'")
    print (" %8.3f %9.4g %9.4g %9.4g %9.4g %10.7f %10.7f %7.2e %7.2e %7.2e %9.4g %9.4g" % (self.t, enep/N, enek/N, virp, enh, enep+enek, enep+enek+enh, vcmx, vcmy, vcmz, self.itemp, self.itemp-kt) )
    for pas in range(nstep) :
        vcmx = 0.
        vcmy = 0.
        vcmz = 0.
        self.t += dt
        # advance one step
        # Nose-Hoover coordinate first:
        self.NHxi += dthh*self.NHchi
        # Nose-Hoover momentum second:
        self.NHchi += dth*fnh
        # Nose-Hoover coordinate third
        self.NHxi += dthh*self.NHchi
        # momenta fourth
        	#Scaling of momenta
        s = exp(-dth*self.NHchi)
        f = (1 - s)/self.NHchi
        self.px = s*self.px + f*(self.fax+self.fdx)
        self.py = s*self.py + f*(self.fay+self.fdy)
        self.pz = s*self.pz + f*(self.faz+self.fdz)
        # positions fifth
        self.rx += dtm*self.px
        self.ry += dtm*self.py
        self.rz += dtm*self.pz
        self.rx -= rint(self.rx)
        self.ry -= rint(self.ry)
        self.rz -= rint(self.rz)
        # compute forces
        epa,epd,vira,vird,self.fax,self.fdx,self.fay,self.fdy,self.faz,self.fdz = calcener(self.rx,self.ry,self.rz,N,self.L)
        enep = epa/self.L**12 + epd/self.L**6
        virp = vira/self.L**12 + vird/self.L**6
        self.fax/= self.L**12
        self.fdx/= self.L**6
        self.fay/= self.L**12
        self.fdy/= self.L**6
        self.faz/= self.L**12
        self.fdz/= self.L**6
        # momenta seventh
 		#Scaling of momenta
        self.px = s*self.px + f*(self.fax+self.fdx)
        self.py = s*self.py + f*(self.fay+self.fdy)
        self.pz = s*self.pz + f*(self.faz+self.fdz)
        vcmx = sum(self.px)
        vcmy = sum(self.py)
        vcmz = sum(self.pz)
        enek = 0.5*sum( self.px*self.px + self.py*self.py + self.pz*self.pz )/(self.mass*self.L**2)
        # Nose-Hoover coordinate seventh
        self.NHxi += dthh*self.NHchi
        # Nose-Hoover momentum eighth
        self.itemp = self.calc_temp(N)
        fnh = (self.itemp/kt - 1.)/tausq
        self.NHchi += dth*fnh
        # Nose-Hoover coordinate nineth
        self.NHxi += dthh*self.NHchi
        enh = G*kt*(.5*tausq*self.NHchi**2 + self.NHxi)
        # computing gdr and single step printout ...
        ekt  += enek
        ept  += enep
        vir  += virp
        enht += enh
        ektsq += enek*enek
        eptsq += enep*enep
        etsq += (enek+enep)*(enek+enep)
        enhtsq += enh*enh
        enhttsq += (enek+enep+enh)*(enek+enep+enh)
        if (pas+1)%freq==0 :
           # compute g(R)
           self.calc_gdr(N)
           # compute and write out running averages (for example temperature)
           itempavg = self.write_avg(self.itemp, itempavg, tt+1, ftempout)
           # save configuration for VMD in xyz format
           self.write_xyz(N)
           # save Phase-Space configuration for correlation functions
           self.write_PStraj(N)
           print (" %8.3f %9.4g %9.4g %9.4g %9.4g %10.7f %10.7f %7.2e %7.2e %7.2e %9.4g %9.4g" % (self.t, enep/N, enek/N, virp, enh, enep+enek, enep+enek+enh, vcmx, vcmy, vcmz, self.itemp, self.itemp-kt) )
           out_data.write(" %8.3f %9.4g %9.4g %9.4g %9.4g %10.7f %10.7f %7.2e %7.2e %7.2e %9.4g %9.4g\n" % (self.t, enep/N, enek/N, virp, enh, enep+enek, enep+enek+enh, vcmx, vcmy, vcmz, self.itemp, self.itemp-kt))
           tt += 1
        # end of md run
        # final configuration
    out_data.close()
    self.write_input(N, self.t, conf_out='conf_in.b')
    return (tt, ekt, ept, vir, enht, ektsq, eptsq, etsq, enhtsq, enhttsq)

  def calc_temp(self, N):
    from numpy import sum
    g=3*N-3
    v2 = sum((self.px/self.L)**2 + (self.py/self.L)**2 + (self.pz/self.L)**2)
    return v2/g

  def write_avg(self, est, AVG, sample, fout):

      iAVG = (AVG*(sample-1) + est)/sample

      with open(fout, 'a') as fp:
          fp.write("%d\t%.4e\t%.4e\n" % (sample, est, iAVG))

      return iAVG

  def write_xyz(self, N):
      from numpy import zeros, rint
      dx = zeros(N)
      dy = zeros(N)
      dz = zeros(N)
      sig=3.4 # in Angstrom for argon
      rout=open('trajectory.xyz','a')
      rout.write('  %d \n' % N )
      rout.write('\n')
      for i in range(N):
          #ar[i] = "Ar"
          dx[i] = self.rx[i]
          dy[i] = self.ry[i]
          dz[i] = self.rz[i]
          dx[i] -= rint(dx[i])
          dy[i] -= rint(dy[i])
          dz[i] -= rint(dz[i])
          dx[i] *= sig*self.L
          dy[i] *= sig*self.L
          dz[i] *= sig*self.L
          rout.write('Ar   %12.5g   %12.5g   %12.5g\n' % (dx[i],dy[i],dz[i]) )
      rout.close()

  def write_PStraj(self, N):
      from numpy import zeros, rint
      dx = zeros(N)
      dy = zeros(N)
      dz = zeros(N)
      vx = zeros(N)
      vy = zeros(N)
      vz = zeros(N)
#      sig=3.4 # in Angstrom for Argon
#      eps = 1.656778224E-21 # in J for Argon
#      mass = 6.633521357E-26 # in kg for Argon
#      vconv = sqrt(eps/mass)*1e-2 # in Angstrom/ps
      with open('PStraj.out','a') as rout:
          for i in range(N):
              dx[i] = self.rx[i]
              dy[i] = self.ry[i]
              dz[i] = self.rz[i]
#              dx[i] -= rint(dx[i])
#              dy[i] -= rint(dy[i])
#              dz[i] -= rint(dz[i])
              dx[i] *= self.L#*sig
              dy[i] *= self.L#*sig
              dz[i] *= self.L#*sig
              vx[i] = self.px[i]/self.L#*vconv
              vy[i] = self.py[i]/self.L#*vconv
              vz[i] = self.pz[i]/self.L#*vconv
              rout.write('%12.5g\t%12.5g\t%12.5g\t%12.5g\t%12.5g\t%12.5g\n' % (dx[i],dy[i],dz[i],vx[i],vy[i],vz[i]) )

  def calc_gdr(self, N ):
    from numpy import sqrt, rint
    for k in range(N-1) :
        j=k+1
        dx = self.rx[k]-self.rx[j:N]
        dy = self.ry[k]-self.ry[j:N]
        dz = self.rz[k]-self.rz[j:N]
        dx[...]-= rint(dx)
        dy[...]-= rint(dy)
        dz[...]-= rint(dz)
        dx[...] = dx*self.L
        dy[...] = dy*self.L
        dz[...] = dz*self.L
        r2 = dx*dx + dy*dy + dz*dz
        # using the mask array "b" for speedup
        b = r2 < self.r2max
        lm  = sqrt(r2[b])
        for elm in lm :
            self.gcount[int(elm/self.ldel)]+=2.  # factor of 2 for gdr normalization

  def write_gdr(self, N, T, rho, gdr_out='gdr.out'):
      from numpy import zeros, pi, savetxt, column_stack
      V = zeros(self.kg)
      r = zeros(self.kg)
      g = zeros(self.kg)
      for lm in range(self.kg) :
          V[lm] = 4./3.*pi*(self.ldel**3)*(3*lm*lm +3*lm + 1);
          g[lm] = self.gcount[lm]/(V[lm]*(N -1)*T*rho);
          r[lm] = (lm+0.5)*self.ldel
      gout = column_stack( (r, g) )
      savetxt(gdr_out, gout , fmt=('%12.7g ','%12.7g'), header="    'r'     'g(r)'" )

  def calc_vcf(self, N, dt, freq, norm=False) :
      from numpy import sqrt, array, correlate, arange, mean, loadtxt, savetxt, c_, split, swapaxes, trapz

      def correlation_FFT(x1, x2, norm=True, mean=False):

        #computing the lenght of the vectors
        n1 = len(x1)
        n2 = len(x2)

        #checking if the vectors has the same lenght (MANDATORY)
        if (n1!=n2):
            print('different lenght vectors!')
            exit()

        #rename the variable
        n = n1

        #statistical analysis on data
        var1 = x1.var()
        var2 = x2.var()
        xx1 = x1
        xx2 = x2

        if mean==True:
            xx1 = x1 - x1.mean()
            xx2 = x2 - x2.mean()

        #computing correlation
        result = correlate(xx1, xx2, mode="full")[-n:]
        result /= arange(n, 0, -1)

        #normalizing the correlation function
        norm1 = sqrt(var1)
        norm2 = sqrt(var2)

        if norm:
            result /= (norm1*norm2)
            return result
        else:
            return result

      vx, vy, vz = loadtxt('PStraj.out', usecols=(3,4,5), unpack=True)
      if (len(vx)%N==0):
          Nt = len(vx)//N
          t = arange(Nt)*dt*freq
      else:
          raise TypeError('Dimension of velocities is not an integer multiple of the number of particles')

      vx, vy, vz = array(split(vx, Nt)), array(split(vy, Nt)), array(split(vz, Nt))
      vx, vy, vz = swapaxes(vx, 0, 1), swapaxes(vy, 0, 1), swapaxes(vz, 0, 1)

      VCTxx = array([correlation_FFT(vx[i], vx[i], norm=norm) for i in range(N)])
      VCTyy = array([correlation_FFT(vy[i], vy[i], norm=norm) for i in range(N)])
      VCTzz = array([correlation_FFT(vz[i], vz[i], norm=norm) for i in range(N)])

      VCTxy = array([correlation_FFT(vx[i], vy[i], norm=norm) for i in range(N)])
      VCTxz = array([correlation_FFT(vx[i], vz[i], norm=norm) for i in range(N)])
      VCTyx = array([correlation_FFT(vy[i], vx[i], norm=norm) for i in range(N)])
      VCTyz = array([correlation_FFT(vy[i], vz[i], norm=norm) for i in range(N)])
      VCTzx = array([correlation_FFT(vz[i], vx[i], norm=norm) for i in range(N)])
      VCTzy = array([correlation_FFT(vz[i], vy[i], norm=norm) for i in range(N)])

      VCT = array([[VCTxx, VCTxy, VCTxz], [VCTyx, VCTyy, VCTyz], [VCTzx, VCTzy, VCTzz]])

      AvVCT = mean(VCT, axis=2)

      VCF = mean([AvVCT[0][0],AvVCT[1][1],AvVCT[2][2]], axis=0)

      savetxt('vcf.out', c_[t, AvVCT[0][0], AvVCT[0][1], AvVCT[0][2], AvVCT[1][0], AvVCT[1][1], AvVCT[1][2], AvVCT[2][0], AvVCT[2][1], AvVCT[2][2]])

      if not norm:
          D = array([trapz(VCF[:n+1], dx=freq*dt) for n in range(Nt)])
          savetxt('vcf_diff.out', c_[t, D])

  def calc_MSD(self, N, dt, freq):
    from numpy import fft, ones, arange, square, append, zeros, loadtxt, array, split, swapaxes, mean, savetxt, c_
    from scipy.optimize import curve_fit

    def autocorr_FFT(x):

       N=len(x)
       F = fft.fft(x, n=2*N)
       PSD = F * F.conjugate()
       res = fft.ifft(PSD)
       res= (res[:N]).real
       n=N*ones(N)-arange(0,N)

       return res/n

    def msd_FFT(r):

      N=len(r)
      D=square(r).sum(axis=1)
      D=append(D,0)
      S2=sum([autocorr_FFT(r[:, i]) for i in range(r.shape[1])])
      Q=2*D.sum()
      S1=zeros(N)

      for m in range(N):
        Q=Q-D[m-1]-D[N-m]
        S1[m]=Q/(N-m)

      return S1-2*S2

    def lin_fit(x, A, B):

        return A + B*x

    def GetM(x, y):

        in_par = [1, 1]

        par, pcov = curve_fit(lin_fit, x, y, in_par)

        return par

    rx, ry, rz = loadtxt('PStraj.out', usecols=(0,1,2), unpack=True)
    if (len(rx)%N==0):
        Nt = len(rx)//N
        t = arange(Nt)*dt*freq
    else:
        raise TypeError('Dimension of positions is not an integer multiple of the number of particles')

    rx, ry, rz = array(split(rx, Nt)), array(split(ry, Nt)), array(split(rz, Nt))
    r = swapaxes(array([rx, ry, rz]), 0, 2)

    MSD = mean(array([msd_FFT(r[i]) for i in range(N)]), axis=0)

    tcut = t[:]
    MSDcut = MSD[:]
    D0, D = GetM(tcut, MSDcut)
    Dfit = D0 + D*t
    savetxt('msd.out', c_[t, MSD])
    savetxt('msd_diff.out', c_[t, Dfit, D/6.*ones(len(t))])
