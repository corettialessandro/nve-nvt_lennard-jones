from numba import jit
from numpy import  zeros, rint

@jit(nopython=True,cache=True)
def calcener(rx, ry, rz, N, L ) :
    """ epsilon=1 , sigma=1 , rcut=3*sigma """
    c6  =  4.0
    cf6 = 24.0
    #c12 = c6
    #cf12= 2.*c12
    # standard cut-off = 3 for LJ systems
    r2cut=9./L**2
    # WCA cut-off in the minimum
    # r2cut=2.**(1./3.)
    dcut=  c6/r2cut**3
    acut= -dcut/r2cut**3
    fdx = zeros(N)
    fax = zeros(N)
    fdy = zeros(N)
    fay = zeros(N)
    fdz = zeros(N)
    faz = zeros(N)
    # zeroing
    enea=0.;
    vipa=0.
    ened=0.;
    vipd=0.
    # double loop on pairs
    for k in range(N-1) :
        for j in range(k+1,N) :
            dx = rx[k]-rx[j]
            dy = ry[k]-ry[j]
            dz = rz[k]-rz[j]
            dx-= rint(dx)
            dy-= rint(dy)
            dz-= rint(dz)
            r2 = dx**2 + dy**2 + dz**2
            if(r2<r2cut) :
                rr2 = 1./r2
                rr6 = rr2*rr2*rr2
                temp=c6*rr6
                ened+= (dcut - temp)
                enea+= (acut + temp*rr6)
                vir6 = cf6*rr6
                vir12=2.*vir6*rr6
                vipd-= vir6
                vipa+= vir12
                vir6 *=rr2
                vir12*=rr2
          # forces
                fdx[k]-= vir6*dx
                fax[k]+=vir12*dx
                fdy[k]-= vir6*dy
                fay[k]+=vir12*dy
                fdz[k]-= vir6*dz
                faz[k]+=vir12*dz
                fdx[j]+= vir6*dx
                fax[j]-=vir12*dx
                fdy[j]+= vir6*dy
                fay[j]-=vir12*dy
                fdz[j]+= vir6*dz
                faz[j]-=vir12*dz

    return ( enea, ened, vipa, vipd, fax, fdx, fay, fdy, faz, fdz )
