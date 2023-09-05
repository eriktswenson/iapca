
### Swenson - "Building Asymmetry into PCA and Application to ENSO" ###
# - Python script requires sst.mnmean.nc and olr.mon.mean.nc
#   as well as apca.py and inout.py
# - necessary packages: numpy, xarray, matplotlib, pandas, sklearn, statsmodels

# Author: Erik Swenson (latest revision Sep 2023)

import numpy as np
import inout as io
import apca as apca
import matplotlib.pyplot as plt


# Do bootstrapping significance test or not
signif = False

# read in SST and OLR and merge together
Y,C,g,idef,ndef1,ndef2,lat1,lon1,lat2,lon2,varattrs1,varattrs2,time = io.readmerge()

# Checks for equations and statements in paper
check = apca.check(Y)
i1 = np.arange(ndef1)
i2 = ndef1+np.arange(ndef2)
apca.subcheck(Y,i1,i2)

n,m = Y.shape
nit = 100 # number of iterations 
var = np.sqrt(n)*Y/g  # with weighting removed (correct units)

# EOF-1 pattern and PC-1 time series
X,beta,C,fve,delta,iconv,xi,v1,vpc1 = apca.apca(Y,nit=nit,seed=0,pca=True,l=1,prnt=False)
Y1 = np.matmul(X,beta)
eof1 = np.reshape(beta[1,:],(1,-1))
pc1 = np.reshape(X[:,1],(-1,1))
fveeof1 = fve[nit-1]

# EOF-2 pattern and PC-2 time series
X2,beta2,C2,fve,delta,iconv,xi,v1,vpc2 = apca.apca(Y-Y1,nit=nit,seed=0,pca=True,l=1,prnt=False)
Y2 = np.matmul(X2,beta2)
eof2 = np.reshape(beta2[1,:],(1,-1))
pc2 = np.reshape(X2[:,1],(-1,1))

# AEOF-1
X,beta,C,fve,delta,iconv,xi,v1,v2 = apca.apca(Y,nit=nit,xi=pc1,a='xi',prnt=False)
fveaeof1 = fve[nit-1]

# AEOF-1 reconstructed solution (Y∗)
Ystar = np.matmul(X,beta)


### NetCDF files for patterns shown in Fig. 1 ###

nlat1 = len(lat1); nlon1 = len(lon1)   # SST spatial dimensions
nlat2 = len(lat2); nlon2 = len(lon2)   # OLR spatial dimensions

varout = np.zeros((6,nlat1*nlon1+nlat2*nlon2)) + np.nan
varout[0,idef] = (vpc1*eof1/g)                     #  EOF-1
varout[1,idef] = -1*(vpc2*eof2/g)                  #  EOF-2 (flip sign)
varout[2,idef] = v1*v2*beta[1,:]/g                 # AEOF-1 P+
varout[3,idef] = v2*beta[2,:]/g                    # AEOF-1 P-
varout[4,idef] = np.sqrt(n)*(beta[0,:]-C)/g[0,:]   # AEOF-1 P0
varout[5,idef] = np.sqrt(n)*C/g[0,:]               # climatological mean
description = 'SST/OLR EOF-1(t=1), EOF-2(t=2), AEOF-1 p+(t=3),p-(t=4),p0(t=5), Clim(t=6)'

varout1 = np.reshape(varout[:,0:(nlat1*nlon1)],(6,nlat1,nlon1))
varname1 = 'sst'; fname1 = varname1+'.DJF.IAPCA.nc'
attrs1 = dict(title='ERSSTv5 tropical Pacific DJF SST',description=description)
io.write1(varout1,fname1,varname1,lon1,lat1,time=time[0:6],varattrs=varattrs1,attrs=attrs1)

varout2 = np.reshape(varout[:,(nlat1*nlon1):(nlat1*nlon1+nlat2*nlon2)],(6,nlat2,nlon2))
varname2 = 'olr'; fname2 = varname2+'.DJF.IAPCA.nc'
attrs2 = dict(title='NOAA tropical Pacific DJF OLR',description=description)
io.write1(varout2,fname2,varname2,lon2,lat2,time=time[0:6],varattrs=varattrs2,attrs=attrs2)


# Look at SST for AEOF-1
io.plot(varout1[2:5,:,:],lon1,lat1,cint=0.3,colors=['red','blue','gray'],leg=['AEOF-1 P+','AEOF-1 P-','AEOF-1 P0'])


### Fig. 2 time series' ###

# normalized APC-1 (Eq. A4) shown on y-axis of Fig.2a
apc1 = np.sum(X[:,1:3],1)
pos = np.where(X[:,1]!=0)[0]; apc1[pos] = np.sqrt(n)*apc1[pos]/(v1*v2)
neg = np.where(X[:,2]!=0)[0]; apc1[neg] = np.sqrt(n)*apc1[neg]/v2

# time series not chosen by H* (gray dots in Fig. 2a)
apcn = np.zeros(n)
apcn[pos] = np.sqrt(n)*np.matmul((Y[pos,:]-beta[0,:]),beta[2,:].T)/v2
apcn[neg] = np.sqrt(n)*np.matmul((Y[neg,:]-beta[0,:]),beta[1,:].T)/(v1*v2)

# alternative calculation (with spatial weighting and relative normalization)
#R,P = apca.apca(var,nit=nit,xi=pc1,a='xi',weight=g,prnt=False)
#print(np.allclose(np.sum(R[:,1:3],1),apc1))
#print(np.allclose(P[0,:],varout[4,idef]))
#print(np.allclose(P[1:3,:],varout[2:4,idef]))

# project Y∗ onto EOF-1-EOF-2 plane and divide by s.d. of PC-1
x1pc = np.sqrt(n)*np.matmul(Ystar-C,eof1.T)/vpc1        # Y∗  proj onto EOF-1
x2pc = np.sqrt(n)*np.matmul(Ystar-C,eof2.T)/vpc1        # Y∗  proj onto EOF-2
int1pc = np.sqrt(n)*np.matmul(beta[0,:]-C,eof1.T)/vpc1  # 𝜷∗0 proj onto EOF-1
int2pc = np.sqrt(n)*np.matmul(beta[0,:]-C,eof2.T)/vpc1  # 𝜷∗0 proj onto EOF-2

# normalize PCs for Fig. 2
pc1 = np.sqrt(n)*pc1/vpc1    # PC-1 on x-axis of Figs. 2a & 2b
pc2 = np.sqrt(n)*pc2/vpc1    # PC-2 on y-axis of Fig. 2b

# correlation between PC-1 and APC-1
cor1 = np.round(np.corrcoef(pc1[:,0],apc1)[0,1],2)


### Fig. 2a ###
xlim = np.array([-1.8,2.8]); xrange = xlim[1]-xlim[0]
ylim = np.array([-2.5,2.5]); yrange = ylim[1]-ylim[0]
plt.scatter(pc1[pos,0],apc1[pos],color='red',s=40)
plt.scatter(pc1[neg,0],apc1[neg],color='blue',s=40)
plt.scatter(pc1[:,0],apcn,color='gray',s=20)
temp = np.array([-3,3])
plt.plot(temp,temp+np.mean(apc1),color='gray',linestyle='dashed')
plt.plot([0,0],[-3,3],color='gray',linestyle='dashed')
plt.plot([-3,3],[0,0],color='gray',linestyle='dashed')
plt.xlabel('PC-1',fontsize=13); plt.ylabel('APC-1',fontsize=13)
plt.text(xlim[0]+0.05*xrange,ylim[1]-0.1*yrange,'(a)',fontsize=15,horizontalalignment='left')
plt.text(xlim[0]+0.05*xrange,ylim[1]-0.225*yrange,'positive',fontsize=13,color='red',horizontalalignment='left')
plt.text(xlim[0]+0.05*xrange,ylim[1]-0.3*yrange,'negative',fontsize=13,color='blue',horizontalalignment='left')
plt.text(xlim[0]+0.65*xrange,ylim[1]-0.1*yrange,'r = '+str(cor1),fontsize=12,horizontalalignment='left')
plt.xlim(xlim[0],xlim[1]); plt.ylim(ylim[0],ylim[1])
plt.savefig("Fig2a.pdf")
plt.show()


### Fig. 2b ###
xlim = np.array([-1.8,2.8]); xrange = xlim[1]-xlim[0]
ylim = np.array([-1.8,1.0]); yrange = ylim[1]-ylim[0]
plt.scatter(pc1,-1*pc2,color='black',s=25)
imin = np.argmin(x1pc[neg])
x = np.array([int1pc,x1pc[neg][imin]])
y = -1*np.array([int2pc,x2pc[neg][imin]])
plt.plot(x,y,color='blue',linewidth=2.5)
imax = np.argmax(x1pc[pos])
x = np.array([int1pc,x1pc[pos][imax]])
y = -1*np.array([int2pc,x2pc[pos][imax]])
plt.plot(x,y,color='red',linewidth=2.5)
plt.plot(xlim,[0,0],color='gray',linestyle='dashed')
plt.plot([0,0],ylim,color='gray',linestyle='dashed')
plt.text(0.75,-1.1,'AEOF-1',fontsize=13,color='red',horizontalalignment='left')
plt.text(0.75,-1.25,'positive',fontsize=13,color='red',horizontalalignment='left')
plt.text(-1,-1.1,'AEOF-1',fontsize=13,color='blue',horizontalalignment='left')
plt.text(-1,-1.25,'negative',fontsize=13,color='blue',horizontalalignment='left')
plt.text(xlim[0]+0.05*xrange,ylim[1]-0.1*yrange,'(b)',fontsize=15,horizontalalignment='left')
plt.xlabel('PC-1',fontsize=13); plt.ylabel('PC-2',fontsize=13)
plt.xlim(xlim[0],xlim[1]); plt.ylim(ylim[0],ylim[1])
plt.savefig("Fig2b.pdf")
plt.show()


### Statistical significance from bootstrapping approach ###
# - randomly shuffle timing of PC-1 100 times and re-compute EOF-1 and AEOF-1
# - estimate of PDF of variance differences when EOF-1 is truly independent of rest of data
# - AEOF-1 increase in variance only due to random chance, i.e. no real asymmetry present
# - contrast with AIC

if signif:
    nit1 = 100
    fves = np.zeros((nit1,2)); iconvs = np.zeros((nit1,2)); corrs = np.zeros((nit1,2))
    for z in range(0,nit1):
        Xa,beta,C,fve,delta,iconv,xi,v1,v2 = apca.apca(Y,nit=nit,seed=z,prnt=False)
        fves[z,0] = fve[nit-1]; iconvs[z,0] = iconv
        corrs[z,0] = np.corrcoef(np.sum(xi,1),pc1[:,0])[0,1]
        xi = pc1*z/nit1 + np.reshape(np.random.normal(0,1,n)*np.sqrt(np.mean(pc1**2))*(nit1-z)/nit1,(n,1))
        Xa,beta,C,fve,delta,iconv,xi,v1,v2 = apca.apca(Y,nit=nit,seed=z,xi=xi,prnt=False)
        fves[z,1] = fve[nit-1]; iconvs[z,1] = iconv
        corrs[z,1] = np.corrcoef(np.sum(xi,1),pc1[:,0])[0,1]
    nr = 9
    fveu = np.unique(np.round(fves[:,0],nr)); nu = len(fveu)
    fvecount = np.zeros(nu)
    for z in range(0,nu):
        fvecount[z] = np.sum(np.round(fves[:,0],nr)==fveu[z])
    fve1 = fveu[np.argmax(fvecount)]
    imost = np.where(np.round(fves[:,0],nr)==fve1)[0]
    plt.scatter(fves[:,0],iconvs[:,0],color='gray')
    plt.scatter(np.mean(fves[imost,0]),np.mean(iconvs[imost,0]),color='red',s=80)
    plt.xlabel('variance explained',fontsize=13); plt.ylabel('number of iterations',fontsize=13)
    plt.title('true solution emerges most frequently with the fastest convergence rates',fontsize=10)
    plt.show()
    isort = np.reshape(np.argsort(np.reshape(abs(corrs),2*nit1)),(10,-1))
    corr1 = np.zeros(10); count1 = np.zeros(10)
    for z in range(0,10):
        temp = np.reshape(abs(corrs),2*nit1)[isort[z,:]]
        corr1[z] = np.mean(temp)
        temp = np.round(np.reshape(fves,2*nit1)[isort[z,:]],nr)
        count1[z] = np.sum(temp==fve1)*100/20
    plt.plot(corr1,count1,linewidth=2,color='gray')
    plt.xlabel('correlation with PC-1',fontsize=13); plt.ylabel('% of time solution found',fontsize=13)
    plt.title('initialization from PC-1 (or in its vicinity) tends to yield the true solution',fontsize=10)
    plt.show()
    import statsmodels.api as sm
    np.random.seed(1)
    ish = np.arange(n)
    Yrest = Y - Y1
    nit1 = 100; nit2 = 100
    fvesh = np.zeros((nit1,2))
    aicsh = np.zeros((nit1,2))
    X1 = np.ones((n,2))
    for z in range(0,nit1):
        np.random.shuffle(ish)
        Ysh = Yrest + Y1[ish,:]
        X,beta,C,fve,delta,iconv,xi,v1,v2 = apca.apca(Ysh,nit=nit2,seed=0,pca=True,l=1,prnt=False)
        fvesh[z,0] = fve[nit-1]
        aicsh[z,0] = sm.OLS(Ysh,X).fit().aic
        xi = np.reshape(X[:,1],(-1,1))
        X,beta,C,fve,delta,iconv,xi,v1,v2 = apca.apca(Ysh,nit=nit2,xi=xi,prnt=False)
        fvesh[z,1] = fve[nit-1]
        aicsh[z,1] = sm.OLS(Ysh,X).fit().aic
    diff = fveaeof1-fveeof1
    diffsh = fvesh[:,1]-fvesh[:,0]
    diffrange = np.round([100*np.min(diffsh),100*np.max(diffsh)],2)
    print('Bootstrapping results (asymmetry due to random chance):')
    print('AEOF-1 increase in variance over EOF-1 is '+str(diffrange[0])+'-'+str(diffrange[1])+'%')
    print('Actual difference is '+str(np.round(100*diff,2))+'% which is statistically ')
    if (diff>=np.sort(diffsh)[95]):
        if (diff>=np.sort(diffsh)[98]):
            print('significant at the 99% confidence interval')
        else:
            print('significant at the 95% confidence interval')
    else:
        print('insignificant at the 95% confidence interval')
    aicdiff = aicsh[:,1]-aicsh[:,0]
    print('AIC is lower for AEOF-1 '+str(np.round(100*np.sum(aicdiff<0)/nit1))+'% of the time')
    if (np.sum(aicdiff<0)/nit1)>0.1:
        print('so it is not a good indicator of significance')


