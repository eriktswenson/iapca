

# - compute Iterative Asymmetric PCA
# - see Swenson - "Building Asymmetry into PCA and Application to ENSO"

# Author: Erik Swenson (latest revision Sep 2023)
 

import numpy as np

def apca(Y,nit=100,seed='null',xi=[0,0],l=2,pca=False,a='null',prnt=True,weight=[0,0]):

    n,m = Y.shape
#   n,l = xi.shape # (if xi is provided)
    name = 'EOF-1-EOF-'+str(l)+' '
    if pca==False:
        l = 2
        name = 'AEOF-1 '
    else:
        if l==1:
            name = 'EOF-1 '
    
    if (seed!='null'):
        np.random.seed(seed)

# if weight is provided, then account for spatial weighting
# and apply relative normalization at the end
    if len(weight)!=2:
        Y = Y*weight/np.sqrt(n)

# random initial time series'
    if len(xi)==2:
        xi = np.reshape(np.random.normal(0,1,n*l),(n,l))

    C = np.mean(Y,0)
    X = np.zeros((n,1+l)); X[:,0] = 1
    xl = np.zeros((n,l)); xl[:,:] = xi[:,:]

    Y1 = np.zeros((n,m)); Ylast = np.zeros((n,m))
    sse = np.zeros(nit); delta = np.zeros(nit)

# iterate nit times
    for z in range(0,nit):
# Eqs. 9 and 10 construct predictor time series'
        if pca:
            X[:,1:(1+l)] = xl[:,:]
        else:
            diagH = np.sum(xl,1)>0
            pos = np.where(diagH)[0]
            neg = np.where(diagH==False)[0]
            X[:,1:(1+l)] = 0
            X[pos,1] = xl[pos,0]
            X[neg,2] = xl[neg,1]
# Eq. 11 calculate regression solution
        XTXinv = np.linalg.inv(np.matmul(X.T,X))
        cov = np.matmul(X.T,Y)
        beta = np.matmul(XTXinv,cov)
# 𝜷'^k
        betap = np.reshape(beta[1:(1+l),:],(l,-1))
# Eq. 12 projections
        xl = np.matmul(Y-beta[0,:],betap.T)
        for zz in range(0,l):
            sd = np.sqrt(np.sum(betap[zz,:]**2))
            xl[:,zz] = xl[:,zz]/sd
        Ylast[:,:] = Y1[:,:]
        Y1 = np.matmul(X,beta)
# save sum-squared error and change in soln.
        sse[z] = np.sum((Y-Y1)**2)
        delta[z] = np.sum((Ylast-Y1)**2)

# fraction of variance explained
    fve = 1-sse/np.sum((Y-C)**2)
# check if convergence occurs
    delta = delta/np.sum((Y1-C)**2)
    dstring = 'does not converge'; iconv = nit
    temp = np.flip(np.where(delta<1e-6)[0]); ntemp = len(temp)
    if ntemp>0:
        for z in range(0,ntemp):
            if (temp[z]==(nit-1-z)):
                iconv = temp[z]+1
                dstring = 'converges in iteration ' + str(iconv)

    if prnt:
        print(name+str(np.round(fve[nit-1]*100,2))+'%:  '+dstring)

# Solve for v1 (Eq. A6)
    if pca==False:
        c1 = (-1/n)*np.sum(X[pos,1]); n1 = len(pos)
        c2 = (-1/n)*np.sum(X[neg,2]); n2 = len(neg)
        if a!='null':
# a is calculated based on initial time series xi
# - in this case, xi ought to be PC-1 with variance 
#   and a factor of 1/np.sqrt(n) contained in it (Eq. A5)
# - then a = ratio of positive to negative variance
            if a=='xi':
                a = np.sum(xi[pos,0]**2)/np.sum(xi[neg,0]**2)
            A1 = n1*c2*c2-a*np.sum((X[neg,2]+c2)**2)
            B1 = 2*a*n1*c1*c2-2*n2*c1*c2
            C1 = np.sum((X[pos,1]+c1)**2)-a*n2*c1*c1
            v1 = np.zeros(2)
            v1[0] = (-1*B1+np.sqrt(B1**2-4*A1*C1))/(2*A1)
            v1[1] = (-1*B1-np.sqrt(B1**2-4*A1*C1))/(2*A1)
            v1 = v1[v1>0][0]
        else:
            v1 = 1
# Solve for v2 (Eq. A7)
        v2 = np.sqrt((1/(v1**2))*np.sum(X[pos,1]**2)+np.sum(X[neg,2]**2)-n*((c1/v1+c2)**2))
    else:
        v1 = 1
        if l==1:
# mean of time series * -1
            c = (-1/n)*np.sum(X[:,1])
# correction for PCA for time series to have zero mean: 
# - shift mean of time series to the intercept
# - does not change np.matmul(X,beta)
            X[:,1] = X[:,1] + c
            beta[0,:] = beta[0,:] - c*beta[1,:]
            v2 = np.sqrt(np.sum(X[:,1]**2))
# - alternatively make no mean correction
#            v2 = np.sqrt(np.sum(X[:,1]**2)-n*(c**2))
        else:
            v2 = 1

# if weight is provided
    if len(weight)!=2:
# re-apply spatial weighting that was originally removed
# and apply relative normalization with v1 and v2
        Y = np.sqrt(n)*Y/weight
# 𝜷∗0 --> p0, 𝜷∗1 --> p+, 𝜷∗2 --> p- 
        beta[0,:] = np.sqrt(n)*(beta[0,:]-C)/weight
        beta[1,:] = v1*v2*beta[1,:]/weight
        X[:,1] = np.sqrt(n)*X[:,1]/(v1*v2)
        if l==2:
            beta[2,:] = v2*beta[2,:]/weight
            X[:,2] = np.sqrt(n)*X[:,2]/v2
        C = np.sqrt(n)*C/weight
# Return last iteration 
        return(X,beta)
    else:
# Return last iteration
        return(X,beta,C,fve,delta,iconv,xi,v1,v2)


def check(Y,nit=100):

    from sklearn.linear_model import LinearRegression

    n,m = Y.shape
    check = np.zeros(17)==1
    C = np.mean(Y,0)

# PCA from SVD
    U,s,V = np.linalg.svd(Y-C)
    pc = np.matmul(U[:,0:2],np.diag(s[0:2]))
    eof = V[0:2,:]
    Y12 = np.matmul(pc,eof) + C

    X,beta,C,fve,delta,iconv,xi,v1,v2 = apca(Y,nit=nit,seed=0,pca=True,l=2)
    Ystar = np.matmul(X,beta); check[0] = np.allclose(Y12,Ystar)
    if check[0]:
        print('Y∗ yields the climatological mean plus the EOF-1-EOF-2 plane')

# EOF-1 pattern and PC-1 time series
    X,beta,C,fve,delta,iconv,xi,v1,vpc1 = apca(Y,nit=nit,seed=0,pca=True,l=1)
    eof1 = np.reshape(beta[1,:],(1,-1)); pc1 = np.reshape(X[:,1],(-1,1))
    isign = np.sign(np.sum(eof[0,:]*eof1[0,:]))
    eof1 = eof1*isign; pc1 = pc1*isign; pc1 = pc1 - np.mean(pc1)
    Y1 = np.matmul(X,beta)
    fveeof1 = fve[nit-1]
    check[1] = np.allclose(eof[0,:],eof1[0,:])
    check[2] = np.allclose(pc[:,0],pc1[:,0])
    check[3] = np.allclose(Y1,np.matmul(pc1,eof1)+C)
    if (check[1]&check[2]&check[3]):
        print('EOF-1 pattern and PC-1 time series recovered')

    totvar = np.sum((Y-C)**2)
    X,beta,C2,fve,delta,iconv,xi,v1,vpc2 = apca(Y-Y1,nit=nit,seed=0,pca=True,l=1,prnt=False)
    eof2 = np.reshape(beta[1,:],(1,-1)); pc2 = np.reshape(X[:,1],(-1,1))
    pc2 = pc2-np.mean(pc2); Y2 = np.matmul(X,beta)
    fveeof2 = 1-np.sum((Y-C-Y2)**2)/totvar
    print('EOF-2 '+str(np.round(100*fveeof2,2))+'% converges in iteration '+str(iconv))

# AEOF-1
    X,beta,C,fve,delta,iconv,xi,v1,v2 = apca(Y,nit=nit,xi=pc1)
    Ystar = np.matmul(X,beta); xstar = np.sum(X[:,1:3],1)
    Ystarp = np.matmul(X[:,1:3],beta[1:3,:])
    Ystarm = np.mean(Ystar,0); Ystarpm = np.mean(Ystarp,0)
    pos = np.where(X[:,1]!=0)[0]; c1 = (-1/n)*np.sum(xstar[pos]); n1 = len(pos)
    neg = np.where(X[:,2]!=0)[0]; c2 = (-1/n)*np.sum(xstar[neg]); n2 = len(neg)
    print('Positive events '+str(n1)+' out of '+str(n)+' years')
    print('Negative events '+str(n2)+' out of '+str(n)+' years')
    fveaeof1 = fve[nit-1]
# variance ratio a for PC-1: compute correct v1 and v2
    a = np.sum(pc1[pos,0]**2)/np.sum(pc1[neg,0]**2)
    X,beta,C,fve,delta,iconv,xi,v1,v2 = apca(Y,nit=nit,xi=pc1,a=a,prnt=False)
# Eqs. 4 and 5
    check[4] = np.allclose(beta[0,:],C-Ystarpm)
    if check[4]:
        print('Eq. 4 holds')

    check[5] = np.allclose(beta[0,:]-C,c1*beta[1,:]+c2*beta[2,:])
    if check[5]:
        print('Eq. A9 holds')
        css1 = c1*np.sqrt(n)/(v1*v2); css2 = c2*np.sqrt(n)/v2
        print('𝜷∗0 = '+str(np.round(css1,2))+'𝜷∗1 + '+str(np.round(css2,2))+'𝜷∗2 (with weighted patterns)')

    mse = np.sum((Y-Ystar)**2)
    check[6] = np.allclose(mse,np.sum((Y-beta[0,:])**2)-np.sum(Ystarp**2))
    if check[6]:
        print('Eq. 5 holds')

# Eq. 15 "the variance of AEOF-1 may be expressed as"
    var = np.sum((Ystar-Ystarm)**2)
    var1 = np.sum(xstar**2)
    var2 = n*np.sum((beta[0,:]-C)**2)
    check[7] = np.allclose(var,var1-var2)
    if check[7]:
        print('Eq. 15 holds')

    U,s,V = np.linalg.svd(Y[pos,:]-beta[0,:])
    pcseg = np.reshape(U[:,0]*s[0],(-1,1)); eofseg = np.reshape(V[0,:],(1,-1))
    check[8] = np.allclose(np.matmul(pcseg,eofseg),Ystar[pos,:]-beta[0,:])
    U,s,V = np.linalg.svd(Y[neg,:]-beta[0,:])
    pcseg = np.reshape(U[:,0]*s[0],(-1,1)); eofseg = np.reshape(V[0,:],(1,-1))
    check[9] = np.allclose(np.matmul(pcseg,eofseg),Ystar[neg,:]-beta[0,:])
    if (check[8]&check[9]):
        print('Y∗−1𝜷∗0, is equivalent to EOF-1 of Y−1𝜷∗0 computed separately for each segment')

    temp = np.matmul(beta[1:3,:],beta[1:3,:].T)
    check[10] = np.allclose(np.array([1,1]),np.diag(temp))
    if check[10]:
        print('Positive and negative AEOF-1 patterns have unit variance')

    pattcorr = np.sum(beta[1,:]*beta[2,:]); angle = np.arccos(pattcorr)*180/np.pi
    print('Pattern correlation is '+str(np.round(pattcorr,3)))
    print('Asymmetry angle is '+str(np.round(angle,1))+' degrees')

    B = LinearRegression(fit_intercept=False).fit(X,Y).coef_.T
    check[11] = np.allclose(beta,B)
    if check[11]:
        print('Eq. 11 holds: 𝜷∗ is the least-squares solution from regressing Y with X∗')

    x1 = np.matmul((Y-beta[0,:]),beta[1,:].T)
    x2 = np.matmul((Y-beta[0,:]),beta[2,:].T)
    check[12] = np.allclose(xstar[pos],x1[pos])&np.allclose(xstar[neg],x2[neg])
    if check[12]:
        print('Eq. 12 holds: APC-1 equivalent to pattern projections')

    fve12 = 1 - np.sum((Ystar-Y12)**2)/np.sum((Y12-C)**2)
    print('AEOF-1 explains '+str(np.round(100*fve12))+'% of the variance of EOF-1-EOF-2')

# Eq. A3
    d = 2*n*c1*c2*(1-pattcorr)
    varalt = np.sum((xstar+c1+c2)**2) + d
    check[13] = np.allclose(var,varalt)
    if check[13]:
        print('Eq. A3 holds')

# Relative normalization
    x1 = x1/(v1*v2); x2 = x2/v2; css = c1/(v1*v2)+c2/v2
    apc1 = np.zeros(n); apc0 = np.zeros(n)
    apc1[pos] = x1[pos]; apc0[pos] = x2[pos]
    apc1[neg] = x2[neg]; apc0[neg] = x1[neg]
    check[14] = np.allclose(a,np.sum((apc1[pos]+css)**2)/np.sum((apc1[neg]+css)**2))
    check[15] = np.allclose(1,np.sum((apc1+css)**2))
    if (check[14]&check[15]):
        print('Eq. A5 holds')

    check[16] = (np.sum(apc1[pos]>apc0[pos])==n1)&(np.sum(apc1[neg]<apc0[neg])==n2)
    if check[16]:
        print('H∗ gives best choice for fixed patterns')

    print(str(np.sum(check))+' out of 17 checks pass')
    print(' ')

    return(check)

def subcheck(Y,i1,i2,nit=100):

    n,m = Y.shape
    check = np.zeros(3)==1
    C = np.mean(Y,0)

# EOF-1 pattern and PC-1 time series
    X,beta,C,fve,delta,iconv,xi,v1,vpc1 = apca(Y,nit=nit,seed=0,pca=True,l=1,prnt=False)
    eof1 = np.reshape(beta[1,:],(1,-1))
    pc1 = np.reshape(X[:,1],(-1,1)); pc1 = pc1 - np.mean(pc1)
    Y1 = np.matmul(X,beta)
    totvar1 = np.sum((Y-C)[:,i1]**2); fve1 = 1-np.sum((Y-Y1)[:,i1]**2)/totvar1
    totvar2 = np.sum((Y-C)[:,i2]**2); fve2 = 1-np.sum((Y-Y1)[:,i2]**2)/totvar2
    print('EOF-1 explains '+str(np.round(100*fve1,1))+'% and '+str(np.round(100*fve2,1))+'% of sub-domains')

# AEOF-1
    X,beta,C,fve,delta,iconv,xi,v1,v2 = apca(Y,nit=nit,xi=pc1,prnt=False)
    Ystar = np.matmul(X,beta)
    fve1 = 1-np.sum((Y-Ystar)[:,i1]**2)/totvar1
    fve2 = 1-np.sum((Y-Ystar)[:,i2]**2)/totvar2
    print('AEOF-1 explains '+str(np.round(100*fve1,1))+'% and '+str(np.round(100*fve2,1))+'% of sub-domains')

# pattern correlations
    denom = np.sqrt(np.sum(beta[1,i1]**2)*np.sum(beta[1,i1]**2))
    pattcorr1 = np.sum(beta[1,i1]*beta[2,i1])/denom
    denom = np.sqrt(np.sum(beta[1,i2]**2)*np.sum(beta[1,i2]**2))
    pattcorr2 = np.sum(beta[1,i2]*beta[2,i2])/denom
    print('Pattern correlations for sub-domains are '+str(np.round(pattcorr1,2))+' and '+str(np.round(pattcorr2,2)))
    print(' ')

