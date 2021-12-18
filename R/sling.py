#!/usr/bin/env python
# coding: utf-8

# In[3]:


import numpy as np

def standard_lasso(x,y,p,l=100):
    n,d = x.shape
    l = l/(2*n)
    w0 = np.mean(y)
    w = np.zeros(d)

    wr = np.zeros(d)
    v = np.zeros((d,d))
    def ai(i):
        return(np.sum(x[:,i]**2)/n)

    def zi(i):
        tot = 0
        for j in range(n):
            yhj = w0 + np.dot(w, x[j,:]) - w[i]*x[j,i]
            tot += (y[j] - yhj)*x[j,i]
        return(tot/n)  
    
    def KKT(i):
        e = 0.01
        g = np.dot(x[:,i],y - (np.dot(x,w)+w0))/(l*n)
        if w[i]>0:
            return ((g>1-e)&(g<1+e))
        elif w[i]<0:
            return ((g>-1-e)&(g<-1+e))
        else:
            return ((g>-1-e)&(g<1+e))

    for i in range(d):
        for j in range(d):
            v[i,j]=np.dot(x[:,i],x[:,j])
    nkkt = d
    itr = 0

    while(nkkt >= (d/p)):
        itr += 1
        for i in range(d):
            wr = w
            if(zi(i)>l):
                w[i] = (zi(i)-l)/ai(i)
            elif(zi(i)< -1*l):
                w[i] = (zi(i)+l)/ai(i)
            else:
                w[i] = 0
        nkkt = 0
        for i in range(d):
            if not KKT(i):
                nkkt += 1 
    obj = np.sum((y - (np.dot(x,w)+w0))**2)/(2*n)+l*(np.dot(w,w)**0.5)
    return w, itr, obj


def Sling(x,y,p,l=100): 
    n,d = x.shape
    l = l/(2*n)
    w0 = np.mean(y)
    w = np.zeros(d)
    wr = np.zeros(d)
    v = np.zeros((d,d))
    for i in range(d):
        for j in range(d):
            v[i,j]=np.dot(x[:,i],x[:,j])

    def ai(i):
        return(np.sum(x[:,i]**2)/n)

    def zi(i):
        tot = 0
        for j in range(n):
            yhj = w0 + np.dot(w, x[j,:]) - w[i]*x[j,i]
            tot += (y[j] - yhj)*x[j,i]
        return(tot/n)
            
    def KKT(i):
        e = 0.01
        g = np.dot(x[:,1],y - (np.dot(x,w)+w0))/(l*n)
        if w[i]>0:
            return ((g>1-e)&(g<1+e))
        elif w[i]<0:
            return ((g>-1-e)&(g<-1+e))
        else:
            return ((g>-1-e)&(g<1+e))
    
    def z_up(i):
         return w[i]-wr[i]+np.dot(v[i,:],v[i,:])*np.dot(w-wr,w-wr)/n+zr[i]

    def z_low(i):
         return w[i]-wr[i]-np.dot(v[i,:],v[i,:])*np.dot(w-wr,w-wr)/n+zr[i]

    def update(i):
        if(zi(i)>l):
            w[i] = (zi(i)-l)/ai(i)
        elif(zi(i)< -1*l):
            w[i] = (zi(i)+l)/ai(i)
        else:
            w[i] = 0

    U = []    
    itr = 0
    nkkt = d

    while (nkkt>=(d/p)):
        itr += 1
        wr = w
        zr = wr + (np.dot(y,x))/n
        err = 1e-5
        prev = np.array([np.inf]*d)
        while np.all(abs(w-prev)>err):
            prev = w
            for i in U:
                if (z_low(i) > l) | (z_up(i) < (-l)):
                    update(i)
        wr = w
        zr = wr + (np.dot(y,x))/n
        prev = np.array([np.inf]*d)
        while np.all(abs(w-prev)>err):
            prev = w
            for i in U:
                if (z_up(i) > l) | (z_low(i) < (-l)):
                    update(i)
                else:
                    w[i] = 0 
        nkkt = 0
        for i in range(d):
            if not KKT(i):
                nkkt+=1
                if i not in U:
                    U.append(i)
    obj = np.sum((y - (np.dot(x,w)+w0))**2)/(2*n)+l*(np.dot(w,w)**0.5)
    return w, itr, obj


# In[ ]:




