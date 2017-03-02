import numpy as np

n_pairs = 10
grid_size = 300

m_M_pairs = np.zeros(shape=(2,n_pairs))
m_M_verify =np.zeros(shape=(2,n_pairs))

for j in range(0,n_pairs): #n_pairs is the number of time pairs
    N = grid_size**2
    x = 0.8
    m = float(np.random.randint(N/100,size=1)) #arbatry N/n to ensure m/N is small


    s=np.random.choice([0,1],size=(N,),p=[(1-m/N),m/N]) # make grid with s_i= 1 with prob m/N
    #print m,np.sum(s)
    s_1=np.squeeze(np.zeros(shape=(1,N)))


    for i in range(0,N,1):
        if s[i,]==1:
            s_1[i,] = np.random.choice([0,1],p=[(1-x),x])  # s_i detected correctly with prob x to give s^1 the underdetected version of s
        else:
            s_1[i,]=0

    M=np.sum(s_1)

    m_M_pairs[0,j]=m
    m_M_pairs[1,j]=float(M)

    #y=m-m*x
    y=M/x-M         #depends how we define y when verifying =! m and =! M

    s_ver=np.squeeze(np.zeros(shape=(1,N)))
    #reverse protocol, over sampling undersampling M to get m
    for i in range(0,N,1):
        if s_1[i,]==1:
            s_ver[i,] = np.random.choice([0,1],p=[((y/N)),1-(y/N)])
        else:
            s_ver[i,]=0

    #print s_ver
    m_test=(np.sum(s_ver))
    M_test=(x*(M+y))
    #print M,y,x

    m_M_verify[0,j]=(m-m_test)/N
    m_M_verify[1,j]=M-M_test


print m_M_pairs,
print m_M_verify

x_ML = (np.sum(m_M_pairs[1,:]))/(np.sum(m_M_pairs[0,:]))

print x_ML
