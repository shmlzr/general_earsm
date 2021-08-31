#===============================================================================
# S U B R O U T I N E S
#===============================================================================
import numpy as np
import os

def getRANSScalar(case, time, var):
    #file = open(case + '/0/cellCentres','r').readlines()
    file = open(case + '/' + str(time) + '/' + var,'r').readlines()
    #lines=0
    tmp = []
    tmp2 = 10**12
    maxIter = -1
    #v= np.zeros([3,70000])
    #cc = False
    j = 0
    for i,line in enumerate(file):
        if 'internalField' in line:
            tmp = i + 1
            tmp2 = i + 3
            #cc = True
            #print(tmp, tmp2)
        elif i==tmp:
            maxLines = int(line.split()[0])
            maxIter  = tmp2 + maxLines
            v = np.zeros([1,maxLines])
            #print maxLines, maxIter
        elif i>=tmp2 and i<maxIter:
            #tmp3=i
            linei = line.split()
            #v[:,j] = [float(linei[0].split('(')[1]), float(linei[1]), float(linei[2].split(')')[0])]
            v[:,j] = [float(linei[0])]
            j += 1

    return v

def getRANSVector(case, time, var):
    #file = open(case + '/0/cellCentres','r').readlines()
    file = open(case + '/' + str(time) + '/' + var,'r').readlines()
    #lines=0
    tmp = []
    tmp2 = 10**12
    maxIter = -1
    #v= np.zeros([3,70000])
    #cc = False
    j = 0
    for i,line in enumerate(file):
        if 'internalField' in line:
            tmp = i + 1
            tmp2 = i + 3
            #cc = True
            #print tmp, tmp2
        elif i==tmp:
            maxLines = int(line.split()[0])
            maxIter  = tmp2 + maxLines
            v = np.zeros([3,maxLines])
           # print maxLines, maxIter
        elif i>=tmp2 and i<maxIter:
            linei = line.split()
            v[:,j] = [float(linei[0].split('(')[1]), float(linei[1]), float(linei[2].split(')')[0])]
            j += 1
    return v

def getRANSTensor(case, time, var):
    file = open(case + '/' + str(time) + '/' + var,'r').readlines()
    #file = open('wavyWall_Re6850_komegaSST_4L_2D/20000/gradU','r').readlines()
    #lines=0
    tmp = []
    tmp2 = 10**12
    maxIter = -1
    #v= np.zeros([3,70000])
    #cc = False
    j = 0
    for i,line in enumerate(file):
        if 'internalField' in line:
            tmp = i + 1
            tmp2 = i + 3
            #cc = True
#            print tmp, tmp2
        elif i==tmp:
            maxLines = int(line.split()[0])
            maxIter  = tmp2 + maxLines
            v = np.zeros([3,3,maxLines])
#            print maxLines, maxIter
        elif i>=tmp2 and i<maxIter:
            #tmp3=i
            linei = line.split()
            v[0,:,j] = [float(linei[0].split('(')[1]), float(linei[1]), float(linei[2])]
            v[1,:,j] = [float(linei[3]), float(linei[4]), float(linei[5])]
            v[2,:,j] = [float(linei[6]), float(linei[7]), float(linei[8].split(')')[0])]

            j += 1

    return v

def getRANSSymmTensor(case, time, var):
    file = open(case + '/' + str(time) + '/' + var,'r').readlines()
    #file = open('wavyWall_Re6850_komegaSST_4L_2D/20000/gradU','r').readlines()
    #lines=0
    tmp = []
    tmp2 = 10**12
    maxIter = -1
    #v= np.zeros([3,70000])
    #cc = False
    j = 0
    for i,line in enumerate(file):
        if 'internalField' in line:
            tmp = i + 1
            tmp2 = i + 3
            #cc = True
#            print tmp, tmp2
        elif i==tmp:
            maxLines = int(line.split()[0])
            maxIter  = tmp2 + maxLines
            v = np.zeros([3,3,maxLines])
#            print maxLines, maxIter
        elif i>=tmp2 and i<maxIter:
            #tmp3=i
            linei = line.split()
            v[0,:,j] = [float(linei[0].split('(')[1]), float(linei[1]), float(linei[2])]
            v[1,:,j] = [float(linei[1]), float(linei[3]), float(linei[4])]
            v[2,:,j] = [float(linei[2]), float(linei[4]), float(linei[5].split(')')[0])]

            j += 1

    return v

def getRANSPlane(mesh, dimension, nx, ny, t='vector',a='yesAve'):
    xRANS = np.zeros([nx,ny])
    yRANS = np.zeros([nx,ny])
    zRANS = np.zeros([nx,ny])
    print((mesh.shape))

    if a=='yesAve':
        if t=='vector':
            if dimension=='3D':
                for i in range(ny):
                    print('hello')
#                    xRANS[:,i] = mesh[0,i*256*128:256+i*256*128]
#                    yRANS[:,i] = mesh[1,i*256*128:256+i*256*128]
#                    zRANS[:,i] = mesh[2,i*256*128:256+i*256*128]
            elif dimension=='2D':
                for i in range(ny):
                    xRANS[:,i] = mesh[0,i*nx:nx+i*nx]
                    yRANS[:,i] = mesh[1,i*nx:nx+i*nx]
                    zRANS[:,i] = mesh[2,i*nx:nx+i*nx]

            return np.array([xRANS, yRANS, zRANS])

        elif t=='scalar':
            for i in range(ny):
                xRANS[:,i] = mesh[0,i*nx:nx+i*nx]

            return np.array([xRANS])

        elif t=='tensor':
            out = np.zeros([3,3,nx,ny])
            for i in range(ny):
                out[:,:,:,i] = mesh[:,:,i*nx:nx+i*nx]

            return out

    elif a=='noAve':
        if t=='vector':
            if dimension=='3D':
                for i in range(96):
                    xRANS[:,i] = mesh[0,i*256*128:256+i*256*128]
                    yRANS[:,i] = mesh[1,i*256*128:256+i*256*128]
                    zRANS[:,i] = mesh[2,i*256*128:256+i*256*128]
            elif dimension=='2D':
                for i in range(96):
                    xRANS[:,i] = mesh[0,i*256:256+i*256]
                    yRANS[:,i] = mesh[1,i*256:256+i*256]
                    zRANS[:,i] = mesh[2,i*256:256+i*256]

            return np.array([xRANS, yRANS, zRANS])

def getRANSField(mesh, dimension, nx, ny, nz, t='vector'):
    xRANS = np.zeros([nx,ny*2,nz])
    yRANS = np.zeros([nx,ny*2,nz])
    zRANS = np.zeros([nx,ny*2,nz])
    print(mesh.shape)

    kk=0

    if t=='vector':
        if dimension=='3D':
            print("3D")
            for j in range(nz):
                for i in range(ny):
                    xRANS[:,i,j] = mesh[0,i*nx+nx*ny*(j+kk*4):nx+i*nx+nx*ny*(j+kk*4)]
                    yRANS[:,i,j] = mesh[1,i*nx+nx*ny*(j+kk*4):nx+i*nx+nx*ny*(j+kk*4)]
                    zRANS[:,i,j] = mesh[2,i*nx+nx*ny*(j+kk*4):nx+i*nx+nx*ny*(j+kk*4)]

        elif dimension=='2D':
            for i in range(ny):
                xRANS[:,i] = mesh[0,i*nx:nx+i*nx]
                yRANS[:,i] = mesh[1,i*nx:nx+i*nx]
                zRANS[:,i] = mesh[2,i*nx:nx+i*nx]

        return np.array([xRANS, yRANS, zRANS])

    elif t=='scalar':
        if dimension=='3D':
            for i in range(ny):
                for j in range(nz):
                    xRANS[:,i,j] = mesh[0,i*nx+nx*ny*j:nx+i*nx+nx*ny*j]
        elif dimension=='2D':
            for i in range(ny):
                xRANS[:,i] = mesh[0,i*nx:nx+i*nx]

        return np.array([xRANS])

    elif t=='tensor':
        out = np.zeros([3,3,nx,ny])
        for i in range(ny):
            out[:,:,:,i] = mesh[:,:,i*nx:nx+i*nx]

        return out


def calcEigenvalues(tau, k):

    if len(tau.shape)==3:
        l=tau.shape[2]

        tauAni = np.zeros([3,3,l])
        for i in range(l):
                tauAni[:,:,i] = tau[:,:,i]/(2.*k[i]) - np.diag([1/3.,1/3.,1/3.])

        eigVal = np.zeros([3,l])
        for i in range(l):
            a,b=np.linalg.eig(tauAni[:,:,i])
            eigVal[:,i]=sorted(a, reverse=True)


    elif len(tau.shape)==4:
        l=tau.shape[2]
        l2=tau.shape[3]
        print(l,l2)
        tauAni = np.zeros([3,3,l,l2])
        for i in range(l):
            for j in range(l2):
                tauAni[:,:,i,j] = tau[:,:,i,j]/(2.*k[i,j]) - np.diag([1/3.,1/3.,1/3.])

        eigVal = np.zeros([3,l,l2])
        for i in range(l):
            for j in range(l2):
                a,b=np.linalg.eig(tauAni[:,:,i,j])
                eigVal[:,i,j]=sorted(a, reverse=True)


    return eigVal


def barycentricMap(eigVal):

    if len(eigVal.shape)==2:
        l=eigVal.shape[1]

        C1c = eigVal[0,:] - eigVal[1,:]
        C2c = 2*(eigVal[1,:] - eigVal[2,:])
        C3c = 3*eigVal[2,:] + 1
        Cc  = np.array([C1c, C2c, C3c])

        locX = np.zeros([l])
        locY = np.zeros([l])
        for i in range(l):
            locX[i] = Cc[0,i] + 0.5*Cc[2,i]
            locY[i] = np.sqrt(3)/2 * Cc[2,i]

    elif len(eigVal.shape)==3:
        l=eigVal.shape[1]
        l2=eigVal.shape[2]

        C1c = eigVal[0,:,:] - eigVal[1,:,:]
        C2c = 2*(eigVal[1,:,:] - eigVal[2,:,:])
        C3c = 3*eigVal[2,:,:] + 1
        Cc  = np.array([C1c, C2c, C3c])

        locX = np.zeros([l,l2])
        locY = np.zeros([l,l2])
        for i in range(l):
            for j in range(l2):
                locX[i,j] = Cc[0,i,j] + 0.5*Cc[2,i,j]
                locY[i,j] = np.sqrt(3)/2 * Cc[2,i,j]


    return np.array([locX, locY])


def calcFracAnisotropy(tau):
    if len(tau.shape)==3:
        l=tau.shape[2]

        eigVal = np.zeros([3,l])
        for i in range(l):
            a,b=np.linalg.eig(tau[:,:,i])
            eigVal[:,i]=sorted(a, reverse=True)
        #eigValMean = eigVal.mean(0)
        FA = np.zeros([l])
        for i in range(l):
#            FA[i] = np.sqrt(3/2.) * np.sqrt( (eigVal[0,i] - eigValMean[i])**2. + (eigVal[1,i] - eigValMean[i])**2. + (eigVal[2,i] - eigValMean[i])**2.) / np.sqrt(eigVal[0,i]**2. + eigVal[1,i]**2. + eigVal[2,i]**2.)
#            FA[i] = np.sqrt(1/2.) * np.sqrt( (eigVal[0,i] - eigVal[1,i])**2. + (eigVal[1,i] - eigVal[2,i])**2. + (eigVal[2,i] - eigVal[0,i])**2.) / np.sqrt(eigVal[0,i]**2. + eigVal[1,i]**2. + eigVal[2,i]**2.)
            FA[i] = np.sqrt(1/2.) * np.sqrt((eigVal[0,i] - eigVal[1,i])**2. + (eigVal[1,i] - eigVal[2,i])**2. + (eigVal[2,i] - eigVal[0,i])**2.) / np.sqrt(eigVal[0,i]**2. + eigVal[1,i]**2. + eigVal[2,i]**2.)


    elif len(tau.shape)==4:
        l=tau.shape[2]
        l2=tau.shape[3]
        print(l,l2)

        eigVal = np.zeros([3,l,l2])
        for i in range(l):
            for j in range(l2):
                a,b=np.linalg.eig(tau[:,:,i,j])
                eigVal[:,i,j]=sorted(a, reverse=True)

        print(l, l2)
        #eigValMean = eigVal.mean(0)
        FA = np.zeros([l,l2])
        for i in range(l):
            for j in range(l2):
#                FA[i,j] = np.sqrt(3/2.) * np.sqrt( (eigVal[0,i,j] - eigValMean[i,j])**2. + (eigVal[1,i,j] - eigValMean[i,j])**2. + (eigVal[2,i,j] - eigValMean[i,j])**2.) / np.sqrt(eigVal[0,i,j]**2. + eigVal[1,i,j]**2. + eigVal[2,i,j]**2.)
                FA[i,j] = np.sqrt(1/2.) * np.sqrt( (eigVal[0,i,j] - eigVal[1,i,j])**2. + (eigVal[2,i,j] - eigVal[1,i,j])**2. + (eigVal[2,i,j] - eigVal[0,i,j])**2.) / np.sqrt(eigVal[0,i,j]**2. + eigVal[1,i,j]**2. + eigVal[2,i,j]**2.)

    return FA

def calcInvariant(eigVal):
    if len(eigVal.shape)==2:
        l=eigVal.shape[1]

        II = np.zeros([l])
        III = np.zeros([l])
        for i in range(l):
            II[i] = 2.*(eigVal[0,i]**2. + eigVal[0,i] * eigVal[1,i] + eigVal[1,i]**2.)
            III[i] = -3*eigVal[0,i]*eigVal[1,i]*(eigVal[0,i] + eigVal[1,i])

    elif len(eigVal.shape)==3:
        l=eigVal.shape[1]
        l2=eigVal.shape[2]

        II = np.zeros([l,l2])
        III = np.zeros([l,l2])
        for i in range(l):
            for j in range(l2):
                II[i,j] = 2.*(eigVal[0,i,j]**2. + eigVal[0,i,j] * eigVal[1,i,j] + eigVal[1,i,j]**2.)
                III[i,j] = -3*eigVal[0,i,j]*eigVal[1,i,j]*(eigVal[0,i,j] + eigVal[1,i,j])

    return II,III


def calcInitialConditions(U, turbulenceIntensity, turbLengthScale, nu, d, D):

    k       = 1.5 * (U*turbulenceIntensity)**2.0
    epsilon = 0.16 * k**1.5 / turbLengthScale
    omega   = 1.8 * np.sqrt(k) / turbLengthScale
    #omega_farfield = U/D
    omega_wall = 10 * 6 * nu / (0.0705 * d**2)
    #omega_wall_wilcox = 6 / (0.0708 * yPlus_wilcox**2)
    nuTilda = np.sqrt(1.5)*U*turbulenceIntensity*turbLengthScale
    nuTilda_NASA = 3*nu
    nut_NASA = 3*nu*(3**3)/(3**3 + 7.1**3)
    Re      = U*D/nu
    tmp = {'k': k, 'epsilon': epsilon, 'omega': omega, 'nuTilda': nuTilda, 
           'omega_wall':omega_wall, 'Re':Re}
    return tmp


def create_case(case_template, case_name, location=''):
    os.system('cp -r ' + case_template + ' ' + location + case_name)


def run_case(case, solver):
    print(solver)
    os.system('cd ' + case + ';' + solver + ' >> log;')
