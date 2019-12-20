# enconding: utf-8
##########################################
#*************************************************************************
#*  Copyright (C) 2016 by Shiwei Zhao                                     *
#*  swzhao@scut.edu.cn                                                    *
#*  Programming Language: Python 2                                        *
#*  This program is free software; it is licensed under the terms of the  *
#*  GNU General Public License v2 or later. See file LICENSE for details. *
#*************************************************************************/

##############################################################################################
##This script was written for the visulization of 3D histograms in the following two
##publications. If you find it useful and use it in your research, please consider citing 
##one of the publications.
##[1] Zhao, S, and Zhou, X.(2017). Effects of Particle Asphericity on the Macro-
##and Micro-Mechanical Behaviors of Granular Assemblies. Granular Matter 19 (2):38.
##[2] Zhao, S, Evans, TM, and Zhou X.(2018). Three-Dimensional Voronoi
##Analysis of Monodisperse Ellipsoids during Triaxial Shear. Powder Technology 323:323-36.
import numpy as np

class VTK3Dhist():
    "A class for 3d histogram exporting to a vtk file for visulization with Paraview."
    def __init__(self,dimension = 8):
        #self.file_in,self.file_outpath,self.file_outname = file_in, file_outpath,file_outname
        self.dimension = dimension
        self.dimension2 = self.dimension*self.dimension
        self.elements_num = 6*self.dimension2
        self.delta = 2.0/self.dimension
        self.xx = 0
        self.yy = 0
        self.zz = 0
    #mapping from cube to sphere
    #A direct mapping from cube to (its inscribed) sphere:
    #http://catlikecoding.com/unity/tutorials/cube-sphere/
    def cube2sp(self,x,y,z):
        m,n = x.shape
        xx = np.zeros([m,n])
        yy = np.zeros([m,n])
        zz = np.zeros([m,n])
        for i in range(m):
            for j in range(n):
                x1 = x[i,j]
                y1 = y[i,j]
                z1 = z[i,j]
                #print 1.-x1**2/2-z1**2/2+x1**2.*z1**2/3
                xx[i,j] = x1*np.sqrt(1.-y1**2/2-z1**2/2+y1**2.*z1**2/3)
                yy[i,j] = y1*np.sqrt(1.-x1**2/2-z1**2/2+x1**2.*z1**2/3)
                zz[i,j] = z1*np.sqrt(1.-y1**2/2-x1**2/2+y1**2.*x1**2/3)
        #return xx,yy,zz
        self.xx = xx
        self.yy = yy
        self.zz = zz

    def generate_cube(self):
        " "
        #we construct a 3-d array for storing info from six faces of the cube.
        #x=0,x=1,y=0,y=1,z=0,z=1
        #cubeelements = zeros(6,self.dimension,self.dimension);%
        #cubenodes = zeros(6,self.dimension+1,self.dimension+1);
        #nodes_num = 6*(self.dimension-1)**2+12*(self.dimension-1)+8;
        #elements_num = 6*self.dimension2
        #update the private member variables
        self.dimension2 = self.dimension*self.dimension
        self.elements_num = 6*self.dimension2
        self.delta = 2.0/self.dimension

        elements = np.zeros([12,self.elements_num])
        #
        face_elements = np.zeros([2,self.dimension2])

        for i in range(self.dimension):
            for j in range(self.dimension):
                face_elements[:,i*self.dimension+j]=np.array([i*self.delta,j*self.delta])

        one = np.ones(self.dimension2);

        face1 = face_elements[0,:];
        face2 = face_elements[1,:];
        #num = self.dimension*self.dimension
        #print face1,face2,one
        temp = np.array([one,face1,face2])
        #print 'temp = ', temp
        indices = [0,0,1,1,2,2]
        for i in range(6):
            I = np.eye(3)
            index =  indices[i]
            tmp =  np.copy(I[:,index])#Caution: The problem is that Numpy basic slicing does not create copies of the actual data
            #print "tmp==",tmp
            I[:, index]=I[:,0]
            I[:,0]=tmp
            #print 'tmp2=',tmp
            #print "i====",i+1,index
            #print I
            #print "I temp=",I.dot(np.array([one,face1+self.delta,face2]))
            I[index,0]=np.mod(i,2)*2
            elements[0:3,i*self.dimension2:(i+1)*self.dimension2]=I.dot(temp)
            elements[3:6,i*self.dimension2:(i+1)*self.dimension2]=I.dot(np.array([one,face1+self.delta,face2]))
            elements[6:9,i*self.dimension2:(i+1)*self.dimension2]=I.dot(np.array([one,face1+self.delta,face2+self.delta]))
            elements[9: ,i*self.dimension2:(i+1)*self.dimension2]=I.dot(np.array([one,face1,face2+self.delta]))

        elements = elements-np.ones([12,self.elements_num])
        cube_x = elements[[0,3,6,9],:]
        cube_y = elements[[1,4,7,10],:]
        cube_z = elements[[2,5,8,11],:]

        # =========================================================================
        # Mapping nodes from cube to sphere.
        # =========================================================================
        #print cube_x
        #print cube_y
        self.cube2sp(cube_x,cube_y,cube_z)
        #self.writeVTK('cube.vtk',self.xx,self.yy,self.zz)
    def transform(self,x,y,z):
        #print x.shape
        NumData=len(x)
        #print NumData
        #transform the coordinates (x,y,z) on a sphere to a cubical surface (cube_x,cube_y,cube_z)
        cube_x=np.zeros(NumData)#coordinate x in the cubical system
        cube_y=np.zeros(NumData)# y
        cube_z=np.zeros(NumData)# z
        for i in range(NumData):
            #find the maximum of |x[i]|,|y[i]| and |z[i]|
            abc = [abs(x[i]),abs(y[i]),abs(z[i])]
            index = [0,1,2]#the first one corresponds to the maximum of a,b,and c
            if abc[0] < abc[1]:#swap the indices of a and b
                index[0],index[1] = index[1],index[0]
                if abc[1] < abc[2]:
                    index[0],index[2] = index[2],index[0]
            elif abc[0] < abc[2]:
                index[0],index[2] = index[2],index[0]
            xyz = [x[i],y[i],z[i]]
            xyz_t = [xyz[j] for j in index]
            xyz2 = [2.0*j**2 for j in xyz_t]
            tmp = [0.0,0.0,0.0]
            xyz_r = [0.,0.,0.]
            #print xyz2
            delta = -np.sqrt((-xyz2[1]+xyz2[2]-3)**2-12*xyz2[1])
            tmp[0] = 1.0
            tmp[1] = np.sqrt((delta+xyz2[1]-xyz2[2]+3)/2)
            tmp[2] = np.sqrt((delta-xyz2[1]+xyz2[2]+3)/2)
            for j in range(3):
                xyz_r[index[j]] = tmp[j]*np.sign(xyz[index[j]])
            #print xyz_r,i
            cube_x[i],cube_y[i],cube_z[i] = xyz_r
        return cube_x,cube_y,cube_z

    def local_index(self,x,delta):#find the element index that x falls in
        return np.fix((x+1)/self.delta)+1 #index of the element array
    def find_j_ab(self,n,x,y,z):
        a1 = (n+0.5*(x+1))*self.dimension**2
        a2 = self.dimension*(self.local_index(z,self.delta)-1)+self.local_index(y,self.delta)
        return int(a1+a2)

    def find_j(self,x,y,z):
        #x_ind = self.local_index(x,self.delta)
        #y_ind = self.local_index(y,self.delta)
        #z_ind = self.local_index(z,self.delta)
        n = 0
        if abs(y) == 1:
            x,y = y,x
            n = 2
        elif abs(z) == 1:
            x,z = z,x
            n = 4
        return self.find_j_ab(n,x,y,z),self.find_j_ab(n,-x,-y,-z)


    def writeVTK(self,filename,xx,yy,zz):
        print("writing Polydata into a VTK file")
        f_out = open(filename,'w')

        #writing the file head
        print>>f_out,"# vtk DataFile Version 2.0"
        print>>f_out,"Sway data processing"
        print>>f_out,"ASCII"
        print>>f_out,"DATASET POLYDATA"

        m,n = xx.shape
        print>>f_out,"POINTS   ", n*4," float"  ###index is from 0 in VTK
        #output points
        for i in range(n):#n elements or polygons
            for j in range(m):
                print>>f_out, xx[j,i],yy[j,i],zz[j,i]
        #output polygons
        print>>f_out,"POLYGONS ", n,n*5  #polygon_num polygon_num*5
        for i in range(n):
            print>>f_out,"4 ",i*4,i*4+1,i*4+2,i*4+3
        f_out.close()
    def drawBar(self,p1,p2,p3,p4,p5):
        #p1~p4: vertices on the top face
        #p5: the origin
        #convert to string
        p1_s = ' ' + str(p1)
        p2_s = ' ' + str(p2)
        p3_s = ' ' + str(p3)
        p4_s = ' ' + str(p4)
        p5_s = ' ' + str(p5)

        line = ""
        #top face
        line += "4 "+ p1_s + p2_s + p3_s + p4_s
        #side faces
        line += " 3 "+ p5_s + p2_s + p1_s  #5 2 1
        line += " 3 "+ p5_s + p1_s + p4_s  #5 1 4
        line += " 3 "+ p5_s + p4_s + p3_s  #5 4 3
        line += " 3 "+ p5_s + p3_s + p2_s  #5 3 2
        return line


    def write3DhistogramVTK(self,filename,coeff):
        print("writing Polydata into a VTK file")
        f_out = open(filename,'w')

        #writing the file head
        print>>f_out,"# vtk DataFile Version 2.0"
        print>>f_out,"Sway data processing"
        print>>f_out,"ASCII"
        print>>f_out,"DATASET POLYDATA"

        m,n = self.xx.shape
        print>>f_out,"POINTS   ", n*4 + 1," float"  ###index is from 0 in VTK; here we include the origin
        #output points
        print>>f_out, 0,0,0  #output the origin
        for i in range(n):#n elements or polygons
            for j in range(m):
                print>>f_out, self.xx[j,i]*coeff[i],self.yy[j,i]*coeff[i],self.zz[j,i]*coeff[i]
        #output polygons
        print>>f_out,"POLYGONS ", n*5,n*21  #polygon_num polygon_num*5
        for i in range(n):
            #print>>f_out,"4 ",i*4,i*4+1,i*4+2,i*4+3
            print>>f_out,self.drawBar(i*4+1,i*4+2,i*4+3,i*4+4,0)
        #output magnitude
        print>>f_out,"CELL_DATA ", n*5  #
        print>>f_out,"SCALARS cell_scalars float 1"
        print>>f_out,"LOOKUP_TABLE default"
        #print values of all cells (polygons)
        for i in range(n):
            print>>f_out,coeff[i],coeff[i],coeff[i],coeff[i],coeff[i]
        f_out.close()


    #The following functions may need customizing for your own purpose.
    def VTK3Dhistogram_basic(self,file_input,file_output):#three-column data
        self.dimension = 8;
        self.generate_cube()#we construct a 3-d array for storing info from six faces of the cube.
        # Reading data
        orient = np.loadtxt(file_input,delimiter='\t', skiprows=1)
        x=orient[:,0]
        y=orient[:,1]
        z=orient[:,2]
        #print x.shape
        NumData=len(x)
        cube_x,cube_y,cube_z = self.transform(x,y,z)

        # Cheking to see an orientation is inside which element in cubical system.

        coeff = np.zeros(self.elements_num)
        for i in range(NumData):
            j1,j2 = self.find_j(cube_x[i],cube_y[i],cube_z[i])
            coeff[j1-1] += 1
            coeff[j2-1] += 1

        coeff= (coeff/NumData)/(4.0*np.pi/self.elements_num)

        self.write3DhistogramVTK(file_output+".vtk",coeff)
    def VTK3DhistogramCN(self,file_input,file_output,shift):#for orientation of Voronoi cell
        #shift = 3; % 0 v_max, shortest axis direction; 3 v_min, longest axis direction
        self.dimension = 10;
        self.generate_cube()#we construct a 3-d array for storing info from six faces of the cube.
        # Reading contacts' orientations.
        orient = np.loadtxt(file_input,delimiter=' ', skiprows=1)
        #v_max and v_min are imported
        x=orient[:,1+shift];
        y=orient[:,2+shift];
        z=orient[:,3+shift];
        #print x.shape
        NumData=len(x)
        cube_x,cube_y,cube_z = self.transform(x,y,z)

        # Cheking to see an orientation is inside which element in cubical system.

        coeff = np.zeros(self.elements_num)
        for i in range(NumData):
            j1,j2 = self.find_j(cube_x[i],cube_y[i],cube_z[i])
            coeff[j1-1] += 1
            coeff[j2-1] += 1


        coeff= (coeff/NumData)/(4.0*np.pi/self.elements_num)
        if shift == 0:#v_max
            suffix = 'shortest.vtk'
        elif shift == 3:#v_min
            suffix = 'longest.vtk'
        else:
            suffix ='.vtk'
        self.write3DhistogramVTK(file_output+suffix,coeff)



    def VTK3Dhistogram(self,file_input,file_outputpath,file_outputname):
        "file_input: path of the input file.\nfile_outpath: directory of the output file.\nfile_outputname: name of the output file."
        shift = 0
        #shift = 3; % 0 v_max, shortest axis direction; 3 v_min, longest axis direction
        self.generate_cube()
        # =========================================================================
        # Reading contacts' orientations.
        data = np.loadtxt(file_input,delimiter=' ', skiprows=2) #No. c_type object1.id object2.id nforce(x,y,z) sforce(x,y,z)
        #find the range of data
        ids = data[:,0]
        m = int(2*len(ids)-np.sum(ids)) #total number of lines with particle-particle contacts
        #data1 = data[:m,3+shift:6+shift]
        fn_v = data[:m,3:6]
        ft_v = data[:m,6:9]
        fn = np.sqrt(fn_v[:,0]**2.0+fn_v[:,1]**2.0+fn_v[:,2]**2.0)
        ft = np.sqrt(ft_v[:,0]**2.0+ft_v[:,1]**2.0+ft_v[:,2]**2.0)
        avg_fn1 = np.average(fn)

        #v_max and v_min are imported
        #contact normal: [x,y,z]
        x=fn_v[:,0]/fn;
        y=fn_v[:,1]/fn;
        z=fn_v[:,2]/fn;
        for i in range(len(fn)):
            if fn[i]>10.0*avg_fn1:#the value larger than 10<Fn> has a pretty low probability of less than 1e-5
                fn[i]=10.0*avg_fn1
        #print "size of data: ",orient.shape
        #concatenate vectors
        #x = np.concatenate((x,-x))
        #y = np.concatenate((y,-y))
        #z = np.concatenate((z,-z))

        ##########
        NumData = len(x)
        cube_x,cube_y,cube_z = self.transform(x,y,z)
        # Checking to see an orientation is inside which element in cubical system.

        cn_coeff = np.zeros(self.elements_num)
        fn_coeff = np.zeros(self.elements_num)
        ft_coeff = np.zeros(self.elements_num)
        #print elements_num
        for i in range(NumData):
            #print mappx[i],mappy[i],mappz[i]
            j1,j2 = self.find_j(cube_x[i],cube_y[i],cube_z[i])
            #contact normal
            cn_coeff[j1-1] += 1     # Number of normal forces pointing within an element.
            cn_coeff[j2-1] += 1
            #normal contact force
            fn_coeff[j1-1] += fn[i]     # normal forces pointing within an element.
            fn_coeff[j2-1] += fn[i]
            #tangential contact force
            ft_coeff[j1-1] += ft[i]     # Number of normal forces pointing within an element.
            ft_coeff[j2-1] += ft[i]

        cn_coeff1= (cn_coeff/NumData/2.0)/(4.0*np.pi/self.elements_num)

        fn_coeff = fn_coeff/cn_coeff
        avg_fn = sum(fn_coeff)/self.elements_num
        fn_coeff /=avg_fn
        ft_coeff = ft_coeff/cn_coeff/avg_fn

        fname = file_outputpath+'/cn'+file_outputname+'.vtk'
        self.write3DhistogramVTK(fname,cn_coeff1)
        fname = file_outputpath+'/fn'+file_outputname+'.vtk'
        self.write3DhistogramVTK(fname,fn_coeff)
        fname = file_outputpath+'/ft'+file_outputname+'.vtk'
        self.write3DhistogramVTK(fname,ft_coeff)


    def VTK3Dhistogram2(self,file_input,file_outputpath,file_outputname):
        shift = 0
        #shift = 3; % 0 v_max, shortest axis direction; 3 v_min, longest axis direction
        self.generate_cube()
        # =========================================================================
        # Reading contacts' orientations.
        data = np.loadtxt(file_input,delimiter=' ', skiprows=2) #No. c_type object1.id object2.id nforce(x,y,z) sforce(x,y,z)
        #find the range of data
        ids = data[:,0]
        m = 2*len(ids)-np.sum(ids) #total number of lines with particle-particle contacts
        #data1 = data[:m,3+shift:6+shift]
        fn_v = data[:m,3:6]
        ft_v = data[:m,6:9]
        fn1 = np.sqrt(fn_v[:,0]**2.0+fn_v[:,1]**2.0+fn_v[:,2]**2.0)
        ft1 = np.sqrt(ft_v[:,0]**2.0+ft_v[:,1]**2.0+ft_v[:,2]**2.0)
        #v_max and v_min are imported
        #contact normal: [x,y,z]
        x = list()
        y = list()
        z = list()
        fn = list()
        ft = list()
        for i in range(len(ft1)):
            #print ft_v[i,0],ft[i]
            if ft1[i] >0:
                x.append(ft_v[i,0]/ft1[i])
                y.append(ft_v[i,1]/ft1[i])
                z.append(ft_v[i,2]/ft1[i])
                fn.append(fn1[i])
                ft.append(ft1[i])
                #print "ft=0"
            #print ft_v[i,0]/ft[i]
        #x=ft_v[:,0]/ft;
        #y=ft_v[:,1]/ft;
        #z=ft_v[:,2]/ft;
        #print "size of data: ",orient.shape
        #concatenate vectors
        #x = np.concatenate((x,-x))
        #y = np.concatenate((y,-y))
        #z = np.concatenate((z,-z))

        #print x.shape
        NumData=len(x)
        #print NumData
        cube_x,cube_y,cube_z = self.transform(x,y,z)
        # Cheking to see an orientation is inside which element in cubical system.

        cn_coeff = np.zeros(elements_num)
        fn_coeff = np.zeros(elements_num)
        ft_coeff = np.zeros(elements_num)
        #print elements_num
        for i in range(NumData):
            j1,j2 = self.find_j(cube_x[i],cube_y[i],cube_z[i])
            #contact normal
            cn_coeff[j1-1] += 1     # Number of normal forces pointing within an element.
            cn_coeff[j2-1] += 1
            #normal contact force
            fn_coeff[j1-1] += fn[i]     # normal forces pointing within an element.
            fn_coeff[j2-1] += fn[i]
            #tangential contact force
            ft_coeff[j1-1] += ft[i]     # Number of normal forces pointing within an element.
            ft_coeff[j2-1] += ft[i]

        cn_coeff1= (cn_coeff/NumData/2.0)/(4.0*np.pi/elements_num)
        avg_fn = np.average(fn)
        fn_coeff = fn_coeff/cn_coeff/avg_fn
        ft_coeff = ft_coeff/cn_coeff/avg_fn

        #fname = file_outputpath+'/cn'+file_outputname+'.vtk'
        #self.write3DhistogramVTK(fname,xx,yy,zz,cn_coeff1)
        #fname = file_outputpath+'/fn'+file_outputname+'.vtk'
        #self.write3DhistogramVTK(fname,xx,yy,zz,fn_coeff)
        fname = file_outputpath+'/ft'+file_outputname+'.vtk'
        self.write3DhistogramVTK(fname,ft_coeff)
