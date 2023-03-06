#
#	ABAQUS LINEAR ADDITIVE (LAYER) MANUFACTURING JOB GENERATION SCRIPT 
#	by D.P. Munro at TU Delft, Delft, The Netherlands, January 2018.
#
if __name__ == '__main__':
#
#	FE discretization
#
    e_x=3								# elements; x-axis
    e_y=3								# y-axis
    e_z=3								# z-axis; build direction 
    L=1e0/e_x	   		  				# element length-scale (mm) 
#
    n_e=e_x*e_y*e_z							# total number of elements
    n_n=(e_x+1)*(e_y+1)*(e_z+1)					# total number of nodes
#
#	Primary Abaqus FEA job file
#
    with open('alm.inp', 'w') as file:
        tmp='***';file.write('%s\n' % tmp)
        tmp='*Heading';file.write('%s\n' % tmp)
        tmp='***';file.write('%s\n' % tmp)
        tmp='*Include, input=alm_nds.inp';file.write('%s\n' % tmp)
        tmp='*Include, input=alm_els.inp';file.write('%s\n' % tmp)
        tmp='*Include, input=alm_sts.inp';file.write('%s\n' % tmp)
        tmp='*Include, input=alm_bcs.inp';file.write('%s\n' % tmp)
        tmp='*Include, input=alm_mat.inp';file.write('%s\n' % tmp)
        tmp='*Include, input=alm_stp.inp';file.write('%s\n' % tmp)
        tmp='***';file.write('%s\n' % tmp)
#
#	Node file
#
    lds=[]
    with open('alm_nds.inp', 'w') as file:
        tmp = '*Node';file.write('%s\n' % tmp)
        n=1
        for k in range(e_z+1):
            for j in range(e_y+1):
                for i in range(e_x+1):
#					Nodal coordinates in the reference configuration:
                    tmp = '\t %d, \t %e, \t %e, \t %e' % (n, 
                        (float(i)*L), (float(j)*L), (float(k)*L))
                    file.write('%s\n' % tmp)
                    if j == e_y and k == e_z:
                        lds.append("%d,%d,1."%(n,2)) 
                    n=n+1
#
#	Element and node connectivity file
#
    with open('alm_els.inp', 'w') as file:
#		Element type:
        tmp = '*Element, type=C3D8'; file.write('%s\n' % tmp)
        c1=1;c2=2;c3=e_x+2;c4=e_x+3;c5=(e_x+1)*(e_y+1)+1;c6=c5+1;c7=c5+e_x+1;c8=c5+e_x+2
        e=1
        for k in range(e_z):
            for j in range(e_y):
                for i in range(e_x):
#					Connectivity:
                    tmp = ' %d,  %d, %d, %d, %d,  %d, %d, %d, %d' % (e,
                        c2, c4, c3, c1, c6, c8, c7, c5)
                    file.write('%s\n' % tmp)
                    c1+=1; c2+=1; c3+=1; c4+=1; c5+=1; c6+=1; c7+=1; c8+=1
                    e=e+1
                c1+=1; c2+=1; c3+=1; c4+=1; c5+=1; c6+=1; c7+=1; c8+=1
            c1+=e_x+1; c2+=e_x+1; c3+=e_x+1; c4+=e_x+1; c5+=e_x+1; 
            c6+=e_x+1; c7+=e_x+1; c8+=e_x+1
#
#	Node-- and element-set file
#
    with open('alm_sts.inp', 'w') as file:
#		Elements in the reference configuration
        tmp = '*Elset,elset=set-ref';file.write('%s\n' % tmp)
        c=1
        for k in range(e_z):
            for j in range(e_y):
                for i in range(e_x):
                    tmp = '%d' % (c);file.write('%s\n' % tmp)
                    c=c+1
#		Nodes in the reference configuration:
        tmp = '*Nset,nset=net-ref';file.write('%s\n' % tmp)
        c=1
        for k in range(e_z+1):
            for j in range(e_y+1):
                for i in range(e_x+1):
                    tmp = '%d' % (c);file.write('%s\n' % tmp)
                    c=c+1
#		Nodes at the baseplate:
        tmp = '*Nset,nset=net-bsp';file.write('%s\n' % tmp)
        c1=1; c2=1
        for j in range(e_y+1):
            for i in range(e_x+1):
                tmp = '\t %d' % (c2);file.write('%s\n' % tmp)
                c2=c2+1
#		Elements per layer:
        for e_k in range(e_z):
#			Elements in layer k:
            tmp = '*Elset,elset=set-lay-%d' % (e_k+1); file.write('%s\n' % tmp)
            for j in range(e_y):
                for i in range(e_x):
                    tmp = '\t %d' % (c1);file.write('%s\n' % tmp)
                    c1=c1+1
            tmp = '*Nset,nset=net-lay-%d' % (e_k+1);file.write('%s\n' % tmp)
            for j in range(e_y+1):
                for i in range(e_x+1):
                    tmp = '\t %d' % (c2);file.write('%s\n' % tmp)
                    c2=c2+1
#
#	Boundary condition file
#
    with open('alm_bcs.inp', 'w') as file:
#		Fixed in all directions at the baseplate:
        tmp = '*Boundary'; file.write('%s\n' % tmp)
        tmp = '\t net-bsp, 1,3'; file.write('%s\n' % tmp)
#
#	Material file
#
    with open('alm_mat.inp', 'w') as file:
        tmp = '*Material, name=mat-alm';file.write('%s\n' % tmp)
        tmp = '*Elastic';file.write('%s\n' % tmp)
        tmp = '%e, %f' % (1.0, 0.3);file.write('%s\n' % tmp)
#       tmp = '*Plastic';file.write('%s\n' % tmp)
#       tmp = '1000e9,0.';file.write('%s\n' % tmp)
#       tmp = '1200e9,0.1';file.write('%s\n' % tmp)
#       tmp = '*Expansion';file.write('%s\n' % tmp)
#       tmp = '%e' % (1e-2);file.write('%s\n' % tmp)
        tmp = '*Solid Section, elset=set-ref, material=mat-alm';file.write('%s\n' % tmp)
#
#	Steps file
#
    with open('alm_stp.inp', 'w') as file:
#		Initial conditions; baseplate at reference temperature:
#		Step 0: Remove the configuration
        tmp = '*Step,nlgeom=no';file.write('%s\n' % tmp)
        tmp = '*Static, solver=iterative cholesky';file.write('%s\n' % tmp)
        tmp = '1.,1.,1.,1.';file.write('%s\n' % tmp)
        tmp = '*Cload';file.write('%s\n' % tmp)
        for l in lds:
            tmp = l
            file.write(tmp+"\n")
        tmp = '*Node file,frequency=1000';file.write('%s\n' % tmp)
        tmp = 'U';file.write('%s\n' % tmp)
        tmp = '*EL file,frequency=1000';file.write('%s\n' % tmp)
        tmp = '';file.write('%s\n' % tmp)
        tmp = '*End Step';file.write('%s\n' % tmp)
#
#       tmp = '*Step,nlgeom=no';file.write('%s\n' % tmp)
#       tmp = '*Static';file.write('%s\n' % tmp)
#       tmp = '0.1,1.,1e-5,1e30';file.write('%s\n' % tmp)
#       tmp = '*Model Change, type=element, remove';file.write('%s\n' % tmp)
#       tmp = 'set-lay-10';file.write('%s\n' % tmp)
#       tmp = 'set-lay-9';file.write('%s\n' % tmp)
#       tmp = 'set-lay-8';file.write('%s\n' % tmp)
#       tmp = 'set-lay-7';file.write('%s\n' % tmp)
#       tmp = 'set-lay-6';file.write('%s\n' % tmp)
#       tmp = 'set-lay-5';file.write('%s\n' % tmp)
#       tmp = 'set-lay-4';file.write('%s\n' % tmp)
#       tmp = 'set-lay-3';file.write('%s\n' % tmp)
#       tmp = 'set-lay-2';file.write('%s\n' % tmp)
#       tmp = '*Temperature';file.write('%s\n' % tmp)
#       tmp = 'net-lay-%d, -2.0' % (1);file.write('%s\n' % tmp)
#       tmp = '*End Step';file.write('%s\n' % tmp)
#
#		Step k: Add layer e_k strain free with immediate contraction
#
#       for e_k in range(1,e_k+1):
#           tmp = '*Step,nlgeom=no';file.write('%s\n' % tmp)
#           tmp = '*Static';file.write('%s\n' % tmp)
#           tmp = '0.1,1.,1e-5,1e30';file.write('%s\n' % tmp)
#           tmp = '*Model Change, type=element, add=strain free';file.write('%s\n' % tmp)
#           tmp = 'set-lay-%d' % (e_k+1);file.write('%s\n' % tmp)
#           tmp = '*End Step';file.write('%s\n' % tmp)
#           tmp = '*Step,nlgeom=no';file.write('%s\n' % tmp)
#           tmp = '*Static';file.write('%s\n' % tmp)
#           tmp = '0.1,1.,1e-5,1e30';file.write('%s\n' % tmp)
#           tmp = '*Temperature';file.write('%s\n' % tmp)
#           tmp = 'net-lay-%d, -2.0' % (e_k+1);file.write('%s\n' % tmp)
#           tmp = '*End Step';file.write('%s\n' % tmp)
#
#
