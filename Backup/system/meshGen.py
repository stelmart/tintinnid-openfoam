#This is a program to create the blockMeshDict file for the tintinnid simulation

from math import *
import sys

#~~~~~~~~~~~~~~~~PARAMETERS~~~~~~~~~~~~~~~~#

r_oral = 25.0                       #radius of oral cavity (must be < r_container). Also the major radius of container torus.
r_torus_multiplier = 8            #minor radius of container torus, as a multiple of l_down (whole numbers only)
h_tin = 100.0                        #height of tintinnid
l_down = 20.0                       #length of cilium during downward stroke. Rounded, depending on mesh.
l_up = 3*l_down/4                        #length of cilium during upward stroke (l_up < l_down for feeding). Rounded, depending on mesh.
num_cilia = 64                     #number of cilia (power of 2)
num_wavelength = 4                #number of wavelengths that fit around the oral cavity (power of 2, less than num_cilia)
frequency = float(sys.argv[1])                    #frequency of cilia in Hz
current_time = 0.0                  #current time in seconds

asymmetry = 0.55                  #0.5 < asymmetry <= 1. A value of 1 indicates that there is no asymmetry in beating patterns.  

phi_divs = 2                      #number of polar mesh divisions between cilia -- used to scale azimuthal divisions as well. Must be >1.
r_divs = 4                        #number of radial divisions from r=0 to r=r_oral.
r_inner_divs = 4                  #number of radial divisions from r_minor=0 to r_minor=l_down.
z_divs = 10                        #number of radial divisions from z=-h_tin to z=0.
theta_divs = 32                   #number of divisions in "azimuthal" angle, theta

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#Don't touch these
convertToMeters = "0.000001"      #sets the length scale wrt 1 meter (e.g. 0.001 means 1 = 1mm)
r_torus = r_torus_multiplier*l_down     
phi_divs_tot = phi_divs*num_cilia
omega_amplitude = pi**2*frequency  #theta(t) = pi/2 * cos(2 pi f t) --> omega(t) = -pi^2 f sin(2 pi f t)
r_minor_divs = r_inner_divs*r_torus_multiplier

meshfile = open('blockMeshDict', 'w')
p = open('../0/p', 'a')
U = open('../0/U', 'a')

#Write headers
meshfile.write("/*--------------------------------*- C++ -*----------------------------------*\ \n");
meshfile.write("| =========                 |                                                 | \n");
meshfile.write("| \\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           | \n");
meshfile.write("|  \\\    /   O peration     | Version:  5.0                                   | \n");
meshfile.write("|   \\\  /    A nd           | Web:      www.OpenFOAM.org                      | \n");
meshfile.write("|    \\\/     M anipulation  |                                                 | \n");
meshfile.write("\*---------------------------------------------------------------------------*/ \n");
meshfile.write("FoamFile \n");
meshfile.write("{ \n");
meshfile.write("    version     2.0; \n");
meshfile.write("    format      ascii; \n");
meshfile.write("    class       dictionary; \n");
meshfile.write("    object      blockMeshDict; \n");
meshfile.write("} \n \n");


#Length scale
meshfile.write("convertToMeters " + convertToMeters + "; \n \n");



#~~~~~~~~~~~~~~~Set vertices~~~~~~~~~~~~~~~~~~~~#
meshfile.write("vertices\n(\n");

#z<=0 vertices
sin_phi = []
cos_phi = []

for i in range(0,phi_divs_tot):
    phi = 2*pi*i/phi_divs_tot
    sin_phi.append(sin(phi))
    cos_phi.append(cos(phi))
    if abs(sin_phi[i]) < 1e-15:
        sin_phi[i] = 0
    if abs(cos_phi[i]) < 1e-15:
        cos_phi[i] = 0


for k in range(0,z_divs+1):
    for j in range(0,r_minor_divs+1):
        for i in range(0,phi_divs_tot):
            meshfile.write("    (" + str((r_oral + j*r_torus/r_minor_divs)*cos_phi[i]) + " " + str((r_oral + j*r_torus/r_minor_divs)*sin_phi[i]) + " " + str(h_tin*(1.0*k/z_divs-1)) + ")\n");


#z>0 torus vertices

sin_theta = []
cos_theta = []

for i in range(0,theta_divs+1):
    theta = pi/2*i/theta_divs
    sin_theta.append(sin(theta))
    cos_theta.append(cos(theta))
    if abs(sin_theta[i]) < 1e-15:
        sin_theta[i] = 0
    if abs(cos_theta[i]) < 1e-15:
        cos_theta[i] = 0

for k in range(1,theta_divs+1):
    for j in range(1,r_minor_divs+1):
        for i in range(0,phi_divs_tot):
            meshfile.write("    (" + str((r_oral + j*r_torus/r_minor_divs*cos_theta[k])*cos_phi[i]) + " " + str((r_oral + j*r_torus/r_minor_divs*cos_theta[k])*sin_phi[i]) + " " + str(j*r_torus/r_minor_divs*sin_theta[k]) + ")\n");


#Upper cylinder, r<=r_oral
for k in range(1,r_divs):
    for j in range(0,r_minor_divs+1):
        for i in range(0,phi_divs_tot):
            meshfile.write("    (" + str(k*r_oral/r_divs*cos_phi[i]) + " " + str(k*r_oral/r_divs*sin_phi[i]) + " " + str(j*r_torus/r_minor_divs) +")\n");


#z-axis, z<=0
for j in range(0,r_minor_divs+1):
    meshfile.write("    (0 0 " + str(j*r_torus/r_minor_divs) +")\n");

    
#End
meshfile.write(");\n\n");





#~~~~~~~~~~~~~~~~Set blocks~~~~~~~~~~~~~~~~~~~~~#



meshfile.write("blocks\n(\n");

#z<0 blocks
for k in range(0,z_divs):
    for j in range(0,r_minor_divs):
        for i in range(0,phi_divs_tot-1):
            meshfile.write("    hex (" + str(i + j*phi_divs_tot + k*phi_divs_tot*(r_minor_divs+1)) + " " + str(i + (j+1)*phi_divs_tot + k*phi_divs_tot*(r_minor_divs+1))  + " " + str(i+1 + (j+1)*phi_divs_tot + k*phi_divs_tot*(r_minor_divs+1)) + " " + str(i+1 + j*phi_divs_tot + k*phi_divs_tot*(r_minor_divs+1)) + " " + str(i + j*phi_divs_tot + (k+1)*phi_divs_tot*(r_minor_divs+1)) + " " + str(i + (j+1)*phi_divs_tot + (k+1)*phi_divs_tot*(r_minor_divs+1)) + " " + str(i+1 + (j+1)*phi_divs_tot + (k+1)*phi_divs_tot*(r_minor_divs+1))  + " " + str(i+1 + j*phi_divs_tot + (k+1)*phi_divs_tot*(r_minor_divs+1)) +") (1 1 1) simpleGrading (1 1 1)\n");

for k in range(0,z_divs):
    for j in range(0,r_minor_divs):
        meshfile.write("    hex (" + str(phi_divs_tot-1 + j*phi_divs_tot + k*phi_divs_tot*(r_minor_divs+1)) + " " + str(phi_divs_tot-1 + (j+1)*phi_divs_tot + k*phi_divs_tot*(r_minor_divs+1))  + " " + str((j+1)*phi_divs_tot + k*phi_divs_tot*(r_minor_divs+1)) + " " + str(j*phi_divs_tot + k*phi_divs_tot*(r_minor_divs+1)) + " " + str(phi_divs_tot-1 + j*phi_divs_tot + (k+1)*phi_divs_tot*(r_minor_divs+1)) + " " + str(phi_divs_tot-1 + (j+1)*phi_divs_tot + (k+1)*phi_divs_tot*(r_minor_divs+1)) + " " + str((j+1)*phi_divs_tot + (k+1)*phi_divs_tot*(r_minor_divs+1))  + " " + str(j*phi_divs_tot + (k+1)*phi_divs_tot*(r_minor_divs+1)) +") (1 1 1) simpleGrading (1 1 1)\n");

        

#Sew z<0 to z>0
i_oral = phi_divs_tot*(r_minor_divs+1)*z_divs
i_init = phi_divs_tot*(r_minor_divs+1)*(z_divs+1)
for j in range(0,r_minor_divs-1):
    for i in range(0,phi_divs_tot-1):
        meshfile.write("    hex (" + str(i_init - r_minor_divs*phi_divs_tot + i + j*phi_divs_tot) + " " + str(i_init - r_minor_divs*phi_divs_tot + i + (j+1)*phi_divs_tot)  + " " + str(i_init - r_minor_divs*phi_divs_tot + i+1 + (j+1)*phi_divs_tot) + " " + str(i_init - r_minor_divs*phi_divs_tot + i+1 + j*phi_divs_tot) + " " + str(i_init + i + j*phi_divs_tot) + " " + str(i_init + i + (j+1)*phi_divs_tot) + " " + str(i_init + i+1 + (j+1)*phi_divs_tot)  + " " + str(i_init + i+1 + j*phi_divs_tot) +") (1 1 1) simpleGrading (1 1 1)\n");

for j in range(0,r_minor_divs-1):
    meshfile.write("    hex (" + str(i_init - r_minor_divs*phi_divs_tot + phi_divs_tot-1 + j*phi_divs_tot) + " " + str(i_init - r_minor_divs*phi_divs_tot + phi_divs_tot-1 + (j+1)*phi_divs_tot)  + " " + str(i_init - r_minor_divs*phi_divs_tot + (j+1)*phi_divs_tot) + " " + str(i_init - r_minor_divs*phi_divs_tot + j*phi_divs_tot) + " " + str(i_init + phi_divs_tot-1 + j*phi_divs_tot) + " " + str(i_init + phi_divs_tot-1 + (j+1)*phi_divs_tot) + " " + str(i_init + (j+1)*phi_divs_tot)  + " " + str(i_init + j*phi_divs_tot) +") (1 1 1) simpleGrading (1 1 1)\n");
        
#z>0 torus

#Calculate which blocks to exclude (cilia BCs). cilia_thetas is a list of which j values to exclude.
cilia_thetas = []
cilia_speeds = [] #dtheta/dt
for i in range(0,num_cilia):
    arg2 = 1.0*i*num_wavelength/num_cilia
    arg1 = (arg2 - floor(arg2)) + frequency*current_time
    cilia_thetas.append(int(round((1.0*theta_divs-2)/2*cos(2*pi*arg1**asymmetry) + (1.0*theta_divs-2)/2 + 1))) #Goes from 1 to theta_divs-1
    if arg1==0:
        cilia_speeds.append(0)
    else:
        cilia_speeds.append(-asymmetry*(arg1**(asymmetry-1))*omega_amplitude*sin(2*pi*arg1**asymmetry))
    if cilia_thetas[i] == theta_divs - 1:
        cilia_thetas[i] = theta_divs - 2

cilia_thetas.append(theta_divs)   #This makes sure the if statement doesn't give an index out of range error

for k in range(0,theta_divs-1):
    for j in range(0,r_minor_divs-1):
        for i in range(0,phi_divs_tot-1):
            if (i%phi_divs == 0 and k == cilia_thetas[i/phi_divs]-1) and ((j+1)*r_torus/r_minor_divs < l_up or ((j+1)*r_torus/r_minor_divs < l_down and cilia_speeds[i/phi_divs] < 0)):
                print("skip")
            else:
                meshfile.write("    hex (" + str(i_init + i + j*phi_divs_tot + k*phi_divs_tot*r_minor_divs) + " " + str(i_init + i + (j+1)*phi_divs_tot + k*phi_divs_tot*r_minor_divs)  + " " + str(i_init + i+1 + (j+1)*phi_divs_tot + k*phi_divs_tot*r_minor_divs) + " " + str(i_init + i+1 + j*phi_divs_tot + k*phi_divs_tot*r_minor_divs) + " " + str(i_init + i + j*phi_divs_tot + (k+1)*phi_divs_tot*r_minor_divs) + " " + str(i_init + i + (j+1)*phi_divs_tot + (k+1)*phi_divs_tot*r_minor_divs) + " " + str(i_init + i+1 + (j+1)*phi_divs_tot + (k+1)*phi_divs_tot*r_minor_divs)  + " " + str(i_init + i+1 + j*phi_divs_tot + (k+1)*phi_divs_tot*r_minor_divs) +") (1 1 1) simpleGrading (1 1 1)\n");

for k in range(0,theta_divs-1):
    for j in range(0,r_minor_divs-1):
        meshfile.write("    hex (" + str(i_init + phi_divs_tot-1 + j*phi_divs_tot + k*phi_divs_tot*r_minor_divs) + " " + str(i_init + phi_divs_tot-1 + (j+1)*phi_divs_tot + k*phi_divs_tot*r_minor_divs)  + " " + str(i_init + (j+1)*phi_divs_tot + k*phi_divs_tot*r_minor_divs) + " " + str(i_init + j*phi_divs_tot + k*phi_divs_tot*r_minor_divs) + " " + str(i_init + phi_divs_tot-1 + j*phi_divs_tot + (k+1)*phi_divs_tot*r_minor_divs) + " " + str(i_init + phi_divs_tot-1 + (j+1)*phi_divs_tot + (k+1)*phi_divs_tot*r_minor_divs) + " " + str(i_init + (j+1)*phi_divs_tot + (k+1)*phi_divs_tot*r_minor_divs)  + " " + str(i_init + j*phi_divs_tot + (k+1)*phi_divs_tot*r_minor_divs) +") (1 1 1) simpleGrading (1 1 1)\n");

        
#Torus wedges
for i in range(0,phi_divs_tot-1):
    meshfile.write("    hex (" + str(i_oral + i) + " " + str(i_oral + phi_divs_tot + i)  + " " + str(i_oral + phi_divs_tot + i+1) + " " + str(i_oral + i+1) + " " + str(i_oral + i) + " " + str(i_init + i) + " " + str(i_init + i+1)  + " " + str(i_oral + i+1) +") (1 1 1) simpleGrading (1 1 1)\n");

meshfile.write("    hex (" + str(i_oral + phi_divs_tot-1) + " " + str(i_oral + 2*phi_divs_tot-1)  + " " + str(i_oral + phi_divs_tot) + " " + str(i_oral) + " " + str(i_oral + phi_divs_tot-1) + " " + str(i_init + phi_divs_tot-1) + " " + str(i_init)  + " " + str(i_oral) +") (1 1 1) simpleGrading (1 1 1)\n");

for k in range(0,theta_divs-1):
    for i in range(0,phi_divs_tot-1):
        if (i%phi_divs == 0 and k == cilia_thetas[i/phi_divs]-1):
            print("skip")
        else:
            meshfile.write("    hex (" + str(i_oral + i) + " " + str(i_init + i + k*phi_divs_tot*r_minor_divs)  + " " + str(i_init + i+1 + k*phi_divs_tot*r_minor_divs) + " " + str(i_oral + i+1) + " " + str(i_oral + i) + " " + str(i_init + i + (k+1)*phi_divs_tot*r_minor_divs) + " " + str(i_init + i+1 + (k+1)*phi_divs_tot*r_minor_divs)  + " " + str(i_oral + i+1) +") (1 1 1) simpleGrading (1 1 1)\n");

for k in range(0,theta_divs-1):
    meshfile.write("    hex (" + str(i_oral + phi_divs_tot-1) + " " + str(i_init + phi_divs_tot-1 + k*phi_divs_tot*r_minor_divs)  + " " + str(i_init + k*phi_divs_tot*r_minor_divs) + " " + str(i_oral) + " " + str(i_oral + phi_divs_tot-1) + " " + str(i_init + phi_divs_tot-1 + (k+1)*phi_divs_tot*r_minor_divs) + " " + str(i_init + (k+1)*phi_divs_tot*r_minor_divs)  + " " + str(i_oral) +") (1 1 1) simpleGrading (1 1 1)\n");


    
#Cylinder
i_cyl = i_init + r_minor_divs*phi_divs_tot*theta_divs

for k in range(0,r_divs-2):
    for j in range(0,r_minor_divs):
        for i in range(0,phi_divs_tot-1):
            meshfile.write("    hex (" + str(i_cyl + i + j*phi_divs_tot + k*phi_divs_tot*(r_minor_divs+1)) + " " + str(i_cyl + i + j*phi_divs_tot + (k+1)*phi_divs_tot*(r_minor_divs+1))  + " " + str(i_cyl + i+1 + j*phi_divs_tot + (k+1)*phi_divs_tot*(r_minor_divs+1)) + " " + str(i_cyl + i+1 + j*phi_divs_tot + k*phi_divs_tot*(r_minor_divs+1)) + " " + str(i_cyl + i + (j+1)*phi_divs_tot + k*phi_divs_tot*(r_minor_divs+1)) + " " + str(i_cyl + i + (j+1)*phi_divs_tot + (k+1)*phi_divs_tot*(r_minor_divs+1)) + " " + str(i_cyl + i+1 + (j+1)*phi_divs_tot + (k+1)*phi_divs_tot*(r_minor_divs+1))  + " " + str(i_cyl + i+1 + (j+1)*phi_divs_tot + k*phi_divs_tot*(r_minor_divs+1)) +") (1 1 1) simpleGrading (1 1 1)\n");

for k in range(0,r_divs-2):
    for j in range(0,r_minor_divs):
        meshfile.write("    hex (" + str(i_cyl + phi_divs_tot-1 + j*phi_divs_tot + k*phi_divs_tot*(r_minor_divs+1)) + " " + str(i_cyl + phi_divs_tot-1 + j*phi_divs_tot + (k+1)*phi_divs_tot*(r_minor_divs+1))  + " " + str(i_cyl + j*phi_divs_tot + (k+1)*phi_divs_tot*(r_minor_divs+1)) + " " + str(i_cyl + j*phi_divs_tot + k*phi_divs_tot*(r_minor_divs+1)) + " " + str(i_cyl + phi_divs_tot-1 + (j+1)*phi_divs_tot + k*phi_divs_tot*(r_minor_divs+1)) + " " + str(i_cyl + phi_divs_tot-1 + (j+1)*phi_divs_tot + (k+1)*phi_divs_tot*(r_minor_divs+1)) + " " + str(i_cyl + (j+1)*phi_divs_tot + (k+1)*phi_divs_tot*(r_minor_divs+1))  + " " + str(i_cyl + (j+1)*phi_divs_tot + k*phi_divs_tot*(r_minor_divs+1)) +") (1 1 1) simpleGrading (1 1 1)\n");
            
#Sew torus to cylinder
i_axis = i_cyl + phi_divs_tot*(r_divs-1)*(r_minor_divs+1)

for j in range(1,r_minor_divs):
    for i in range(0,phi_divs_tot-1):
        meshfile.write("    hex (" + str(i_axis - phi_divs_tot*(r_minor_divs+1) + i + j*phi_divs_tot) + " " + str(i_cyl - phi_divs_tot*r_minor_divs + i + (j-1)*phi_divs_tot)  + " " + str(i_cyl - phi_divs_tot*r_minor_divs + i+1 + (j-1)*phi_divs_tot) + " " + str(i_axis - phi_divs_tot*(r_minor_divs+1) + i+1 + j*phi_divs_tot) + " " + str(i_axis - phi_divs_tot*(r_minor_divs+1) + i + (j+1)*phi_divs_tot) + " " + str(i_cyl - phi_divs_tot*r_minor_divs + i + j*phi_divs_tot) + " " + str(i_cyl - phi_divs_tot*r_minor_divs + i+1 + j*phi_divs_tot)  + " " + str(i_axis - phi_divs_tot*(r_minor_divs+1) + i+1 + (j+1)*phi_divs_tot) +") (1 1 1) simpleGrading (1 1 1)\n");

for j in range(1,r_minor_divs):
    meshfile.write("    hex (" + str(i_axis - phi_divs_tot*(r_minor_divs+1) + phi_divs_tot-1 + j*phi_divs_tot) + " " + str(i_cyl - phi_divs_tot*r_minor_divs + phi_divs_tot-1 + (j-1)*phi_divs_tot)  + " " + str(i_cyl - phi_divs_tot*r_minor_divs + (j-1)*phi_divs_tot) + " " + str(i_axis - phi_divs_tot*(r_minor_divs+1) + j*phi_divs_tot) + " " + str(i_axis - phi_divs_tot*(r_minor_divs+1) + phi_divs_tot-1 + (j+1)*phi_divs_tot) + " " + str(i_cyl - phi_divs_tot*r_minor_divs + phi_divs_tot-1 + j*phi_divs_tot) + " " + str(i_cyl - phi_divs_tot*r_minor_divs + j*phi_divs_tot)  + " " + str(i_axis - phi_divs_tot*(r_minor_divs+1) + (j+1)*phi_divs_tot) +") (1 1 1) simpleGrading (1 1 1)\n");

for i in range(0,phi_divs_tot-1):
    meshfile.write("    hex (" + str(i_axis - phi_divs_tot*(r_minor_divs+1) + i) + " " + str(i_oral + i)  + " " + str(i_oral + i+1) + " " + str(i_axis - phi_divs_tot*(r_minor_divs+1) + i+1) + " " + str(i_axis - phi_divs_tot*r_minor_divs + i) + " " + str(i_cyl - phi_divs_tot*r_minor_divs + i) + " " + str(i_cyl - phi_divs_tot*r_minor_divs + i+1)  + " " + str(i_axis - phi_divs_tot*r_minor_divs + i+1) +") (1 1 1) simpleGrading (1 1 1)\n");

meshfile.write("    hex (" + str(i_axis - phi_divs_tot*(r_minor_divs+1) + phi_divs_tot-1) + " " + str(i_oral + phi_divs_tot-1)  + " " + str(i_oral) + " " + str(i_axis - phi_divs_tot*(r_minor_divs+1)) + " " + str(i_axis - phi_divs_tot*r_minor_divs + phi_divs_tot-1) + " " + str(i_cyl - phi_divs_tot*r_minor_divs + phi_divs_tot-1) + " " + str(i_cyl - phi_divs_tot*r_minor_divs)  + " " + str(i_axis - phi_divs_tot*r_minor_divs) +") (1 1 1) simpleGrading (1 1 1)\n");


#Cylinder wedges (z-axis)
for j in range(0,r_minor_divs):
    for i in range(0,phi_divs_tot-1):
        meshfile.write("    hex (" + str(i_axis + j) + " " + str(i_cyl + i + j*phi_divs_tot)  + " " + str(i_cyl + i+1 + j*phi_divs_tot) + " " + str(i_axis + j) + " " + str(i_axis + j+1) + " " + str(i_cyl + i + (j+1)*phi_divs_tot) + " " + str(i_cyl + i+1 + (j+1)*phi_divs_tot)  + " " + str(i_axis + j+1) +") (1 1 1) simpleGrading (1 1 1)\n");

for j in range(0,r_minor_divs):
    meshfile.write("    hex (" + str(i_axis + j) + " " + str(i_cyl + phi_divs_tot-1 + j*phi_divs_tot)  + " " + str(i_cyl + j*phi_divs_tot) + " " + str(i_axis + j) + " " + str(i_axis + j+1) + " " + str(i_cyl + phi_divs_tot-1 + (j+1)*phi_divs_tot) + " " + str(i_cyl + (j+1)*phi_divs_tot)  + " " + str(i_axis + j+1) +") (1 1 1) simpleGrading (1 1 1)\n");

        
#End
meshfile.write(");\n\n");

#~~~~~~~~~~~~~~~~~Set edges~~~~~~~~~~~~~~~~~~~~~#
meshfile.write("edges\n(\n");

meshfile.write(");\n\n");

###############################################################################

#~~~~~~~~~~~~~~~~Set boundary~~~~~~~~~~~~~~~~~~~#
p.write("boundaryField\n{\n");
p.write("    fixedWalls\n    {\n        type            zeroGradient;\n    }\n\n");
p.write("    openBoundaries\n    {\n        type            fixedValue;\n        value           uniform 0;\n    }\n\n");

U.write("boundaryField\n{\n");
U.write("    fixedWalls\n    {\n        type            noSlip;\n    }\n\n");
U.write("    openBoundaries\n    {\n        type            zeroGradient;\n    }\n\n");



meshfile.write("boundary\n(\n");

#~~~~~~~~~~~~~~fixedWalls
meshfile.write("    fixedWalls\n    {\n        type wall;\n        faces\n        (\n");

#Tintinnid shaft
for k in range(0,z_divs):
    for i in range(0,phi_divs_tot-1):
        meshfile.write("            (" + str(i + k*phi_divs_tot*(r_minor_divs+1)) + " " + str(i+1 + k*phi_divs_tot*(r_minor_divs+1)) + " " + str(i+1 + (k+1)*phi_divs_tot*(r_minor_divs+1)) + " " + str(i + (k+1)*phi_divs_tot*(r_minor_divs+1)) + ")\n");

for k in range(0,z_divs):
    meshfile.write("            (" + str(phi_divs_tot-1 + k*phi_divs_tot*(r_minor_divs+1)) + " " + str(k*phi_divs_tot*(r_minor_divs+1)) + " " + str((k+1)*phi_divs_tot*(r_minor_divs+1)) + " " + str(phi_divs_tot-1 + (k+1)*phi_divs_tot*(r_minor_divs+1)) + ")\n");

#Tintinnid top (wedges, main, sew to r=r_oral)
for i in range(0,phi_divs_tot-1):
    meshfile.write("            (" + str(i_axis) + " " + str(i_cyl + i) + " " + str(i_cyl + i+1) + " " + str(i_axis) + ")\n");

meshfile.write("            (" + str(i_axis) + " " + str(i_cyl + phi_divs_tot-1) + " " + str(i_cyl) + " " + str(i_axis) + ")\n");


for k in range(0,r_divs-2):
    for i in range(0,phi_divs_tot-1):
        meshfile.write("            (" + str(i_cyl + i + k*phi_divs_tot*(r_minor_divs+1)) + " " + str(i_cyl + i + (k+1)*phi_divs_tot*(r_minor_divs+1)) + " " + str(i_cyl + i+1 + (k+1)*phi_divs_tot*(r_minor_divs+1)) + " " + str(i_cyl + i+1 + k*phi_divs_tot*(r_minor_divs+1)) + ")\n");

for k in range(0,r_divs-2):
    meshfile.write("            (" + str(i_cyl + phi_divs_tot-1 + k*phi_divs_tot*(r_minor_divs+1)) + " " + str(i_cyl + phi_divs_tot-1 + (k+1)*phi_divs_tot*(r_minor_divs+1)) + " " + str(i_cyl + (k+1)*phi_divs_tot*(r_minor_divs+1)) + " " + str(i_cyl + k*phi_divs_tot*(r_minor_divs+1)) + ")\n");


for i in range(0,phi_divs_tot-1):
    meshfile.write("            (" + str(i_axis - phi_divs_tot*(r_minor_divs+1) + i) + " " + str(i_oral + i) + " " + str(i_oral + i+1) + " " + str(i_axis - phi_divs_tot*(r_minor_divs+1) + i+1) + ")\n");

meshfile.write("            (" + str(i_axis - phi_divs_tot*(r_minor_divs+1) + phi_divs_tot-1) + " " + str(i_oral + phi_divs_tot-1) + " " + str(i_oral) + " " + str(i_axis - phi_divs_tot*(r_minor_divs+1)) + ")\n");
    
meshfile.write("        );\n\n");
meshfile.write("    }\n");

#~~~~~~~~~~~~~~openBoundaries

#Template
#meshfile.write("            (" + str() + " " + str() + " " + str() + " " + str() + ")\n");

meshfile.write("    openBoundaries\n    {\n        type wall;\n        faces\n        (\n");

#Open bottom
for j in range(0,r_minor_divs):
    for i in range(0,phi_divs_tot-1):
        meshfile.write("            (" + str(i + j*phi_divs_tot) + " " + str(i + (j+1)*phi_divs_tot) + " " + str(i+1 + (j+1)*phi_divs_tot) + " " + str(i+1 + j*phi_divs_tot) + ")\n");

for j in range(0,r_minor_divs):
    meshfile.write("            (" + str(phi_divs_tot-1 + j*phi_divs_tot) + " " + str(phi_divs_tot-1 + (j+1)*phi_divs_tot) + " " + str((j+1)*phi_divs_tot) + " " + str(j*phi_divs_tot) + ")\n");

#Open cylinder
for k in range(0,z_divs):
    for i in range(0,phi_divs_tot-1):
        meshfile.write("            (" + str(phi_divs_tot*r_minor_divs + i + k*phi_divs_tot*(r_minor_divs+1)) + " " + str(phi_divs_tot*r_minor_divs + i + (k+1)*phi_divs_tot*(r_minor_divs+1)) + " " + str(phi_divs_tot*r_minor_divs + i+1 + (k+1)*phi_divs_tot*(r_minor_divs+1)) + " " + str(phi_divs_tot*r_minor_divs + i+1 + k*phi_divs_tot*(r_minor_divs+1)) + ")\n");

for k in range(0,z_divs):
    meshfile.write("            (" + str(phi_divs_tot*r_minor_divs + phi_divs_tot-1 + k*phi_divs_tot*(r_minor_divs+1)) + " " + str(phi_divs_tot*r_minor_divs + phi_divs_tot-1 + (k+1)*phi_divs_tot*(r_minor_divs+1)) + " " + str(phi_divs_tot*r_minor_divs + (k+1)*phi_divs_tot*(r_minor_divs+1)) + " " + str(phi_divs_tot*r_minor_divs + k*phi_divs_tot*(r_minor_divs+1)) + ")\n");

#Open torus
for k in range(0,theta_divs):
    for i in range(0,phi_divs_tot-1):
        meshfile.write("            (" + str(i_init - phi_divs_tot + i + k*phi_divs_tot*r_minor_divs) + " " + str(i_init - phi_divs_tot + i + (k+1)*phi_divs_tot*r_minor_divs) + " " + str(i_init - phi_divs_tot + i+1 + (k+1)*phi_divs_tot*r_minor_divs) + " " + str(i_init - phi_divs_tot + i+1 + k*phi_divs_tot*r_minor_divs) + ")\n");

for k in range(0,theta_divs):
    meshfile.write("            (" + str(i_init - phi_divs_tot + phi_divs_tot-1 + k*phi_divs_tot*r_minor_divs) + " " + str(i_init - phi_divs_tot + phi_divs_tot-1 + (k+1)*phi_divs_tot*r_minor_divs) + " " + str(i_init - phi_divs_tot + (k+1)*phi_divs_tot*r_minor_divs) + " " + str(i_init - phi_divs_tot + k*phi_divs_tot*r_minor_divs) + ")\n");

#Open top circle (wedges, main, sew to r=r_oral)

for i in range(0,phi_divs_tot-1):
    meshfile.write("            (" + str(i_axis + r_minor_divs) + " " + str(i_cyl + phi_divs_tot*r_minor_divs + i) + " " + str(i_cyl + phi_divs_tot*r_minor_divs + i+1) + " " + str(i_axis + r_minor_divs) + ")\n");

meshfile.write("            (" + str(i_axis + r_minor_divs) + " " + str(i_cyl + phi_divs_tot*r_minor_divs + phi_divs_tot-1) + " " + str(i_cyl + phi_divs_tot*r_minor_divs) + " " + str(i_axis + r_minor_divs) + ")\n");


for k in range(0,r_divs-2):
    for i in range(0,phi_divs_tot-1):
        meshfile.write("            (" + str(i_cyl + phi_divs_tot*r_minor_divs + i + k*phi_divs_tot*(r_minor_divs+1)) + " " + str(i_cyl + phi_divs_tot*r_minor_divs + i + (k+1)*phi_divs_tot*(r_minor_divs+1)) + " " + str(i_cyl + phi_divs_tot*r_minor_divs + i+1 + (k+1)*phi_divs_tot*(r_minor_divs+1)) + " " + str(i_cyl + phi_divs_tot*r_minor_divs + i+1 + k*phi_divs_tot*(r_minor_divs+1)) + ")\n");

for k in range(0,r_divs-2):
    meshfile.write("            (" + str(i_cyl + phi_divs_tot*r_minor_divs + phi_divs_tot-1 + k*phi_divs_tot*(r_minor_divs+1)) + " " + str(i_cyl + phi_divs_tot*r_minor_divs + phi_divs_tot-1 + (k+1)*phi_divs_tot*(r_minor_divs+1)) + " " + str(i_cyl + phi_divs_tot*r_minor_divs + (k+1)*phi_divs_tot*(r_minor_divs+1)) + " " + str(i_cyl + phi_divs_tot*r_minor_divs + k*phi_divs_tot*(r_minor_divs+1)) + ")\n");


for i in range(0,phi_divs_tot-1):
    meshfile.write("            (" + str(i_axis - phi_divs_tot + i) + " " + str(i_cyl - phi_divs_tot + i) + " " + str(i_cyl - phi_divs_tot + i+1) + " " + str(i_axis - phi_divs_tot + i+1) + ")\n");

meshfile.write("            (" + str(i_axis - 1) + " " + str(i_cyl - 1) + " " + str(i_cyl - phi_divs_tot) + " " + str(i_axis - phi_divs_tot) + ")\n");


meshfile.write("        );\n\n");
meshfile.write("    }\n");


#~~~~~~~~~~~~~~Weird boundaries

#Cilia wedges

r_minor_j = l_down/r_inner_divs/2
for k in range(0,theta_divs-1):
    for i in range(0,phi_divs_tot):
        if i % phi_divs == 0 and k == cilia_thetas[i/phi_divs]-1:
            #Front face
            meshfile.write("    ciliaWedge" + str(i/phi_divs) + "f\n    {\n        type wall;\n        faces\n        (\n");
            meshfile.write("            (" + str(i_oral + i) + " " + str(i_init + i + k*phi_divs_tot*r_minor_divs) + " " + str(i_init + i + (k+1)*phi_divs_tot*r_minor_divs) + " " + str(i_oral + i) + ")\n");
            meshfile.write("        );\n\n");
            meshfile.write("    }\n");
                
            p.write("    ciliaWedge" + str(i/phi_divs) + "f\n    {\n        type            zeroGradient;\n    }\n\n");

            theta_temp = pi/2*(k+1.5)/theta_divs
            phi_temp = 2*pi*i/phi_divs_tot
            U.write("    ciliaWedge" + str(i/phi_divs) + "f\n    {\n        type            fixedValue;\n        value           uniform (" + str(-r_minor_j*cilia_speeds[i/phi_divs]*sin(theta_temp)*cos(phi_temp)) + " " + str(-r_minor_j*cilia_speeds[i/phi_divs]*sin(theta_temp)*sin(phi_temp)) + " " + str(r_minor_j*cilia_speeds[i/phi_divs]*cos(theta_temp)) + ");\n    }\n\n");

            #Back face
            meshfile.write("    ciliaWedge" + str(i/phi_divs) + "back\n    {\n        type wall;\n        faces\n        (\n");
            meshfile.write("            (" + str(i_oral + i+1) + " " + str(i_init + i+1 + k*phi_divs_tot*r_minor_divs) + " " + str(i_init + i+1 + (k+1)*phi_divs_tot*r_minor_divs) + " " + str(i_oral + i+1) + ")\n");
            meshfile.write("        );\n\n");
            meshfile.write("    }\n");
                
            p.write("    ciliaWedge" + str(i/phi_divs) + "back\n    {\n        type            zeroGradient;\n    }\n\n");

            theta_temp = pi/2*(k+1.5)/theta_divs
            phi_temp = 2*pi*(i+1)/phi_divs_tot
            U.write("    ciliaWedge" + str(i/phi_divs) + "back\n    {\n        type            fixedValue;\n        value           uniform (" + str(-r_minor_j*cilia_speeds[i/phi_divs]*sin(theta_temp)*cos(phi_temp)) + " " + str(-r_minor_j*cilia_speeds[i/phi_divs]*sin(theta_temp)*sin(phi_temp)) + " " + str(r_minor_j*cilia_speeds[i/phi_divs]*cos(theta_temp)) + ");\n    }\n\n");

            #Bottom face
            meshfile.write("    ciliaWedge" + str(i/phi_divs) + "bot\n    {\n        type wall;\n        faces\n        (\n");
            meshfile.write("            (" + str(i_oral + i) + " " + str(i_init + i + k*phi_divs_tot*r_minor_divs) + " " + str(i_init + i+1 + k*phi_divs_tot*r_minor_divs) + " " + str(i_oral + i+1) + ")\n");
            meshfile.write("        );\n\n");
            meshfile.write("    }\n");
                
            p.write("    ciliaWedge" + str(i/phi_divs) + "bot\n    {\n        type            zeroGradient;\n    }\n\n");

            theta_temp = pi/2*(k+1)/theta_divs
            phi_temp = 2*pi*(i+0.5)/phi_divs_tot
            U.write("    ciliaWedge" + str(i/phi_divs) + "bot\n    {\n        type            fixedValue;\n        value           uniform (" + str(-r_minor_j*cilia_speeds[i/phi_divs]*sin(theta_temp)*cos(phi_temp)) + " " + str(-r_minor_j*cilia_speeds[i/phi_divs]*sin(theta_temp)*sin(phi_temp)) + " " + str(r_minor_j*cilia_speeds[i/phi_divs]*cos(theta_temp)) + ");\n    }\n\n");

            #Top face
            meshfile.write("    ciliaWedge" + str(i/phi_divs) + "t\n    {\n        type wall;\n        faces\n        (\n");
            meshfile.write("            (" + str(i_oral + i) + " " + str(i_init + i + (k+1)*phi_divs_tot*r_minor_divs) + " " + str(i_init + i+1 + (k+1)*phi_divs_tot*r_minor_divs) + " " + str(i_oral + i+1) + ")\n");
            meshfile.write("        );\n\n");
            meshfile.write("    }\n");
                
            p.write("    ciliaWedge" + str(i/phi_divs) + "t\n    {\n        type            zeroGradient;\n    }\n\n");

            theta_temp = pi/2*(k+2)/theta_divs
            phi_temp = 2*pi*(i+0.5)/phi_divs_tot
            U.write("    ciliaWedge" + str(i/phi_divs) + "t\n    {\n        type            fixedValue;\n        value           uniform (" + str(-r_minor_j*cilia_speeds[i/phi_divs]*sin(theta_temp)*cos(phi_temp)) + " " + str(-r_minor_j*cilia_speeds[i/phi_divs]*sin(theta_temp)*sin(phi_temp)) + " " + str(r_minor_j*cilia_speeds[i/phi_divs]*cos(theta_temp)) + ");\n    }\n\n");


            
#Cilia body

for j in range(0,r_inner_divs-1):
    r_minor_j = (j+1.5)*l_down/r_inner_divs
    for k in range(0,theta_divs-1):
        for i in range(0,phi_divs_tot):
            if i % phi_divs == 0 and k == cilia_thetas[i/phi_divs]-1 and ((j+1)*r_torus/r_minor_divs < l_up or cilia_speeds[i/phi_divs] < 0):
                #Front face
                meshfile.write("    cilia" + str(i/phi_divs) + "_" + str(j) + "f\n    {\n        type wall;\n        faces\n        (\n");
                meshfile.write("            (" + str(i_init + i + j*phi_divs_tot + k*phi_divs_tot*r_minor_divs) + " " + str(i_init + i + (j+1)*phi_divs_tot + k*phi_divs_tot*r_minor_divs) + " " + str(i_init + i + (j+1)*phi_divs_tot + (k+1)*phi_divs_tot*r_minor_divs) + " " + str(i_init + i + j*phi_divs_tot + (k+1)*phi_divs_tot*r_minor_divs) + ")\n");
                meshfile.write("        );\n\n");
                meshfile.write("    }\n");
                
                p.write("    cilia" + str(i/phi_divs) + "_" + str(j) + "f\n    {\n        type            zeroGradient;\n    }\n\n");

                theta_temp = pi/2*(k+1.5)/theta_divs
                phi_temp = 2*pi*i/phi_divs_tot
                U.write("    cilia" + str(i/phi_divs) + "_" + str(j) + "f\n    {\n        type            fixedValue;\n        value           uniform (" + str(-r_minor_j*cilia_speeds[i/phi_divs]*sin(theta_temp)*cos(phi_temp)) + " " + str(-r_minor_j*cilia_speeds[i/phi_divs]*sin(theta_temp)*sin(phi_temp)) + " " + str(r_minor_j*cilia_speeds[i/phi_divs]*cos(theta_temp)) + ");\n    }\n\n");

                #Back face
                meshfile.write("    cilia" + str(i/phi_divs) + "_" + str(j) + "back\n    {\n        type wall;\n        faces\n        (\n");
                meshfile.write("            (" + str(i_init + i+1 + j*phi_divs_tot + k*phi_divs_tot*r_minor_divs) + " " + str(i_init + i+1 + (j+1)*phi_divs_tot + k*phi_divs_tot*r_minor_divs) + " " + str(i_init + i+1 + (j+1)*phi_divs_tot + (k+1)*phi_divs_tot*r_minor_divs) + " " + str(i_init + i+1 + j*phi_divs_tot + (k+1)*phi_divs_tot*r_minor_divs) + ")\n");
                meshfile.write("        );\n\n");
                meshfile.write("    }\n");
                
                p.write("    cilia" + str(i/phi_divs) + "_" + str(j) + "back\n    {\n        type            zeroGradient;\n    }\n\n");

                theta_temp = pi/2*(k+1.5)/theta_divs
                phi_temp = 2*pi*(i+1)/phi_divs_tot
                U.write("    cilia" + str(i/phi_divs) + "_" + str(j) + "back\n    {\n        type            fixedValue;\n        value           uniform (" + str(-r_minor_j*cilia_speeds[i/phi_divs]*sin(theta_temp)*cos(phi_temp)) + " " + str(-r_minor_j*cilia_speeds[i/phi_divs]*sin(theta_temp)*sin(phi_temp)) + " " + str(r_minor_j*cilia_speeds[i/phi_divs]*cos(theta_temp)) + ");\n    }\n\n");

                #Bottom face
                meshfile.write("    cilia" + str(i/phi_divs) + "_" + str(j) + "bot\n    {\n        type wall;\n        faces\n        (\n");
                meshfile.write("            (" + str(i_init + i + j*phi_divs_tot + k*phi_divs_tot*r_minor_divs) + " " + str(i_init + i + (j+1)*phi_divs_tot + k*phi_divs_tot*r_minor_divs) + " " + str(i_init + i+1 + (j+1)*phi_divs_tot + k*phi_divs_tot*r_minor_divs) + " " + str(i_init + i+1 + j*phi_divs_tot + k*phi_divs_tot*r_minor_divs) + ")\n");
                meshfile.write("        );\n\n");
                meshfile.write("    }\n");
                
                p.write("    cilia" + str(i/phi_divs) + "_" + str(j) + "bot\n    {\n        type            zeroGradient;\n    }\n\n");

                theta_temp = pi/2*(k+1)/theta_divs
                phi_temp = 2*pi*(i+0.5)/phi_divs_tot
                U.write("    cilia" + str(i/phi_divs) + "_" + str(j) + "bot\n    {\n        type            fixedValue;\n        value           uniform (" + str(-r_minor_j*cilia_speeds[i/phi_divs]*sin(theta_temp)*cos(phi_temp)) + " " + str(-r_minor_j*cilia_speeds[i/phi_divs]*sin(theta_temp)*sin(phi_temp)) + " " + str(r_minor_j*cilia_speeds[i/phi_divs]*cos(theta_temp)) + ");\n    }\n\n");

                #Top face
                meshfile.write("    cilia" + str(i/phi_divs) + "_" + str(j) + "t\n    {\n        type wall;\n        faces\n        (\n");
                meshfile.write("            (" + str(i_init + i + j*phi_divs_tot + (k+1)*phi_divs_tot*r_minor_divs) + " " + str(i_init + i + (j+1)*phi_divs_tot + (k+1)*phi_divs_tot*r_minor_divs) + " " + str(i_init + i+1 + (j+1)*phi_divs_tot + (k+1)*phi_divs_tot*r_minor_divs) + " " + str(i_init + i+1 + j*phi_divs_tot + (k+1)*phi_divs_tot*r_minor_divs) + ")\n");
                meshfile.write("        );\n\n");
                meshfile.write("    }\n");
                
                p.write("    cilia" + str(i/phi_divs) + "_" + str(j) + "t\n    {\n        type            zeroGradient;\n    }\n\n");

                theta_temp = pi/2*(k+2)/theta_divs
                phi_temp = 2*pi*(i+0.5)/phi_divs_tot
                U.write("    cilia" + str(i/phi_divs) + "_" + str(j) + "t\n    {\n        type            fixedValue;\n        value           uniform (" + str(-r_minor_j*cilia_speeds[i/phi_divs]*sin(theta_temp)*cos(phi_temp)) + " " + str(-r_minor_j*cilia_speeds[i/phi_divs]*sin(theta_temp)*sin(phi_temp)) + " " + str(r_minor_j*cilia_speeds[i/phi_divs]*cos(theta_temp)) + ");\n    }\n\n");

                #Cilia caps
                if j == r_inner_divs-2 or (cilia_speeds[i/phi_divs] >= 0 and (j+2)*r_torus/r_minor_divs >= l_up):
                    meshfile.write("    ciliaCap" + str(i/phi_divs) + "\n    {\n        type wall;\n        faces\n        (\n");
                    meshfile.write("            (" + str(i_init + i + (j+1)*phi_divs_tot + k*phi_divs_tot*r_minor_divs) + " " + str(i_init + i + (j+1)*phi_divs_tot + (k+1)*phi_divs_tot*r_minor_divs) + " " + str(i_init + i+1 + (j+1)*phi_divs_tot + (k+1)*phi_divs_tot*r_minor_divs) + " " + str(i_init + i+1 + (j+1)*phi_divs_tot + k*phi_divs_tot*r_minor_divs) + ")\n");
                    meshfile.write("        );\n\n");
                    meshfile.write("    }\n");
                
                    p.write("    ciliaCap" + str(i/phi_divs) + "\n    {\n        type            zeroGradient;\n    }\n\n");

                    theta_temp = pi/2*(k+1.5)/theta_divs
                    phi_temp = 2*pi*(i+0.5)/phi_divs_tot
                    r_minor_cap = (j+2)*l_down/r_inner_divs
                    U.write("    ciliaCap" + str(i/phi_divs) + "\n    {\n        type            fixedValue;\n        value           uniform (" + str(-r_minor_cap*cilia_speeds[i/phi_divs]*sin(theta_temp)*cos(phi_temp)) + " " + str(-r_minor_cap*cilia_speeds[i/phi_divs]*sin(theta_temp)*sin(phi_temp)) + " " + str(r_minor_cap*cilia_speeds[i/phi_divs]*cos(theta_temp)) + ");\n    }\n\n");



#End

meshfile.write(");\n\n");



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


p.write("}\n");
U.write("}\n");

meshfile.write("mergePatchPairs\n(\n");

meshfile.write(");\n\n");
