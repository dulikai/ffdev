
import numpy as np

#
# 2010.10
# dulikai
# shandong univ.
# Jinan
#

# the implimentation is as :
# want to get A-X-E, the coordinate of E is need...
# imaging ..
# A' --- X' --- E'
# and A' is the origin point of the coordinate
# A'X' = AX, X'E' = XE
# SO align A'X'E' in the x axies, their coordinate are known
# the rotation matrix XX' is also known
# the EE' should also known based on simple principle.

# fh = open("abc.txt", "r")
# fout = open("out.txt", "w")

# line = fh.readline()
# list = line.split()
# cA = np.array([float(list[1]), float(list[2]), float(list[3])])

# line = fh.readline()
# list = line.split()
# cX = np.array([float(list[1]), float(list[2]), float(list[3])])


### read coord file content to atomlist variable
def read_tm_coord(filename="coord"):
	fin = open(filename, "r")
	atomlist = []
	i = 0
	
	line=fin.readline()
	while line!="":
		if line.find("$coord") == -1:
			line=fin.readline()
		else:
			break;
		# print line
		# print "first while\n"
	line=fin.readline()		
	while line!="":
		if line.find("$") == -1:
			item = line.split()
			atom = [float(item[0]), float(item[1]), float(item[2]), item[3]]
			atomlist.append(atom)
			i = i+1
			# print line
		else:
			break
		line=fin.readline()
		
	return atomlist

### read point_charges file content to atomlist variable
def read_tm_pc(filename="point_charges"):
	fin = open(filename, "r")
	pclist = []
	i = 0
	
	line=fin.readline()
	while line!="":
		if line.find("$point_charges") == -1:
			line=fin.readline()
		else:
			break;
		# print line
		# print "first while\n"
	line=fin.readline()		
	while line!="":
		if line.find("$") == -1:
			item = line.split()
			pc = [float(item[0]), float(item[1]), float(item[2]), float(item[3])]
			pclist.append(pc)
			i = i+1
			# print line
		else:
			break
		line=fin.readline()
		
	return pclist


### output to xyz format coord file content (atomlist)
def print_tm2xyz(atomlist, filename="tm.xyz"):
	
	if atomlist == []:
		print "empty list is found. no work were done\n"
		exit()
	fout = open(filename, "w")
	
	natom = len(atomlist)
	
	print >> fout, "%d\n" %natom
	
	for atom in atomlist:
		emt = atom[3]
		coordx = float(atom[0])*0.529177
		coordy = float(atom[1])*0.529177
		coordz = float(atom[2])*0.529177
		
		print >> fout, "%-s %.14f %.14f %.14f" %(emt, coordx, coordy, coordz)
		
	return filename
	

### output to tm format coord file content (atomlist)	
def print_tm2tm(atomlist, filename="coord.tm"):
	
	if atomlist == []:
		print "empty list is found. no work were done\n"
		exit()
	fout = open(filename, "w")
	
	natom = len(atomlist)
	
	print >> fout, "$coord" 
	
	for atom in atomlist:
		emt = atom[3]
		coordx = float(atom[0])
		coordy = float(atom[1])
		coordz = float(atom[2])
		
		print >> fout, "%.14f %.14f %.14f %-s" %(coordx, coordy, coordz, emt)
		
	print >> fout, "$end\n" 
	return filename
	

	
### output to xyz format coord file content (atomlist)
def print_pc2xyz(pclist, filename="pc.xyz"):
	
	if atomlist == []:
		print "empty list is found. no work were done\n"
		exit()
	fout = open(filename, "w")
	
	npc = len(pclist)
	
	print >> fout, "%d\n" %npc
	
	for pc in pclist:
		emt = "X"
		coordx = float(pc[0])*0.529177
		coordy = float(pc[1])*0.529177
		coordz = float(pc[2])*0.529177
		
		print >> fout, "%-s %.14f %.14f %.14f" %(emt, coordx, coordy, coordz)
		
	return filename
	

### output to tm format coord file content (atomlist)	
def print_pc2pc(pclist, filename="point_charges.tm"):
	
	if atomlist == []:
		print "empty list is found. no work were done\n"
		exit()
	fout = open(filename, "w")

	print >> fout, "$point_charges" 
	
	for pc in pclist:
		emt = float(pc[3])
		coordx = float(pc[0])
		coordy = float(pc[1])
		coordz = float(pc[2])
		
		print >> fout, "%.14f %.14f %.14f %.10f" %(coordx, coordy, coordz, emt)
		
	print >> fout, "$end\n" 
	return filename
	

	
	
	
def pto4v(point):
	p4 = np.append(point, 1.0).reshape(4,1)
	
	return p4
	
	
# FUNCTION:
# obtain 4v movable matrix
# Formula: R*A(old) = E(new), so R=A^-1; T*O(old)=zero(new), so T=-O
# use R*T to obtain make the 4v matrix for rot. and trans.
# variable:
# origin is the coordinate of the new orgin in old axis system
# Pa, Pb, Pc ; the three point to define the new axis
def get_rt_mat(origin, Pa, Pb, Pc):
	
	# define the vector for new axis
	# the x y z axis
	vab = np.subtract(Pb, Pa)
	vac = np.subtract(Pc, Pa)
	vz = np.cross(vab, vac)
	# the origin 
	vtrans = np.negative(origin)	
	
	# the normalized new x y z vector in old axis represent
	ex = vab/np.sqrt(np.dot(vab,vab))
	ez = vz/np.sqrt(np.dot(vz,vz))	
	vy = np.cross(vz, ex)
	ey = vy/np.sqrt(np.dot(vy, vy))
	print "the normalized x y z vector:"
	print ex
	print ey
	print ez
	
	# Make the rot. mat.
	matA = np.array([ex,ey,ez])	
	print "mat. A before transpose (xyz in row represent:"
	print matA
	matA = matA.transpose()
	print "mat. A after transpose, the coordinate axis (xyz in col):"
	print matA
	matR = np.linalg.inv(matA)
	print "after inv, the real rotation matrix (in row):"
	print matR
	
	# try to make 4v rot. mat.
	ex = np.append(matR[:,0], 0)
	ey = np.append(matR[:,1], 0)
	ez = np.append(matR[:,2], 0)
	trans = np.array([0.0, 0.0, 0.0, 1.0])	
	matR4v = np.array([ex,ey,ez,trans])	
 	matR4v= matR4v.transpose()
	print "final 4v rot. mat.:"
	print matR4v
	
	# try to make 4v trans. mat.
	matT4v = np.array([[1.0, 0.0, 0.0,vtrans[0]], [0.0, 1.0, 0.0, vtrans[1]], [0.0, 0.0, 1.0, vtrans[2]], [0.0, 0.0, 0.0, 1.0]])
	print "final 4v trans. mat.:"
	print matT4v
	
	# try to get the final 4v mat. ( R * T mat.)
	move4v = np.dot(matR4v, matT4v)
	print "final 4v R*T mat.:"
	print move4v
	
	# for checking by users
	a = np.dot(move4v, pto4v(Pa))
	b = np.dot(move4v, pto4v(Pb))
	c = np.dot(move4v, pto4v(Pc))
	print "the input Pa Pb Pc in new axis (shown in order: new old):"
	print Pa, a.reshape(1,4)
	print Pb, b.reshape(1,4)
	print Pc, c.reshape(1,4)
	
	return move4v

	
def get_new_point(point, move4v):
	pointx = pto4v(point)
	new_pointx = np.dot(move4v,pointx)
	
	point3 = new_pointx[0:3]
	# print "old_p4, new_p4, point3:"
	# print pointx
	# print new_pointx
	# print point3
	
	return point3
	

	
	
def move_point_sets(points, move4v):
	
	new_points = []
	
	if points == []:
		print "empty list is found. no work were done\n"
		exit()
	
	for coord in points:
		point = coord
		item = get_new_point(point, move4v)	
		new_points.append(item)
	
	return new_points
		
		
		
# the file coord new coordinates
def gen_tm_coord(atomlist, move4v):
	points = []
	new_atomlist = []
	for atom in atomlist:
		item = np.array([float(atom[0]), float(atom[1]), float(atom[2])])
		points.append(item)

	new_points = move_point_sets(points, move4v)
	
	for atom, coord in zip(atomlist, new_points):
		new_item = [coord[0], coord[1], coord[2], atom[3]]
		new_atomlist.append(new_item)		
	
	return new_atomlist
	
# the file point_charges new coordinates
def gen_tm_pc(pclist, move4v):
	points = []
	new_pclist = []
	for pc in pclist:
		item = np.array([float(pc[0]), float(pc[1]), float(pc[2])])
		points.append(item)

	new_points = move_point_sets(points, move4v)
	
	for pc, coord in zip(pclist, new_points):
		new_item = [coord[0], coord[1], coord[2], pc[3]]
		new_pclist.append(new_item)		
	
	return new_pclist

	
	
def get_axis_points(atomlist, table):
	
	P = []
	isign = 1
	if len(table) < 3:
		print "fatal error: please provide 3 point to define the coordinate sys. [in order: xyz"
	for i in table:
		if i < 0:
			i = -i
			isign = -1
		i = i-1
		a = np.array([atomlist[i][0]*isign, atomlist[i][1]*isign, atomlist[i][2]*isign])
		P.append(a)
		
	return P
	
	
	
	
# Main Program	

atomlist = read_tm_coord()

table = [41, 9, 37]
P = get_axis_points(atomlist, table)

print P

Pa,Pb,Pc = P
move4v = get_rt_mat(Pa, Pa, Pb, Pc)

# get_new_point(Pa, move4v)	

new_atomlist = gen_tm_coord(atomlist, move4v)

print_tm2xyz(atomlist, "tmold.xyz")
print_tm2xyz(new_atomlist)
print_tm2tm(new_atomlist)

# Point Charges
pclist = read_tm_pc()
new_pclist = gen_tm_pc(pclist, move4v)
print_pc2xyz(pclist, "pcold.xyz")
print_pc2xyz(new_pclist)
print_pc2pc(new_pclist)




def get_one_ep(cA, cX, len_XE=1.0):


	vAX = cX - cA
	len_AX = np.sqrt(vAX.dot(vAX))

	ciX = np.array([len_AX, 0.0, 0.0])
	ciE = np.array([len_AX+len_XE, 0.0, 0.0])

	rot_viXX = cX - ciX
	len_viXX = np.sqrt(np.dot(rot_viXX, rot_viXX))
	rot_viEE = len_AX/(len_AX+len_XE) * rot_viXX

	cE = cA + ciE + rot_viEE 

	# cX_new = ciX + rot_viXX

	return cE
	
	
# print >> fout, cA, cX, cE



# fout.close()
# np.savetxt("a.txt", cE)

# dist = 1.0
# Vxa = Pa - Px
# scaler = np.sqrt(Vxa.dot(Vxa))
# Vxe_z = (scaler*dist-Vxa[0]-Vxa[1])/Vxa[2]

# print Vxe_z, scaler
# Vxe = ([1, 1, Vxe_z])


# Vxe = scaler*Vxe
# print Vxe

# Pe = Vxe + Px

# print Pa,Px, Vxa, Vxe, Pe

	
# fh.close()

