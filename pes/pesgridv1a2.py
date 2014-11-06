#! python
from matrix import *
from copy import *
import os

# pes Grid, an extended version.
# scan bond distance, angle, dihedral is implemented.
# in this version, X --- Y is the unmovable atoms. this maybe better for pes scan.
# version records:
# ====
# 2012
# VERsion 1.0 alpha 1: 1a1
# This version process one coord file to generate model.
# ====
# UPDATE:
# 2013-01-14 10 a.m.
# version 1.0 alpha 2: 1a2
# this version could read gaussian 09 BSSE fragment type gjf file
# ====
#
# A-X --- Y-M
# dist(XY); angle(AXY)1; angle(XYM)2; dihedral(AXYM); 
# Note: A X Y M coulde be a set of atoms
# 
# if known A-X --- Y-M
# Then generate new Y in the AxY plane, 
# axis x is XB, z is pe axb, y is cross.
# You must provide an xyz file.
class pesgrid():
	""" pes grid data """
	def __init__(self):
		""" initialize data """
		print "Scan and Generate models: version 1.0 alpha 2(1a2)"
		print "A X Y M could be a set of atoms."
		
		# init. geom. file and format is required.
		self.filename = "input"
		self.filetype = "xyz"
		self.ctrfile = "control.txt"
		# jobtype
		self.jobtype = "E"
		# model
		self.model = []
		self.model2 = []
		# molecule data in xyz file is stored.
		# natom, atom, title, etc.
		self.init_mol = {}
		self.modify_mol = {} # modified based on init mol.
		self.work_mol = {}
		self.tmp_mol = {}
		# modify the following lines
		self.plane = [] # coordinates of origin, pa, pb, pc; and their label id
		self.theta1_map = {'min':90, 'max':270, 'stepsize':5, 'grid':[]} # unit degree
		self.theta2_map = {'min':90, 'max':180, 'stepsize':10, 'grid':[]} # unit degree
		self.radical_map = {'min':2.5, 'max':5.5, 'stepsize':0.1, 'grid':[]} # unit A
		self.phi_map = {'min':90, 'max':180, 'stepsize':10, 'grid':[]} # unit degree
		# A-X --- Y-M
		self.aatom = [6]
		self.xatom = [12]
		self.yatom = [14]
		self.matom = [13]		
		# Donnor / Acceptor
		self.donnor = ['1-12']
		self.acceptor = ['13-16']
		
		self.keypoint = []
		
		self.rd_cmd()
		self.unpack()
		self.map = {'radical': self.radical_map,'theta1': self.theta1_map,'theta2': self.theta2_map,'phi': self.phi_map }

	def rd_cmd(self):
		"""
		two kind mode.
		one for stream; one for file data
		"""
		# working directory
		line = raw_input("Enter the woring directory(press enter to use default: [default: .]\n>")
		mydir = line.strip()
		if mydir == "":
			mydir = '.'
		os.chdir(mydir)
		print os.getcwd()
		# CMD LINE
		line = raw_input("Input control file name (simply enter for CUI):\n>")
		str = line.strip()
		if str == "":
			self.rd_cmd_stream()
			self.ctrfile = ""
		else:
			self.ctrfile = str
			self.rd_cmd_control()
		return
	
	def rd_cmd_stream(self):
		"""
		read in control data from keyboard.
		jump here if control file name is EMPTY
		"""
		print "the coordinate file should contain a dimer. (Donnor and Acceptor)"
		# init. geom. file data and format.
		line = raw_input("Enter the inital file of the coordinate \n (xyz/gjf format only):\n> ")
		rec = line.strip().split(".")
		self.filename = "".join(rec[0:-1])
		self.filetype = rec[-1]
		# Donnor/Acceptor data
		line = raw_input("Input DONNOR list. recommend (1-12); (format: 1-2,3-6,7,9):\n>")
		self.donnor = line.strip().split(",")
		line = raw_input("Input ACCEPTOR list. recommend (1-12); (format: 1-2,3-6,7,9):\n>")
		self.acceptor = line.strip().split(",")
		# AXYM number.
		line = raw_input("Input A number (format: 1,2,3,4):\n> ")
		self.aatom = line.strip().split(",")
		line = raw_input("Input X number (format: 1,2,3,4):\n> ")
		self.xatom = line.strip().split(",")
		line = raw_input("Input Y number (format: 1,2,3,4):\n> ")
		self.yatom = line.strip().split(",")
		line = raw_input("Input M number (format: 1,2,3,4):\n> ")
		self.matom = line.strip().split(",")
		# Grid map
		self.__rd_cmd_stream_map()		
		# Job type.
		self.__rd_cmd_stream_jobtype()		
		return
	def __rd_cmd_stream_jobtype(self):
		"""
		read jobtype
		"""
		print "@ Job type selection @ "
		print "radical(A);\n theta1(B);\n theta2(C);\n phi(D);\n radical+theta1(E)\n"
		line = raw_input("select from A B C D E:\n>")
		str = line.strip()
		self.jobtype = str
		return	
	def __rd_cmd_stream_map(self):
		"""
		read map data.
		"""
		print "@ Grid Map @ "
		print "format: min,max,stepsize. "
		# Distance/Radical Map
		line = raw_input("Radical or Distance Map (i.e. 2.5,5.5,0.1):\n>")
		rec = line.strip().split(",")
		if len(rec) < 3:
			print "Omit this record. No Radical Map information"
		else:
			min, max, stepsize = float(rec[0]), float(rec[1]), float(rec[2])
			self.radical_map = {'min':min, 'max':max, 'stepsize':stepsize, 'grid':[]} # unit A
		# Theta1 Map	
		line = raw_input("Theta1 Map (i.e. 90,180,10):\n>")
		rec = line.strip().split(",")
		if len(rec) < 3:
			print "Omit this record. No Theta1 Map information"
		else:
			min, max, stepsize = float(rec[0]), float(rec[1]), float(rec[2])
			self.theta1_map = {'min':min, 'max':max, 'stepsize':stepsize, 'grid':[]} # unit degree
		# Theta2 Map
		line = raw_input("Theta2 Map (i.e. 90,180,10):\n>")
		rec = line.strip().split(",")
		if len(rec) < 3:
			print "Omit this record. No Theta2 Map information"
		else:
			min, max, stepsize = float(rec[0]), float(rec[1]), float(rec[2])
			self.theta2_map = {'min':min, 'max':max, 'stepsize':stepsize, 'grid':[]} # unit degree
		# PHI Map
		line = raw_input("PHI Map (i.e. 90,180,10):\n>")
		rec = line.strip().split(",")
		if len(rec) < 3:
			print "Omit this record. No PHI Map information"
		else:
			min, max, stepsize = float(rec[0]), float(rec[1]), float(rec[2])
			self.phi_map = {'min':min, 'max':max, 'stepsize':stepsize, 'grid':[]} # unit degree
		return
		
	# initialize from control.txt file, if it exist.
	def rd_cmd_control(self):
		""" read control file """
		ctrfile = self.ctrfile
		fp = open(ctrfile, 'r')
		# filename
		line = fp.readline()
		line = fp.readline()
		rec = line.strip().split(".")
		self.filename = "".join(rec[0:-1])
		self.filetype = rec[-1]
		# donnor
		line = fp.readline()
		line = fp.readline()
		self.donnor = line.strip().split(",")
		# acceptor
		line = fp.readline()
		line = fp.readline()
		self.acceptor = line.strip().split(",")
		# AXYM
		line = fp.readline()
		line = fp.readline()
		self.aatom = line.strip().split(",")
		line = fp.readline()
		self.xatom = line.strip().split(",")
		line = fp.readline()
		self.yatom = line.strip().split(",")
		line = fp.readline()
		self.matom = line.strip().split(",")
		# Grid Map
		self.__rd_cmd_control_map(fp)
		# Job type
		line = fp.readline()
		line = fp.readline()
		self.jobtype = line.strip()
		fp.close()
		return
	def __rd_cmd_control_map(self, fp):
		"""
		read control map
		"""
		line = fp.readline()
		# read map data
		# Distance/Radical Map
		line = fp.readline()
		rec = line.strip().split(",")
		if len(rec) < 3:
			print "Omit this record. No Radical Map information"
		else:
			min, max, stepsize = float(rec[0]), float(rec[1]), float(rec[2])
			self.radical_map = {'min':min, 'max':max, 'stepsize':stepsize, 'grid':[]} # unit A
		# Theta1 Map	
		line = fp.readline()
		rec = line.strip().split(",")
		if len(rec) < 3:
			print "Omit this record. No Theta1 Map information"
		else:
			min, max, stepsize = float(rec[0]), float(rec[1]), float(rec[2])
			self.theta1_map = {'min':min, 'max':max, 'stepsize':stepsize, 'grid':[]} # unit degree
		# Theta2 Map
		line = fp.readline()
		rec = line.strip().split(",")
		if len(rec) < 3:
			print "Omit this record. No Theta2 Map information"
		else:
			min, max, stepsize = float(rec[0]), float(rec[1]), float(rec[2])
			self.theta2_map = {'min':min, 'max':max, 'stepsize':stepsize, 'grid':[]} # unit degree
		# PHI Map
		line = fp.readline()
		rec = line.strip().split(",")
		if len(rec) < 3:
			print "Omit this record. No PHI Map information"
		else:
			min, max, stepsize = float(rec[0]), float(rec[1]), float(rec[2])
			self.phi_map = {'min':min, 'max':max, 'stepsize':stepsize, 'grid':[]} # unit degree
		return
	
	def __unpack_list(self,mylist):
		""" unpack the list like 1-5 to 1,2,3,4,5 automatically """
		tmplist = []
		for emt in mylist:
			tmp = str(emt).split('-')
			if len(tmp) == 1:
				tmplist.append(int(tmp[0]))
			elif len(tmp) == 2:
				for i in range(int(tmp[0]), int(tmp[1])+1, 1):
					tmplist.append(int(i))		
			else:
				print "ERROR OCCUR IN donnor/acceptor list"
		return	tmplist	
	def __str2int_xyam(self):
		""" xyam str 2 int """
		xatom = []
		yatom = []
		aatom = []
		matom = []
		for i in self.xatom:
			xatom.append(int(i))
		for i in self.yatom:
			yatom.append(int(i))
		for i in self.aatom:
			aatom.append(int(i))
		for i in self.matom:
			matom.append(int(i))
		self.xatom = xatom
		self.yatom = yatom
		self.aatom = aatom
		self.matom = matom
		return
	def unpack(self):
		""" unpack lists """
		self.donnor = self.__unpack_list(self.donnor)
		self.acceptor = self.__unpack_list(self.acceptor)	
		self.__str2int_xyam()
		return

	def rd_init_geom(self):
		"""
		read in init. coordinate information
		"""
		filetype = self.filetype
		if filetype == "xyz":
			self.rd_xyz()
		elif filetype == "gjf":
			self.rd_gjf()
		else:
			print "only xyz/gjf format is supported"
			exit()
		return
		
	def rd_xyz(self):
		""" file => data grid """
		filename=self.filename+'.'+self.filetype
		fp = open(filename, "r")
		line = fp.readline()
		natom = int(line)
		self.natom = natom
		self.init_mol['natom'] = natom
		self.init_mol['atom'] = []
		self.init_mol['title'] = ""
		# jump one line
		line = fp.readline()
		# rd atoms
		for i in range(natom):
			iflag = -1
			line = fp.readline()
			items = line.split()
			atomname = items[0]
			coord = np.array([ float(items[1]), float(items[2]), float(items[3]) ])
			atom = {'name': atomname, 'coord': coord, 'iflag': iflag}
			self.init_mol['atom'].append(atom)
			
		self.modify_mol = deepcopy(self.init_mol)
		self.work_mol = deepcopy(self.init_mol)
		return 		
		
	def __check_gjf_frg(self, line):
		"""
			check gjf fragment type 03 or 09 version
		"""
		iflag = 0
		if line.find('=') == -1:
			myline = line
			items = myline.split()
			atomname = items[0]
			coord = np.array([ float(items[1]), float(items[2]), float(items[3]) ])
		else:
			myline = line.replace('(',' ').replace(')',' ').replace('=',' ')
			items = myline.split()
			atomname = items[0]
			coord = np.array([ float(items[3]), float(items[4]), float(items[5]) ])
		rec = {'name': atomname, 'coord': coord, 'iflag': iflag}
		return rec
	def rd_gjf(self):
		""" read gjf file to get coordinate information """
		gfile = self.filename+'.'+self.filetype
		fp = open(gfile, 'r')
		# read header
		line = 'STARTER'
		while line.strip() != "":
			line = fp.readline()
		# read title; charge/spin; then empty line
		line = fp.readline()
		line = fp.readline()
		line = fp.readline()		
		# read coordinate geom data
		self.init_mol['atom'] = []
		self.init_mol['title'] = ""
		natom = 0
		# rd atoms
		line = fp.readline()
		while line.strip() != "":
			# iflag = 0
			# items = line.split()
			# atomname = items[0]
			# coord = np.array([ float(items[1]), float(items[2]), float(items[3]) ])
			# rec = {'name': atomname, 'coord': coord, 'iflag': iflag}
			rec = self.__check_gjf_frg(line)
			self.init_mol['atom'].append(rec)
			line = fp.readline()
			natom = natom + 1
			
		self.init_mol['natom'] = natom
		self.modify_mol = deepcopy(self.init_mol)
		self.work_mol = deepcopy(self.init_mol)
		# other data
		return		
		
	def wrt_xyz_once(self, type='init'):
		""" wrt one mol data => file """
		fp = open(type+'.xyz', 'w')
		if type == 'init':
			mol = self.init_mol
		elif type == 'tmp':
			mol = self.tmp_mol
		elif type == 'work':
			mol = self.work_mol
		else:
			print "Nothing done to print:" + type
			exit(0)
		natom = mol['natom']
		atom = mol['atom']
		print >>fp, "%d\n" % natom 
		for i in range(natom):
			record = atom[i]
			atomname = record['name']
			coord = record['coord']
			iflag = record['iflag']			
			print >>fp, "%5s%12.6f%12.6f%12.6f%5d" % (atomname, coord[0], coord[1], coord[2], iflag)
		return		
	def wrt_xyz_model(self, type="1d"):
		""" write more than one mol data ==>file """
		fp = open('model.xyz', 'w')
		if type == '1d':
			model = self.model
		elif type == '2d':
			model = self.model2
		else:
			print "No Other model can be output. Message: wrt_xyz_model"
			exit(0)
		for mol in model:
			natom = mol['natom']
			atom = mol['atom']
			title = mol['title']
			print >>fp, "%d" % (natom)
			print >>fp, "%s" % title
			for i in range(natom):
				record = atom[i]
				atomname = record['name']
				coord = record['coord']
				iflag = record['iflag']			
				print >>fp, "%5s%12.6f%12.6f%12.6f%5d" % (atomname, coord[0], coord[1], coord[2], iflag)
		return
			
		
	def get_cluster_coord(self, varlist):
		""" get the center of cluster atoms """
		nvar = len(varlist)
		atom = self.work_mol['atom']
		new_coord = np.array([0.0, 0.0, 0.0])
		for i in varlist:
			coord = atom[i-1]['coord']
			new_coord = np.add(new_coord, coord)
		new_coord = np.divide(new_coord, nvar)		
		return new_coord

	# define the new coordinate system
	def set_plane(self, type='donnor'):
		""" reset the atom fromm the plane """
		mol = self.work_mol
		atom = mol['atom']
		donnor = ['donnor', 'theta1', 'radical', 'phi']
		# for donnor or acceptor
		if type in set (donnor):
			# origin pa pb pc, a is the central atom. x axis is pa--->pb		
			# pd is another key point value
			pa = self.get_cluster_coord(self.yatom)
			pb = self.get_cluster_coord(self.xatom)
			pc = self.get_cluster_coord(self.aatom)
			pd = self.get_cluster_coord(self.matom)
			origin = self.get_cluster_coord(self.xatom)
			self.plane = [origin, pa, pb, pc, pd, 'donnor']
		elif type == 'acceptor' or type == 'theta2':
			# origin pa pb pc			
			pa = self.get_cluster_coord(self.xatom)
			pb = self.get_cluster_coord(self.yatom)
			pc = self.get_cluster_coord(self.matom)
			pd = self.get_cluster_coord(self.aatom)
			origin = self.get_cluster_coord(self.yatom)
			self.plane = [origin, pa, pb, pc, pd, 'acceptor']
		else:
			print "only acceptor/donnor type is possible."
			exit(0)
		return	
	def modify_plane(self):
		""" check if pa pb pc in a line, if it is, move c about 0.1 A to reset plane """
		origin,pa,pb,pc,pd,cid=self.plane
		vab = np.subtract(pb, pa)
		vac = np.subtract(pc, pa)
		dab = sqrt(np.dot(vab,vab))
		dac = sqrt(np.dot(vac, vac))
		costheta = np.dot(vab, vac)/(dab*dac)
		if costheta + 1.0 < 1.0e-6:
			pc = self.plane[3]
			pc = np.array([pc[0]+0.1, pc[1]+0.1, pc[2]])
			self.plane[3] = pc
			print "Waring: new PLANE were reset, AUTOmatically. Maybe this is not what you want. check initial geometry."
			print "Message: modify_plane"
		else:
			print "Ok, three atoms were not on a line. Good. continuing..."
		return
		
	def unpack_keypoint(self):
		""" set the key point in their related order """
		pa = self.get_cluster_coord(self.aatom)
		pb = self.get_cluster_coord(self.xatom)
		pc = self.get_cluster_coord(self.yatom)
		pd = self.get_cluster_coord(self.matom)
		self.keypoint = [pa,pb,pc,pd]
		return
		
			
# at first, we need to get the rotated and translated geometry.
# at this coordinate system, more manipulations can be easily done.
# define X --- Y as the x axis. then, A or M is another point to form the plane.
# Here, if the X Y A/M is on a line, A/M will automatically adjust about 0.1 A in both x/y direction
	def new_init_geom(self):
		""" comment """
		origin,pa,pb,pc,pd,cid = self.plane
		# transformation matrix
		move4v = get_rt_mat(origin,pa,pb,pc)
		mol = self.work_mol
		natom = mol['natom']
		atom = mol['atom']
		new_atom = []
		for i in range(natom):	
			record = atom[i]
			point = record['coord']
			# need to reshape to 3x1 vector. default is 1x3 vector
			coord = get_new_point(point, move4v).reshape(3)			
			new_record = {'name': record['name'], 'coord': coord, 'iflag': record['iflag']}
			new_atom.append(new_record)
			self.tmp_mol = {'natom':natom, 'atom':new_atom, 'title':''}	
		return
		
	def reset_work_mol(self, molid=""):
		""" reset work_mol use tmp_mol """
		if molid !="":
			self.work_mol = deepcopy(molid)
		else:
			self.work_mol = deepcopy(self.tmp_mol)
		return
		
	def wrt_map(self, type='theta1'):
		""" write output map info """
		map = self.map[type]
		fp = open('map.txt', 'w')
		min = map['min']
		max = map['max']
		stepsize = map['stepsize']
		nstepsize = int(abs(max-min)/abs(stepsize))
		print >>fp, "min: %12.6f\nmax: %12.6f\nstepsize: %12.6f\nnumber of step: %10d" % (min, max, stepsize, nstepsize)
		print >>fp, "%12s%12s%12s%12s%36s" % ('radius', 'theta1', 'theta2', 'phi', 'coordinate value')
		for grid in map['grid']:
			rad = grid['radius']
			theta1 = grid['theta1']
			theta2 = grid['theta2']
			phi = grid['phi']
			coord = grid['coord']
			print >>fp, "%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f"  % (rad, theta1, theta2, phi, coord[0], coord[1], coord[2])
		fp.close()
		return
		
# the grid map is about the definition of the variable in pes scan.	
# work_mol and plane should be reset for each kind of pes scan.
# the self.plane, the range variable should has been set.
	def gen_theta1_map(self):
		""" generate the theta1 angle variables. """
		theta = self.theta1_map
		min = theta['min']
		max = theta['max']
		stepsize = theta['stepsize']
		plane = self.plane
		c2 = plane[2]
		c3 = plane[3]
		v23 = np.subtract(c3,c2)
		r0 = sqrt(np.dot(v23,v23))
		nstepsize = int(abs(max-min)/abs(stepsize))
		for i in range(nstepsize+1):
			rad = r0
			angle = min+i*stepsize
			th = 180.0-angle				
			thr = pi*(th/180.0)
			x = r0*cos(thr)
			y = r0*sin(thr)
			zd = 0.0
			coord = np.array([x, y, zd])
			grid = {'index': [i], 'coord': coord, 'radius':rad, 'theta1':angle, 'theta2':-1, 'phi':-1}
			theta['grid'].append(grid)
		self.theta1_map = theta
		# since acceptor/donnor list is started from 1, so default 0 starter should consider.
		for i in self.acceptor:
			k = i - 1
			iflag = 0
			self.work_mol['atom'][k]['iflag'] = iflag
		for i in self.donnor:
			k = i - 1
			iflag = 1
			if (i in set(self.xatom)):
				iflag = 0
			self.work_mol['atom'][k]['iflag'] = iflag
		return	
	def gen_theta2_map(self):
		""" generate the theta angle 2 variables. """
		theta = self.theta2_map
		min = theta['min']
		max = theta['max']
		stepsize = theta['stepsize']
		plane = self.plane
		c2 = plane[2]
		c3 = plane[3]
		v23 = np.subtract(c3,c2)
		r0 = sqrt(np.dot(v23,v23))
		nstepsize = int(abs(max-min)/abs(stepsize))
		for i in range(nstepsize+1):
			rad = r0
			angle = min+i*stepsize
			th = 180.0-angle				
			thr = pi*(th/180.0)
			x = r0*cos(thr)
			y = r0*sin(thr)
			zd = 0.0
			coord = np.array([x, y, zd])
			grid = {'index': [i], 'coord': coord, 'radius':rad, 'theta1':-1, 'theta2':angle, 'phi':-1}
			theta['grid'].append(grid)
		self.theta2_map = theta
		# since acceptor/donnor list is started from 1, so default 0 starter should consider.
		for i in self.donnor:
			k = i - 1
			iflag = 0
			self.work_mol['atom'][k]['iflag'] = iflag
		for i in self.acceptor:
			k = i - 1
			iflag = 1
			if (i in set(self.yatom)):
				iflag = 0
			self.work_mol['atom'][k]['iflag'] = iflag
		
		return	
	# suppose X --- Y is on the axis x
	def gen_radical_map(self):
		""" generate the radical grid map """
		radical = self.radical_map
		min = radical['min']
		max = radical['max']
		stepsize = radical['stepsize']
		cX = self.keypoint[1]
		cY = self.keypoint[2]
		nstepsize = int(abs(max-min)/abs(stepsize))
		# suppose x --- y is on the x axis.
		if cX[0] - cY[0] < 0:
			cid = 'donnor'  # origin is on the acceptor atoms, move donnor
			iflag_a = 0
			iflag_d = 1
		else:
			cid = 'acceptor'	# origin is on the donnor atoms, move acceptor		
			iflag_a = 1
			iflag_d = 0
		for i in range(nstepsize+1):
			rad = (min + i * stepsize) 
			x = -rad 	 # move the atom not at origin <-- X --- Y(O)
			y = 0.0
			z = 0.0
			coord = np.array([x, y, z])
			grid = {'index': [i], 'coord': coord, 'radius':rad, 'theta1':-1, 'theta2':-1, 'phi':-1}
			radical['grid'].append(grid)
		self.radical_map = radical

		for i in self.donnor:
			k = i - 1
			self.work_mol['atom'][k]['iflag'] = iflag_d
		for i in self.acceptor:
			k = i - 1
			self.work_mol['atom'][k]['iflag'] = iflag_a
		return
	def gen_phi_map(self):
		""" 
			the dihedral phi, is A-X --- Y-M dihedral.
			maybe this is interesting to study
			SINCE in this process, we maintain the theta1 theta2, and the distance between X/Y
			SO the dihedral is possible to solve exactly.
			THIS is: move the atom A/M(refer as Z) to definited position.
			DEFINE Z(x,y,z)
			IN our coordinate system, YXA or XYM. 
			TO satisfy the above constraint, x is known to be fixed; y/z is movable.
			THEY(x,y) satisfy: z^2+y^2 = x^2; z/y = tan[theta()] KNOWN x is constant.
		"""
		phi = self.phi_map
		min = phi['min']
		max = phi['max']
		stepsize = phi['stepsize']
		cX = self.keypoint[1]
		cY = self.keypoint[2]
		
		# suppose x --- y is on the x axis.
		if cX[0] - cY[0] < 0:
			cid = 'donnor'  # origin is on the acceptor atoms, move donnor   X-A
			cz = self.keypoint[0]
			czz = cX
			iflag_a = 0
			iflag_d = 1
		else:
			cid = 'acceptor'	# origin is on the donnor atoms, move acceptor	Y-M
			cz = self.keypoint[3]
			czz = cY
			iflag_a = 1
			iflag_d = 0
		nstepsize = int(abs(max-min)/abs(stepsize))		
		v23 = np.subtract(cz,czz)
		r0 = sqrt(np.dot(v23,v23))
		for i in range(nstepsize+1):
			rad = r0
			angle = min+i*stepsize
			thr = angle/180.0*pi
			x = cz[0]  # not changed
			y = rad * cos(thr)
			zd = rad * sin(thr)
			coord = np.array([x, y, zd])
			grid = {'index': [i], 'coord': coord, 'radius':rad, 'theta1':-1, 'theta2':-1, 'phi':angle}
			phi['grid'].append(grid)
		self.phi_map = phi
		for i in self.donnor:
			k = i - 1
			self.work_mol['atom'][k]['iflag'] = iflag_d
		for i in self.acceptor:
			k = i - 1
			self.work_mol['atom'][k]['iflag'] = iflag_a
		for i in self.xatom:
			k = i - 1
			self.work_mol['atom'][k]['iflag'] = 0
		for i in self.yatom:
			k = i - 1
			self.work_mol['atom'][k]['iflag'] = 0
		return
	
	def gen_map(self, type='theta1'):
		""" generate various map data : theta1 theta2 radical phi"""
		if type == 'theta1':
			self.theta1_map['grid'] = []
			self.gen_theta1_map()
		elif type == 'theta2':
			self.theta2_map['grid'] = []
			self.gen_theta2_map()
		elif type == 'radical':
			self.radical_map['grid'] = []
			self.gen_radical_map()
		elif type == 'phi':
			self.phi_map['grid'] = []
			self.gen_phi_map()
		else:
			print "Unavaible map Type. Exit..."
			exit(0)
		return
		
	# POINT ROTATE around arbitary vector
	# P is the point, A is arbitary 
	# P' = Pcostheta + (A cross P)sintheta + A(AdotP)(1 - costheta)
	# http://www.cppblog.com/wicbnu/archive/2009/08/13/93215.html
	def point_rotate_ar_vector(self, point, vector, theta):
		""" point rotate around arbitary vector by theta degree angle """
		d = sqrt(np.dot(vector, vector))
		v = np.divide(vector, d)
		p1 = np.multiply(point, cos(theta))
		vcap = np.cross(v, point)
		p2 = np.multiply(vcap, sin(theta))
		dap = np.dot(v, point)*(1-cos(theta))
		p3 = np.multiply(v, dap)
		p = np.add(p1, p2)
		p = np.add(p, p3)
		return p			
		
	def __get_axis_angle(self, oldcoord, grid):
		""" get the rotated angle between oldcoord ---> grid;
			the oldcoord x grid vector were also given """
		a = sqrt(np.dot(grid, grid))
		b = sqrt(np.dot(oldcoord, oldcoord))
		costheta = np.dot(grid, oldcoord)/(a*b)
		costheta = round(costheta, 14)
		theta = acos(costheta)
		vdir = np.cross(oldcoord,grid)		
		return theta, vdir				
		
	def __gen_theta1_geom_mol(self, grid):
		""" move one molecule coordinate, ie. one frame """
		#work_mol
		mol = self.work_mol
		natom = mol['natom']
		atom = mol['atom']	
		title = mol['title']
		point = grid['coord']
		index = grid['index']
		title = "%s %s" % (title, index)
		oldcoord = self.keypoint[0]
		theta, vdir = self.__get_axis_angle(oldcoord, point)
		new_atom = []		
		cmin = 1.0e-4
		for i in range(natom):
			record = atom[i]
			if record['iflag'] != 0 and theta > cmin:
				# A-X --- Y-M, the M and X vector
				coord = self.point_rotate_ar_vector(record['coord'], vdir, theta)
			else:
				coord = record['coord']
			tmp_record = {'name': record['name'], 'coord': coord, 'iflag': record['iflag'], 'index': grid['index'] }
			new_atom.append(tmp_record)
		self.tmp_mol = {'natom':natom, 'atom':new_atom, 'title':title}	
		return
	def __gen_theta2_geom_mol(self, grid):
		""" move one molecule coordinate, ie. one frame """
		#work_mol
		mol = self.work_mol
		natom = mol['natom']
		atom = mol['atom']	
		title = mol['title']
		point = grid['coord']
		index = grid['index']
		title = "%s %s" % (title, index)
		oldcoord = self.keypoint[3]
		theta, vdir = self.__get_axis_angle(oldcoord, point)
		new_atom = []		
		cmin = 1.0e-4
		for i in range(natom):
			record = atom[i]
			if record['iflag'] != 0 and theta > cmin:
				# A-X --- Y-M, the M and X vector
				coord = self.point_rotate_ar_vector(record['coord'], vdir, theta)
			else:
				coord = record['coord']
			tmp_record = {'name': record['name'], 'coord': coord, 'iflag': record['iflag'], 'index': grid['index'] }
			new_atom.append(tmp_record)
		self.tmp_mol = {'natom':natom, 'atom':new_atom, 'title':title}	
		return			
	def __gen_radical_geom_mol(self, grid):
		""" move one molecule coordinate, ie. one frame """
		#work_mol
		mol = self.work_mol
		natom = mol['natom']
		atom = mol['atom']	
		title = mol['title']
		point = grid['coord']
		index = grid['index']
		title = "%s %s" % (title, index)
		cX = self.keypoint[1]
		cY = self.keypoint[2]
		oldcoord = cY
		if cX[0] - cY[0] < 0.0:
			oldcoord = cX
		vec_move = np.subtract(point, oldcoord)		
		new_atom = []		
		for i in range(natom):
			record = atom[i]
			if record['iflag'] != 0:
				# A-X --- Y-M, the M and X vector
				coord = np.add(record['coord'], vec_move)
			else:
				coord = record['coord']
			tmp_record = {'name': record['name'], 'coord': coord, 'iflag': record['iflag'], 'index': grid['index'] }
			new_atom.append(tmp_record)
		self.tmp_mol = {'natom':natom, 'atom':new_atom, 'title':title}	
		return	
	def __gen_phi_geom_mol(self, grid):
		""" move one molecule coordinate, ie. one frame """
		#work_mol
		mol = self.work_mol
		natom = mol['natom']
		atom = mol['atom']	
		title = mol['title']
		point = grid['coord']
		index = grid['index']
		title = "%s %s" % (title, index)
		cX = self.keypoint[1]
		cY = self.keypoint[2]
		oldcoord = self.keypoint[3]
		if cX[0] - cY[0] < 0.0:
			oldcoord = self.keypoint[0]
		xaxis = np.array([oldcoord[0], 0.0, 0.0])
		oldcoord = np.subtract(oldcoord, xaxis)
		point = np.subtract(point, xaxis)
		theta, vdir = self.__get_axis_angle(oldcoord, point)
		new_atom = []		
		cmin = 1.0e-4
		for i in range(natom):
			record = atom[i]
			if record['iflag'] != 0 and theta > cmin:
				# A-X --- Y-M, the M and X vector
				coord = self.point_rotate_ar_vector(record['coord'], vdir, theta)
			else:
				coord = record['coord']
			tmp_record = {'name': record['name'], 'coord': coord, 'iflag': record['iflag'], 'index': grid['index'] }
			new_atom.append(tmp_record)
		self.tmp_mol = {'natom':natom, 'atom':new_atom, 'title':title}	
		return	
		
	def gen_theta1_geom(self):
		""" for a series of model """
		model = []
		map = self.theta1_map
		# each map
		for grid in map['grid']:
			self.__gen_theta1_geom_mol(grid)
			model.append(self.tmp_mol)
		self.model = model
		return 
	def gen_theta2_geom(self):
		""" rotate around theta1 or theta2 angle """
		model = []
		map = self.theta2_map
		# for each map
		for grid in map['grid']:
			self.__gen_theta2_geom_mol(grid)
			model.append(self.tmp_mol)
		self.model = model
		return 
	def gen_radical_geom(self):
		""" move along radical """
		model = []
		map = self.map['radical']
		# for each map
		for grid in map['grid']:
			self.__gen_radical_geom_mol	(grid)
			model.append(self.tmp_mol)
		self.model = model
		return 
	def gen_phi_geom(self):
		""" rotate around phi dihedral """
		model = []
		map = self.phi_map
		# for each map
		for grid in map['grid']:
			self.__gen_phi_geom_mol(grid)
			model.append(self.tmp_mol)
		self.model = model
		return 
	
	def gen_geom(self, type='theta1'):
		""" generate various geometery parameter"""
		# reset model list
		self.model = []
		if type == 'theta1':
			self.gen_theta1_geom()
		elif type == 'theta2':
			self.gen_theta2_geom()
		elif type == 'radical':
			self.gen_radical_geom()
		elif type == 'phi':
			self.gen_phi_geom()
		else:
			print "Unavaible geom Type in gen_geom. Exit..."
			exit(0)
		print "One model Done."
		return
		
	# In fact, it needn't to reset plane, but it doesn't matter to set or not...
	# One should init_new_geom.
	# IN FOR CIRCLE...
	# self.set_plane(type='theta1')
	# self.modify_plane()
	# self.new_init_geom()
	# self.wrt_xyz_once(type='tmp')
	# print "The initial xyz coordinate in new axis were write out in tmp.xyz:"
	# print self.plane
	def scan_theta1_radical(self):
		""" scan theta1 and radical variable """
		model2 = []
		print "LOOP for raical; THEN theta"
		# set working coordinate system.
		self.set_plane(type='radical')
		self.modify_plane()
		self.new_init_geom()
		self.wrt_xyz_once(type='tmp')
		print "The initial xyz coordinate in new axis were write out in tmp.xyz:"
		print self.plane
		self.reset_work_mol()
		self.unpack_keypoint()
		self.gen_map(type='radical')
		self.gen_geom(type='radical')
		mdl = deepcopy(self.model)
		# for each radical/theta1 pairs
		for mol in mdl:
			# reset work mol
			self.reset_work_mol(molid=mol)
			# In fact, it needn't to reset plane, but it doesn't matter to set or not...
			self.unpack_keypoint()
			self.gen_map(type='theta1')			
			self.gen_geom(type='theta1')
			model2.extend(self.model)
		self.model = model2
		return
		
	# if one want to scan a distance, seperated by two point.
	# Here is the case: coordinates of dimmer is supposed to be given.
	def scan_distance(self):
		""" scan in the radical direction between two point. """
		# set working coordinate system.
		self.set_plane(type='donnor')
		self.modify_plane()
		self.new_init_geom()
		self.wrt_xyz_once(type='tmp')
		print "The initial xyz coordinate in new axis has been written out in tmp.xyz:"
		print self.plane
		self.reset_work_mol()
		self.unpack_keypoint()
		self.gen_map(type='radical')
		self.gen_geom(type='radical')
		return
		
	# if one want to scan an angle, by three point.
	# Here is the case: .
	def scan_angle(self, jobtype='donnor'):
		""" 
		scan in the angle of three point. 
		There are two type: donnor, acceptor.
		"""
		# set working coordinate system.
		self.set_plane(type=jobtype)
		self.modify_plane()
		self.new_init_geom()
		self.wrt_xyz_once(type='tmp')
		print "The initial xyz coordinate in new axis were write out in tmp.xyz:"
		# print self.plane
		self.reset_work_mol()
		self.unpack_keypoint()
		if jobtype == 'donnor':
			atype = 'theta1'
		elif jobtype == 'acceptor':
			atype = 'theta2'
		self.gen_map(type=atype)
		self.gen_geom(type=atype)
		return

	# if one want to scan an dihedral, by four point.
	# Here is the case: .
	def scan_dihedral(self, jobtype='donnor'):
		""" 
		scan in the dihedral of four point. 
		There are two type: donnor, acceptor.
		"""
		# set working coordinate system.
		self.set_plane(type=jobtype)
		self.modify_plane()
		self.new_init_geom()
		self.wrt_xyz_once(type='tmp')
		print "The initial xyz coordinate in new axis were write out in tmp.xyz:"
		# print self.plane
		self.reset_work_mol()
		self.unpack_keypoint()
		self.gen_map(type="phi")
		self.gen_geom(type="phi")
		return

	def calculator(self):
		"""
		calculator wrap
		"""
		jobtype = self.jobtype
		print "@ Job type selection @ "
		print "radical(A);\n theta1(B);\n theta2(C);\n phi(D);\n radical+theta1(E)\n"
		print "CURRENT JOBTYPE IS : %10s" % jobtype
		if jobtype.lower() == "a": # radical
			self.scan_distance()
		elif jobtype.lower() == "b":
			self.scan_angle(jobtype='donnor')
		elif jobtype.lower() == "c":
			self.scan_angle(jobtype='acceptor')
		elif jobtype.lower() == "d":
			self.scan_dihedral()
		elif jobtype.lower() == "e":
			self.scan_theta1_radical()
		return
		
# data structure of mol.
# -natom -atom -title
# -atom:list
# --name --coord --iflag --index

# calculation...
pes = pesgrid()
pes.rd_init_geom()
pes.wrt_xyz_once(type='init')
# pes.scan_theta1_radical()
# pes.scan_distance()
pes.calculator()
pes.wrt_xyz_model()









		