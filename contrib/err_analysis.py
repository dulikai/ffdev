# python
import os
import math

#
# author: Liu Fang
# 2014.11.4
# @ Jinan, Shandong, China
# @ Shandong Univ.
#

print "ENERGY and ERR ploter. version 1.0"

class ePlotPro():
	"""
	process data for further plotting purpose.
	"""
	def __init__(self):
		"""
		init. value 
		"""
		
		self.grid = {'dft_method':{},'bench_method':{},'ngrid':0,'distance':[], 'angle':[], 'value':[], 'energy':[], 'nvalue':0}
		self.cui = {}
		self.cui['gridfile1'] = 'dft4.dat'
		self.cui['gridfile2'] = 'mp2.dat'
		self.cui['errfile']='err-analyze.txt'
		# self.rd_cmd()
		return

	def rd_dft_energy_grid_bsse(self):
		"""
		read in gaussian output file.
		this has been processed by shell.
		This is the bsse calc. analysis
		so, five value per record

		"""
		filename = self.cui['gridfile1']
		nrec = 5
		au2kcal = 627.5 # kcal/mol
		fp = open(filename, 'r')
		line = fp.readline()
		# print line
		ngrid = line.strip()
		while line.strip() != '':
			energy = []
			line = fp.readline()
			method = line.strip()
			for i in range(int(ngrid)):
			
				line = fp.readline()
				value = line.split()
				rec = {'total':value[0], 'bsse_a':value[1],'bsse_b':value[2],'a':value[3], 'b':value[4]}
				rec['bsse_diff'] = (float(rec['total'])-float(rec['bsse_a'])-float(rec['bsse_b']))*au2kcal
				rec['diff'] = (float(rec['total'])-float(rec['a'])-float(rec['b']))*au2kcal
				energy.append(rec)
			self.grid['dft_method'][method]=energy
			
			line = fp.readline()
		# print self.grid['dft_method']
		self.grid['ngrid'] = ngrid	
		fp.close()
		return	

	def rd_benchmark_energy_grid_bsse(self):
		"""
		read the benchmark energy(usually is mp2 or ccsdt)
		in gaussian output file.
		this has been processed by shell.
		This is the bsse calc. analysis
		so, five value per record

		"""
		nrec = 5
		au2kcal = 627.5 # kcal/mol
		energy = []
		filename = self.cui['gridfile2']
		fp = open(filename, 'r')
		line = fp.readline()
		ngrid = line.strip()
		line = fp.readline()
		method = line.strip()
		line = fp.readline()
		for i in range(int(ngrid)):
			value = line.split()
			# for i in len(value):
				# value[i]= value[i].upper().replace('D', 'E'))
			rec = {'total':value[0], 'bsse_a':value[1],'bsse_b':value[2],'a':value[3], 'b':value[4]}
			rec['bsse_diff'] = (float(rec['total'])-float(rec['bsse_a'])-float(rec['bsse_b']))*au2kcal
			rec['diff'] = (float(rec['total'])-float(rec['a'])-float(rec['b']))*au2kcal
			energy.append(rec)
			line = fp.readline()
		self.grid['bench_method']['method']=method
		self.grid['bench_method']['energy']=energy
		# print self.grid['bench_method']
		fp.close()
		return
		
	def cal_dft_bench_err(self):
		"""
		root mean squared error RMSE= sqrt((dft_energy-bench_energy)**2/n)
		mean unsigned/absolute error MUE= abs(dft_energy-bench_energy)/n
		mean signed error MSE= (dft_energy-bench_energy)/n
		"""
		# print '2'
		
		dft_energy = self.grid['dft_method']
		# print dft_energy
		bench_energy = self.grid['bench_method']['energy']
		print '2'
		# print bench_energy
		ngrid = self.grid['ngrid']
		print ngrid

		
		for rec in dft_energy:
			DIF,dif = 0.0,0.0
			ABS,abs = 0.0,0.0
			SQU,squ = 0.0,0.0
			# print rec
			for i in range(int(ngrid)):
				dif=float(dft_energy[rec][i]['bsse_diff'])-float(bench_energy[i]['bsse_diff'])
				# print dif
				if dif < 0:
					abs = -dif
				else:
					abs = dif
				# abs=math.abs(float(dif))
				squ=float(dif)**2
				DIF = DIF+dif
				ABS = ABS+abs
				SQU = SQU+squ
			MSE = DIF/float(ngrid)
			MUE = ABS/float(ngrid)
			RMSE = math.sqrt(SQU/float(ngrid))
			err = {'err':[MSE,MUE,RMSE]}
			dft_energy[rec].append(err)
			# print dft_energy[rec]
			# dft_energy[rec]['err']=[MSE,MUE,RMSE]
			# dft_energy[rec].append(err)
			
			# print dft_energy[rec]
			# rec.append(err)
			# print rec
			self.grid['dft_metohd']=rec
		return
		
	def wrt_dft_bench_err_file(self):
		method = self.grid['dft_method']
		filename = self.cui['errfile']
		fp = open(filename, 'w')
		print >>fp, "%10s%10s%10s%10s" %('','MSE','MUE','RMSE')
		# print method
		for rec in method:
			# print method[rec][-1]['err'][0],method[rec][-1]['err'][1],method[rec][-1]['err'][2]
			print >>fp, "%10s%12.5f%12.5f%12.5f" %(rec,method[rec][-1]['err'][0],method[rec][-1]['err'][1],method[rec][-1]['err'][2])
			# print >>fp, "%5s%12.5f%12.5f%12.5f" %(rec,float(method[rec]['err'][0]),float(method[rec]['err'][1]),float(method[rec]['err'][2]))
		fp.close()
		return
			
		
	

if __name__ == "__main__":		
    dat = ePlotPro()
    # dat.set_info()
    dat.rd_dft_energy_grid_bsse()
    dat.rd_benchmark_energy_grid_bsse()
    dat.cal_dft_bench_err()
    dat.wrt_dft_bench_err_file()
