import numpy as np
import math
from sympy import *
class Node:
	def __init__(self, elem, coords, restraint=None):
		self.elem = elem
		self.x = coords[0]
		self.y = coords[1]

		self.dx = symbols(f'Delta_{elem}x')
		self.dy = symbols(f'Delta_{elem}y')
		self.dt = symbols(f'theta_{elem}')

		self.fx = symbols(f'F_{elem}x')
		self.fy = symbols(f'F_{elem}y')
		self.moment= symbols(f'M_{elem}')

		self.set_restraint(restraint)
	def add_load(self, load=(0,0,0)):
		self.fx = load[0]
		self.fy = load[1]
		self.moment = load[2]
	def get_coords(self):
		return (self.x, self.y)
	def get_disp(self):
		return np.array([self.dx, self.dy, self.dt])
	def get_name(self):
		return self.elem
	def get_Q(self):
		return np.array([self.fx, self.fy, self.moment])

	def set_vals(self, ans):
		if hasattr(self.fx, 'subs'):
			self.fx = self.fx.subs(ans)
		if hasattr(self.fy, 'subs'):
			self.fy = self.fy.subs(ans)
		if hasattr(self.moment, 'subs'):
			self.moment = self.moment.subs(ans)
		if hasattr(self.dx, 'subs'):
			self.dx = self.dx.subs(ans)
		if hasattr(self.dy, 'subs'):
			self.dy = self.dy.subs(ans)
		if hasattr(self.dt, 'subs'):
			self.dt = self.dt.subs(ans)
	def set_restraint(self, restraint):
		self.restraint = restraint
		#### make restraint types
		if restraint == None:
			self.fx = 0
			self.fy = 0
			self.moment= 0
		elif restraint == 'pin':
			self.dx = 0
			self.dy = 0

			self.moment = 0
		elif restraint == 'fixed':
			self.dx = 0
			self.dy = 0
			self.dt = 0

		elif restraint == 'rolling_y':
			self.dx = 0

			self.fy = 0
			self.moment = 0
		elif restraint == 'rolling_x':
			self.dy = 0

			self.fx = 0
			self.moment = 0
		elif restraint == 'sliding_y':
			self.dt = 0
			self.dx = 0

			self.fy = 0
		elif restraint == 'sliding_x':
			self.dt = 0
			self.dy = 0

			self.fx = 0
		else:
			print('This is not one of the supported support types.\nPlease choose one of the following or manually constrain the node:\n (pin, fixed, roling_y, rolling_x, sliding_y, sliding_x)')



class Beam:
	def __init__(self,beam, nodes, E, A, I):
		self.E = E
		self.A = A
		self.I = I
		self.nodeA = nodes[0]
		self.nodeB = nodes[1]
		self.L = self.get_length()
		self.compute_elem_stiffness()
	def compute_elem_stiffness(self):
		L = self.L
		EAL = self.E * self.A / self.L
		EI = self.E * self.I

		coordsA = self.nodeA.get_coords()
		coordsB = self.nodeB.get_coords()
		if (coordsB[0] - coordsA[0]) == 0:
			if (coordsB[1] - coordsA[1]) >0:
				phi = pi / 2
			elif (coordsB[1] - coordsA[1]) < 0:
				phi = pi / -2
			else:
				raise ValueError('Node A is in the same Location as Node B. Please specify different Nodes.')
		else:
			phi = atan((coordsB[1] - coordsA[1]) / (coordsB[0] - coordsA[0]))

		kprime = np.zeros((6,6)).astype(object)
		kprime[0,:] = [EAL, 0,0, -1*EAL,0,0]
		kprime[1,:] = [0, 12*EI/L**3, 6*EI/L**2, 0, -12*EI/L**3, 6*EI/L**2]
		kprime[2,:] = [0, 6*EI/L**2, 4*EI/L, 0, -6*EI/L**2, 2*EI/L]
		kprime[3,:] = [-1*EAL, 0,0,EAL,0,0]
		kprime[4,:] = [0, -12 * EI/L**3, -6*EI/L**2, 0, 12*EI/L**3, -1*6*EI/L**2]
		kprime[5,:] = [0, 6*EI/L**2, 2*EI/L, 0, -6*EI/L**2, 4*EI/L]

		T = np.zeros((6,6)).astype(object)
		basis = np.array([[cos(phi), sin(phi)], [-1 * sin(phi), cos(phi)] ])
		T[:2,:2] = basis
		T[3:5:,3:5] = basis
		T[2,2] = 1
		T[5,5] = 1
		self.T =T
		self.kprime = kprime
		self.k = T.transpose().dot(kprime).dot(T)
	def local_forces(self):
		Q = np.append(self.nodeA.get_Q(), self.nodeB.get_Q())
		D = np.append(self.nodeA.get_disp(), self.nodeB.get_disp())
		localDelta = self.T.dot(D)
		internalForces = self.kprime.dot(localDelta)
		self.internal = internalForces


	def get_length(self):
		coordsA = self.nodeA.get_coords()
		coordsB = self.nodeB.get_coords()
		L = sqrt((coordsB[1] - coordsA[1])**2 + (coordsB[0] - coordsA[0])**2)
		return L
	def get_k(self):
		return self.k
	def get_node_name(self):
		return (self.nodeA.get_name(), self.nodeB.get_name())


class Structure:
	def __init__(self, units=None):
		self.constructed = False
		self.solved = False
		if units==None:
			self.force = None
			self.dist = None
			self.moment = None
		else:
			self.force = units['f']
			self.dist = units['l']
			self.moment = self.force * self.dist
		self.nodes = {}
		self.beams = {}


	def add_beam(self, name, nodes, E, A, I):
		nodes = list(nodes)
		nodes[0] = self.get_node(str(nodes[0]))
		nodes[1] = self.get_node(str(nodes[1]))

		beam = Beam(name, nodes, E, A, I)
		self.beams[name] = beam
	def add_load(self, nodeName, load=(0,0,0)):
		node = self.get_node(nodeName)
		node.add_load(load)
	def add_node(self, name, coords, restraint=None):
		node = Node(name, coords, restraint)
		self.nodes[name] = node
	def combine_k(self):
		size = 3 * len(self.nodes.keys())
		K = np.zeros((size,size)).astype(object)
		for beamName in self.beams.keys():
			K = self.insert(K, beamName)
		self.K = K
		self.constructed = True

	def get_beam(self, name):
		return self.beams[name]
	def get_node(self, name):
		return self.nodes[name]


	def insert(self, K, beamName):
		beam = self.beams[beamName]
		nodes = beam.get_node_name()
		first  = (int(nodes[0]) - 1) * 3
		second = (int(nodes[1]) - 1) * 3
		Kelem = beam.get_k()

		K[first:first+3, first:first+3] += Kelem[:3, :3]
		K[second:second+3, first:first+3] += Kelem[3:, :3]

		K[second:second+3, second:second+3] += Kelem[3:, 3:]
		K[first:first+3,second:second+3] += Kelem[:3, 3:]
		return K

	def get_K(self):
		if self.constructed:
			return self.K
		else:
			self.combine_k()
			return self.K
	def solve(self, exclude=[]):
		if self.solved:
			raise Exception("This structure has already been solved.")
		size = 3 * len(self.nodes.keys())
		D = np.zeros((size)).astype(object)
		Q = np.zeros((size)).astype(object)
		K = self.get_K().astype(object)
		for nodeNo in self.nodes.keys():
			node = self.get_node(nodeNo)
			nodeNo = (int(nodeNo)-1)*3
			D[nodeNo:nodeNo+3] = node.get_disp()
			Q[nodeNo:nodeNo+3] = node.get_Q()
		eqs = []
		for i, row in enumerate(np.dot(K, D)):
			if Q[i] != row:
				eqs.append(Eq(Q[i], row))
		ans = solve(eqs, exclude=exclude)
		if type(ans) == list:
			if len(ans) == 1:
				ans = ans[0]
			else:
				raise ValueError('More than one solution. Please ensure the system is properly constrained')
		self.solved = True
		print(ans)
		for node in self.nodes.values():
			node.set_vals(ans)
		for beam in self.beams.values():
			beam.local_forces()



