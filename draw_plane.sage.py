

# This file was *autogenerated* from the file draw_plane.sage
from sage.all_cmdline import *   # import sage library

_sage_const_1 = Integer(1); _sage_const_2 = Integer(2); _sage_const_0 = Integer(0); _sage_const_3 = Integer(3); _sage_const_6 = Integer(6)# a program to draw a plane graphs (e.g. generated by plantri)
# author: Manfred Scheucher 2023

from sys import *
from copy import *


def cyclic_rotations(L):
	for i in range(len(L)):
		yield L[i:]+L[:i]


def add_edge(E_plus,u,v):
	if (u,v) > (v,u): (u,v) = (v,u)
	assert((u,v) not in E_plus)
	E_plus.append((u,v))


def triangulate(f,E_plus):
	for v in f:
		if f.count(v) > _sage_const_1 : # if one vertex occurs multiple times
			iv = f.index(v)
			f = f[iv:]+f[:iv] # w.l.o.g. rotate face so that v is first node
			x = f[-_sage_const_1 ] 
			y = f[+_sage_const_1 ]
			add_edge(E_plus,x,y) # add edge from -1 to 1 to remove 0 from face
			triangulate(f[_sage_const_1 :],E_plus)   # and triangulate the rest
			return

	# otherwise the vertices induce an outerplane graph
	# hence we can pick a degree 2 vertex 
	# and triangulate the face with its adjacent edges
	E_induced = {e for e in E_plus if len(set(e)&set(f)) == _sage_const_2 }
	u = [v for v in V if len({e for e in E_induced if v in e}) == _sage_const_2 ][_sage_const_0 ]

	iu = f.index(u)
	f = f[iu:]+f[:iu]
	for v in f[_sage_const_2 :-_sage_const_1 ]:
		add_edge(E_plus,u,v)


ct = _sage_const_0 
for g in open(argv[_sage_const_1 ]): # generate plane graphs with "plantri 5 -a -c1 -p > plane5.ascii"
	g = g.replace("\n","")

	ct += _sage_const_1 
	print("plane graph #",ct,":",g)

	n,rotations = g.split(" ")
	n = int(n)
	rotations = rotations.split(',')

	V = [chr(ord('a')+i) for i in range(n)]
	print("V",len(V),V)

	rot = {V[i]: rotations[i] for i in range(n)}
	print("rot",rot)

	E = []
	for v in V:
		for w in rot[v]: 
			if v<w:
				E.append((v,w))

	F = []
	for v in V:
		for w in rot[v]: 
			face = [w]
			prv = v
			cur = w

			while _sage_const_1 :
				nxt = rot[cur][rot[cur].index(prv)-_sage_const_1 ]
				prv,cur = cur,nxt
				if prv == v and cur == w: break
				face.append(cur)

			if face == min(cyclic_rotations(face)): # any k-face is found k times (each of its vertices can be staring note)
				F.append(face)

	print("E",len(E),E)
	print("F",len(F),F)
	assert(len(V)-len(E)+len(F)==_sage_const_2 ) # check euler formula
	
	
	# build triangulation
	E_plus = deepcopy(E) # copy

	for f in F:
		triangulate(f,E_plus)

	print("E_plus",len(E_plus),E_plus)
	assert(len(E) <= len(E_plus))
	assert(len(E_plus)==_sage_const_3 *len(V)-_sage_const_6 )

	G_plus = Graph(E_plus)
	G_plus.is_planar(set_pos=_sage_const_1 )
	
	G = Graph(E)
	G.set_pos(G_plus.get_pos())

	fp = argv[_sage_const_1 ]+"_"+str(ct)+".png"
	G.plot().save(fp)
	print("wrote to file",fp)
	print()

