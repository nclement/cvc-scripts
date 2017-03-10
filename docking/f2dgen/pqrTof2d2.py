
from QSlimLib import qslimlib
import numpy, sys, os

def decimate(vert, tri, normals, targetTri):
    #print 'decimating', percent
    if len(vert)==0 or len(tri)==0:
        return [], [], []

    model = qslimlib.QSlimModel(verts, tri, norms=normals)

    nverts = model.get_vert_count()
    nfaces = model.get_face_count()
    newverts = numpy.zeros((nverts, 3)).astype('f')
    newfaces = numpy.zeros((nfaces, 3)).astype('i')
    newnorms = numpy.zeros((nverts, 3)).astype('f')

    model.slim_to_target(targetTri)
    numFaces = model.num_valid_faces()
    model.outmodel(newverts, newfaces, outcolors=None, outnorms=newnorms)
    numVertices = max(newfaces.flatten())+1
    # build lookup table of vertices used in faces
    d = {}
    for t in newfaces[:numFaces]:
        d[t[0]] = 1
        d[t[1]] = 1
        d[t[2]] = 1
    vl = {}
    decimVerts = []
    decimFaces = []
    decimNorms = []
    nvert = 0
    vertices=newverts[:numVertices]
    norms=newnorms[:numVertices]
    for t in newfaces[:numFaces]:
        for i in (0,1,2):
                if not vl.has_key(t[i]):
                    vl[t[i]] = nvert
                    decimVerts.append(vertices[t[i]])
                    decimNorms.append(norms[t[i]])
                    nvert += 1
        decimFaces.append( (vl[t[0]], vl[t[1]], vl[t[2]] ) ) 

    return decimVerts, decimFaces, decimNorms


if __name__=='__main__':
    if len(sys.argv) != 4:  # the program name, receptor PQR, and skin thickness
        # stop the program and print an error message
        sys.exit("Usage: pqrTof2d.py <receptor.pqr> <receptor.rawn> <skin_thickness in A>")

    #protsys = sys.argv[1]
    protsys = sys.argv[1]
    meshrawn = sys.argv[2]
    skinThickness = sys.argv[3]

    #if not os.path.exists(protsys+'_r_b.pqr'):
        #raise ValueError, "File " + protsys + '_r_b.pqr not found'
    #receptor = protsys+'_r_b.pqr'
    receptor = sys.argv[1]

    #if not os.path.exists(protsys+'_l_b.pqr'):
        #raise ValueError, "File " + protsys+'_l_b.pqr not found'
    #ligand = protsys+'_l_b.pqr'



    # read molecule
    from MolKit import Read
    mols = Read(receptor)
 
    # get coords and radii
    coords = mols.allAtoms.coords
    rads = mols.allAtoms.radius

    # compute MSMS
#    from mslib import MSMS
#    msms = MSMS(coords=coords, radii=rads)
#    msms.compute(density=6.0)

    # retrieve triangles
#    vf, vi, tri = msms.getTriangles()

    # decimate the mesh
#    verts = vf[:, :3]
#   faces = tri[:, :3]
#    normals = vf[:, 3:6]

    meshfile = open(meshrawn)
    line = meshfile.readline()
    vfnumber = line.split()

    verts = []
    faces = []
    normals = []
  
    for i in range(0, int(vfnumber[0])):
        line = meshfile.readline()
        linesp = line.split()
        points = [float(linesp[0]), float(linesp[1]), float(linesp[2])]
        verts.append(points)
        points = [float(linesp[3]), float(linesp[4]), float(linesp[5])]
        normals.append(points)
	
    for i in range(0, int(vfnumber[1])):
        line = meshfile.readline()
        linesp = line.split()
        points = [int(linesp[0]), int(linesp[1]), int(linesp[2])]
        faces.append(points)
		
		



    # I think 2*nbAtoms triangles corresponds to nbAtoms vertices
    targetTri = 2*len(mols.allAtoms) 
    dv, ft, dn = decimate(verts, faces, normals, targetTri)

    print "decimated to %d tri, %d verts"%(len(ft), len(dv))
    skinCenters = numpy.array(dv) + numpy.array(dn)*float(skinThickness)

    #write a f2d file
    f = open(receptor[:-4]+'_'+skinThickness+'.f2d', 'w')
    #f.write('# file written by pqrTof2d.py\n')
    for i,a in enumerate(mols.allAtoms):
        # write real atoms as 'I'
        res = a.parent
        c = a.coords
        f.write('ATOM   %5d %4s %3s %4s  %9.3f %9.3f  %9.3f %7.3f %5.3f I\n'%(
            i, a.name, res.type, res.number, c[0], c[1], c[2], a.charge, a.radius))

    resNum = int(res.number)+1
    # write skin atoms as 'E'
    i += 1
    for j,c in enumerate(skinCenters):
        #f.write('HETATM %5d %4s %3s %4s  %9.3f %9.3f  %9.3f %7.3f %5.3f E\n'%(
            #i+j, 'O', 'HOH', str(resNum), c[0], c[1], c[2], 0.0, float(skinThickness)))
        f.write('HETATM %5d %4s %3s %4s  %9.3f %9.3f  %9.3f %7.3f %5.3f E\n'%(
            i+j, 'O', 'HOH', str(resNum), c[0], c[1], c[2], 0.0, 1.4))
    f.close()
