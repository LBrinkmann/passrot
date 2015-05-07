from collections import namedtuple
import numpy as np

def check_consistency(scata,scatb):
    if scata.nqabs != scatb.nqabs:
        raise NameError('Files are not identical in nqabs. First = '+str(scata.nqabs)+'. Second = '+str(scatb.nqabs))
    for qabs_a,qabs_b in zip(scata.qabs.abs,scatb.qabs.abs):
        if abs(qabs_a-qabs_b)>0.001:
            raise NameError('Files are not identical in qabs. First = '+str(qabs_a)+'. Second = '+str(qabs_b))
    for index_a,index_b in zip(scata.qabs.index,scatb.qabs.index):
        if index_a != index_b:
            raise NameError('Files are not identical in indexofqabs. First = '+str(index_a)+'. Second = '+str(index_b))
    for nq_a,nq_b in zip(scata.qabs.nq,scatb.qabs.nq):
        if nq_a != nq_b:
            raise NameError('Files are not identical in number of q vecs for qabs='+str(qabs_a)+'. First = '+str(nq_a)+'. Second = '+str(nq_b))
    for q_a,q_b in zip(scata.q,scatb.q):
        if max(np.subtract(q_a,q_b))>0.001:
            raise NameError('Files are not identical in q vecs. First = '+str(q_a)+'. Second = '+str(q_b))


class ProteinRef:
    """Scattering Amplitude in the Protein Reference System"""
    def __init__(self):
        self.nq = 0
        self.q = []
        self.nqabs = 0
        self.qabs = namedtuple('qabs','abs,index,nq')
        self.qabs.abs = []
        self.qabs.index = []
        self.qabs.nq = []
        self.sa = namedtuple('sf','S,S2,S4,ReS2,ImS2')
        self.sa.S = []
        self.sa.S2 = []
        self.sa.S4 = []
        self.sa.ReS2 = []
        self.sa.ImS2 = []

    def read_envelope(self,f_name):
        self.f_xvg = f_name
        with open(self.f_xvg,'r') as f:
            f_in = f.readlines()
        if 'f_in' not in locals():
            raise NameError('File '+self.f_xvg+' not found.')
        else:
            print "Read file "+self.f_xvg+" !"
        qabslast = 999999999999
        nqvec = 1
    	for line1 in f_in:
            line=line1.split()
            if len(line)>0 and line[0]=='#':
                pass
            elif len(line)>0:
                floats=map(float,line[0:5])
                qnow=np.array(floats[0:3])
                qabs=np.sqrt(qnow.dot(qnow))
                if abs(qabs-qabslast)>0.001:
                    if self.nqabs > 0:
                        self.qabs.nq.append(nqvec)
                        self.qabs.index.append(self.qabs.index[self.nqabs-1] + nqvec)
                    else:
                        self.qabs.index.append(0)
                    self.nqabs += 1
                    self.qabs.abs.append(qabs)
                    nqvec = 0
                nqvec += 1
                qabslast=qabs
                self.nq += 1
                self.q.append(qnow)
                self.sa.S.append(complex(*floats[3:5]))
                self.sa.S2.append((complex(*floats[3:5])*complex(*floats[3:5]).conjugate()).real)
                self.sa.S4.append(0)
                self.sa.ReS2.append(0)
                self.sa.ImS2.append(0)
        self.qabs.nq.append(nqvec)
       
    def read_file(self,f_name):
        self.f_xvg = f_name
        with open(self.f_xvg,'r') as f:
            f_in = f.readlines()
        if 'f_in' not in locals():
            raise NameError('File '+self.f_xvg+' not found.')
        else:
            print "Read file "+self.f_xvg+" !"

    	for line1 in f_in:
            line=line1.split()
            if len(line)>0 and line[0]=='#':
                if line[1]=='nabs':
                    self.nqabs = int(line[3])
                elif line[1]=='indexofqabs':
                    f_index = int(line[3])
                elif line[1]=='qabs':
                    f_qabs = float(line[3])
                    self.qabs.abs.append(f_qabs)
                elif line[1]=='nofqvec':
                    f_nqvec = int(line[3])
                    self.qabs.nq.append(f_nqvec)
                    if f_index > 0:
                        self.qabs.index.append(self.qabs.index[f_index-1] + self.qabs.nq[f_index-1])
                    else:
                        self.qabs.index.append(0)
            elif len(line)>0:
                floats=map(float,line[0:9])
                qnow=np.array(floats[0:3])
                self.q.append(qnow)
                self.sa.S.append(complex(*floats[3:5]))
                self.sa.S2.append(floats[5])
                self.sa.S4.append(floats[6])
                self.sa.ReS2.append(floats[7])
                self.sa.ImS2.append(floats[8])
                qabs=np.sqrt(qnow.dot(qnow))
                if abs(qabs-f_qabs)>0.001:
                    raise NameError('qabs in input file does not correspond to absolute q value. qabs_input = %10.5f, qabs = %10.5f' % (qabs,f_qabs))
        self.nq = np.sum(self.qabs.nq)
