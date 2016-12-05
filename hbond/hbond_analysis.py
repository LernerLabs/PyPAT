#!/usr/bin/env python

"""
Michael Lerner's hbond analysis

Right now, just handles pasting together ptraj output.
"""

import copy,pprint,os,sys

class Atom:
    def __init__(self,atom_number=None,resi=None,atom_name=None):
        """
        If we are involved in a solvent hydrogen bond, we expect
        that atom_number, resi and atom_name will all be 'solv'
        """
        self.atom_number = atom_number
        self.resi = resi
        self.atom_name = atom_name
    def __str__(self):
        if len(str(self.atom_number)) > 4:
            raise Exception ('Atom Number Too Long: %s'%self.atom_number)
        spaces = (4-len(str(self.atom_number))) * ' '
        number_str = spaces + '(%s)' % self.atom_number
        return "%3s %4s %s" %(self.resi,
                              self.atom_name,
                              number_str)
    __repr__ = __str__
    def __eq__(self,other):
        #print "Comparing atoms"
        #print "resi", (self.resi == other.resi)
        #print "name", (self.atom_name == other.atom_name)
        #print "number", (self.atom_number == other.atom_number)
        return (self.resi == other.resi) and (self.atom_name == other.atom_name) and (self.atom_number == other.atom_number)
    def __ne__(self,other):
        return not (self == other)
    
    
class HBond:

    def init_from_str(self,s):
        '''
        str looks like what is returned by self._atom_str:
         41    O  (659) -->  94    H (1455) -  94    N (1454)
        '''
        d_resi,d_name,d_number,arrow,ah_resi,ah_name,ah_number,dash,a_resi,a_name,a_number = s.split()
        d_number,ah_number,a_number = [int(i.replace('(','').replace(')','')) for i in d_number,ah_number,a_number]
        d_resi,a_resi,ah_resi = [int(i) for i in d_resi,a_resi,ah_resi]
        self.donor = Atom(d_number,d_resi,d_name)
        self.acceptorh = Atom(ah_number,ah_resi,ah_name)
        self.acceptor = Atom(a_number,a_resi,a_name)
    def __init__(self,line):
        '''
        Initialize ourself from a line that looks like this:
                DONOR         ACCEPTORH      ACCEPTOR
          atom# :res@atom   atom# :res@atom atom# :res@atom %occupied  distance       angle              lifetime maxocc
        |  2546 :160@OA23|  1018   :63@HG   1017   :63@OG  |  99.50  2.641 ( 0.10)  20.89 ( 9.75)    100.0 ( 47.0)    147 |@@@@@@@@@@@@@*@@@@@|
        |  2545 :160@OA22|   705   :44@HH22  703   :44@NH2 |  98.51  2.756 (10.09)  17.97 (19.79)     99.0 (127.0)    126 |*@@@@@@@@@@*@@@@@@@|

        The numbers in parens are standard deviations.

        maxocc is the maximum number of consecutive frames where the H-bond exists.

        Here is a note from cheetham (http://amber.scripps.edu/Questions/mail/322.html)::

            The maxocc is the maximum number of consecutive frames that the
            interaction is found in the trajectory (i.e. 39 consecutive frames).
            
            The lifetime is the average time an interaction occurred...
            
            For example, assume that each space below represents 1ps and a star
            
            means it is occupied:
            
                     10        20         30        40        50
                 *****     *****     **********          *****|
            
            The occupancy would be 5 + 5 + 10 + 5 / 50 or 50%
            The maxocc would be 10
            The lifetime would be 5 + 5 + 10 + 5 / 4 = 6.25 ps (assuming 1 ps between
            frames; the time per frame can be specified on the hbond command line)
            
        So, we will need to adjust the lifetimes.  We have 5ps per frame, so our lifetime
        will need to be multiplied by 5.

        Adding things is not finished yet, but it works for
          - occ
          - distance
          - angle
          - graph
        '''

        # count tells us how many frames have been added together.
        self.count = 1
        if line is None:
            self.donor = Atom()
            self.acceptorh = Atom()
            self.acceptor = Atom()
            self.occ = self.dist = self.dist_stdev = self.angle = self.angle_stdev = 0.0
            self.lifetime = self.lifetime_stdev = 0.0
            self.maxocc = 0
            self.graph =  'XXXXXXXXXXXXXXXXXXX'
            self.graph =  '                   '
            self.graph =  '          '
            return
        
        line = line.strip()
        try:
            leading_junk,donor,acceptor,stats,graph,trailing_junk = line.split('|')
        except ValueError:
            print "Could not hbond",line
            raise
        self.donor = self._ptraj_hbond_chunk_to_atom(donor)
        self.acceptorh = self._ptraj_hbond_chunk_to_atom(' '.join(acceptor.split()[:2]))
        self.acceptor = self._ptraj_hbond_chunk_to_atom(' '.join(acceptor.split()[2:]))
        occ,dist = stats.split('(')[0].strip().split()
        dist_stdev = stats.split('(')[1].split(')')[0].strip()
        angle = stats.split(')')[1].split('(')[0].strip()
        angle_stdev = stats.split('(')[2].split(')')[0].strip()
        lifetime = stats.split('(')[-2].split()[-1].strip()
        lifetime_stdev = stats.split('(')[-1].split(')')[0].strip()
        maxocc = stats.split()[-1]
        self.occ,self.dist,self.dist_stdev,self.angle,self.angle_stdev = [float(i) for i in occ,dist,dist_stdev,angle,angle_stdev]
        self.lifetime,self.lifetime_stdev = [float(i)*5 for i in lifetime,lifetime_stdev]
        self.maxocc = int(maxocc)*5
        #print self.occ,self.dist,self.dist_stdev,self.angle,self.angle_stdev,self.lifetime,self.lifetime_stdev,self.maxocc
        #print graph,'->',graph.strip().strip('|')
        #self.graph = graph.strip().strip('|')
        self.graph = graph
    def _ptraj_hbond_chunk_to_atom(self,chunk):
        ''' chunk is something like "  2546 :160@OA23  " '''
        if chunk.strip() in ('solvent acceptor',
                             'solvent donor',
                             '',):
            if 'acceptor' in chunk:
                resn = 'acc'
            elif 'donor' in chunk:
                resn = 'don'
            else:
                resn = ''
            return Atom('solv','slv',resn)
        else:
            return Atom(int(chunk.split(':')[0].strip()),
                        int(chunk.split(':')[1].split('@')[0].strip()),
                        chunk.split(':')[1].split('@')[1].strip(),
                        )
        
    def __str__(self):
        return self._atom_str() + ' ' + self._occ_graph_str()
    def _atom_str(self):
        return "%s --> %s - %s"%(self.donor,
                                 self.acceptorh,
                                 self.acceptor,
                                 )
    def _occ_graph_str(self):
        return "occ:%6.2f(%2s) |%s|"%(self.occ,
                                      self.count,
                                      self.graph,)
    __repr__ = __str__
    def __add__(self,other):
        '''
        bleh.  i used to do things like
        result = copy.deepcopy(self)
        result.occ = (self.occ + other.occ)/2.0
        result.dist = (self.dist + other.dist)/2.0
        result.angle = (self.angle + other.angle)/2.0
        result.graph = self.graph + other.graph

        but, really, i need to add an attribute that
        tells me how long each hbond is for, so that when
        i add the fifth one in, it does not divide by two.

        for now, i will zero out all of the things just to
        make sure people understand that they are not to be
        believed.
        '''
        if type(self) != type(other):
            raise Exception('no can do, hombre')
        sep = '|'
        if other.acceptor.atom_number is None:
            result = copy.deepcopy(self)
            result.dist = result.angle = 0.0
            result.graph = self.graph + sep + other.graph
            result.count = self.count + other.count
        elif self.acceptor.atom_number is None or (self._atom_str() == other._atom_str()):
            result = copy.deepcopy(other)
            result.dist = result.angle = 0.0
            result.graph = self.graph + sep + other.graph
            result.count = self.count + other.count
        else:
            raise Exception('Can only add hbonds with the same donors and acceptors\n%s != %s'%(self._atom_str(),other._atom_str()))
        #
        # Add various parts now
        #
        result.occ = (self.count * self.occ + other.count * other.occ)/(self.count + other.count)
        return result

def hbond_lines(lines):
    reading = False
    for line in lines:
        if line.strip() == '  atom# :res@atom   atom# :res@atom atom# :res@atom %occupied  distance       angle              lifetime maxocc'.strip():
            reading = True
        if not reading or line.strip().startswith('atom') or not line.replace('-','').strip():
            continue
        yield line
def hbonds(f):
    return [HBond(line) for line in hbond_lines(f)]

def test_file_parsing():
    pprint.pprint(hbonds(file('/Users/mglerner/work/Dynamics-DHFR/MD_Files/ptrajtestfiles/hbond_1rx1_ns2.out')))

def is_relevant(hbond,criteria):
    """
    Tells us if a hbond is relevant
    """
    if 'm20' in criteria:
        if ((9 <= hbond.donor.resi <= 24) or (9 <= hbond.acceptor.resi <= 24) or (9 <= hbond.acceptorh.resi <= 24)):
            #print 'm20',hbond.donor.resi,hbond.acceptor.resi,hbond.acceptorh.resi
            return True
    if 'nap' in criteria:
        if (hbond.donor.resi == 160) or (hbond.acceptor.resi == 160) or (hbond.acceptorh.resi == 160):
            #print "nap",hbond.donor.resi,hbond.acceptor.resi,hbond.acceptorh.resi
            return True
    if 'newpocket' in criteria:
        resis = 137,153,155,30,33,111
        if (hbond.donor.resi in resis) or (hbond.acceptor.resi in resis):
            return True
    if 'other' in criteria:
        return not (is_relevant(hbond,'m20') or is_relevant(hbond,'nap'))
    return False
    #return (9 <= hbond.donor.resi <= 24) or (9 <= hbond.acceptor.resi <= 24) or (9 <= hbond.acceptorh.resi <= 24) or (hbond.donor.resi == 160) or (hbond.acceptor.resi == 160) or (hbond.acceptorh.resi == 160)

def print_relevant_hbonds(fname,criterias):
    """
    Reads in one file and prints out the relevant stuff.
    This doesn't do any adding up of hbonds, etc.  It's
    mostly designed to work with my solvent hbonds.
    """
    hbonds = [HBond(line) for line in hbond_lines(file(fname))]

    for criteria in criterias:
        output = []
        for hbond in hbonds:
            if is_relevant(hbond,criteria):
                output.append((hbond.occ,hbond._atom_str()+' '+hbond._occ_graph_str()))
        output.sort()
        output.reverse()
        output = [o[1] for o in output]
        print '\n'.join(output)
    

def sum_hbonds(struct,criteria):
    print struct,criteria

    
    def add_hbonds(all_hbond_names,hbonds):
        """
        for every name in all_hbond_names that is not found in
        hbonds, add a blank hbond.
        """
        for hbond_name in all_hbond_names:
            found = False
            for hbond in hbonds:
                if hbond._atom_str() == hbond_name:
                    found = True
                    break
            if not found:
                hbond = HBond(None)
                hbond.init_from_str(hbond_name)
                hbonds.append(hbond)

    hbond_output_dir = '/Users/mglerner/work/Dynamics-DHFR/MD_Files/HBontOutputFiles'
    all_hbonds = {}
    all_hbond_names = {}
    for i in range(1,11):
        all_hbonds[i] = [HBond(line) for line in hbond_lines(file(os.path.join(hbond_output_dir,'hbond_%s_NS%s.0.out'%(struct,i))))]
        for j in all_hbonds[i]:
            all_hbond_names[j._atom_str()] = None
    for i in all_hbonds:
        add_hbonds(all_hbond_names,all_hbonds[i])

    combined_hbonds = {}
    for hbond in all_hbonds[1]:
        combined_hbonds[hbond._atom_str()] = hbond
    for i in range(2,11):
        for hbond in all_hbonds[i]:
            combined_hbonds[hbond._atom_str()] = combined_hbonds[hbond._atom_str()] + hbond
    output = []
    for k,v in combined_hbonds.iteritems():
        if is_relevant(v,criteria):
            #print k,v._occ_graph_str()
            output.append((v.occ,k+' '+v._occ_graph_str()))
    output.sort()
    output.reverse()
    output = [o[1] for o in output]
    print '\n'.join(output)
    


def test_hbond_constructors_and_overrides():
    hb1 = HBond('|  2546 :160@OA23|  1018   :63@HG   1017   :63@OG  |  99.50  2.641 ( 0.10)  20.89 ( 9.75)    100.0 ( 47.0)    147 |@@@@@@@@@@@@@*@@@@@|')
    print hb1

    hb2 = HBond('|  2545 :160@OA22|   705   :44@HH22  703   :44@NH2 |  98.51  2.756 (10.09)  17.97 (19.79)     99.0 (127.0)    126 |*@@@@@@@@@@*@@@@@@@|')
    print hb2

    #print hb1+hb2

    hb3 = HBond('|  2546 :160@OA23|  1018   :63@HG   1017   :63@OG  |  20.50  2.641 ( 0.10)  20.89 ( 9.75)    20.0 ( 17.0)     47 |------ooo--ooo@@@@@|')
    print hb3

    print hb1+hb3

        
