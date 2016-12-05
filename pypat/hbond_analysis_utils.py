#!/usr/bin/env python

"""
Michael Lerner's hbond analysis, modified by Steve Spronk

Right now, just handles pasting together ptraj output.
"""

import copy,pprint,os,sys
from scipy import sqrt
from string import ascii_letters
from hbond_tool_utils import *

class Atom:
    def __init__(self, atom_name = None, resi_name = None, resi_num = None):
        """
	Solvent atoms will have atom OW or HW and resi WAT.
        """
        self.atom_name = atom_name
	self.resi_name = resi_name
	self.resi_num = resi_num
    def __eq__(self, other):
        return self.atom_name == other.atom_name and \
	       self.resi_name == other.resi_name and \
	       self.resi_num == other.resi_num
    def __ne__(self, other):
        return not (self == other)
    
class HBond:
    """
    Class to provide a mechanism for handing data contained in the output
    from ptraj
    """

    # ---------------
    # Initializations
    # ---------------


    def __init__(self, line = None, segment_size = 1000, resi_map = None):
        '''
        Initialize ourself from a line that looks like this:
                DONOR         ACCEPTORH      ACCEPTOR
          atom# :res@atom   atom# :res@atom atom# :res@atom %occupied  distance       angle              lifetime maxocc
        |  2546 :160@OA23|  1018   :63@HG   1017   :63@OG  |  99.50  2.641 ( 0.10)  20.89 ( 9.75)    100.0 ( 47.0)    147 |@@@*@@@@@|
        |  2545 :160@OA22|   705   :44@HH22  703   :44@NH2 |  98.51  2.756 (10.09)  17.97 (19.79)     99.0 (127.0)    126 |*@@@*@@@@|
        | solvent donor  |   127    :9@HD21  126    :9@ND2 |   2.00  3.193 ( 0.00)  46.59 ( 0.01)      2.0 (  0.0)      1 |    .    |
        |  5612 :361@OG  |         solvent acceptor        |   2.00  2.915 ( 0.00)  11.31 ( 0.00)      2.0 (  0.0)      1 |    .    |

        The numbers in parentheses are standard deviations.

        Here is a note from cheatham (http://amber.scripps.edu/Questions/mail/322.html)::

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
            
        Adding hbonds only works for some attributes (occupancy, distance, distance
	standard deviation, angle, angle standard deviation, and graph). 

        But because we have split the trajectory into segments, the lifetime and maxocc 
        are not truly a reflection of the H-bonds across the whole trajectory.
        Therefore the manipulation of lifetime and maxocc data are not implemented
        in the current version of hbond_analysis.
        '''

        # num_frames tells us how many frames have been added together.
        self.num_frames = segment_size
        if line is None:
            self.donor = Atom()
            self.acceptorh = Atom()
            self.acceptor = Atom()
            self.occ_pct = self.occ_num = self.dist = self.dist_stdev = self.angle = self.angle_stdev = 0.0
            self.graph =  '          '
            return
        
        line = line.strip()
        try:
            leading_junk, donor, acceptor, stats, graph, trailing_junk = line.split('|')
        except ValueError:
            print "Could not hbond", line
            raise

	# Parse line:

        self.donor = self._ptraj_hbond_chunk_to_atom(donor, resi_map)
        self.acceptorh = self._ptraj_hbond_chunk_to_atom(' '.join(acceptor.split()[:2]), resi_map)
        self.acceptor = self._ptraj_hbond_chunk_to_atom(' '.join(acceptor.split()[2:]), resi_map)

        occ_pct,dist = stats.split('(')[0].strip().split()
        dist_stdev = stats.split('(')[1].split(')')[0].strip()
        angle = stats.split(')')[1].split('(')[0].strip()
        angle_stdev = stats.split('(')[2].split(')')[0].strip()

	# Make necessary type adjustments and calculations

        self.occ_pct,self.dist,self.dist_stdev,self.angle,self.angle_stdev = [float(i) for i in occ_pct,dist,dist_stdev,angle,angle_stdev]
        self.graph = graph
	self.occ_num = int(round(self.occ_pct / 100.0 * self.num_frames))
	if self.occ_num < 2:
	    self.dist_stdev = self.angle_stdev = 0.0
	self.straight_from_ptraj = True

    def _ptraj_hbond_chunk_to_atom(self, chunk, resi_map = None):
        ''' chunk is something like "  2546 :160@OA23  " '''
        if chunk.strip() in ('solvent donor', ''):
            return Atom(atom_name = 'OW', resi_name = 'Wat', resi_num = 999999)
	elif chunk.strip() == 'solvent acceptor':
	    return Atom(atom_name = 'HW', resi_name = 'Wat', resi_num = 999999)
        else:
	    resi_name = chunk.split(':')[1].split('@')[0].strip()
	    if resi_map != None:

		resi_name = resi_map[int(resi_name)]
		try:
		    resi_num = int(resi_name)           # no aa code
		except ValueError:
		    if resi_name[1] in ascii_letters:   # 3-letter aa code
			resi_num = int(resi_name[3:])
		    else:                               # 1-letter aa code
			resi_num = int(resi_name[1:])
	    else:
		resi_num = int(resi_name)
	    atom_name = chunk.split(':')[1].split('@')[1].strip()
            return Atom(atom_name, resi_name, resi_num)

    def init_from_atomstr(self, s, segment_size = 1000):
        '''
        atomstr looks like what is returned by self._atom_str:
            102 NH1--HH11 ... O     88
          Tyr71    OH--HH ... OG1   Asp228
        '''
        a_resi, a_atom, ah_atom, dots, d_atom, d_resi = s.replace('--',' ').split()

	if a_resi == 'Wat':
	    a_resi_num = 999999
	else:
	    try:
		a_resi_num = int(a_resi)        # no aa code
	    except ValueError:
		if a_resi[1] in ascii_letters:  # 3-letter aa code
		    a_resi_num = int(a_resi[3:])
		else:                           # 1-letter aa code
		    a_resi_num = int(a_resi[1:])
	
	if d_resi == 'Wat':               # Same for donor atom
	    d_resi_num = 999999
        else:
            try:
                d_resi_num = int(d_resi)
            except ValueError:
                if d_resi[1] in ascii_letters:
                    d_resi_num = int(d_resi[3:])
                else:
                    d_resi_num = int(d_resi[1:])

        self.donor     = Atom( d_atom, d_resi, d_resi_num)
        self.acceptor  = Atom( a_atom, a_resi, a_resi_num)
        self.acceptorh = Atom(ah_atom, a_resi, a_resi_num) # H is always in same residue as heavy atom it's bonded to

	self.num_frames = segment_size
	self.occ_num = 0
	self.occ_pct = self.dist = self.dist_stdev = self.angle = self.angle_stdev = 0.0
	self.graph =  '          '
	self.straight_from_ptraj = True

    def init_from_str(self, s):
        """
        str looks like what is output by __str__ :
         Lys142   NZ--HZ3 ... OE1  Glu134  27.80( 2500) |--... .*-.|oo---x- - |
	or what is output by _attr_str:
         Lys142   NZ--HZ3 ... OE1  Glu134  27.80( 2500) 2.850(0.17) 29.14(15.66) |--... .*-.|oo---x- - |
        """
        atom_str_len = 34
        hbond_name = s[:atom_str_len].strip()
        hbond_attr = s[atom_str_len:].strip()

        # Take care of Atoms first

        self.init_from_atomstr(hbond_name)

        # Then take care of attributes

        try:
	    attr_list = hbond_attr.split(')')

	    # Attributes from __str__
            self.occ_pct = float(attr_list[0].split('(')[0])
	    self.num_frames = int(attr_list[0].split('(')[1])
	    self.graph = attr_list[-1].strip()[1:-1]  # The [1:-1] takes care of leading and trailing '|'

	    # If present, attributes from _attr_str
	    attr_list = attr_list[1:-1]
	    if attr_list != []:
                self.dist = float(attr_list[0].split('(')[0])
                self.dist_stdev = float(attr_list[0].split('(')[1])
                self.angle = float(attr_list[1].split('(')[0])
                self.angle_stdev = float(attr_list[1].split('(')[1])
        except:
            print "String could not be converted to hbond:", s
            raise
	
	self.occ_num = int(round(self.occ_pct / 100.0 * self.num_frames))
	self.straight_from_ptraj = False

    # ---------------
    # Representations
    # ---------------

    def __str__(self):
        return self._atom_str() + ' ' + self._occ_graph_str()

    def _atom_str(self):
	"""
	Returns the atoms identifying the Hbond as a formatted string.
	Examples:
            102 NH1--HH11 ... O    88 
          Tyr71    OH--HH ... OG1  Asp228
	"""
	spaces = (7 - len(self.acceptor.resi_name)) * ' '
	bond_string = spaces + self.acceptor.resi_name
        acceptor_str = "%s--%s"%(self.acceptor.atom_name,
                                 self.acceptorh.atom_name)
        spaces = (10 - len(acceptor_str)) * ' '
        bond_string += spaces + acceptor_str + " ... "
        spaces = (5 - len(self.donor.atom_name)) * ' '
        bond_string += self.donor.atom_name + spaces
	spaces = (7 - len(self.donor.resi_name)) * ' '
        return bond_string + self.donor.resi_name + spaces

    def _attr_str(self):
	"""
	Returns the attributes in a formatted string.	
	"""
        return "%6.2f(%5s)%6.3f(%4.2f)%6.2f(%5.2f) |%s|"\
                %(self.occ_pct, self.num_frames, self.dist, self.dist_stdev,
                  self.angle, self.angle_stdev, self.graph,)

    def _occ_graph_str(self):
        """
        Returns the occupancy, count, and graph in a formatted string.
        """
        return "%6.2f(%5s) |%s|"%(self.occ_pct, self.num_frames, self.graph)

    __repr__ = __str__

    # ----------
    # Operations
    # ----------

    def __add__(self,other):
        """
        Combines the statistics of two hbonds.  The new number of frames, number of
	occupied frames, occupancy percentage, distance, angle, distance standard 
	deviation, angle standard devation, and graph are all accurately calculated.  

	A note on the standard deviation calculations:  ptraj calculates sigma as the 
	standard deviation (which has N in the denominator).  This is not strictly 
	correct, as this formula only holds true if we know all of the data.  However, 
	we know that our data only contains a sampling from the actual ensemble, so
	it we should use the estimated population standard deviation (S) of the
	statistics, which has N-1 in the denominator of the calculation.  
        """
        if type(self) != type(other):
            raise Exception('Cannot add hbond to non-hbond %s object: %s'%(type(other),other))

	if self._atom_str() != other._atom_str():
	    raise Exception('Can only add hbonds with the same donors and acceptors\n' \
                            '%s != %s'%(self._atom_str(),other._atom_str()))

	result = HBond()

	result.donor = Atom(self.donor.atom_name, self.donor.resi_name, self.donor.resi_num)
	result.acceptor = Atom(self.acceptor.atom_name, self.acceptor.resi_name, self.acceptor.resi_num)
        result.acceptorh = Atom(self.acceptorh.atom_name, self.acceptor.resi_name, self.acceptor.resi_num)
	 
	result.num_frames = self.num_frames + other.num_frames
	sep = '|'
	result.graph = self.graph + sep + other.graph
	result.occ_num = self.occ_num + other.occ_num
	result.occ_pct = result.occ_num * 100.0 / result.num_frames
	result.straight_from_ptraj = False

	if result.occ_num > 0:

            result.dist  = (self.occ_num * self.dist  + other.occ_num * other.dist ) / result.occ_num
            result.angle = (self.occ_num * self.angle + other.occ_num * other.angle) / result.occ_num

            # It's relatively complicated to calculate the new standard deviation.  See my Notebook 3,
            # pp. 72-4 for the derivation.  We must make a distinction on whether or not the data is
	    # straight from the ptraj files, because when we are looking at the data from ptraj 
	    # (straight_from_ptraj = True) the std. dev. is actually sigma as opposed to S, the estimated 
	    # standard deviation of the population.  In practice, these values are close (for relatively 
	    # large N), but I want to be precise with my statistics.

	    if result.occ_num == 1:
	        result.dist_stdev = result.angle_stdev = 0.0

	    else:
	        dist_sumsq = angle_sumsq = 0.0
	        if self.straight_from_ptraj:
	            dist_sumsq  += self.dist_stdev      * self.dist_stdev      * self.occ_num + \
			           self.dist            * self.dist            * self.occ_num
                    angle_sumsq += self.angle_stdev     * self.angle_stdev     * self.occ_num + \
			           self.angle           * self.angle           * self.occ_num
	        else:
                    dist_sumsq  += self.dist_stdev      * self.dist_stdev      * (self.occ_num - 1) + \
                                   self.dist            * self.dist            * self.occ_num
                    angle_sumsq += self.angle_stdev     * self.angle_stdev     * (self.occ_num - 1) + \
                                   self.angle           * self.angle           * self.occ_num
                if other.straight_from_ptraj:                                                                                                  
                    dist_sumsq  += other.dist_stdev     * other.dist_stdev     * other.occ_num + \
                                   other.dist           * other.dist           * other.occ_num
                    angle_sumsq += other.angle_stdev    * other.angle_stdev    * other.occ_num + \
                                   other.angle          * other.angle          * other.occ_num
                else:                                                                     
                    dist_sumsq  += other.dist_stdev     * other.dist_stdev     * (other.occ_num - 1) + \
                                   other.dist           * other.dist           * other.occ_num
                    angle_sumsq += other.angle_stdev    * other.angle_stdev    * (other.occ_num - 1) + \
                                   other.angle          * other.angle          * other.occ_num

	        result.dist_stdev     = sqrt((dist_sumsq  - result.occ_num*result.dist    *result.dist    ) / (result.occ_num - 1))
                result.angle_stdev    = sqrt((angle_sumsq - result.occ_num*result.angle   *result.angle   ) / (result.occ_num - 1))

	#else:
	#    result.dist = result.dist_stdev = result.angle = result.angle_stdev = 0.0

        return result

    def compress_graph(self):
	"""
	Compresses the graph of a trajectory into one half the size.
	Each pair of characters is replaced by a single character
	that is representative of the percentage of occupancy for
	the union of the two segments.  Unfortunately, the actual
	occupancy percentage of the union can not be absolutely
	determined from the two symbols of the graph, so the new
	graph may not be precise.  See my Notebook 3, pp. 78-79
	for a detailed analysis of how I determined how two 
	symbols should be combined.
	"""
	graph_sections = self.graph.split('|')
	new_graph = ''
	for graph_num in range(len(graph_sections)):
	    for i in range(0, 10, 2):
		pair = graph_sections[graph_num][i:i+2]
		if pair[0] == pair[1]:
		    new_graph += pair[0]

		elif pair == ' .' or pair == '. ':
		    new_graph += '.'
                elif pair == ' -' or pair == '- ':
                    new_graph += '.'
                elif pair == ' o' or pair == 'o ':
                    new_graph += '-'
                elif pair == ' x' or pair == 'x ':
                    new_graph += '-'
                elif pair == ' *' or pair == '* ':
                    new_graph += 'o'
                elif pair == ' @' or pair == '@ ':
                    new_graph += 'o'
                elif pair == '.-' or pair == '-.':
                    new_graph += '-'
                elif pair == '.o' or pair == 'o.':
                    new_graph += '-'
                elif pair == '.x' or pair == 'x.':
                    new_graph += 'o'
                elif pair == '.*' or pair == '*.':
                    new_graph += 'o'
                elif pair == '.@' or pair == '@.':
                    new_graph += 'o'
                elif pair == '-o' or pair == 'o-':
                    new_graph += 'o'
                elif pair == '-x' or pair == 'x-':
                    new_graph += 'o'
                elif pair == '-*' or pair == '*-':
                    new_graph += 'o'
                elif pair == '-@' or pair == '@-':
                    new_graph += 'x'
                elif pair == 'ox' or pair == 'xo':
                    new_graph += 'o'
                elif pair == 'o*' or pair == '*o':
                    new_graph += 'x'
                elif pair == 'o@' or pair == '@o':
                    new_graph += 'x'
                elif pair == 'x*' or pair == '*x':
                    new_graph += 'x'
                elif pair == 'x@' or pair == '@x':
                    new_graph += '*'
		elif pair == '*@' or pair == '@*':
                    new_graph += '*'
	    
	    if graph_num % 2 == 1:
		new_graph += '|'

	if new_graph[-1] == '|':
	    self.graph = new_graph[:-1]
	else:
	    self.graph = new_graph

    # ------ End class HBond ----

def hbond_lines(lines):
    reading = False
    for line in lines:
        if line.strip() == '  atom# :res@atom   atom# :res@atom atom# :res@atom %occupied  distance       angle              lifetime maxocc'.strip():
            reading = True
        if not reading or line.strip().startswith('atom') or not line.replace('-','').strip():
            continue
        yield line
def hbonds_from_ptraj(f, segment_size = 1000, resi_map = None):
    return [HBond(line, segment_size, resi_map) for line in hbond_lines(f)]

def is_resinum_of_interest(hbond, criteria = ['all']):
    """
    Tells us if a hbond has a residue number among those we want to view
    """
    if 'all' in criteria:
        return True
    if hbond.donor.resi_num in criteria or hbond.acceptor.resi_num in criteria:
	return True
    else:
	return False

def is_atom_of_interest(hbond, criteria = ['all']):
    """
    Tells us if an hbond has an atom type among those we want to view
    """
    if 'all' in criteria:
	return True
    if 'protein_only' in criteria:
	if hbond.donor.atom_name == 'OW' or hbond.acceptor.atom_name == 'OW':
	    return False
	else:
	    return True
    if 'bb_only' in criteria:
	if hbond.donor.atom_name == 'O' and hbond.acceptor.atom_name == 'N':
	    return True
    if 'not_bb' in criteria:
	if hbond.donor.atom_name != 'O' or hbond.acceptor.atom_name != 'N':
	    return True
    if hbond.donor.atom_name in criteria or \
       hbond.acceptor.atom_name in criteria or \
       hbond.acceptorh.atom_name in criteria:
	return True
    else:
	return False

def combine_hbonds(hbond_files, segment_size = 1000,
		resi_map = None, output_file = None,
		resi_criteria = ['all'], atom_criteria = ['all'],
		occ_thresh = 0.0, occ_graph_only = False,
		hbond_data_dir = None):
    """
    Reads through a set of files that have been output by ptraj and compiles
    all the data.  

    hbond_files:  the hbond_files output from ptraj to be combined.
    segment_size:  the number of frames included in each segment of the 
	trajectory. (default: 1000)
    resi_map:  a dictionary mapping the name of each residue onto the residue
	number.  If 'None,' the residue name will simply be the number. 
	(default: None)
    output_file:  the name of the output file.  If None, the results will be 
        written to stdout. (default: None)
    resi_criteria:  a list containing residue number criteria to include in the
	output. (default: ['all'])
    atom_criteria:  a list containing atom name criteria to include in the 
	output. (default: ['all'])
    occ_thresh:  the minimum occupancy threshold that the hbonds must have 
	to be reported. (default: 0.0)
    occ_graph_only:  if True, only the atom string, occupancy, and graph of 
	each hbond will be written to output. (default: False)
    hbond_data_dir:  the directory that contains the hbond data files.  If
	'None,' the file names will be used without modification, and the
	output will be written to the current directory. (default: None) 
    """

    # Do error checking of file names

    files_to_remove = []
    for each_file in hbond_files:
        if hbond_data_dir != None:
            full_file = os.path.join(hbond_data_dir, each_file)
        else:
            full_file = each_file
        if not os.path.exists(full_file):
            print 'Warning:  File ' + full_file + ' does not exist.\n' + \
                  '  Will be ignored.'
            files_to_remove.append(each_file)
    for each_file in files_to_remove:
        hbond_files.remove(each_file)
    if len(hbond_files) == 0:
        sys.exit('ERROR:  No input files provided.\n')

    # Create list of hbonds in each file, and a master hbond dict

    hbonds_from_file = {}  # {filename: list of hbond objects}
    combined_hbonds = {}   # {hbond string: hbond object}
 
    for each_file in hbond_files:
	if hbond_data_dir != None:
	    hbond_file = os.path.join(hbond_data_dir, each_file)
	else:
	    hbond_file = each_file
	try:
	    hbond_f = file(hbond_file)
	except:
	    sys.exit('ERROR:  Could not open ' + hbond_file + '.\n')
        hbonds_from_file[each_file] = hbonds_from_ptraj(hbond_f, segment_size, resi_map)
        for hbond in hbonds_from_file[each_file]:
            combined_hbonds[hbond._atom_str()] = None

    # Run through the master hbond dict, and find out the missing hbonds
    # in each file.  If any are missing, create an hbond with no occupancy.

    for each_file in hbond_files:
        for hbond_str in combined_hbonds:
            found = False
            for hbond in hbonds_from_file[each_file]:
                if hbond._atom_str() == hbond_str:
                    found = True
                    break
            if not found:
                hbond = HBond()
                hbond.init_from_atomstr(hbond_str, segment_size)
                hbonds_from_file[each_file].append(hbond)

    # Do the addition of the hbonds from each file to create the
    # final combined hbond object.

    for hbond in hbonds_from_file[hbond_files[0]]:
        combined_hbonds[hbond._atom_str()] = hbond
    for each_file in hbond_files[1:]:
        for hbond in hbonds_from_file[each_file]:
            combined_hbonds[hbond._atom_str()] = combined_hbonds[hbond._atom_str()] + hbond
	    
    # Write output to file or stdout

    output = []
    for hbond in combined_hbonds.values():
        if is_resinum_of_interest(hbond, resi_criteria) and \
	   is_atom_of_interest(hbond, atom_criteria) and \
	   hbond.occ_pct > occ_thresh:
            if not occ_graph_only:
		output.append((hbond.occ_pct, hbond._atom_str() + ' ' + hbond._attr_str()))
	    else:
		output.append((hbond.occ_pct, str(hbond)))
    output.sort()
    output.reverse()
    output = [o[1] for o in output]
    output_str = '\n'.join(output)

    if hbond_data_dir == None:
	output_dir = '.'
    else:
	output_dir = hbond_data_dir

    if output_file == None:
	print output_str
    else:
	try:
	    output_file = os.path.join(output_dir, output_file)
	    output_f = file(output_file, 'w')
	except IOError:
	    print 'Warning:  Could not open ' + output_file + '.\n'
	    print output_str
	else:
	    output_f.write(output_str + '\n')
	    output_f.close()

def subset_hbonds(hbond_file, output_file = None,
		resi_criteria = ['all'], atom_criteria = ['all'], 
		occ_thresh = 0.0, occ_graph_only = False, 
		sort = 'occ_pct', compress = False,
		hbond_data_dir = None):
    """
    Following combination of hbonds by combine_hbond(), this function can be
    used to write to stdout or a file only a subset of all the data present.

    hbond_file:  the hbond file with data to be analyzed.
    output_file:  the name of the output file.  If None, the results will be
        written to stdout. (default: None)
    resi_criteria:  a list containing residue number criteria to include in the
        output. (default: ['all'])
    atom_criteria:  a list containing atom name criteria to include in the
        output. (default: ['all'])
    occ_thresh:  the minimum occupancy threshold that the hbonds must have
        to be reported. (default: 0.0)
    occ_graph_only:  if True, only the atom string, occupancy, and graph of
        each hbond will be written to output. (default: False)
    sort:  one of 'occ_pct', 'donor', 'acceptor', 'dist', or 'angle' that
        indicates how to sort the output. (default: occ_pct)
    compress:  if True, the graphs will be compressed by compress_graph().
        (default: False)
    hbond_data_dir:  the directory that contains the hbond data files.  If
        'None,' the file names will be used without modification, and the
        output will be written to the current directory. (default: None)
    """
    # Do error checking of file names.

    if not hbond_file:
	sys.exit('ERROR:  No input file provided.\n')
    if type(hbond_file) is type([]):
	if len(hbond_file) > 1:
	    print 'Warning:  More than 1 input file provided.\n' + \
		  '  Will only use first one: ' + hbond_file[0] 
	hbond_file = hbond_file[0]
    if hbond_data_dir != None:
        full_file = os.path.join(hbond_data_dir, hbond_file)
    else:
        full_file = hbond_file
    try:
	hbond_f = file(full_file)
    except IOError:
	sys.exit('ERROR:  Could not open ' + full_file + '.\n')

    # Create list of hbonds in the input file, check to see if they
    # satisfy the necessary criteria for output.

    hbond_list = []
    for line in hbond_f:
        hbond = HBond()
        hbond.init_from_str(line)
        hbond_list.append(hbond)

    output = []
    for hbond in hbond_list:
        if is_resinum_of_interest(hbond, resi_criteria) and \
	   is_atom_of_interest(hbond, atom_criteria) and \
	   hbond.occ_pct > occ_thresh:
            if compress:
                hbond.compress_graph()

            if occ_graph_only:
                hbond_str = str(hbond)
            else:
                hbond_str = hbond._atom_str() + ' ' + hbond._attr_str()

	    if sort not in 'occ_pct acceptor donor dist angle'.split():
                print 'Warning:  Unknown sorting method: ' + sort + '.\n'  + \
                      '  Will sort by occupancy percentage.'
                sort = 'occ_pct'

            if sort == 'occ_pct':
                output.append((hbond.occ_pct, 
			       hbond.acceptor.resi_num,
			       hbond_str))
            elif sort == 'acceptor':
                output.append((hbond.acceptor.resi_num, 
			       hbond.acceptor.atom_name, 
                               hbond.donor.resi_num,
			       hbond.donor.atom_name,
                               hbond_str))
            elif sort == 'donor':
                output.append((hbond.donor.resi_num, 
			       hbond.donor.atom_name,
                               hbond.acceptor.resi_num,
			       hbond.acceptor.atom_name,
                               hbond_str))
            elif sort == 'dist':
                output.append((hbond.dist, hbond_str))
            else:                          # sort must be 'angle'
                output.append((hbond.angle, hbond_str))

    # Write output

    output.sort()
    if sort == 'occ_pct':
        output.reverse()
    output = [o[-1] for o in output]
    output_str = '\n'.join(output)

    if hbond_data_dir == None:
        output_dir = '.'
    else:
        output_dir = hbond_data_dir

    if output_file == None:
        print output_str
    else:
        try:
            output_file = os.path.join(output_dir, output_file)
            output_f = file(output_file, 'w')
        except IOError:
            print 'Warning:  Could not open ' + output_file + '.\n'
            print output_str
        else:
            output_f.write(output_str + '\n')
            output_f.close()

def compare_hbonds(hbond_files, identifiers = [], output_file = None,
                resi_criteria = ['all'], atom_criteria = ['all'],
                occ_thresh = 0.0, occ_graph_only = False,
                sort = 'occ_diff', compress = False,
                hbond_data_dir = None):
    """
    Following combination of hbonds by combine_hbond() for distinct
    trajectories, this function can be used to present the data as a 
    side-by-side comparison of hbond occupancies.

    hbond_files:  the hbond files with data to be analyzed.
    identifiers:  the list of names associated with each hbond_file.  If
        the list is empty, each file will simply be assigned a number.
        (default: [])
    output_file:  the name of the output file.  If None, the results will be
        written to stdout. (default: None)
    resi_criteria:  a list containing residue number criteria to include in the
        output. (default: ['all'])
    atom_criteria:  a list containing atom name criteria to include in the
        output. (default: ['all'])
    occ_thresh:  the minimum occupancy threshold that the hbonds must have
        to be reported. (default: 0.0)
    occ_graph_only:  if True, only the atom string, occupancy, and graph of
        each hbond will be written to output. (default: False)
    sort:  one of 'occ_diff', 'occ_pct', 'donor', or 'acceptor' that
        indicates how to sort the output. (default: occ_diff)
    compress:  if True, the graphs will be compressed by compress_graph().
        (default: False)
    hbond_data_dir:  the directory that contains the hbond data files.  If
        'None,' the file names will be used without modification, and the
        output will be written to the current directory. (default: None)
    """
    # Set up identifier strings
    
    for i in range(len(hbond_files)):
	if i >= len(identifiers):
	    identifiers.append(str(i + 1))
    max_id_length = max(len(id) for id in identifiers)
    for i in range(len(identifiers)):
	num_spaces = max_id_length - len(identifiers[i])
	identifiers[i] = num_spaces * ' ' + identifiers[i]

    # Do error checking on file names

    files_to_remove = []
    for each_file in hbond_files:
	if hbond_data_dir != None:
	    full_file = os.path.join(hbond_data_dir, each_file)
	else:
	    full_file = each_file
	if not os.path.exists(full_file):
	    print 'Warning:  File ' + full_file + ' does not exist.\n' + \
		  '  Will be ignored.'
	    files_to_remove.append(each_file)
    for each_file in files_to_remove:
	i = hbond_files.index(each_file)
	identifiers.remove(identifiers[i])
	hbond_files.remove(each_file)
    if len(hbond_files) == 0:
	sys.exit('ERROR:  No input files provided.\n')

    if hbond_data_dir != None:
	for i in range(len(hbond_files)):
	    hbond_files[i] = os.path.join(hbond_data_dir, hbond_files[i])

    # Create dictionaries for each file indicating their hbonds
	
    hb_dict_list = []     # One dictionary per hbond input file
    combined_hbonds = {}  # {hbond_string: None} just keeps cumulative track
    for each_file in hbond_files:
	hb_dict = {}      # {hbond_string: hbond object}
	for line in file(each_file):
	    hbond = HBond()
	    hbond.init_from_str(line)
	    if is_resinum_of_interest(hbond, resi_criteria) and \
	       is_atom_of_interest(hbond, atom_criteria):
		if compress:
		    hbond.compress_graph()
		hb_dict[hbond._atom_str()] = hbond
		combined_hbonds[hbond._atom_str()] = None
	hb_dict_list.append(hb_dict)

    # Run through the master list of all hbonds.  If a given
    # dictionary doesn't have an entry for one, create one with
    # zero occupancy. 

    for hb_dict in hb_dict_list:
        for hbond_str in combined_hbonds:
            found = False
            for hbond_str_dict in hb_dict:
                if hbond_str_dict == hbond_str:
                    found = True
                    break
            if not found:
                hbond = HBond()
                hbond.init_from_atomstr(hbond_str)
                hb_dict[hbond_str] = hbond

    # Compile and sort relevant data

    if sort not in 'occ_diff occ_pct donor acceptor'.split():
	print 'Warning:  Unknown sorting method: ' + sort + '.\n' + \
	      '  Will use occ_diff to sort.'
	sort = 'occ_diff'

    output = []
    for hbond_str in combined_hbonds:
	hb_list = [ hb_dict[hbond_str] for hb_dict in hb_dict_list ]

	max_occ = max(hbond.occ_pct for hbond in hb_list)
	min_occ = min(hbond.occ_pct for hbond in hb_list)
	occ_diff = max_occ - min_occ

        if sort == 'occ_diff' and occ_diff > occ_thresh:
            output.append((occ_diff, 
			   hb_list[0].acceptor.resi_num,
			   hb_list))
	elif sort == 'occ_pct' and max_occ > occ_thresh:
	    output.append((max_occ,
			   hb_list[0].acceptor.resi_num,
			   hb_list))
        elif sort == 'donor' and occ_diff > occ_thresh:
            output.append((hb_list[0].donor.resi_num, 
			   hb_list[0].donor.atom_name,
			   hb_list[0].acceptor.resi_num,
			   hb_list[0].acceptor.atom_name,
			   hb_list))
        elif sort == 'acceptor' and occ_diff > occ_thresh:
            output.append((hb_list[0].acceptor.resi_num, 
			   hb_list[0].acceptor.atom_name,
			   hb_list[0].donor.resi_num,
			   hb_list[0].donor.atom_name,
			   hb_list))

    output.sort()
    if sort == 'occ_diff' or sort == 'occ_pct':
        output.reverse()
    output = [o[-1] for o in output]

    # Write output

    output_str = ''
    for each_hbond in output:
	for i in range(len(each_hbond)):
	    hbond = each_hbond[i]
            if occ_graph_only:
		output_str += identifiers[i] + ': ' + str(hbond) + '\n'
            else:
		output_str += identifiers[i] + ': ' + \
                              hbond._atom_str() + ' ' + hbond._attr_str() + '\n'
        output_str += '\n'

    if hbond_data_dir == None:
	output_dir = '.'
    else:
	output_dir = hbond_data_dir

    if output_file == None:
        print output_str[:-2]  # Removes the last 2 newlines
    else:
        try:
            output_file = os.path.join(output_dir, output_file)
            output_f = file(output_file, 'w')
        except IOError:
            print 'Warning:  Could not open ' + output_file + '.\n'
            print output_str
        else:
            output_f.write(output_str[:-1])
            output_f.close()


