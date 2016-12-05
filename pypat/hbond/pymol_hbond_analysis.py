#!/usr/bin/env python

"""

This version uses get_distance, etc. and takes around 10 hours per nanosecond.
That means that it'll take around 100 hours for my 10ns simulation.  If I split
it up onto 8 processors, that'll finish in about 12 or 13 hours.  So, about a
day to do the full hbond analysis.

TODO
----

Perhaps build in a PyMOL interface that will highlight bridging waters
in a trajectory that you're currently viewing.

 - One way of doing this would be to just write out a new trajectory
   that included the appropriate waters in the appropriate states and
   then have the user load that in (or possibly load them in ourselves).

 - Another way would be to hide all of the waters that never make
   hydrogen bonds and then use dist mode=2.


REAL DOCUMENTATION
------------------

There are two main parts to this.

1) Use PyMOL to figure out what the bridging interactions are
   at each snapshot and spit them out to files.  That's typically
   done via a driver script like this::

       #!/usr/bin/env python
       import sys
       PYPAT_CODE_DIR = '/some/path/'

       sys.path.append(PYPAT_CODE_DIR)
       from pypat.hbond import pymol_hbond_analysis
       r = range(2501,3001,500)
       starts_and_stops = zip(r[:-1],r[1:])
       for (start,stop) in starts_and_stops:
           print
           print 'DOING',start,stop
           print
           pymol_hbond_analysis.find_bridging_waters_in_trajectory('1ra1',start,stop)


2) Parse those results, invert them so they're in terms of bridging
   interactions rather than individual waters, spit out the results.
   That's typically done like this (and is done in the __main__
   loop here)::

    import glob
    d = '/some/directory/'
    fnames = glob.glob(os.path.join(d,'1rx1_hbonds_*_*.txt'))
    a = get_hbond_trajectories('1rx1',
                               timestep=5,
                               fnames=fnames,
                               combination_method='loose',
                               min_required_dwell_time=3,
                               looseness=2,
                               dist_cutoff=3.5,
                               )
    including_resis =  137,153,155,30,33,111
    including_resis = None
    print a.get_trajectory_string(minocc=0.20,
                                  numchunks=50,
                                  including_resis=including_resis,
                                  )

   

"""


try:
    enumerate
except NameError:
    def enumerate(thing):
        result = []
        idx = 0
        for t in thing:
            result.append((idx,t))
            idx += 1
        return result
        
        
try:
    sum
except NameError:
    def sum(thing):
        result = 0
        for t in thing:
            result += t
        return result



import sys,os


class TrajectoryFormatter:
    '''
    From ptraj action.c
          if ( m > 0.95 * pt )
            fprintf(fpout, "@");
          else if ( m > 0.80 * pt )
            fprintf(fpout, "*");
          else if ( m > 0.60 * pt )
            fprintf(fpout, "x");
          else if ( m > 0.40 * pt )
            fprintf(fpout, "o");
          else if ( m > 0.20 * pt )
            fprintf(fpout, "-");
          else if ( m > 0.05 * pt )
            fprintf(fpout, ".");
          else 
            fprintf(fpout, " ");
    '''
    def __init__(self):
        pass
    def equal_split(self,L,N):
        """
        Split L into N parts, returning a list containing those parts.  The last part may be smaller than the others.
        """
        part = int((len(L) + N - 1)/ N)
        _L = []
        for i in range(0,N):
            _L.append( L[part*i : part*i + part] )
        return [i for i in _L if i]

    def format(self,mintime,maxtime,timestep,numchunks,trajectory):
        time_chunks = self.equal_split(list(range(mintime,maxtime+timestep,timestep)),numchunks)
        result = ''
        for r in time_chunks:
            occ = 0.
            for t in r:
                if trajectory[t]: occ += 1
            occ = occ/len(r)

            if occ > 0.95:
                result += "@"
            elif occ > 0.80:
                result += "*"
            elif occ > 0.60:
                result += "x"
            elif occ > 0.40:
                result += "o"
            elif occ > 0.20:
                result += "-"
            elif occ > 0.05:
                result += "."
            else:
                result += " "
        return result
    def get_description_string(self):
        
        return """
Trajectory formatting codes:
1.0-0.95 0.95-0.80 0.80-0.60 0.60-0.40 0.40-0.20 0.20-0.05 0.05-0.0
@@@@@@@@ ********* xxxxxxxxx ooooooooo --------- .........
"""

default_trajectory_formatter = TrajectoryFormatter()
            
        

class BridgingWaterTrajectoryAnalyzer:
    """

    """
    def __init__(self,structure,timestep=5):
        """

        You are expected to run code like this:
        
        a = BridgingWaterTrajectoryAnalyzer(structure,timestep)
        for fname in fnames:
            a.read_in_file(fname)
        a.finalize()
        return a

        If you don't call finalize, we will still be able to output some
        statistics, but our trajectories will just be lists of
        SingleSnapshotBridgingWater instances, so you'll miss out on a
        lot of the output power.

        Parameters
        ----------

        structure: protein structure, e.g. 1RX1

        timestep: the number of picoseconds each frame represents.
        """
        self.struct = structure
        self.timestep = timestep
        self.mintime = 1000000000
        self.maxtime = 0
        #
        # self.trajectories is a dictionary that maps SingleSnapshotBridgingWater objects
        # to a series of times.  We'll set it up as a defaultdict so that you can ask
        # for trajectories[SSBW][time] and it'll return True or False.  True if we've
        # explicitly added things and false if we haven't.
        #
        """
        I think this might actually work better if self.trajectories[x] returned
        a trajectory, where x is the hash of a SSBW.  we then call
        self.trajectories[x][time] = x.wat_resi.

        so, that means that a trajectory has to have a __getitem__ and __setitem__.
        """
        self.trajectories = {}
    def __str__(self):
        results = []
        for bwt in self.trajectories:
            result = '%-103s'%bwt + ': ' + bwt.get_trajectory_string()

            results.append(result)
        return '\n'.join(results)

    def get_trajectory_string(self,minocc,numchunks,including_resis,sort_by):
        """
        Parameters
        ----------

        minocc: the minimum occupancy time required to return a string

        numchunks: the number of chunks into which we will divide the output.
                   If this is 10, you'll get 10 points in the output string, etc.
                   If this is None, you'll get the whole thing.

        including_resis: we will print out all trajectories that have at least one
                         leg involving at least one of these residues.

        sort_by: The attribute by which the trajectories will be sorted.  Usual
                 attributes are occ and average_dwell_time.  This should be a
                 string, as we will use it in a getattr() call.

        """
        global default_trajectory_formatter
        bwts = [(getattr(bwt,sort_by),bwt) for bwt in self.trajectories]
        print("BWTA getting %s trajectory strings with minocc %s, numchunks %s, including_resis %s"%(len(bwts),
                                                                                                     minocc,
                                                                                                     numchunks,
                                                                                                     including_resis,
                                                                                                     ))
        print(default_trajectory_formatter.get_description_string())
        bwts.sort()
        bwts.reverse()
        bwts = [thing[1] for thing in bwts]
        return '\n'.join(['%-103s :: '%bwt + bwt.get_trajectory_string(numchunks) for bwt in bwts if bwt.occ >= minocc and bwt.contains_resi(including_resis)])
            
    def as_edges(self,bridging_water_info,edge_length=2,is_valid=None,dist_cutoff=9999999.0,angle_cutoff=0.0):
        """
        Returns a bridging_water_info object as one SingleSnapshotBridgingWater per edge.

        e.g. ['7464', ('53', 'PRO', 'O', 'donatingto',2.2,'angle',160.0),
                      ('52', 'ARG', 'HE', 'acceptingfrom',2.4,'angle',165.0),
                      ('52', 'ARG', 'HH22', 'acceptingfrom',3.5,'angle',164.0)]

             will be returned as
               [SSBW(7464,('53', 'PRO', 'O',  'donatingto',   2.2,'angle',160.0),('52', 'ARG', 'HE',   'acceptingfrom',2.4,'angle',165.0),),
                SSBW(7464,('53', 'PRO', 'O',  'donatingto',   2.2,'angle',160.0),('52', 'ARG', 'HH22', 'acceptingfrom',3.5,'angle',164.0)),
                SSBW(7464,('52', 'ARG', 'HE', 'acceptingfrom',2.4,'angle',165.0),('52', 'ARG', 'HH22', 'acceptingfrom',3.5,'angle',164.0))]

        This is all subject to the is_valid function.  For instance, an is_valid function
        that says that you can't be bridging between two hydrogens on the same residue would
        reject the last entry.

        Parameters
        ----------

        is_valid: a function to tell us if a bridging water is valid.  For example, we might
                  want to reject bridges between two hydrogens on the same residue.  If this
                  is None, we will use our default is_valid function, defined below.

        dist_cutoff: if you're using our is_valid function, it will reject edges where either
                     bond is longer than dist_cutoff.

        angle_cutoff: if you're using our is_valid function, it will reject edges where either
                      angle is less than angle_cutoff

        edge_length: length of the edges that we care about.  Default is two. I suppose it's
                     not really an 'edge' if it's bigger than two, but you know what I mean.
                     It's the list of things that are connected to this water.
        """
        if edge_length != 2:
            raise NotImplementedError('Only edge lengths of two are supported')
        results = []
        if is_valid is None:
            """
            Anything with None in it is invalid, and is just a default arg to the
            Hbond constructor

            We do not allow bridges between two hydrogens that are part of the same
            residue.
            """
            def is_valid(atom1,atom2,dist_cutoff=dist_cutoff,angle_cutoff=angle_cutoff):
                """

                an atom looks like ('111', 'TYR', 'O', 'donatingto',3.2, 'angle',160.0)
                """
                #
                # Reject default args (None)
                #
                if None in atom1:
                    return False
                if None in atom2:
                    return False
                #
                # Reject bridges between two hydros that are part of the same residue
                #
                if atom1[0] == atom2[0]: # same residue
                    if 0: print("same residue hydros",atom1,atom2)
                    if atom1[2].startswith('H') and atom2[2].startswith('H'):
                        #print "rejecting",atom1,atom2
                        return False
                #
                # Bridges must be between two different residues
                #
                if atom1[0] == atom2[0]:
                    if 0: print("same residue",atom1,atom2)
                    return False
                
                #
                # All edges must be <= the distance cutoff
                #
                if (atom1[4] > dist_cutoff) or (atom2[4] > dist_cutoff):
                    if 0:
                        print("bad dist",atom1,atom2)
                    return False
                #
                # All angles must be >= the angle cutoff
                #
                if (atom1[6] < angle_cutoff) or (atom2[6] < angle_cutoff):
                    if 0:
                        print("bad angle",atom1,atom2)
                    return False
                return True
        wat_resi,atoms = bridging_water_info[0],bridging_water_info[1:]
        for i,atom1 in enumerate(atoms):
            for atom2 in atoms[i+1:]:
                if is_valid(atom1,atom2):
                    results.append(SingleSnapshotBridgingWater(wat_resi,atom1,atom2))
        return results

    def finalize(self,combination_method,min_required_dwell_time,looseness):
        """
        Turn our trajectories into proper trajectories so that we can print out nicely.

        It also calculates a lot of statistics.

        Please look at BridgingWaterTrajectory.finalize() for more documentation.
        """
        sys.stdout.write("Now finalizing %s trajectories\n"%len(self.trajectories))
        sys.stdout.flush()
        for bwt in self.trajectories:
            bwt.mintime = self.mintime
            bwt.maxtime = self.maxtime
            bwt.timestep = self.timestep
            bwt.finalize(combination_method=combination_method,
                         min_required_dwell_time=min_required_dwell_time,
                         looseness=looseness)
            

    def read_in_file(self,fname,times=None,dist_cutoff=999999.,angle_cutoff=0):
        """

        Parameters
        ----------

        fname: The name of the file we should read in.
               The contents of this file will be eval()'d, and we expect
               to be able to turn the result of that into a dictionary.
               The keys in that dictionary are PyMOL object names and the
               values are lists of bridging waters found in those objects.

        times: A dictionary mapping the PyMOL object names to times.
               If times is None, we will determine the times from the
               filename as follows:

               1rx1_hbonds_635_640.txt means that times is range(635,640)
               with the special exception that, due to how PyMOL treats
               trajectories, there's never a time == 0, so 0_x will be
               range(1,x).

               In that case, we'll map the times to the lexigraphical
               ordering of the keys in fname.

               Times will be in picoseconds, and will be multiplied by
               self.timestep in order to convert from snapshot number
               to picoseconds.

               This function takes care of some bookkeeping as well by
               ensuring that self.mintime and self.maxtime are correct.
               Any other functions that touch self.trajectories should
               make sure to do this!

        dist_cutoff: if you're using our is_valid function, it will
                     reject edges where either bond is longer than
                     dist_cutoff.

        angle_cutoff: if you're using our is_valid function, it will
                      reject edges where either angle is smaller than
                      angle_cutoff.

        """
        f = file(fname)
        try:
            x = eval(f.read())
        except SyntaxError:
            print("Could not read file %s .. possibly still being written to?"%fname)
            return
        f.close()
        if times is None:
            just_fname = os.path.splitext(os.path.split(fname)[-1])[0]
            parts = just_fname.split('_')
            start,stop = int(parts[-2]),int(parts[-1])
            if start == 0:
                start = 1
            times = {}
            for (o,t) in zip(sorted(x.keys()),list(range(start,stop))):
                times[o] = t*self.timestep
            if 0:
                print("Time mapping",times)
        for obj_name,bridging_water_infos in x.items():
            t = times[obj_name]
            if t < self.mintime:
                self.mintime = t
            if t > self.maxtime:
                self.maxtime = t
            for bridging_water_info in bridging_water_infos:
                wat_resi = int(bridging_water_info[0])
                for bw in self.as_edges(bridging_water_info,dist_cutoff=dist_cutoff,angle_cutoff=angle_cutoff):
                    self.add_to_trajectory(bw,wat_resi,t)
    def add_to_trajectory(self,bw,wat_resi,t):
        if bw not in self.trajectories:
            il = bw.get_nonempty_interaction_list(include_distances=True)
            if 0:
                print(il)
            if len(il) != 2:
                raise NotImplementedError('We only support bridges with two edges at this time %s'%il)
            if 0:
                print(il[0])
            resi1,resn1,atomname1,interaction1,dist1 = il[0]
            resi2,resn2,atomname2,interaction2,dist2 = il[1]
            #
            # TODO: FIXME: we should pay more attention to the distance.
            # We should store it with the time in self.trajectories.
            #
            if 0:
                print("interaction1",interaction1)
                print("interaction2",interaction2)
            bwt = BridgingWaterTrajectory((resi1,resn1,atomname1,interaction1,dist1,),
                                          (resi2,resn2,atomname2,interaction2,dist2,),
                                          )

            self.trajectories[bwt] = bwt
            if 0:
                print("added",(resi1,resn1,atomname1,interaction1,dist1,),(resi2,resn2,atomname2,interaction2,dist2,))
            
        # It's worth knowing that I've overloaded __setitem__ so that it will
        # actually just append wat_resi to the list of times.x
        try:
            self.trajectories[bw][t] = wat_resi
        except KeyError:
            print("trouble adding",bw)
            raise

        
class BridgingWaterTrajectory:
    """
        I think this might actually work better if self.trajectories[x] returned
        a trajectory, where x is the hash of a SSBW.  we then call
        self.trajectories[x][time] = x.wat_resi.

        so, that means that a trajectory has to have a __getitem__ and __setitem__.

    """

    def __init__(self, xxx_todo_changeme, xxx_todo_changeme1,
                 ):
        #
        # If we every use more trajectory formatters,
        # we'll make it an argument to __init__.
        #
        (resi1,resn1,atomname1,interaction1,dist1) = xxx_todo_changeme
        (resi2,resn2,atomname2,interaction2,dist2) = xxx_todo_changeme1
        global default_trajectory_formatter
        from collections import defaultdict

        self.resi1 = int(resi1)
        self.resn1 = resn1
        self.atomname1 = atomname1
        self.interaction1 = interaction1
        self.resi2 = int(resi2)
        self.resn2 = resn2
        self.atomname2 = atomname2
        self.interaction2 = interaction2

        self._times = defaultdict(list)

        self.mintime = None
        self.maxtime = None
        self.timestep = None
        self.occ = None

        self.tf = default_trajectory_formatter

    def get_nonempty_interaction_list(self,include_distances=False):
        #
        # We never include the wat_resi, because we're just saying
        # what the interactions are.
        #
        # This is necessary for comparing to SSBTs
        # so, in the general case, we will not include distances.
        # However, it's also useful for building up other lists of
        # info, so we will allow the possibility of including
        # distances.
        #
        if include_distances:
            return [(self.resi1,self.resn1,self.atomname1,self.interaction1,dist1),
                    (self.resi2,self.resn2,self.atomname2,self.interaction2,dist2),
                    ]
        else:
            return [(self.resi1,self.resn1,self.atomname1,self.interaction1,),
                    (self.resi2,self.resn2,self.atomname2,self.interaction2,),
                    ]

    def __getitem__(self,time):
        return self._times[time]
    def __setitem__(self,time,wat_resi):
        if 0:
            print("setting",time,wat_resi)
        self._times[time].append(wat_resi)
        return self._times[time]
    def __hash__(self):
        """
        SUPER IMPORTANT NOTE:
        It is very important to make sure that this hash function
        is the same as the hash function used for an SSBW.  Otherwise,
        the indexing into .trajectories won't work at all.

        We don't include distance for the same reason that we don't include
        wat_resi.
        """
        inter = sorted([(self.resi1,self.resn1,self.atomname1,self.interaction1,),
                 (self.resi2,self.resn2,self.atomname2,self.interaction2,),
                 ])
        if 0:
            print("Will hash..",tuple(inter),"to",hash(tuple(inter)),"in bwt")
        return hash(tuple(inter))
    def __repr__(self):
        try:
            verbose = False
            if verbose:
                result = '<BI %7.1fps occ: %5.1f adt: %5.1fps %3i [%s] (%3s %3s %-4s):%13s, (%3s %3s %-4s):%13s,>'%(self.maxtime - self.mintime + self.timestep,
                                                                                                                    self.occ * 100,
                                                                                                                    self.average_dwell_time,
                                                                                                                    self.num_waters_seen,
                                                                                                                        self.__wat_resis,
                                                                                                                    self.resn1,self.resi1,self.atomname1,
                                                                                                                    self.interaction1,
                                                                                                                    self.resn2,self.resi2,self.atomname2,
                                                                                                                    self.interaction2,
                                                                                                                    )
            else:
                result = '<BI %7.1fps occ: %5.1f adt: %5.1fps %3i (%3s %3s %-4s):%13s, (%3s %3s %-4s):%13s,>'%(self.maxtime - self.mintime + self.timestep,
                                                                                                                self.occ * 100,
                                                                                                                self.average_dwell_time,
                                                                                                                self.num_waters_seen,
                                                                                                                #     self.__wat_resis,
                                                                                                                self.resn1,self.resi1,self.atomname1,
                                                                                                                self.interaction1,
                                                                                                                self.resn2,self.resi2,self.atomname2,
                                                                                                                self.interaction2,
                                                                                                                )
                
        except TypeError:
            print("Args")
            print((self.maxtime, self.mintime, self.timestep,
                   self.occ,100,
                   self.average_dwell_time,
                   self.num_waters_seen,
                   self.resn1,self.resi1,self.atomname1,
                   self.interaction1,
                   self.resn2,self.resi2,self.atomname2,
                   self.interaction2,))
            raise
        return result
    def contains_resi(self,resis):
        """
        Returns true of anything in resis is contained in any of our legs.

        Parameters
        ----------

        resis: a single resi, a list of resis, or None.  If None, we will
               return True for everything.
        """
        if resis is None:
            return True
        if isinstance(resis, type('')):
            resis = [int(i) for i in resis.split()]
        elif isinstance(resis, type(1)):
            resis = [resis,]
        if len(self.get_nonempty_interaction_list()) != 2:
            raise NotImplementedError("Only bwts of length two are supported")
        for resi in resis:
            if resi == self.resi1:
                return True
            if resi == self.resi2:
                return True
        return False
    def get_trajectory_string(self,numchunks):
        """
        TODO: FIXME: verify that this is correct

        """
        if numchunks is None:
            result = ''
            for t in range(self.mintime,self.maxtime+self.timestep,self.timestep):
                if self[t]:
                    #result += '.'
                    result += '%s,'%self[t]
                else:
                    result += ' '
        else:
            result = self.tf.format(self.mintime,
                                    self.maxtime,
                                    self.timestep,
                                    numchunks,
                                    self,
                                    )
                                    
        return result
    def finalize(self,combination_method,min_required_dwell_time,looseness):
        """
        Calculate occupancy, dwell times.

        Parameters
        ----------

        combination_method: We see several trajectories where the occupancy by a given
                            water molecule will look like ..... . . .....
                            If you choose method 'strict' you'll see that as the same
                            water molecule occupying the bridge 4 different times.
                            The average dwell time there would be 4+1+1+4/4.
            
                            If you choose method 'combine' we will simply keep track
                            of how much time a given water molecule spends in the
                            bridge for the entire simulation.  That fixes up the
                            above case, but has some trouble.  For instance, of
                            water molecule 3254 occupies the bridge for 4ps, then
                            leaves for 1ns, then comes back for 8ps, we'll record
                            a dwell time of 12ps.
            
                            The difference between these methods is very significant.
                            In one particular case, for instance, I see it change
                            the dwell time from 465ps to 37.4ps.
            
                            If you choose method 'loose' we will try to be a little
                            smarter about combining the small trajectories.  We will
                            allow a skip of one between them.  That is,.. ....  ...
                            will be turned into .......  ... which is, I think,
                            a little better.

        looseness: If we're using method 'loose', this tells us how many
                   snapshots can be missing and still let us combine the
                   trajectory.
                
        min_required_dwell_time: measured in timesteps, not picoseconds.
        
        """
        from collections import defaultdict
        #
        # Occupancy
        #
        occ = 0.
        r = list(range(self.mintime,self.maxtime+self.timestep,self.timestep))
        for t in r:
            if self._times[t]:
                occ += 1
        self.occ = occ/len(r)

        #
        # Dwell times
        #
        #combination_method = 'strict'
        if combination_method == 'strict':
            raise NotImplementedError("'strict' method is no longer supported.  Please look through the source and re-enable it if you want.")
            counts = []
            last_wat = None
            count = 1
            for t in r:
                if self._times[t] == last_wat:
                    count += 1
                else:
                    if last_wat not in (False,None,[]):
                        counts.append(count)
                    count = 1
                last_wat = self._times[t]
            if last_wat != None:
                counts.append(count)
            self.dwell_times = [c * self.timestep for c in counts]
        elif combination_method == 'combine':
            raise NotImplementedError("'combine' method is no longer supported.  Please look through the source and re-enable it if you want.")
            counts = defaultdict(int)
            for t in r:
                wat_resi = self._times[t]
                if wat_resi not in (False,None,[]):
                    counts[wat_resi] += 1
            self.dwell_times = [i * self.timestep for i in list(counts.values())]
        elif combination_method == 'loose':
            """
            Writing this algorithm myself is actually a little tricky.  Instead,
            I'm just going to use strings.  For each water that we see, I'll
            make a trajectory string that contains 'x' when the water is there and
            ' ' when it's not.  Then, I can use string.replace('x x', 'xxx').
            """
            self.dwell_times = []
            #
            # I used to say this:
            #wat_resis = [i for i in set(self._times.values()) if i not in (None,False,[])]
            # but that doesn't work anymore.  We've made self._times
            # a list that contains all of the wat_resis satisfying this
            # particular bridge at a given time.
            #
            wat_resis = set()
            for i in list(self._times.values()):
                for j in i:
                    wat_resis.add(j)
            for wat_resi in wat_resis:
                traj_str = ''
                for t in r:
                    if wat_resi in self._times[t]:
                        traj_str += 'x'
                    else:
                        traj_str += ' '
                #print wat_resi,traj_str
                """
                A Note About Regular Expressions:

                This is not good enough:

                for i in range(1,looseness+1):
                    traj_str = traj_str.replace('x'+' '*i+'x','x'+'x'*i+'x')
                    print wat_resi,traj_str,'replacing',i

                because the regular expression engine returns only non-overlapping
                matches.  That means that we'll see something like this:

                In [71]: 'x x xx  x'.replace('x x','xxx')
                Out[71]: 'xxx xx  x'
                
                In order to take care of overlapping matches, we have to apply the
                pattern multiple times.  Question: how many times?  I think that
                two is enough, because the pattern and replacement is simple enough
                that we really only care about its edges.

                Who really cares about that, though?  I'll just keep replacing until
                I stop finding what I'm looking for.
                """
                for i in range(1,looseness+1):
                    pattern     = 'x'+' '*i+'x'
                    replacement = 'x'+'x'*i+'x'
                    while traj_str.find(pattern) >= 0:
                        traj_str = traj_str.replace(pattern,replacement)
                        #print wat_resi,traj_str,'replacing',i
                for traj in traj_str.split():
                    self.dwell_times.append(len(traj) * self.timestep)
        self.dwell_times = [i for i in self.dwell_times if i >= self.timestep*min_required_dwell_time]
        self.num_waters_seen = len(self.dwell_times)
        self.average_dwell_time = 0
        if self.dwell_times:
            self.average_dwell_time = sum(self.dwell_times)/len(self.dwell_times)
        self.__wat_resis=wat_resis #TODO: FIXME: make sure this actually works
        #
        # TODO: clarify this: there's a question about exactly what I mean by "num_waters_seen"
        # as is, it's the number of .. well, make sure it makes sense later, with 2 molecules at once, etc.
        #
        
        #
        # This is a little flakey .. if we're there for 10 steps, gone for one
        # and back for 10, that will show up as seen twice.
        #
            
        
    
    
class SingleSnapshotBridgingWater:

    def __init__(self,
                 wat_resi, xxx_todo_changeme2, xxx_todo_changeme3, xxx_todo_changeme4=(None,None,None,None,None,None,None,), xxx_todo_changeme5=(None,None,None,None,None,None,None,), xxx_todo_changeme6=(None,None,None,None,None,None,None,), xxx_todo_changeme7=(None,None,None,None,None,None,None,), xxx_todo_changeme8=(None,None,None,None,None,None,None,), xxx_todo_changeme9=(None,None,None,None,None,None,None,),
                 ):
        """

        Parameters
        ----------
        resi1,atomname1,interaction1,dist1: one of the residues, and in interaction
                                            in {donatingto,acceptingfrom}.
                                            153,donatingto would mean that the water is
                                            donating an hbond to residue 153 with distance dist1.

        At the moment, there's no reason to record the angles, so we won't.
        
        """
        (resi1,resn1,atomname1,interaction1,dist1,_angle1,angle1) = xxx_todo_changeme2
        (resi2,resn2,atomname2,interaction2,dist2,_angle2,angle2) = xxx_todo_changeme3
        (resi3,resn3,atomname3,interaction3,dist3,_angle3,angle3) = xxx_todo_changeme4
        (resi4,resn4,atomname4,interaction4,dist4,_angle4,angle4) = xxx_todo_changeme5
        (resi5,resn5,atomname5,interaction5,dist5,_angle5,angle5) = xxx_todo_changeme6
        (resi6,resn6,atomname6,interaction6,dist6,_angle6,angle6) = xxx_todo_changeme7
        (resi7,resn7,atomname7,interaction7,dist7,_angle7,angle7) = xxx_todo_changeme8
        (resi8,resn8,atomname8,interaction8,dist8,_angle8,angle8) = xxx_todo_changeme9
        self.wat_resi = wat_resi
        if resi1 is not None: resi1 = int(resi1)
        if resi2 is not None: resi2 = int(resi2)
        if resi3 is not None: resi3 = int(resi3)
        if resi4 is not None: resi4 = int(resi4)
        if resi5 is not None: resi5 = int(resi5)
        if resi6 is not None: resi6 = int(resi6)
        if resi7 is not None: resi7 = int(resi7)
        if resi8 is not None: resi8 = int(resi8)
        self.interactions = {(resi1,resn1,atomname1,dist1):interaction1,
                             (resi2,resn2,atomname2,dist2):interaction2,
                             (resi3,resn3,atomname3,dist3):interaction3,
                             (resi4,resn4,atomname4,dist4):interaction4,
                             (resi5,resn5,atomname5,dist5):interaction5,
                             (resi6,resn6,atomname6,dist6):interaction6,
                             (resi7,resn7,atomname7,dist7):interaction7,
                             (resi8,resn8,atomname8,dist8):interaction8,
                            }
    def get_nonempty_interaction_list(self,include_distances=False):
        #
        # We never include the wat_resi, because we're just saying
        # what the interactions are.
        #
        # This is necessary for comparing to BWTs
        # so, in the general case, we will not include distances.
        # However, it's also useful for building up other lists of
        # info, so we will allow the possibility of including
        # distances.
        #
        result = []
        for k,v in self.interactions.items():
            if None in k:
                continue
            if v is None:
                continue
            if include_distances:
                result.append((k[0],k[1],k[2],v,k[3]))
            else:
                result.append((k[0],k[1],k[2],v))
        return result

    def __repr__(self):
        """

        string representation.  For purposes of principal, we will not include
        the bridging water here.
        """
        #result = '(%s,{'%self.wat_resi
        result = '<SSBW %s '%self.wat_resi
        for k,v in self.interactions.items():
            if None in k:
                continue
            if v is None:
                continue
            result += "%s:'%s',"%(k,v)
        #result +='})'
        result +='>'
        return result
        #return str((self.wat_resi,self.interactions))
    def __hash__(self):
        """
        """
        #
        # if we get rid of self.wat_resi, we may not be able to
        # correctly keep track of dwell times for particular waters.
        #
        inter = sorted(self.get_nonempty_interaction_list())
        #return hash((self.wat_resi,tuple(inter)
        #             ))
        if 0:
            print("Will hash",inter)
        return hash(tuple(inter))

    def __eq__(self,other):
        if 0:
            print("Comparing",self,other,sorted(self.get_nonempty_interaction_list()) == sorted(other.get_nonempty_interaction_list()))
        return sorted(self.get_nonempty_interaction_list()) == sorted(other.get_nonempty_interaction_list())
        #return self.interactions == other.interactions
    def __ne__(self,other):
        return sorted(self.get_nonempty_interaction_list) != sorted(other.get_nonempt_interaction_list)
        #return self.interactions != other.interactions


class SingleSnapshotHbondEmitter:

    def __init__(self,states,hbond_dist_cutoff=3.5,hbond_angle_cutoff=None):
        '''

        Parmeters
        ---------

        states: the number of states we will analyze.
        
        hbond_dist_cutoff: HeavyAtom-Hydro distance, so we default to a generous 3.5.

        hbond_angle_cutoff: At the moment, we will mimic ptraj and not have an angle cutoff.

        Our lists of donors and acceptors are taken from http://amber.scripps.edu/tutorials/basic/tutorial3/files/analyse_hbond.ptraj
        '''
        self.states = states
        self.hbond_dist_cutoff = hbond_dist_cutoff
        #if hbond_angle_cutoff is not None:
        #    raise NotImplementedError('hbond_angle_cutoff not implemented yet (we are like ptraj default here).')
        self.hbond_angle_cutoff = hbond_angle_cutoff

        self.logfilename = 'HbondEmitterLogFile.txt'
    def log(self,message):
        f = file(self.logfilename,'a')
        f.write(message)
        f.close()
        sys.stdout.write(message)
        sys.stdout.flush()

    def create_selections(self):

        from pymol import cmd,stored
        from . import hbond_definitions
        hbond_definitions.do_standard_selections()

        if 0:
            #Once we do the per-state selections, we may be able
            # to resurrect this.
            cmd.select('prot_acceptors','prot_acceptors and %s'%sel)
            cmd.select('prot_donors','prot_donors and %s'%sel)
    def analyze_and_emit(self,outf):
        #
        # Note, I think it's a little clearer this way.  OTOH, we could
        # quite easily change the format so that analyze() takes in a
        # file handle and writes each bridge out to it on one line,
        # something like (state,bridge)
        #
        results = self.analyze()
        outf.write(str(results))
        
    def analyze(self):
        """
        go through the whole thing at once

        You might be wondering why I have all of these safety checks in here.
        It turns out that they've all been necessary at one time or another.

        In particular, I often find that the last water in a trajectory will
        trigger the PROBLEM WITH errors, and I'll often see things like
        FAILED H2DIST 3 1, neither of which should ever happen.
        """
        from pymol import cmd,stored
        self.create_selections()
        bridges = {}
        for state in range(1,self.states+1):
            bridges[state] = {}
            stored.wat_list = []
            #
            # Here, we have a question about how to make this as fast as
            # possible. After talking to Warren, it looks like the best thing
            # to do is to make a new object (tmp_wat_papd .. temporary waters +
            # protein acceptors + protein donors) that contains only the atoms
            # from a particular state that we are interested in. This sort of
            # object creation should be faster than passing the "state" argument
            # to get_distance and get_angle.
            #
            # We'll be a little bit generous in creating that object: we'll grab
            # everything within dist_cutoff of any relevant atoms, rather than
            # just searching for the heavy-heavy distances.
            #
            cmd.create('tmp_wat_papd','(prot_donors or prot_acceptors or neighbor prot_donors) or (byres (resn WAT+HOH) within %s of (prot_donors or prot_acceptors or neighbor prot_donors))'%(self.hbond_dist_cutoff),state,1)
            cmd.select('tmp_wat_near','byres (tmp_wat_papd and (resn WAT+HOH)) within %s of (tmp_wat_papd in (prot_donors or prot_acceptors or neighbor prot_donors))'%(self.hbond_dist_cutoff+0.1))
            cmd.iterate('tmp_wat_near and elem o','stored.wat_list.append(resi)')
            num_wat_resis = len(stored.wat_list)
            self.log('\nI will loop over %s nearby waters for state %s\n'%(num_wat_resis,state))
            for (wat_idx,w_resi) in enumerate(stored.wat_list):
                pct_done = int(100*(wat_idx+1)/num_wat_resis)
                if divmod(pct_done,10)[-1] == 0: self.log('X')
                else: self.log('.')

                bridges[state][w_resi] = []
                stored.pa_list = []
                cmd.iterate('(tmp_wat_papd in prot_acceptors) within %s of (tmp_wat_papd and resi %s and elem h)'%(self.hbond_dist_cutoff,w_resi),'stored.pa_list.append((resi,resn,name))')
                for (resi,resn,name) in stored.pa_list:
                    # Here, we're looking for distances to the water hydrogens,
                    # which is why we use resi <w_resi> and name H1|H2.
                    #
                    # Also, the select(left, (not left) and (neighbor left))
                    # business is because our hbond definition file lists
                    # the hydrogens, but we want heavy-atom distances.

                    # Calculate the distance/angle w.r.t the water's H1
                    left_hydro = cmd.select('left_hydro',"resi %s and name H1 and tmp_wat_papd"%w_resi)
                    if left_hydro != 1: self.log("PROBLEM WITH resi %s and name H1 and tmp_wat_papd in state %s"%(w_resi,state))

                    left_heavy = cmd.select('left_heavy','(not left_hydro) and (neighbor left_hydro) and tmp_wat_papd')
                    if left_heavy != 1: self.log("PROBLEM WITH new left sele resi %s and name H1 and tmp_wat_papd in state %s"%(w_resi,state))

                    right = cmd.select('right',"resi %s and name %s and tmp_wat_papd"%(resi,name))
                    if right != 1: self.log("PROBLEM WITH resi %s and name %s and tmp_wat_papd in state %s"%(resi,name,state))

                    try:
                        h1dist = cmd.get_distance('left_heavy','right') # state not included because we're in tmp_wat_papd
                    except:
                        self.log("FAILED H1DIST %s %s"%(left_heavy,right))
                        cmd.delete('adist')
                        h1dist = cmd.distance('adist','left_heavy','right')
                        self.log("DIST SAYS %s"%h1dist)
                    try:
                        h1angle = cmd.get_angle('left_heavy','left_hydro','right')
                    except:
                        self.log('FAILED H1ANGLE %s %s %s'%(left_heavy,left_hydro,right))
                        cmd.delete('anangle')
                        h1angle = cmd.angle('anangle','left_heavy','left_hydro','right')

                    # Calculate the distance/angle w.r.t the water's H2
                    left_hydro = cmd.select('left_hydro',"resi %s and name H2 and tmp_wat_papd"%w_resi)
                    if left_hydro != 1: self.log("PROBLEM WITH resi %s and name H2 and tmp_wat_papd in state %s(%s)"%(w_resi,state,left))

                    left_heavy = cmd.select('left_heavy','(not left_hydro) and (neighbor left_hydro) and tmp_wat_papd')
                    if left_heavy != 1: self.log("PROBLEM WITH new left sele resi %s and name H2 and tmp_wat_papd in state %s"%(w_resi,state))

                    right = cmd.select('right',"resi %s and name %s and tmp_wat_papd"%(resi,name))
                    if right != 1: self.log("PROBLEM WITH resi %s and name %s and tmp_wat_papd in state %s(%s)"%(resi,name,state,right))

                    try:
                        h2dist = cmd.get_distance('left_heavy','right') # state not included because we're in tmp_wat_papd
                    except:
                        self.log("FAILED H2DIST %s %s"%(left,right))
                        cmd.delete('adist')
                        h2dist = cmd.distance('adist','left_heavy','right')
                        self.log("DIST SAYS %s"%h2dist)
                    try:
                        h2angle = cmd.get_angle('left_heavy','left_hydro','right')
                    except:
                        self.log('FAILED H2ANGLE %s %s %s'%(left_heavy,left_hydro,right))
                        cmd.delete('anangle')
                        h2angle = cmd.angle('anangle','left_heavy','left_hydro','right')

                    # Given the choice, we chose the one with the smallest distance.
                    if h1dist <= h2dist:
                        if (h1dist <= self.hbond_dist_cutoff) and (h1angle >= self.hbond_angle_cutoff):
                            bridges[state][w_resi].append((resi,resn,name,'donatingto',h1dist,'angle',h1angle))
                        elif (h2dist <= self.hbond_dist_cutoff) and (h2angle >= self.hbond_dist_cutoff):
                            bridges[state][w_resi].append((resi,resn,name,'donatingto',h2dist,'angle',h2angle))
                    else:
                        if (h2dist <= self.hbond_dist_cutoff) and (h2angle >= self.hbond_dist_cutoff):
                            bridges[state][w_resi].append((resi,resn,name,'donatingto',h2dist,'angle',h2angle))
                        elif (h1dist <= self.hbond_dist_cutoff) and (h1angle >= self.hbond_angle_cutoff):
                            bridges[state][w_resi].append((resi,resn,name,'donatingto',h1dist,'angle',h1angle))
                    # old dist-only code:
                    #if (h2dist <= self.hbond_dist_cutoff) and (h2angle >= self.hbond_dist_cutoff):
                    #    bridges[state][w_resi].append((resi,resn,name,'donatingto',min(h1dist,h2dist)))
                stored.pd_list = []
                cmd.iterate('(tmp_wat_papd in prot_donors) within %s of (tmp_wat_papd and resi %s and elem o)'%(self.hbond_dist_cutoff,w_resi),'stored.pd_list.append((resi,resn,name))')
                for (resi,resn,name) in stored.pd_list:
                    left = cmd.select('left',"resi %s and name o and tmp_wat_papd"%w_resi)
                    if left != 1: self.log("PROBLEM WITH resi %s and name o and tmp_wat_papd in state %s(%s)"%(w_resi,state,left))

                    right_hydro = cmd.select('right_hydro',"resi %s and name %s and tmp_wat_papd"%(resi,name))
                    if right_hydro != 1: self.log("PROBLEM WITH resi %s and name %s and tmp_wat_papd in state %s(%s)"%(resi,name,state,right))

                    right_heavy = cmd.select('right_heavy','(not right_hydro) and (neighbor right_hydro) and tmp_wat_papd')
                    if right_heavy != 1: self.log("PROBLEM WITH new right sele resi %s and name %s and tmp_wat_papd in state %s(%s)"%(resi,name,state,right))
                    
                    try:
                        hdist = cmd.get_distance('left','right_heavy') # state not included because we're in tmp_wat_papd
                    except:
                        self.log("FAILED HDIST %s %s"%(left,right))
                        cmd.delete('adist')
                        hdist = cmd.distance('adist','left','right_heavy')
                        self.log("DIST SAYS %s"%hdist)
                    try:
                        hangle = cmd.get_angle('left','right_hydro','right_heavy')
                    except:
                        self.log('FAILED HANGLE %s %s %s'%(left,right_hydro,right_heavy))
                        cmd.delete('anangle')
                        hangle = cmd.angle('anangle','left','right_hydro','right_heavy')
                    if (hdist <= self.hbond_dist_cutoff) and (hangle >= self.hbond_angle_cutoff):
                        bridges[state][w_resi].append((resi,resn,name,'acceptingfrom',hdist,'angle',hangle))
        real_bridges = {}
        for state in bridges:
            for wat_resi,bridge in bridges[state].items():
                if len(bridge) >= 2:
                    if state not in real_bridges: real_bridges[state] = []
                    real_bridges[state].append([wat_resi,] + bridge)
        return real_bridges


def find_bridging_waters_in_trajectory(name,start,stop,hbond_dist_cutoff,hbond_angle_cutoff,overwrite=False):
    """
    Read in a trajectory file and find bridging hbonds.
    
    Parameters
    ----------

    name: We will read in <name>.trj and <name>.top from the
          current working directory.

    start,stop: The states of the trajectory are like indices in
                a list.  We will process list[start:stop].  You
                should know that there is no state 0.

    overwrite: If this is False, we will not overwrite existing files.
               Rather, we will print a message and do nothing.

    Output
    ------

    The results will be written to <name>_hbonds_<start>_<stop>.txt.
    If that file exists, it will be overwritten.  That file will
    contain a Python dictionary of all of the bridging waters.  The
    keys in that dictionary will be state numbers, relative to the
    start:stop range.  That is, if it's the first state we actually
    read in after all of the skipping, it'll be called state 1.
    That's why the emitter doesn't have to know what the absolute
    starts and stops are.
    """
    fname = '%s_hbond_%s_%s.txt'%(name,start,stop)
    if os.path.exists(fname) and not overwrite:
        print(fname,"already exists.  Skipping.")
        return
    print("Opening",fname,"for writing")
    f = file(fname,'w')

    from pymol import cmd,stored
    cmd.delete('all')
    print("ready to load top")
    cmd.load(name+'.top')
    print("ready to load trj")
    cmd.load_traj(name+'.trj',start=start,stop=stop)
    #
    # Yeah, tell me about it. But, PyMOL will often load up water molecules from
    # an AMBER trajectory so that the hydrogens are bonded to each other.
    #
    cmd.unbond('hydro','hydro')
    states = cmd.count_states(name)
    e = SingleSnapshotHbondEmitter(states=states,hbond_dist_cutoff=hbond_dist_cutoff,hbond_angle_cutoff=hbond_angle_cutoff)
    e.analyze_and_emit(f)
    print("Closing",fname)
    f.close()

if __name__ == '__pymol__':
    pymol.cmd.extend('fbwt',find_bridging_waters_in_trajectory)

def get_hbond_trajectories(structure,timestep,fnames,combination_method,min_required_dwell_time,looseness,dist_cutoff,angle_cutoff):
    """
    Takes in structure,timestep,fnames, returns a BridgingWaterTrajectoryAnalyzer.

    This assumes that times can be automatically determined from
    filenames.
    """
    print("Bridging Water Trajectory calculated with combination_method %s, looseness %s, min_required_dwell_time %s, dist_cutoff %s, angle_cutoff %s"%(combination_method,
                                                                                                                                                        looseness,
                                                                                                                                                        min_required_dwell_time,
                                                                                                                                                        dist_cutoff,
                                                                                                                                                        angle_cutoff,
                                                                                                                                                        ))
    a = BridgingWaterTrajectoryAnalyzer(structure,timestep)
    #
    # There are two places where we filter things out.
    # 1) When we're reading in the files, we filter by static information
    #    like dist_cutoff and angle_cutoff
    # 2) When we're finalizing, we filter by trajectory information
    #    like min_required_dwell_time.
    #
    for fname in fnames:
        a.read_in_file(fname,dist_cutoff=dist_cutoff,angle_cutoff=angle_cutoff)
    a.finalize(combination_method=combination_method,
               min_required_dwell_time=min_required_dwell_time,
               looseness=looseness)
    return a

if __name__ == '__main__':
    #
    # Main code moved to script in example directory.
    # Put some test functions here.
    #
    pass
