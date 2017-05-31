#!/cluster/home2/mglerner/anaconda3/bin/python
#!/usr/bin/env python

import sys,os,re

#
# Each section that we care about looks like this:
#
exampleSection = """
 NSTEP =     0  TIME(PS) =    0.000  TEMP(K) =   457.87  PRESS =      0.00
 Etot   = -119068.4556  EKtot   =   33623.6366  EPtot      = -152692.0922
 BOND   =     130.2398  ANGLE   =     448.8823  DIHED      =     849.5638
 1-4 NB =     515.8140  1-4 EEL =    6193.3324  VDWAALS    =   15714.9288
 EELEC  = -176544.8532  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 Ewald error estimate:   0.6201E-04
 ------------------------------------------------------------------------------

===============================================================================
                      NMR restraints for step      0
 Energy (this step): Bond =    0.000   Angle =     0.000   Torsion =     0.000
 Energy (tot. run) : Bond =    0.000   Angle =     0.000   Torsion =     0.000

 DEVIATIONS:    Target=(r2+r3)/2                 Target = closer of r2/r3
            This step         Entire run        This step        Entire run
           ave.    rms       ave.    rms       ave.    rms      ave.     rms
 Bond      0.000   0.000     0.000   0.000     0.000   0.000     0.000   0.000
 Angle     0.000   0.000     0.000   0.000     0.000   0.000     0.000   0.000
 Torsion   0.000   0.000     0.000   0.000     0.000   0.000     0.000   0.000
===============================================================================
"""


class StepInfo(dict):
    #
    # see the comments in __getitem__ for info about allKeys.
    # we also use it to cheat .. we know that 'EAMBER (non-constraint)' will
    # show up eventually, but maybe not in the first step, so we'll seed our
    # keys list.
    #
    allKeys = {'EAMBER (non-constraint)':True}
    def __init__(self, inString):
        #
        # The useful lines are all of the form
        # A = B  C=D E = F
        #
        # That is, it's always name=value pairs, but
        #  - there may or may not be spaces around either
        #    side of the =.
        #  - the name may be one or more words
        #
        # We assume that we've only been passed useful lines.
        #
        # So, we need to be a little clever.  We concatenate the
        # lines into one big string.  Then we split it based on
        # equals signs.  Everything in between the equals signs
        # then gets split up based on whitespace.  So, the name-value
        # pairs are [everything but the first thing in group i-1]-
        # [first thing in group i].  We have to take care to get
        # the first and last groups correct.
        # 
        inString = ' '.join(inString.split()) # get rid of spurious whitespace
        parts = inString.split('=')
        names = [parts[0].strip()]
        values = []
        for part in parts[1:-1]:
            try:
                values.append(float(part.split()[0]))
            except ValueError:
                values.append(0.0)
            names.append(' '.join(part.split()[1:]))
        values.append(float(parts[-1]))
        if len(names) != len(values):
            sys.exit('oops')
        for name,value in zip(names,values):
            self[name] = value
    def __setitem__(self,name,value):
        self.allKeys[name] = True
        return super(StepInfo,self).__setitem__(name,value)
    def __getitem__(self,item):
        #
        # If it's a "valid key" (meaning that some StepInfo has seen it before),
        # and we haven't seen it, we'll return the empty string.
        #
        if item in self.allKeys:
            val = ''
            try:
                val = super(StepInfo,self).__getitem__(item)
            except KeyError:
                pass
        else:
            val = super(StepInfo,self).__getitem__(item)
        return val
    def keys(self):
        return list(self.allKeys.keys())
    def iteritems(self):
        #
        # we need to make sure that this includes the fake keys in self.allKeys
        #
        for k in list(self.keys()):
            yield k,self[k]
        

def steps(fileobj):
    #
    # When we get to the RMS fluctuations, we're done.
    #
    step = []
    for line in fileobj:
        if line.startswith('      R M S  F L U C T U A T I O N S'):
            print("BREAKING")
            break
        if line.startswith(' NSTEP'):
            if step:
                yield ''.join(step)
            step = [line]
        else:
            #
            # Only include lines that we know should be included
            #
            goodStarts = (' NSTEP', ' Etot', ' BOND', ' 1-4 NB', ' EELEC', ' EKCMT','                                                Density',  ' EAMBER (non-constraint)',)
            #goodStarts = (' NSTEP', ' Etot', ' BOND', ' 1-4 NB', ' EELEC', ' EKCMT','                                                Density',)
            #
            # Only some lines have EAMBER .. StepInfo deals with that.
            #
            for start in goodStarts:
                if line.startswith(start):
                    step.append(line)
    if step:
        yield(''.join(step))


def getSteps(filename):
    f = file(filename)
    for step in steps(f):
        yield StepInfo(step)
    f.close()

def getFilename(s,extension = '.txt',dir=''):
    """everything that's not a letter or number turns into an underscore"""
    return os.path.join(dir,re.sub('\W','_',s) + extension)

usage = """
parseSanderOut.py -f file.out -d dirname

will put the data in file.out into nice files in dirname.
The directory called dirname must not exist when this script is run.
Among those files are:

  data/allout.txt      tab-delimited text file with all information
                       (suitable for gnumeric, koffice or excel)
  data/Etot.txt, etc.  individual files with different types of sander
                       output (two columns per file, one is time(ps))
  results.html         html file showing you lots graphs of your data
  images/*             postscript and gif graphs of your data

So, for the most part, you probably just want to point your favorite
webbrowser at dirname/results.html.

NOTE: this script assumes that you have gnuplot and convert installed
and in your path.
"""

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser(usage=usage)
    parser.add_option('-f','--file',
                      dest='datafileName',
                      default='sander.out',
                      help='The name of the sander output file to parse [default %default]',
                      )
    parser.add_option('-d','--dir',
                      dest='outdir',
                      default='data',
                      help='The name of the output directory (must not exist) [default %default]',
                      )

    options,args = parser.parse_args()
    datafileName, dir = options.datafileName,options.outdir
    try:
        f = file(datafileName)
        f.close()
    except:
        print("****************************************")
        print()
        print("Could not read %s." % datafileName)
        print()
        print("****************************************")
        print()
        sys.exit(usage)
    if os.path.exists(dir):
        print("****************************************")
        print()
        print("%s exists.  Please (re)move it or specify" % dir)
        print("a different destination directory.")
        print()
        print("****************************************")
        print()
        sys.exit(usage)

    os.system('mkdir %s' % dir)
    #
    # We write out alldata.txt, which contains all of the data.
    # We also write out <thing>.txt where <thing> is all of the
    # other types of data, e.g. EELEC.data.
    #
    allout = file(os.path.join(dir,'alldata.txt'),'w')
    # fileMap maps data type (EELEC) to file object and filename.
    # e.g. fileMap = {'EELEC':{'filename':'EELEC.txt','file':<fileobj>},...}
    #
    fileMap = {} 
    #
    # We'll set the keys with the first step and reuse them
    # for all subsequent setps.  This guarantees that we have a consistent
    # number of columns.
    #
    first = True
    nonTimeKeys = None
    for step in getSteps(datafileName):
        #
        # We want TIME(PS) to be the first column
        #
        if nonTimeKeys is None:
            nonTimeKeys = list(step.keys())
            try:
                nonTimeKeys.remove('TIME(PS)')
            except:
                sys.exit('Could not find time(ps) in ' + str(nonTimeKeys))
        if first:
            first = False
            #
            # Set up fileMap
            #
            for key in nonTimeKeys:
                filename = getFilename(key)
                fileMap[key] = {'filename':filename,
                                'file':file(os.path.join(dir,filename),'w'),}
            #
            # Write the headers
            #
            allout.write('# TIME(PS)\t' + '\t'.join(nonTimeKeys) + '\n')
            for key in list(fileMap.keys()):
                fileMap[key]['file'].write('# TIME(PS)\t' + key + '\n')
            
        #
        # line keeps track of the line for allout.txt
        #
        line = '%s\t' % step['TIME(PS)']
        for k in nonTimeKeys:
            if k != 'TIME(PS)':
                line += '%s\t' % step[k]
                fileMap[k]['file'].write('%s\t%s\n' % (step['TIME(PS)'],
                                                       step[k]))
        line += '\n'
        allout.write(line)
    #
    # Close all of the files
    #
    allout.close()
    for k in list(fileMap.keys()):
        fileMap[k]['file'].close()

    #
    # Now we use gnuplot to make a bunch of images
    #
    imgDir = os.path.join(dir,'images')
    os.system('mkdir -p ' + imgDir)
    gpltCmd = 'set terminal postscript\n'
    for k in list(fileMap.keys()):
        imgName = getFilename(k,'.ps',imgDir)
        gpltCmd += "set output '%s'\n" % imgName
        gpltCmd += "plot '%s'\n" % os.path.join(dir,fileMap[k]['filename'])
    gpltFilename = os.path.join(dir,'commands.gnuplot')
    gpltFile = file(gpltFilename,'w')
    gpltFile.write(gpltCmd)
    gpltFile.close()
    os.system('/usr/bin/env gnuplot %s' % gpltFilename)
    #
    # gnuplot can't always make gifs, so we make postscript files and convert
    # them into gifs.
    #
    thumbPart = '-size 200x200 -geometry 200x200'
    for k in list(fileMap.keys()):
        ps = getFilename(k,'.ps',imgDir)
        gif = getFilename(k,'.gif',imgDir)
        thumb = getFilename(k,'_thumb.gif',imgDir)
        os.system('/usr/bin/env convert -rotate 90 %s %s' % (ps,gif))
        os.system('/usr/bin/env convert -rotate 90 %s %s %s' % (thumbPart,
                                                                ps,
                                                                thumb))

    #
    # This will be a really simple html file
    #
    numCols = 4
    htmlFile = file(os.path.join(dir,'results.html'),'w')
    htmlFile.write("<html><head></head><body>")
    htmlFile.write('<table><tr><td colspan="%s">These are the thumbnails .. click on them for the big picture(s)</td>' % numCols)
    count = 0
    for k in list(fileMap.keys()):
        if count % numCols == 0:
            htmlFile.write('</tr><tr>')
        count += 1
        bigName = getFilename(k,'.gif','images')
        thumbName = getFilename(k,'_thumb.gif','images')
        htmlFile.write('<td>%s<br>'%(k))
        htmlFile.write('<a href="#%s"><img src="%s"></a>'%(k,thumbName))
        htmlFile.write('</td>\n')
    htmlFile.write('</tr></table>\n')
    for k in list(fileMap.keys()):
        bigName = getFilename(k,'.gif','images')
        htmlFile.write('<br>')
        htmlFile.write('<a name="%s"><a href="%s"><img src="%s"></a></a>\n' % (k,getFilename(k,'.txt'), bigName))
    htmlFile.write("</body></html>\n")
    htmlFile.close()
