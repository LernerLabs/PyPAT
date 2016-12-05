#!/usr/bin/env python

import sys,os,glob

if __name__ == '__main__':
    from pypat import tool_utils
    import glob
    from optparse import OptionParser

    usage = tool_utils.usage + """
Please also make sure that the imagemagick utility convert
is installed and in your path.

This will spit out an html file that will show you your images.  If you
need to look at the images on another machine, tar up the html file and
output-dir/images together and move that to the other machine.

The html file will be called <struct>BigAnimatedMovies.html.
"""
    parser = OptionParser(option_class=tool_utils.MyOption,usage=tool_utils.usage)

    tool_utils.add_standard_options(parser)
    parser.add_option('--plot-types',dest="plottypes",
                      default='ca avg max min abs straight mainheavy allheavy sidechainhbond hbond'.split(),
                      type="strlist",
                      help="Comma-separated list of plot types.  [default: %default]",
                      )
    parser.add_option('--no-slow-movies',dest="slowmovies",
                      default=True,
                      action="store_false",
                      help="Set this if you do not want to generate the movies that have 0.5s spacing between the frames.",
                      )
    parser.add_option('--movie-link',dest="movielink",
                      default='fast',
                      help="'fast' if you want the thumbnails to link to the fast images, anything else for the slow ones. [default: %default]",
                      )
    options,args = parser.parse_args()



    html_txt = '''
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
<html>
  <head>
    <title>Correlated Dynamics Movies</title>
  </head>
  
  <body>
    <h1>Correlated Dynamics Movies</h1>
    
    These are animated gifs of a sliding 1NS window throughout the MD simulation.  So, e.g., %s 2.2NS is the correlated dynamics of the %s structure from 1.7NS to 1.7NS (the center of the window is 1.2NS).
    <ul>
      <li><b>straight</b> every atom is shown</li>
      <li><b>mainheavy</b> mainchain heavy (nonhydrogen) atoms</li>
      <li><b>allheavy</b> all heavy (nonhydrogen) atoms</li>
      <li><b>hbond</b> NOSP vs. hydrogen atoms</li>
      <li><b>sidechainhbond</b> sidechain NOSP vs. sidechain hydrogen atoms</li>
      <li><b>ca</b> alpha carbons plotted against eachother, one per residue.  This is the standard plot in the literature.</li>
      <li><b>abs</b> the largest absolute value for each residue-residue pair</li>
      <li><b>min</b> the minimum value for each residue-residue pair</li>
      <li><b>max</b> the maximum value for each residue-residue pair</li>
    </ul>

    The full-sized animated gifs, which you can see by clicking on the versions shown here,
    will play more slowly than the thumbnails.  If you wish to see the speeded-up version,
    remove the "0.5" from the filename.  E.g. look at animated_%s_resi_mainheavy_correl.gif
    instead of animated_0.5_%s_resi_mainheavy_correl.gif.
    
    <table>
      <tr>
	<td><b>%s</b> (closed loop starting structure)</td>
      </tr>
      ''' % (options.structurename,
             options.structurename,
             options.structurename,
             options.structurename,
             options.structurename,)




    for plot_type in options.plottypes:
        # New naming conventions mean that these are zero-padded and will be in the correct order.
        #filenames =  glob.glob(os.path.join(options.outputdir,'images',options.structurename+" NS??? resi "+plot_type+" correl*",)) +  glob.glob(os.path.join(options.outputdir,'images',options.structurename+" NS???? resi "+plot_type+" correl*",)) +  glob.glob(os.path.join(options.outputdir,'images',options.structurename+" NS????? resi "+plot_type+" correl*",)) #does anyone do any 100ns+ simulations??
        filenames =  glob.glob(os.path.join(options.outputdir,'images',options.structurename+" NS* resi "+plot_type+" correl*",))
        if not filenames:
            filenames =  glob.glob(os.path.join(options.outputdir,'images',options.structurename+"* ns*"+plot_type+" correl*",))
        if not filenames:
            print "COULD NOT FIND files for",plot_type
            continue
        filenames.sort()
        prog = 'convert'
        print filenames
        #sys.exit()
        print [os.path.join(options.outputdir,'images',"animated_"+options.structurename+"_resi_"+plot_type+"_correl.gif",),]
        args = ['-loop','0',] + filenames + [os.path.join(options.outputdir,'images',"animated_"+options.structurename+"_resi_"+plot_type+"_correl.gif",),]
        tool_utils.run(prog,args,verbose=True)

        prog,args = 'convert',('-resize','256x',
                               os.path.join(options.outputdir,'images',"animated_"+options.structurename+"_resi_"+plot_type+"_correl.gif",),
                               os.path.join(options.outputdir,'images',"animated_"+options.structurename+"_resi_"+plot_type+"_correl_thumb.gif",),
                               )
        tool_utils.run(prog,args,verbose=True)
        if options.slowmovies:
            prog,args = 'convert',('-delay','50',
                                   os.path.join(options.outputdir,'images',"animated_"+options.structurename+"_resi_"+plot_type+"_correl.gif",),
                                   os.path.join(options.outputdir,'images',"animated_0.5_"+options.structurename+"_resi_"+plot_type+"_correl.gif",),
                                   )
        tool_utils.run(prog,args,verbose=True)

        if (options.movielink == 'fast') or not options.slowmovies:
            html_txt += '''      <tr>
            <td><a href="images/animated_%s_resi_%s_correl.gif"><img src="images/animated_%s_resi_%s_correl_thumb.gif/"/></a></td>
          </tr>\n'''%(
                      options.structurename,
                      plot_type,

                      options.structurename,
                      plot_type,
                      )
        else:
            html_txt += '''      <tr>
            <td><a href="images/animated_0.5_%s_resi_%s_correl.gif"><img src="images/animated_%s_resi_%s_correl_thumb.gif/"/></a></td>
          </tr>\n'''%(
                      options.structurename,
                      plot_type,

                      options.structurename,
                      plot_type,
                      )
            
    html_txt += '''    </table>
    <hr>
  </body>
</html>
'''
    f = file(os.path.join(options.outputdir,options.structurename+'BigAnimatedMovies.html'),'w')
    f.write(html_txt)
    f.close()

