#!/usr/bin/env python

'''
make all of our correl and covar plots with scipy.
'''

import sys,os

if sys.platform == 'darwin':
    PYPAT_CODE_DIR = '/Users/mglerner/nethome/work/PyPAT/code'
elif sys.platform == 'linux2':
    #PYPAT_CODE_DIR = '/users/mlerner/work/src/PyPAT/code'
    PYPAT_CODE_DIR = '/v/bigbox2/home/mglerner/work/PyPAT/code'
    
sys.path.append(PYPAT_CODE_DIR)

if __name__ == '__main__':
    from pypat.plotting import make_correl_plots_for_movie, make_one_correl_plot
    from pypat import tool_utils
    from optparse import OptionParser
    import pylab

    usage = tool_utils.usage + """
This will spit out an html file that will show you your images.  If you
need to look at the images on another machine, tar up the html file and
output-dir/images together and move that to the other machine.
"""
    parser = OptionParser(option_class=tool_utils.MyOption,usage=usage)

    tool_utils.add_standard_options(parser)
    tool_utils.add_window_options(parser)

    parser.add_option("--cmap",dest="cmap",default="Normal",
                      help="Color map to use when making the plots. Our custom cmaps are Normal and Scaled.  Standard matplotlib cmaps %s are also supported."%[m for m in pylab.cm.datad.keys() if not m.endswith("_r")] + "[default: %default]",
                      )
                      
    parser.add_option('--plot-types',dest="plottypes",
                      default='ca avg max min abs straight mainheavy allheavy sidechainhbond hbond'.split(),
                      type="strlist",
                      help="Comma-separated list of plot types.  [default: %default]")

    parser.add_option('--single-time',dest='singletime',
                      default=None,
                      help="If you specify a single time here, we will ignore the windowing options and display an interactive plot of the time you specify.  The time should be, e.g., NS0.5",
                      )
    parser.add_option('--mark-resis',dest='markresis',
                      default=[],
                      type="zerobasedintlist",
                      help="A list of residues to mark on the plots.  [default: %default]")
    parser.add_option('--highlight',dest='highlight',
                      default=0.2,
                      type='float',
                      help="How strongly to highlight the marked residues.  Note that --highlight-mode tells us how exactly we will do the highlighting.  0.1 and 0.2 are decent values if you want to use this feature for most plots, although you'll need something stronger for the absolute value plots. [default: %default]")
    parser.add_option('--highlight-mode',dest='highlightmode',
                      default='positive',
                      help="When highlight-mode is 'negative' we put a white block down on top of the marked residues, the opacity of which is controled by --highlight.  When it's 'positive', we put that white block down on squares of residues that are *not* highlighted instead.  When it's 'supernegative', we do just like 'negative' except that the block will be twice as opaque where the highlighted rows and columns intersect.  Positive and supernegative seem to be more useful than negative.  [default: %default]")
    parser.add_option('--skip-resis',dest='skipresis',
                      default=[],
                      #
                      # Note to programmers: setup_ca expects this to come in as 1-based.
                      #
                      type="onebasedintlist",
                      help="A list of residues that will be skipped in the plots.  [default: %default]")
    parser.add_option('--no-ticks',dest='ticks',
                      default=True,
                      action='store_false',
                      help="do not include tick marks on the axes",
                      )
    parser.add_option('--dpi',dest='dpi',
                      default=200,
                      type='int',
                      help='dpi for figures [default: %default]',
                      )
    parser.add_option('--title',dest='title',
                      default=None,
                      help='If no title is specified, one will be automatically generated. Note that the title is part of the filename that we save. [default: %default]',
                      )

    options,args = parser.parse_args()
    desired=tool_utils.get_desired(options)

    #
    # If you want to make all of the plots for a given time and you want to
    # write them out to files, use this:
    if options.singletime is None:
        all_times = [i[-1] for i in desired]
        make_correl_plots_for_movie(structures=[options.structurename,],
                                    all_times=all_times,
                                    cmaps=[options.cmap,],
                                    detail_levels=['fine',],
                                    plot_types=options.plottypes,
                                    overwrite=True,
                                    image_dir=os.path.join(options.outputdir,'images'),
                                    dat_dir=options.outputdir,
                                    ref_pdb_fname=os.path.join(options.outputdir,options.structurename+'_ref.pdb.1'),
                                    mark_resis=options.markresis,
                                    highlight=options.highlight,
                                    highlight_mode=options.highlightmode,
                                    skip_resis=options.skipresis,
                                    ticks=options.ticks,
                                    dpi=options.dpi,
                                    title=options.title,
                                    )
    #
    # If you just want to make one plot and you want it to show up in a window,
    # use this:
    #
    else:
        if len(options.plottypes) != 1:
            sys.exit("In single-plot mode, you must specify exactly one plot type.")
        print "title",options.title
        make_one_correl_plot(struct=options.structurename,
                             times=options.singletime,
                             cmap=options.cmap,
                             detail_level='fine',
                             plot_type=options.plottypes[0],
                             save_fig=False,
                             overwrite=True,
                             dat_dir=options.outputdir,
                             out_dir=os.path.join(options.outputdir,'images'),
                             ref_pdb_fname=os.path.join(options.outputdir,options.structurename+'_ref.pdb.1'),
                             mark_resis=options.markresis,
                             highlight=options.highlight,
                             highlight_mode=options.highlightmode,
                             skip_resis=options.skipresis,
                             ticks=options.ticks,
                             dpi=options.dpi,
                             title=options.title,
                             )
    
