#!/usr/bin/env python

import sys, os

def parse_residue_list_file(filename = None):
    '''
    Receives a file name and parses that file to create a dictionary
    of residue lists associated with the lists' names, which can then
    be used in function parse_residue_list().  
    
    Although it looks like 'None' is the default for the filename 
    keyword, the function will look for a file named 'residue_lists' 
    as the default.  The reason for the None is that the warning 
    message output by this function if the file was not found is 
    suppressed if no name was actually passed to the function. 

    The 'residue_lists' file should contain lines with the following
    pattern:
    NAME : RESIDUE_STRING
    where NAME is the name to be associated with the residue list
    and RESIDUE_STRING is a comma-separated string that can be 
    parsed to create a list of residue names.  The residue string 
    may contain any or all of the following:
	* individual residue numbers
	* sequences of numbers, with the first and last numbers 
	    separated by a '-'
	* strings associated with residue lists that have
	    already been defined in the file
    '''
    if not filename:
	filename = 'residue_lists'
	warning_str = ''
    else:
	warning_str = 'Warning:  Could not open file ' + filename + \
		   '.\n  No standard residue lists are available.'

    try:
	f = file(filename)
    except:
	print warning_str
	return {}

    standard_residue_lists = {}

    for line in f:
	if not line.strip() or line.strip().startswith('#'):
	    continue
	try:
	    resi_list = []
	    key, string_to_parse = line.split(':')
	    key = key.strip()
	    for piece in string_to_parse.strip().split(','):
		if standard_residue_lists.has_key(piece):
		    resi_list.extend(standard_residue_lists[piece])
		elif '-' in piece:
		    first, last = piece.split('-')
		    resi_list.extend(range(int(first), int(last) + 1))
		else:
		    resi_list.append(int(piece))
	except:
	    print 'Warning:  Ignoring line in file ' + filename + '\n' + \
		  '  because of improper syntax:\n' + line
	resi_list = sorted(list(set(resi_list)))
	standard_residue_lists[key] = resi_list

    return standard_residue_lists


def parse_residue_string(string_to_parse, offset = 0, filename = None):
    '''
    Receives a comma-separated string and returns a list of residues.  
    The string may contain any or all of the following:
	* individual residue numbers
	* a range of numbers, separated by a '-'
	* strings associated with valid residue lists in the 
	    standard file 'residue_lists' (a different file name
	    can be passed to this function using the filename keyword)
    A list of individual residue numbers will be printed out.

    Option offset can be used to create, e.g., a zero-based list
    for indexing.
    '''
    standard_residue_lists = parse_residue_list_file(filename)

    piece_list = string_to_parse.strip().split(',')
    if 'all' in piece_list:
	return ['all']

    try:
        resi_list = []
        for piece in piece_list:
            if standard_residue_lists.has_key(piece):
                for resi in standard_residue_lists[piece]:
		    resi_list.append(resi - offset)
            elif '-' in piece:
                first, last = piece.split('-')
                resi_list.extend(range(int(first) - offset, int(last) + 1 - offset))
            else:
                resi_list.append(int(piece) - offset)
    except:
	print 'Warning:  Residue list string ' + string_to_parse + \
	      ' has improper syntax.'
	return []

    resi_list = sorted(list(set(resi_list)))

    output_str = 'Residue list ' + string_to_parse + ' represents ' + \
		  str(len(resi_list)) + ' residues:\n'
    for resi in resi_list:
	output_str += str(resi) + ','
    output_str = output_str[:-1] + '\n'
    print output_str

    return resi_list


def parse_residue_list(string_to_parse, offset = 0, filename = None):
    '''
    Receives a string and returns one or two lists of residues.  If
    the string contains a ':', the strings on both sides of the colon
    are parsed by parse_residue_string and the two lists are returned.
    If no colon is present, the string itself is parsed and a single
    list is returned.   
    '''
    string_pieces = string_to_parse.strip().split(':')

    if len(string_pieces) == 1:
	return parse_residue_string(string_pieces[0], offset, filename)
    elif len(string_pieces) == 2:
	return parse_residue_string(string_pieces[0], offset, filename), \
	       parse_residue_string(string_pieces[1], offset, filename)
    else:
	print 'Warning:  Residue list string ' + string_to_parse + \
	      ' has improper syntax.'
	return []


def parse_atom_list(string_to_parse):
    '''
    Receives a string and splits by commas.  Used to parse the list of 
    atoms given to various hbond programs.
    '''
    piece_list = string_to_parse.strip().split(',')
    return piece_list
   
