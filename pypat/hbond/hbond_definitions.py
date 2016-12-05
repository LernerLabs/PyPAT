#!/usr/bin/env python

############################
#
# A note to users:
# Feel free to edit this file to add your own
# hbond definitions. Python allows you to use single-
# and double-quotes pretty much interchangeably. In this
# file, however, you should use only double-quotes to quote
# things. Otherwise, it's too easy to mess things up because
# nucleic acids often use single-quotes as part of the atom names
# and we use those in our selections. In particular, I may run an
# automated procedure that adds in the "*" convention as well as the
# "'" convention, and those will break if you use "'" for anything else.
#
############################

from pymol import cmd

#
# RNA Defs from Mark Ditzler
#
def select_rna_acceptors():
    cmd.select("prot_acceptors","prot_acceptors or (resn RG and name O2')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RG and name O3')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RG and name O4')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RG and name O5')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RG and name O1P)")
    cmd.select("prot_acceptors","prot_acceptors or (resn RG and name O2P)")
    cmd.select("prot_acceptors","prot_acceptors or (resn RG and name O6)")
    cmd.select("prot_acceptors","prot_acceptors or (resn RG and name N7)")
    cmd.select("prot_acceptors","prot_acceptors or (resn RG and name N3)")

    cmd.select("prot_acceptors","prot_acceptors or (resn RG3 and name O2')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RG3 and name O3')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RG3 and name O4')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RG3 and name O5')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RG3 and name O1P)")
    cmd.select("prot_acceptors","prot_acceptors or (resn RG3 and name O2P)")
    cmd.select("prot_acceptors","prot_acceptors or (resn RG3 and name O6)")
    cmd.select("prot_acceptors","prot_acceptors or (resn RG3 and name N7)")
    cmd.select("prot_acceptors","prot_acceptors or (resn RG3 and name N3)")

    cmd.select("prot_acceptors","prot_acceptors or (resn RG5 and name O2')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RG5 and name O3')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RG5 and name O4')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RG5 and name O5')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RG5 and name O1P)")
    cmd.select("prot_acceptors","prot_acceptors or (resn RG5 and name O2P)")
    cmd.select("prot_acceptors","prot_acceptors or (resn RG5 and name O6)")
    cmd.select("prot_acceptors","prot_acceptors or (resn RG5 and name N7)")
    cmd.select("prot_acceptors","prot_acceptors or (resn RG5 and name N3)")

    cmd.select("prot_acceptors","prot_acceptors or (resn RC and name O2')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RC and name O3')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RC and name O4')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RC and name O5')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RC and name O1P)")
    cmd.select("prot_acceptors","prot_acceptors or (resn RC and name O2P)")
    cmd.select("prot_acceptors","prot_acceptors or (resn RC and name O2)")
    cmd.select("prot_acceptors","prot_acceptors or (resn RC and name N3)")

    cmd.select("prot_acceptors","prot_acceptors or (resn RC3 and name O2')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RC3 and name O3')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RC3 and name O4')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RC3 and name O5')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RC3 and name O1P)")
    cmd.select("prot_acceptors","prot_acceptors or (resn RC3 and name O2P)")
    cmd.select("prot_acceptors","prot_acceptors or (resn RC3 and name O2)")
    cmd.select("prot_acceptors","prot_acceptors or (resn RC3 and name N3)")

    cmd.select("prot_acceptors","prot_acceptors or (resn RC5 and name O2')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RC5 and name O3')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RC5 and name O4')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RC5 and name O5')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RC5 and name O1P)")
    cmd.select("prot_acceptors","prot_acceptors or (resn RC5 and name O2P)")
    cmd.select("prot_acceptors","prot_acceptors or (resn RC5 and name O2)")
    cmd.select("prot_acceptors","prot_acceptors or (resn RC5 and name N3)")

    cmd.select("prot_acceptors","prot_acceptors or (resn RA and name O2')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RA and name O3')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RA and name O4')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RA and name O5')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RA and name O1P)")
    cmd.select("prot_acceptors","prot_acceptors or (resn RA and name O2P)")
    cmd.select("prot_acceptors","prot_acceptors or (resn RA and name N7)")
    cmd.select("prot_acceptors","prot_acceptors or (resn RA and name N3)")
    cmd.select("prot_acceptors","prot_acceptors or (resn RA and name N1)")

    cmd.select("prot_acceptors","prot_acceptors or (resn RA3 and name O2')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RA3 and name O3')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RA3 and name O4')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RA3 and name O5')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RA3 and name O1P)")
    cmd.select("prot_acceptors","prot_acceptors or (resn RA3 and name O2P)")
    cmd.select("prot_acceptors","prot_acceptors or (resn RA3 and name N7)")
    cmd.select("prot_acceptors","prot_acceptors or (resn RA3 and name N3)")
    cmd.select("prot_acceptors","prot_acceptors or (resn RA3 and name N1)")

    cmd.select("prot_acceptors","prot_acceptors or (resn RA5 and name O2')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RA5 and name O3')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RA5 and name O4')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RA5 and name O5')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RA5 and name O1P)")
    cmd.select("prot_acceptors","prot_acceptors or (resn RA5 and name O2P)")
    cmd.select("prot_acceptors","prot_acceptors or (resn RA5 and name N7)")
    cmd.select("prot_acceptors","prot_acceptors or (resn RA5 and name N3)")
    cmd.select("prot_acceptors","prot_acceptors or (resn RA5 and name N1)")

    cmd.select("prot_acceptors","prot_acceptors or (resn RU and name O2')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RU and name O3')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RU and name O4')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RU and name O5')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RU and name O1P)")
    cmd.select("prot_acceptors","prot_acceptors or (resn RU and name O2P)")
    cmd.select("prot_acceptors","prot_acceptors or (resn RU and name O2)")
    cmd.select("prot_acceptors","prot_acceptors or (resn RU and name O4)")

    cmd.select("prot_acceptors","prot_acceptors or (resn RU3 and name O2')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RU3 and name O3')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RU3 and name O4')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RU3 and name O5')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RU3 and name O1P)")
    cmd.select("prot_acceptors","prot_acceptors or (resn RU3 and name O2P)")
    cmd.select("prot_acceptors","prot_acceptors or (resn RU3 and name O2)")
    cmd.select("prot_acceptors","prot_acceptors or (resn RU3 and name O4)")

    cmd.select("prot_acceptors","prot_acceptors or (resn RU5 and name O2')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RU5 and name O3')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RU5 and name O4')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RU5 and name O5')")
    cmd.select("prot_acceptors","prot_acceptors or (resn RU5 and name O1P)")
    cmd.select("prot_acceptors","prot_acceptors or (resn RU5 and name O2P)")
    cmd.select("prot_acceptors","prot_acceptors or (resn RU5 and name O2)")
    cmd.select("prot_acceptors","prot_acceptors or (resn RU5 and name O4)")
def select_rna_donors():
    """
    We only care about the hydrogen here, but we include the attached
    heavy atom in a comment for completeness.
    """
    #acceptor mask :RG@N1 :RG@H1
    cmd.select("prot_donors","prot_donors or (resn RG and name H1)")
    #acceptor mask :RG@N2 :RG@H21
    cmd.select("prot_donors","prot_donors or (resn RG and name H21)")
    #acceptor mask :RG@N2 :RG@H22
    cmd.select("prot_donors","prot_donors or (resn RG and name H22)")
    #acceptor mask :RG@O2' :RG@HO'2
    cmd.select("prot_donors","prot_donors or (resn RG and name HO'2)")

    #acceptor mask :RG3@N1 :RG3@H1
    cmd.select("prot_donors","prot_donors or (resn RG3 and name H1)")
    #acceptor mask :RG3@N2 :RG3@H21
    cmd.select("prot_donors","prot_donors or (resn RG3 and name H21)")
    #acceptor mask :RG3@N2 :RG3@H22
    cmd.select("prot_donors","prot_donors or (resn RG3 and name H22)")
    #acceptor mask :RG3@O2' :RG3@HO'2
    cmd.select("prot_donors","prot_donors or (resn RG3 and name HO'2)")
    #acceptor mask :RG3@O3' :RG3@H3T
    cmd.select("prot_donors","prot_donors or (resn RG3 and name H3T)")

    #acceptor mask :RG5@N1 :RG5@H1
    cmd.select("prot_donors","prot_donors or (resn RG5 and name H1)")
    #acceptor mask :RG5@N2 :RG5@H21
    cmd.select("prot_donors","prot_donors or (resn RG5 and name H21)")
    #acceptor mask :RG5@N2 :RG5@H22
    cmd.select("prot_donors","prot_donors or (resn RG5 and name H22)")
    #acceptor mask :RG5@O2' :RG5@HO'2
    cmd.select("prot_donors","prot_donors or (resn RG5 and name HO'2)")
    #acceptor mask :RG5@O5' :RG5@H5T
    cmd.select("prot_donors","prot_donors or (resn RG5 and name H5T)")

    #acceptor mask :RC@N4 :RC@H41
    cmd.select("prot_donors","prot_donors or (resn RC and name H41)")
    #acceptor mask :RC@N4 :RC@H42
    cmd.select("prot_donors","prot_donors or (resn RC and name H42)")
    #acceptor mask :RC@O2' :RC@HO'2
    cmd.select("prot_donors","prot_donors or (resn RC and name HO'2)")

    #acceptor mask :RC3@N4 :RC3@H41
    cmd.select("prot_donors","prot_donors or (resn RC3 and name H41)")
    #acceptor mask :RC3@N4 :RC3@H42
    cmd.select("prot_donors","prot_donors or (resn RC3 and name H42)")
    #acceptor mask :RC3@O2' :RC3@HO'2
    cmd.select("prot_donors","prot_donors or (resn RC3 and name HO'2)")
    #acceptor mask :RC3@O3' :RC3@H3T
    cmd.select("prot_donors","prot_donors or (resn RC3 and name H3T)")

    #acceptor mask :RC5@N4 :RC5@H41
    cmd.select("prot_donors","prot_donors or (resn RC5 and name H41)")
    #acceptor mask :RC5@N4 :RC5@H42
    cmd.select("prot_donors","prot_donors or (resn RC5 and name H42)")
    #acceptor mask :RC5@O2' :RC5@HO'2
    cmd.select("prot_donors","prot_donors or (resn RC5 and name HO'2)")
    #acceptor mask :RC5@O5' :RC5@H5T
    cmd.select("prot_donors","prot_donors or (resn RC5 and name H5T)")

    #acceptor mask :RA@N6 :RA@H61
    cmd.select("prot_donors","prot_donors or (resn RA and name H61)")
    #acceptor mask :RA@N6 :RA@H62
    cmd.select("prot_donors","prot_donors or (resn RA and name H62)")
    #acceptor mask :RA@O2' :RA@HO'2
    cmd.select("prot_donors","prot_donors or (resn RA and name HO'2)")

    #acceptor mask :RA3@N6 :RA3@H61
    cmd.select("prot_donors","prot_donors or (resn RA3 and name H61)")
    #acceptor mask :RA3@N6 :RA3@H62
    cmd.select("prot_donors","prot_donors or (resn RA3 and name H62)")
    #acceptor mask :RA3@O2' :RA3@HO'2
    cmd.select("prot_donors","prot_donors or (resn RA and name HO'2)")
    #acceptor mask :RA3@O3' :RA3@H3T
    cmd.select("prot_donors","prot_donors or (resn RA3 and name H3T)")

    #acceptor mask :RA5@N6 :RA5@H61
    cmd.select("prot_donors","prot_donors or (resn RA5 and name H61)")
    #acceptor mask :RA5@N6 :RA5@H62
    cmd.select("prot_donors","prot_donors or (resn RA5 and name H62)")
    #acceptor mask :RA5@O2' :RA5@HO'2
    cmd.select("prot_donors","prot_donors or (resn RA and name HO'2)")
    #acceptor mask :RA5@O5' :RA5@H5T
    cmd.select("prot_donors","prot_donors or (resn RA5 and name H5T)")

    #acceptor mask :RU@N3 :RU@H3
    cmd.select("prot_donors","prot_donors or (resn RU and name H3)")
    #acceptor mask :RU@O2' :RU@HO'2
    cmd.select("prot_donors","prot_donors or (resn RU and name HO'2)")

    #acceptor mask :RU3@N3 :RU3@H3
    cmd.select("prot_donors","prot_donors or (resn RU3 and name H3)")
    #acceptor mask :RU3@O2' :RU3@HO'2
    cmd.select("prot_donors","prot_donors or (resn RU3 and name HO'2)")
    #acceptor mask :RU3@O3' :RU3@H3T
    cmd.select("prot_donors","prot_donors or (resn RU3 and name H3T)")

    #acceptor mask :RU5@N3 :RU5@H3
    cmd.select("prot_donors","prot_donors or (resn RU5 and name H3)")
    #acceptor mask :RU5@O2' :RU5@HO'2
    cmd.select("prot_donors","prot_donors or (resn RU5 and name HO'2)")
    #acceptor mask :RU5@O5' :RU5@H5T
    cmd.select("prot_donors","prot_donors or (resn RU5 and name H5T)")

def select_standard_prot_donors_and_acceptors():
    """
    This is not standard because
    1) it knows about our terminal residues
    2) it includes our ligands
    """
    cmd.select("prot_acceptors","resn GLN and name OE1")
    cmd.select("prot_acceptors","prot_acceptors or (resn GLN and name OE1)")
    cmd.select("prot_acceptors","prot_acceptors or (resn GLN and name NE2)")
    cmd.select("prot_acceptors","prot_acceptors or (resn ASN and name OD1)")
    cmd.select("prot_acceptors","prot_acceptors or (resn ASN and name ND2)")
    cmd.select("prot_acceptors","prot_acceptors or (resn TYR and name OH)")
    cmd.select("prot_acceptors","prot_acceptors or (resn ASP and name OD1)")
    cmd.select("prot_acceptors","prot_acceptors or (resn ASP and name OD2)")
    cmd.select("prot_acceptors","prot_acceptors or (resn GLU and name OE1)")
    cmd.select("prot_acceptors","prot_acceptors or (resn GLU and name OE2)")
    cmd.select("prot_acceptors","prot_acceptors or (resn SER and name OG)")
    cmd.select("prot_acceptors","prot_acceptors or (resn THR and name OG1)")
    cmd.select("prot_acceptors","prot_acceptors or (resn HIS and name ND1)")
    cmd.select("prot_acceptors","prot_acceptors or (resn HIE and name ND1)")
    cmd.select("prot_acceptors","prot_acceptors or (resn HID and name NE2)")
    
    cmd.select("prot_donors","resn ASN and name HD21")
    cmd.select("prot_donors","prot_donors or (resn ASN and name 1HD2)")
    #acceptor mask  :ASN@ND2 :ASN@HD21
    cmd.select("prot_donors","prot_donors or (resn ASN and name HD21)")
    cmd.select("prot_donors","prot_donors or (resn ASN and name 1HD2)")
    #acceptor mask  :ASN@ND2 :ASN@HD22
    cmd.select("prot_donors","prot_donors or (resn ASN and name HD22)")
    cmd.select("prot_donors","prot_donors or (resn ASN and name 2HD2)")
    #acceptor mask  :TYR@OH  :TYR@HH
    cmd.select("prot_donors","prot_donors or (resn TYR and name HH)")
    #acceptor mask  :GLN@NE2 :GLN@HE21
    cmd.select("prot_donors","prot_donors or (resn GLN and name HE21)")
    cmd.select("prot_donors","prot_donors or (resn GLN and name 1HE2)")
    #acceptor mask  :GLN@NE2 :GLN@HE22
    cmd.select("prot_donors","prot_donors or (resn GLN and name HE22)")
    cmd.select("prot_donors","prot_donors or (resn GLN and name 2HE2)")
    #acceptor mask  :TRP@NE1 :TRP@HE1
    cmd.select("prot_donors","prot_donors or (resn TRP and name HE1)")
    #acceptor mask  :LYS@NZ  :LYS@HZ1
    cmd.select("prot_donors","prot_donors or (resn LYS and name HZ1)")
    #acceptor mask  :LYS@NZ  :LYS@HZ2
    cmd.select("prot_donors","prot_donors or (resn LYS and name HZ2)")
    #acceptor mask  :LYS@NZ  :LYS@HZ3
    cmd.select("prot_donors","prot_donors or (resn LYS and name HZ3)")
    #acceptor mask  :SER@OG  :SER@HG
    cmd.select("prot_donors","prot_donors or (resn SER and name HG)")
    #acceptor mask  :THR@OG1 :THR@HG1
    cmd.select("prot_donors","prot_donors or (resn THR and name HG1)")
    #acceptor mask  :ARG@NH2 :ARG@HH21
    cmd.select("prot_donors","prot_donors or (resn ARG and name HH21)")
    cmd.select("prot_donors","prot_donors or (resn ARG and name 1HH2)")
    #acceptor mask  :ARG@NH2 :ARG@HH22
    cmd.select("prot_donors","prot_donors or (resn ARG and name HH22)")
    cmd.select("prot_donors","prot_donors or (resn ARG and name 2HH2)")
    #acceptor mask  :ARG@NH1 :ARG@HH11
    cmd.select("prot_donors","prot_donors or (resn ARG and name HH11)")
    #acceptor mask  :ARG@NH1 :ARG@HH12
    cmd.select("prot_donors","prot_donors or (resn ARG and name HH12)")
    #acceptor mask  :ARG@NE  :ARG@HE
    cmd.select("prot_donors","prot_donors or (resn ARG and name HE)")
    #acceptor mask  :HIS@NE2 :HIS@HE2
    cmd.select("prot_donors","prot_donors or (resn HIS and name HE2)")
    #acceptor mask  :HIE@NE2 :HIE@HE2
    cmd.select("prot_donors","prot_donors or (resn HIE and name HE2)")
    #acceptor mask  :HID@ND1 :HID@HD1
    cmd.select("prot_donors","prot_donors or (resn HID and name HD1)")
    #acceptor mask  :HIP@ND1,NE2 :HIP@HE2,HD1
    cmd.select("prot_donors","prot_donors or (resn HIP and name HIP@HE2,HD1)")
    #-- Backbone donors and acceptors for this particular molecule
    #   N-H for prolines do not exist so are not in the mask.
    #  
    
    
    #donor mask @O
    cmd.select("prot_acceptors","prot_acceptors or (name o and not resn WAT+HOH)")
    #   In our case, prolines are residues 21, 25, 39, 53, 55, 
    #   66, 89, 105, 126, 130.  We would say 1-159, but we exclude
    #   prolines and terminii.
    #
    #acceptor mask :2-20,22-24,26-38,40-52,54,56-65,67-88,90-104,106-125,127-129,131-158@N :1-158@H
    cmd.select("prot_donors","prot_donors or (name H and resi 2-158)")
    #Terminal residues have different atom names
    #donor mask @OXT
    cmd.select("prot_acceptors","prot_acceptors or name OXT")
    #acceptor mask :1@N :1@H1
    #acceptor mask :1@N :1@H2
    #acceptor mask :1@N :1@H3
    cmd.select("prot_donors","prot_donors or (resi 1 and name H1+H2+H3)")
    
    
def select_nap_donors_and_acceptors():
    """
    Ligand specific selections for NADPH (NAP)
    """
    #-- NADPH
    #acceptor mask :NAP@N6A  :NAP@H61
    cmd.select("prot_donors","prot_donors or (resn NAP and name H61)")
    #acceptor mask :NAP@N6A  :NAP@H62
    cmd.select("prot_donors","prot_donors or (resn NAP and name H62)")
    #acceptor mask :NAP@O'A3 :NAP@HOA3
    cmd.select("prot_donors","prot_donors or (resn NAP and name HOA3)")
    #acceptor mask :NAP@O'N3 :NAP@HON3
    cmd.select("prot_donors","prot_donors or (resn NAP and name HON3)")
    #acceptor mask :NAP@O'N2 :NAP@HON2
    cmd.select("prot_donors","prot_donors or (resn NAP and name HON2)")
    #acceptor mask :NAP@N7N  :NAP@H72
    cmd.select("prot_donors","prot_donors or (resn NAP and name H72)")
    #acceptor mask :NAP@N7N  :NAP@H71
    cmd.select("prot_donors","prot_donors or (resn NAP and name H71)")
    cmd.select("prot_acceptors","prot_acceptors or (resn NAP and name N1A)")
    cmd.select("prot_acceptors","prot_acceptors or (resn NAP and name N3A)")
    cmd.select("prot_acceptors","prot_acceptors or (resn NAP and name N7A)")
    cmd.select("prot_acceptors","prot_acceptors or (resn NAP and name OA23)")
    cmd.select("prot_acceptors","prot_acceptors or (resn NAP and name OA22)")
    cmd.select("prot_acceptors","prot_acceptors or (resn NAP and name OA24)")
    cmd.select("prot_acceptors","prot_acceptors or (resn NAP and name O'A2)")
    cmd.select("prot_acceptors","prot_acceptors or (resn NAP and name O'A3)")
    cmd.select("prot_acceptors","prot_acceptors or (resn NAP and name O'A4)")
    cmd.select("prot_acceptors","prot_acceptors or (resn NAP and name O'A5)")
    cmd.select("prot_acceptors","prot_acceptors or (resn NAP and name OPA1)")
    cmd.select("prot_acceptors","prot_acceptors or (resn NAP and name OPA2)")
    cmd.select("prot_acceptors","prot_acceptors or (resn NAP and name OPN1)")
    cmd.select("prot_acceptors","prot_acceptors or (resn NAP and name O3P)")
    cmd.select("prot_acceptors","prot_acceptors or (resn NAP and name OPN2)")
    cmd.select("prot_acceptors","prot_acceptors or (resn NAP and name O'N5)")
    cmd.select("prot_acceptors","prot_acceptors or (resn NAP and name O'N4)")
    cmd.select("prot_acceptors","prot_acceptors or (resn NAP and name O'N3)")
    cmd.select("prot_acceptors","prot_acceptors or (resn NAP and name O'N2)")
    cmd.select("prot_acceptors","prot_acceptors or (resn NAP and name O7N)")
    

    
################################

def do_standard_selections():
    select_standard_prot_donors_and_acceptors()
    select_rna_donors()
    select_rna_acceptors()

cmd.extend("do_standard_selections",do_standard_selections)
