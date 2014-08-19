#! /usr/bin/env python
from prody import *
import hamiltonian
import argparse

def main():
    #no log messages:
    prody.confProDy(verbosity='none')

    # parse command line arguments
    parser = argparse.ArgumentParser(description='Calculate MDENM energies from a pdb \
                                    will calculate energy using modes from pdb\
                                    and then from reference--> crystal should\
                                    be the reference')
    parser.add_argument('--pdb', required=True,
                        help='Molecule we want to examine. It will also be used as topology')
    parser.add_argument('--reference', required=True,
                        help='Rerence pdb to which we will rmsd everything')
    args = parser.parse_args()


    #Load the structures
    pdb = parsePDB(args.pdb)
    calphas = pdb.select('calpha')
    ref = prody.parsePDB(args.reference)
    ref_alpha = ref.select('calpha')
    #Make sure we are in same reference set
    t = calcTransformation(calphas,ref_alpha)
    t.apply(calphas)

    native = ref_alpha.copy()
    pred = calphas.copy()

    h = hamiltonian.EDENMHamiltonian( native.getCoords() )
    Forw_E_ED = h.evaluate_energy( pred.getCoords())

    h = hamiltonian.EDENMHamiltonian( pred.getCoords() )
    Back_E_ED = h.evaluate_energy( native.getCoords())

    number_of_residues = ref.numResidues()
    rmsdED = calcRMSD(calphas,target=ref_alpha)

    print "%s %.2f %.2f %.2f " % (args.pdb, rmsdED, Forw_E_ED / number_of_residues, Back_E_ED / number_of_residues),


if __name__ == '__main__':
    main()

