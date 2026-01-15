antechamber -i tfsi.pdb -fi pdb -o tfsi_gaff2.mol2 -fo mol2 -c bcc -nc -1 -at gaff2
parmchk2 -i tfsi_gaff2.mol2 -f mol2 -o tfsi.frcmod
tleap -f tleap.in