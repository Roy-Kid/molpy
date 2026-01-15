# %%
#import *
import numpy as np
import os
import nglview as ngl
#from ase.io import read, write
WDIR='~/GroPolBul/tutorial/'
#%mkdir {WDIR}
%cd {WDIR}

#Desired polymer electrolyte system
#%mkdir PEO_CH3
%cd PEO_CH3

# %% [markdown]
# # Step-1: Parameterizing monomer

# %%
#Use shor polymer with repeat group, head and tail group
#repeat group: COCCOCCOCCOCCOC
#head and tail group: C
#Use antechaber (ambertools) to optimize and prameterise with AM1-BCC charges 
#antechamber options help:
#Usage: antechamber -i     input file name
#                   -fi    input file format
#                   -o     output file name
#                   -fo    output file format
#                   -c     charge method
#                   -nc    net molecular charge (int)
#                   -rn    residue name

!antechamber -i PEO_initial.pdb -fi pdb -o PEO.ac -fo ac -at gaff -an y -c bcc -nc 0 -rn PEO 


# %%
#Visualize the molecule
mon=read('sqm.pdb')
vi=ngl.show_ase(mon);vi.add_label(radius=2,color='black',label_type='atomindex')
vi


# %% [markdown]
# # Step-2: Defining the HEAD, CHAIN, and TAIL of monomer

# %%
#Defining CHAIN HEAD and TAIL in the monomer directory
ac=open('PEO.ac',mode='r') #Reading .ac file
[next(ac) for _ in range(2)] #Skipping first two lines of text
l=ac.readlines() #Reading lines

#Breaking the parts of monomer to CHAIN, HEAD and TAIL
#Atom index where head and tail of monomer; Check from above ngl view of mol
#Change these values accordingly
head_id=1
tail_id=7

head_omit=[3, 23, 24, 25] #Atoms to omit near head
tail_omit=[8, 32, 33, 34] #Atoms to omit near tail


# %%
#Write the atom indexes of CHAIN, HEAD and TAIL
chain=open('PEO.chain','w+');head=open('PEO.head','w+');tail=open('PEO.tail','w+')

chain.write('HEAD_NAME '+str(l[head_id].split()[2])+'\n')
tail.write('HEAD_NAME '+str(l[head_id].split()[2])+'\n')
chain.write('TAIL_NAME '+str(l[tail_id].split()[2])+'\n')
head.write('TAIL_NAME '+str(l[tail_id].split()[2])+'\n')

for i in range(len(head_omit)):
	chain.write('OMIT_NAME '+str(l[int(head_omit[i])].split()[2])+'\n')
	tail.write('OMIT_NAME '+str(l[int(head_omit[i])].split()[2])+'\n')
chain.write('PRE_HEAD_TYPE '+str(l[tail_id].split()[9])+'\n')
tail.write('PRE_HEAD_TYPE '+str(l[tail_id].split()[9])+'\n')
tail.write('CHARGE 0')

for i in range(len(tail_omit)):
	chain.write('OMIT_NAME '+str(l[int(tail_omit[i])].split()[2])+'\n')
	head.write('OMIT_NAME '+str(l[int(tail_omit[i])].split()[2])+'\n')
chain.write('POST_TAIL_TYPE '+str(l[head_id].split()[9])+'\n')
chain.write('CHARGE 0')
head.write('POST_TAIL_TYPE '+str(l[head_id].split()[9])+'\n')
head.write('CHARGE 0')

chain.close();head.close();tail.close()

#Use prepgen to prepare the CHAIN, HEAD and TAIL res files
#Adds dummy atoms at the desired positions
!prepgen -i PEO.ac -o PEO.prepi -f prepi -m PEO.chain -rn PEO -rf PEO.res 
!prepgen -i PEO.ac -o HPT.prepi -f prepi -m PEO.head -rn HPT -rf HPT.res 
!prepgen -i PEO.ac -o TPT.prepi -f prepi -m PEO.tail -rn TPT -rf TPT.res 

# %% [markdown]
# # Step-3: Build and parameterize single polymer chain

# %%
#Build the desired homopolymer using above prepi files
#%mkdir polymer
n_mono_repeat=5
n_mono_pol=25

repeat=" ".join(['PEO'] * int(int(n_mono_pol/n_mono_repeat)-2))
print('HPT '+str(repeat)+' TPT')

#Write the tleap input file to combine the preparatory files and build polymer chain
tleap=open('polymer/PEO_tleap.in','w+')
tleap.write('''source leaprc.gaff
loadamberprep PEO.prepi
loadamberprep HPT.prepi
loadamberprep TPT.prepi
mol = sequence {HPT '''+str(repeat)+''' TPT}
savepdb mol polymer/PEO_'''+str(n_mono_pol)+'''mer.pdb
saveamberparm mol polymer/PEO_'''+str(n_mono_pol)+'''mer.prmtop polymer/PEO_'''+str(n_mono_pol)+'''mer.inpcrd
quit''')
tleap.close()

!tleap -s -f polymer/PEO_tleap.in > polymer/PEO_tleap.out

# %%
#Converting AMBER to GROMACS using "intermol".
!python ../../convert.py --amb_in polymer/PEO_{n_mono_pol}mer.inpcrd polymer/PEO_{n_mono_pol}mer.prmtop --gromacs 

# %%
#Visualize the molecule
mol=read('polymer/PEO_25mer.pdb')
vi=ngl.show_ase(mol);#vi.add_label(radius=2,color='black',label_type='atomindex')
vi


# %%



