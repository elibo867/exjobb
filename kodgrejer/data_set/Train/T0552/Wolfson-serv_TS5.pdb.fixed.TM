
 *****************************************************************************
 *                                 TM-SCORE                                  *
 * A scoring function to assess the similarity of protein structures         *
 * Based on statistics:                                                      *
 *       0.0 < TM-score < 0.17, random structural similarity                 *
 *       0.5 < TM-score < 1.00, in about the same fold                       *
 * Reference: Yang Zhang and Jeffrey Skolnick, Proteins 2004 57: 702-710     *
 * For comments, please email to: zhng@umich.edu                             *
 *****************************************************************************

Structure1: input/casp  Length=   98
Structure2: input/casp  Length=   98 (by which all scores are normalized)
Number of residues in common=   98
RMSD of  the common residues=   16.023

TM-score    = 0.1866  (d0= 3.61)
MaxSub-score= 0.1245  (d0= 3.50)
GDT-TS-score= 0.1888 %(d<1)=0.0816 %(d<2)=0.1122 %(d<4)=0.2143 %(d<8)=0.3469
GDT-HA-score= 0.1122 %(d<0.5)=0.0408 %(d<1)=0.0816 %(d<2)=0.1122 %(d<4)=0.2143

 -------- rotation matrix to rotate Chain-1 to Chain-2 ------
 i          t(i)         u(i,1)         u(i,2)         u(i,3)
 1     20.0708453294   0.0013471591   0.7085462481  -0.7056630921
 2     64.8878774104   0.5102866660  -0.6073607826  -0.6088681288
 3     56.1433817518  -0.8600033161  -0.3592702243  -0.3623799143

Superposition in the TM-score: Length(d<5.0)= 21  RMSD=  2.91
(":" denotes the residue pairs of distance < 5.0 Angstrom)
PFSIETMPVPKKLKVGETAEIRCQLHRDTKYFIRYFQPDGAGTLKMSDGTVLLPNDLYPLPGETFRLYYTSASTDQQTVDVYFQDSFGQLQQLTFSFN
 :::::::::  :   ::::::::::        :                                                               
PFSIETMPVPKKLKVGETAEIRCQLHRDTKYFIRYFQPDGAGTLKMSDGTVLLPNDLYPLPGETFRLYYTSASTDQQTVDVYFQDSFGQLQQLTFSFN
12345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678

