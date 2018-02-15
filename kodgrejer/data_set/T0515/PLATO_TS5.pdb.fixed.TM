
 *****************************************************************************
 *                                 TM-SCORE                                  *
 * A scoring function to assess the similarity of protein structures         *
 * Based on statistics:                                                      *
 *       0.0 < TM-score < 0.17, random structural similarity                 *
 *       0.5 < TM-score < 1.00, in about the same fold                       *
 * Reference: Yang Zhang and Jeffrey Skolnick, Proteins 2004 57: 702-710     *
 * For comments, please email to: zhng@umich.edu                             *
 *****************************************************************************

Structure1: input/casp  Length=  348
Structure2: input/casp  Length=  348 (by which all scores are normalized)
Number of residues in common=  348
RMSD of  the common residues=  162.856

TM-score    = 0.1634  (d0= 6.79)
MaxSub-score= 0.0604  (d0= 3.50)
GDT-TS-score= 0.0884 %(d<1)=0.0345 %(d<2)=0.0546 %(d<4)=0.0862 %(d<8)=0.1782
GDT-HA-score= 0.0496 %(d<0.5)=0.0230 %(d<1)=0.0345 %(d<2)=0.0546 %(d<4)=0.0862

 -------- rotation matrix to rotate Chain-1 to Chain-2 ------
 i          t(i)         u(i,1)         u(i,2)         u(i,3)
 1     50.1425495900   0.6412343255   0.5514754255   0.5335666734
 2     13.1249053719   0.3338368908   0.4255966873  -0.8410828675
 3     45.6857760363  -0.6909207409   0.7174554445   0.0888043633

Superposition in the TM-score: Length(d<5.0)= 29  RMSD=  3.44
(":" denotes the residue pairs of distance < 5.0 Angstrom)
MIETPYYLIDKAKLTRNMERIAHVREKSGAKALLALKCFATWSVFDLMRDYMDGTTSSSLFEVRLGRERFGKETHAYSVAYGDNEIDEVVSHADKIIFNSISQLERFADKAAGIARGLRLNPQRLGEWDVPKVERVMDRINGFMIHNNCENKDFGLFDRMLGEIEERFGALIARVDWVSLGGGIHFTGDDYPVDAFSARLRAFSDRYGVQIYLEPGEASITKSTTLEVTVLDTLYNKNLAIVDSSIEAHMLDLLIYRETAKVLPNEGSHSYMICGKSCLAGDVFGEFRFAEELKVGDRISFQDAAGYTMVKKNWFNGVKMPAIAIRELDGSVRTVREFTYADYEQSLS
                                ::      ::                                 :                ::: :::::::::::::::::::                         ::                                                                                                                                                                                                              
MIETPYYLIDKAKLTRNMERIAHVREKSGAKALLALKCFATWSVFDLMRDYMDGTTSSSLFEVRLGRERFGKETHAYSVAYGDNEIDEVVSHADKIIFNSISQLERFADKAAGIARGLRLNPQRLGEWDVPKVERVMDRINGFMIHNNCENKDFGLFDRMLGEIEERFGALIARVDWVSLGGGIHFTGDDYPVDAFSARLRAFSDRYGVQIYLEPGEASITKSTTLEVTVLDTLYNKNLAIVDSSIEAHMLDLLIYRETAKVLPNEGSHSYMICGKSCLAGDVFGEFRFAEELKVGDRISFQDAAGYTMVKKNWFNGVKMPAIAIRELDGSVRTVREFTYADYEQSLS
123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678

