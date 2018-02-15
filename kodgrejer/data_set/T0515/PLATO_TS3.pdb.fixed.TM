
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
RMSD of  the common residues=  165.441

TM-score    = 0.1342  (d0= 6.79)
MaxSub-score= 0.0519  (d0= 3.50)
GDT-TS-score= 0.0747 %(d<1)=0.0345 %(d<2)=0.0489 %(d<4)=0.0776 %(d<8)=0.1379
GDT-HA-score= 0.0460 %(d<0.5)=0.0230 %(d<1)=0.0345 %(d<2)=0.0489 %(d<4)=0.0776

 -------- rotation matrix to rotate Chain-1 to Chain-2 ------
 i          t(i)         u(i,1)         u(i,2)         u(i,3)
 1     59.8175651653   0.9328808890   0.3589323809  -0.0300132124
 2      9.3392485132   0.3600166711  -0.9317494125   0.0472337686
 3     72.7705375059  -0.0110110640  -0.0548687368  -0.9984328611

Superposition in the TM-score: Length(d<5.0)= 27  RMSD=  3.54
(":" denotes the residue pairs of distance < 5.0 Angstrom)
MIETPYYLIDKAKLTRNMERIAHVREKSGAKALLALKCFATWSVFDLMRDYMDGTTSSSLFEVRLGRERFGKETHAYSVAYGDNEIDEVVSHADKIIFNSISQLERFADKAAGIARGLRLNPQRLGEWDVPKVERVMDRINGFMIHNNCENKDFGLFDRMLGEIEERFGALIARVDWVSLGGGIHFTGDDYPVDAFSARLRAFSDRYGVQIYLEPGEASITKSTTLEVTVLDTLYNKNLAIVDSSIEAHMLDLLIYRETAKVLPNEGSHSYMICGKSCLAGDVFGEFRFAEELKVGDRISFQDAAGYTMVKKNWFNGVKMPAIAIRELDGSVRTVREFTYADYEQSLS
                     :  :: :::::                                                                     ::: ::::::::::: : : : ::                                                                                                                                                                                                                               
MIETPYYLIDKAKLTRNMERIAHVREKSGAKALLALKCFATWSVFDLMRDYMDGTTSSSLFEVRLGRERFGKETHAYSVAYGDNEIDEVVSHADKIIFNSISQLERFADKAAGIARGLRLNPQRLGEWDVPKVERVMDRINGFMIHNNCENKDFGLFDRMLGEIEERFGALIARVDWVSLGGGIHFTGDDYPVDAFSARLRAFSDRYGVQIYLEPGEASITKSTTLEVTVLDTLYNKNLAIVDSSIEAHMLDLLIYRETAKVLPNEGSHSYMICGKSCLAGDVFGEFRFAEELKVGDRISFQDAAGYTMVKKNWFNGVKMPAIAIRELDGSVRTVREFTYADYEQSLS
123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678

