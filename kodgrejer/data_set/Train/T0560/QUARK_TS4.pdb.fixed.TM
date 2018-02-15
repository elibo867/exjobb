
 *****************************************************************************
 *                                 TM-SCORE                                  *
 * A scoring function to assess the similarity of protein structures         *
 * Based on statistics:                                                      *
 *       0.0 < TM-score < 0.17, random structural similarity                 *
 *       0.5 < TM-score < 1.00, in about the same fold                       *
 * Reference: Yang Zhang and Jeffrey Skolnick, Proteins 2004 57: 702-710     *
 * For comments, please email to: zhng@umich.edu                             *
 *****************************************************************************

Structure1: input/casp  Length=   64
Structure2: input/casp  Length=   64 (by which all scores are normalized)
Number of residues in common=   64
RMSD of  the common residues=    2.993

TM-score    = 0.6626  (d0= 2.74)
MaxSub-score= 0.6717  (d0= 3.50)
GDT-TS-score= 0.7383 %(d<1)=0.4063 %(d<2)=0.6719 %(d<4)=0.8750 %(d<8)=1.0000
GDT-HA-score= 0.5508 %(d<0.5)=0.2500 %(d<1)=0.4063 %(d<2)=0.6719 %(d<4)=0.8750

 -------- rotation matrix to rotate Chain-1 to Chain-2 ------
 i          t(i)         u(i,1)         u(i,2)         u(i,3)
 1     10.7994831221   0.2856047951   0.2056715956  -0.9360176792
 2    -77.0320945298   0.8660468555   0.3628246281   0.3439783906
 3     36.2730092382   0.4103568508  -0.9088770456  -0.0744967851

Superposition in the TM-score: Length(d<5.0)= 56  RMSD=  2.02
(":" denotes the residue pairs of distance < 5.0 Angstrom)
KKIVGANAGKVWHALNEADGISIPELARKVNLSVESTALAVGWLARENKVVIERKNGLIEIYNE
:::::::::::::::::::::::::::::::::::::::::::::::: :::     ::::  :
KKIVGANAGKVWHALNEADGISIPELARKVNLSVESTALAVGWLARENKVVIERKNGLIEIYNE
1234567890123456789012345678901234567890123456789012345678901234

