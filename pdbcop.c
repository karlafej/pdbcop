#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define TOL_D 0.2;
#define TOL_B 3.5;
/*#define M_PI acos(-1.0)*/


FILE *finput1, *finput2; /* pointer na vstupni soubor */
char *fname1, *fname2; 
FILE *f_shiftoutput, *f_Boutput, *fasta;

#define WARN1 "    \0"
#define WARN2 "**  \0"
#define WARN3 "*** \0"
#define WARN4 "****\0"
#define WARN5 "!!!!\0"



int LONGLIST=0;
int SEQUENCE=0;
int DIFFER=0;
int ANISO=0;
/* int WATER_CONTACTS=0; */

char record[250], record2[250];
long int Aparam, Oparam;        /* Number of records with "alternate" and/or    BYLO LONG INT? nejde %4d
                                   "Occupancy=0.50 atoms */
double Bmin, Bmax, Omin, Omax;  /* Minimal and maximal B value and
                                   Minimal and maximal occupancy respectively*/
double Bav;                     /* Average B value */
int atcount;                    /* Number of atoms */
int atBcount=0;
int atAcount=0;
int watercount=0;
int newres=1;


double told, tolb;

struct inmemtype{               /* Structure for record storage in a memory */
  char *PDBstring;
  struct inmemtype *nextr;
};

struct inmemtype 
  *alternate_first,             /* The first record of "alternate" list */  
  *alternate_last,              /* The last record of "alternate" list */
  *Occup_first,                 /* The first record of "Occup=50.00" list */
  *Occup_last,                  /* The last record of "Occup=50.00" list */
  *Zccup_first,                 /* The first record of "Occup<=0.1" list */
  *Zccup_last,                  /* The last record of "Occup<=0.1" list */
  *locator;                     /* Just for moving along the list */

char Chstring_old[1];  /* Storing of the chain identifier*/
typedef struct bchain {        /* Structure for storage of B values of individual chains */
  char chain[1];
  double Bmax, Bmin, Bav;
  int atcount, rescount;
}BCHAIN;

BCHAIN **chain_array= NULL;
int chain_nu=-1;

typedef struct changeat {        /* Structure for storage of changed atoms */
  char Atiden[27], warning[5];
  double dd,dx,dy,dz,dB,B,O;
}CHANGEAT;

CHANGEAT **shiftat_array= NULL;   /*shifted atoms */
int shiftat_nu=-1;

CHANGEAT **difB_array=NULL;    /* large B changes */
int difB_nu=-1;

CHANGEAT **difAl_array=NULL;   /* both shifts and B changes */
int difAl_nu=-1;

typedef struct atoms {          /* structure for storing input PDB file */
    char *PDBstring;
    char resname[4];
    int Resnum;
    char iCode;
    double B;
}ATOMS;

ATOMS **atomA_array = NULL;


char atomA[10], atomB[10];

typedef struct anisoadp {
  char *ADPstring;
  double U11,U22,U33,U12,U13,U23;
}ANISOADP;

ANISOADP **adpA_array = NULL;
int adpAcount=0;

int ResBprev, ResShiftprev;  /* prev residues - export of PDB files with biggest changes*/

typedef struct bextreme
{
  char Atiden[27];
  double B;
}BEXTREME;

BEXTREME lowB_array[10];
BEXTREME highB_array[10];

typedef struct nonconsh {   /* nonstandard atoms */
char atom[3]; 
char resname[4];
int checked;  
} NONCONSH;

NONCONSH **nonconsh_array=NULL;
int nonstd_nu=-1;

typedef struct nonstandard {    /*nonstandard residues */
char chain; 
int  Atomcount, Resnum, checked;  
char resname[4]; 
} NONSTD;

NONSTD **nonstd_array=NULL;
int nonconsh_nu=-1; 
/*char resname[3]; */             /* proc [3]? kde se to pouziva? */
int Atomcount=0;    /*number of atoms in a nonstandard residue*/

typedef struct seq {
  char chain;
  char *seqstr;
  int chain_size;
} SEQ;

SEQ **sequence_array=NULL;
int sequence_nu=-1;

/* eigenvalues 
Joachim Kopp
Numerical diagonalization of hermitian 3x3 matrices
arXiv.org preprint: physics/0610206
Int. J. Mod. Phys. C19 (2008) 523-548
*/

// Constants
#define M_SQRT3    1.73205080756887729352744634151   // sqrt(3)

// Macros
#define SQR(x)      ((x)*(x))                        // x^2 



void printhelp(char *comm)            /*Prints a help screen... */
{
  fprintf(stderr,"\n Use: %s [-l] [-d]  file_final (file2) [shift_tol[]] [delta_B_tol[]] \n\n",comm);
  fprintf(stderr," Options: -l  ... lists all atoms with occupancy different\n");
  fprintf(stderr,"                  from 0.00 and 1.00\n");
  fprintf(stderr,"          -d  ... compares atom positions and B-factors in two PDB files \n");
  fprintf(stderr,"                  shift_tol: minimal shift to be printed\n");
  fprintf(stderr,"                  delta_B: minimal difference in B to be printed\n");
  fprintf(stderr,"                  (Alternate locations, occupancies and B-factors are printed\n");
  fprintf(stderr,"                  only for the first file.)\n");
  fprintf(stderr,"          -s  ... prints out sequence\n\n");
  
  /*printf("           -w  ... zatim nic - blizke kontakty vod\n");*/
};

void initparam(void)            /* Intializes the program parameters */
{  
  Aparam=Oparam=0;
  Bmin=100.00;
  Bmax=0.00;
  Omin=1.00;
  Omax=0.00;
  Bav=0.00;
  atcount=0;
  alternate_first=alternate_last=Occup_first=Occup_last=Zccup_first=Zccup_last=NULL;
  Chstring_old[0]='\0';
  ResBprev=0;
  ResShiftprev=0;


};

void getrecordA()            /* Reads a record from the file defined */
{
  int i, CONTINUE;

  i=0; 
  CONTINUE=1;                   /* It's just for case when there is no newline
                                   at the end of file */
 
  while ((!feof(finput1))&&(CONTINUE)){
    record [i++]=fgetc(finput1);
    if (record[i-1]=='\n') CONTINUE=0;
  };
  record[i-1]=0;
      
};

void getrecordB()
{
  int i, CONTINUE;

  i=0; 
  CONTINUE=1;                   /* It's just for case when there is no newline
                                   at the end of file */
 
  while ((!feof(finput2))&&(CONTINUE)){
    record2 [i++]=fgetc(finput2);
    if (record2[i-1]=='\n') CONTINUE=0;
  };
  record2[i-1]=0;
  /*locator->PDBstring=(char *) malloc(strlen(record2));*/
};

void getatomsA()   /*stores an array of atom records */
{ 
  int i;
  int Resnum;
  double B, U;
  char Resstring[4], Bstring[6], Ustring[7];
  while (!feof(finput1)) {
    getrecordA();

    if (((record[0]=='A')&&(record[1]=='T')&&(record[2]=='O')&&(record[3]=='M'))||((record[0]=='H')&&(record[1]=='E')&&(record[2]=='T')&&(record[3]=='A'))) {
       if((atomA_array = (ATOMS **)realloc(atomA_array, (atAcount + 1) * sizeof(ATOMS *)))==NULL) {
        fprintf(stderr,"Could not allocate memory\n" );
        return;
       }
       if((atomA_array[atAcount]=(ATOMS *)malloc(sizeof(ATOMS)))==NULL) {
        fprintf(stderr,"Could not allocate memory\n" );
        return;
       }
       if((atomA_array[atAcount]->PDBstring=(char *) malloc(strlen(record)))==NULL) {
        fprintf(stderr,"Could not allocate memory\n" );
        return;
       }
      for (i=0;i<4;++i) Resstring[i]=record[i+22];                         
      Resstring[i]='\0';
      Resnum=atof(Resstring); 
      atomA_array[atAcount]->Resnum=Resnum;
      for (i=0;i<6;++i) Bstring[i]=record[i+60];  /* The B factor value is read... */
      Bstring[i]='\0';
      B=atof(Bstring); 
      atomA_array[atAcount]->B=B;
      for (i=0;i<3;++i) (atomA_array[atAcount]->resname[i])=record[i+17];                         
      atomA_array[atAcount]->resname[i]='\0';
      atomA_array[atAcount]->iCode=record[26];
      atomA_array[atAcount++]->PDBstring = strdup(record);
  }
  else if ((record[0]=='A')&&(record[1]=='N')&&(record[2]=='I')&&(record[3]=='S')&&(record[4]=='O')&&(record[5]=='U')) {  /* anisotropic ANISOADP */
    
      if((adpA_array = (ANISOADP **)realloc(adpA_array, (adpAcount + 1) * sizeof(ANISOADP *)))==NULL) {
        fprintf(stderr,"Could not allocate memory\n" );
        return;
      }
        if((adpA_array[adpAcount]=(ANISOADP *)malloc(sizeof(ANISOADP)))==NULL) {
        fprintf(stderr,"Could not allocate memory\n" );
        return;
       }
     if((adpA_array[adpAcount]->ADPstring=(char *) malloc(strlen(record)))==NULL) {
        fprintf(stderr,"Could not allocate memory\n" );
        return;
       }
      for (i=0;i<7;++i) Ustring[i]=record[i+28];  
      Ustring[i]='\0';
      U=atof(Ustring); 
      U=U/10000;
      adpA_array[adpAcount]->U11=U;
      for (i=0;i<7;++i) Ustring[i]=record[i+35];  
      U=atof(Ustring); 
      U=U/10000;
      adpA_array[adpAcount]->U22=U;
      for (i=0;i<7;++i) Ustring[i]=record[i+42];  
      U=atof(Ustring); 
      U=U/10000;
      adpA_array[adpAcount]->U33=U;
      for (i=0;i<7;++i) Ustring[i]=record[i+49];  
      U=atof(Ustring); 
      U=U/10000;
      adpA_array[adpAcount]->U12=U;
      for (i=0;i<7;++i) Ustring[i]=record[i+56]; 
      U=atof(Ustring); 
      U=U/10000;
      adpA_array[adpAcount]->U13=U;
      for (i=0;i<7;++i) Ustring[i]=record[i+63];  
      U=atof(Ustring); 
      U=U/10000;
      adpA_array[adpAcount]->U23=U;
      
      adpA_array[adpAcount++]->ADPstring = strdup(record);
      ANISO=1; 

    };

  };
};


void analyserecord(void)         /* Self-explanatory (I think) */
{  
  int i,k,j,l;
  char Bstring[6], Ostring[6];   /* To read the B and occupancy values 
                                    from record */
  double B, Bcor, O, Ocor;       /* To operate with B and occupancy values
                                    as with a number */
  int Resnumprev;

  B=0, O=0;

  char atomname[2];
  int Resnumtmp=0; 
  char iCodetmp='\0';                  /* for work wint nonstandard residues and atoms */
  char resnametmp[4]={'\0'};
  char chaintmp[2];

    for (i=0;i<10;i++) {
      lowB_array[i].B=500;
      highB_array[i].B=-500;
    };

  for (k=0; k<atAcount; k++) {
    
    if (atomA_array[k]->PDBstring[16]!=' ') {       /* The atom has alternate location! */
      Aparam++;

      locator=(struct inmemtype *) malloc(sizeof(struct inmemtype));
      locator->PDBstring=atomA_array[k]->PDBstring;
      locator->nextr=NULL;
                                 /* The record is taken to the memory */

      if (alternate_first==NULL) { 
                                 /*Is this a first "alternate" atom? */ 
        alternate_first=locator;
        alternate_last=locator;
      }else alternate_last->nextr=locator;
                                 /* If no, then put it at the end 
                                    of the list... */
      alternate_last=locator;
    };

    if (Chstring_old[0]!=(atomA_array[k]->PDBstring[21])) {  /* new chain? initialization of parameters */
      chain_nu++;
      chain_array = (BCHAIN **)realloc(chain_array, (chain_nu+1)*sizeof(BCHAIN *));
      chain_array[chain_nu]=(BCHAIN *)malloc(sizeof(BCHAIN));
      chain_array[chain_nu]->chain[0]=atomA_array[k]->PDBstring[21];
      chain_array[chain_nu]->Bmin=100.00;
      chain_array[chain_nu]->Bmax=0.00;
      chain_array[chain_nu]->Bav=0.00;
      chain_array[chain_nu]->atcount=0;
      chain_array[chain_nu]->rescount=0;
      Chstring_old[0]=atomA_array[k]->PDBstring[21];
      Resnumprev=0;
    };
    
    B=atomA_array[k]->B;             /* ...converted from a string 
                                    into a number... */
    
    chain_array[chain_nu]->atcount++;
    chain_array[chain_nu]->Bav=chain_array[chain_nu]->Bav+B;
    if(chain_array[chain_nu]->Bmin>B) chain_array[chain_nu]->Bmin=B;
    if(chain_array[chain_nu]->Bmax<B) chain_array[chain_nu]->Bmax=B;
    if((atomA_array[k]->Resnum)!=Resnumprev) {
      chain_array[chain_nu]->rescount++;
      Resnumprev=atomA_array[k]->Resnum;
    }
    
    Bav=Bav+B;                   /* Sum of B values -- counting of average B */
    atcount++;
    
    if (B<Bmin) Bmin=B;          /* Is it a new B minimum? */
    if (B>Bmax) Bmax=B;          /* Or is it a new B maximum? */


     for (i=0;i<10;i++) {
        if(B < lowB_array[i].B) {
            for (j=9; j>i; j--) {
              lowB_array[j].B=lowB_array[j-1].B;
              strcpy(lowB_array[j].Atiden,lowB_array[j-1].Atiden);
            }
            lowB_array[i].B=atomA_array[k]->B;
            strncpy(lowB_array[i].Atiden,atomA_array[k]->PDBstring,26);
            lowB_array[i].Atiden[26]='\0';
            break;
        };
      };

      for (i=0;i<10;i++) {
        if(B>highB_array[i].B) {
          for (j=9; j>i; j--) {
              highB_array[j].B=highB_array[j-1].B;
              strcpy(highB_array[j].Atiden,highB_array[j-1].Atiden);
            }
          highB_array[i].B=B;
          strncpy(highB_array[i].Atiden,atomA_array[k]->PDBstring,26);
          highB_array[i].Atiden[26]='\0';
          break;
        };
      };

    for (i=0;i<6;++i) Ostring[i]=atomA_array[k]->PDBstring[i+54];
                                 /* The occupancy value is read... */
    Ostring[i]='\0';
    O=atof(Ostring);             /* ...converted from a string 
                                    into a number... */
    if (O<Omin) Omin=O;          /* Is it a new occupancy minimum? */
    if (O>Omax) Omax=O;          /* Or is it a new occupancy maximum? */


    if ((O!=0.0)&&(O!=1.0)){          /* When a occupancy value differs from
                                         0.00 and 1.00... */
      Oparam++;
      locator=(struct inmemtype *) malloc(sizeof(struct inmemtype));
      locator->PDBstring=atomA_array[k]->PDBstring;
                                 /* ...put the record to the memory 
                                    (in the same way as with "alternate" 
                                    records ;-) */
      locator->nextr=NULL;
      if (Occup_first==NULL) {
        Occup_first=locator;
        Occup_last=locator;
      }else Occup_last->nextr=locator;
      Occup_last=locator;
    };

    if ((O<=0.1)){          /* When a occupancy value is .le. 0.1 */                                   
      
      locator=(struct inmemtype *) malloc(sizeof(struct inmemtype));
      locator->PDBstring=atomA_array[k]->PDBstring;
                                 /* ...put the record to the memory 
                                    (in the same way as with "alternate" 
                                    records ;-) */
      locator->nextr=NULL;
      if (Zccup_first==NULL) {
        Zccup_first=locator;
        Zccup_last=locator;
      }else Zccup_last->nextr=locator;
      Zccup_last=locator;
    };

    if ((atomA_array[k]->PDBstring[13]=='O')&&(atomA_array[k]->PDBstring[17]=='H')&&(atomA_array[k]->PDBstring[18]=='O')&&(atomA_array[k]->PDBstring[19]=='H')) {
      watercount++;   /*number of water molecules */
    }

  if((((atomA_array[k]->Resnum)!=Resnumtmp)||((atomA_array[k]->iCode)!=iCodetmp))||(strcmp((atomA_array[k]->resname),resnametmp)!=0)) { /*new residue */
        Atomcount=0;
        strncpy(resnametmp,(atomA_array[k]->resname),4);
        iCodetmp=(atomA_array[k]->iCode);
        Resnumtmp=(atomA_array[k]->Resnum);
        chaintmp[0]=(atomA_array[k]->PDBstring[21]);
        chaintmp[1]='\0';
       
  if ((strcmp(resnametmp,"GLY")!=0)&&    /*is it a non-standard residue?*/
      (strcmp(resnametmp,"ALA")!=0)&&
      (strcmp(resnametmp,"VAL")!=0)&&
      (strcmp(resnametmp,"LEU")!=0)&&
      (strcmp(resnametmp,"ILE")!=0)&&
      (strcmp(resnametmp,"SER")!=0)&&
      (strcmp(resnametmp,"TYR")!=0)&&
      (strcmp(resnametmp,"THR")!=0)&&
      (strcmp(resnametmp,"PHE")!=0)&&
      (strcmp(resnametmp,"TRP")!=0)&&
      (strcmp(resnametmp,"HIS")!=0)&&
      (strcmp(resnametmp,"LYS")!=0)&&
      (strcmp(resnametmp,"ARG")!=0)&&
      (strcmp(resnametmp,"GLN")!=0)&&
      (strcmp(resnametmp,"GLU")!=0)&&
      (strcmp(resnametmp,"ASN")!=0)&&
      (strcmp(resnametmp,"ASP")!=0)&&
      (strcmp(resnametmp,"PRO")!=0)&&
      (strcmp(resnametmp,"CYS")!=0)&&
      (strcmp(resnametmp,"MET")!=0)&&
      (strcmp(resnametmp,"HOH")!=0)) {
   
    nonstd_nu++;
    if((nonstd_array = (NONSTD **)realloc(nonstd_array, (nonstd_nu+1)*sizeof(NONSTD *)))==NULL) {
       fprintf(stderr,"Could not allocate memory\n" );
        return;
    }
    if((nonstd_array[nonstd_nu]=(NONSTD *)malloc(sizeof(NONSTD)))==NULL) {
      fprintf(stderr,"Could not allocate memory\n" );
      return;
    }

    strncpy((nonstd_array[nonstd_nu]->resname),resnametmp,strlen(resnametmp)+1);   /*store a non-standard residue*/
    nonstd_array[nonstd_nu]->resname[3]='\0';
    nonstd_array[nonstd_nu]->Resnum=atomA_array[k]->Resnum;
    nonstd_array[nonstd_nu]->chain=(atomA_array[k]->PDBstring[21]);
    nonstd_array[nonstd_nu]->checked=0;
   
   
    for(l=k;l<atAcount;l++) {
    if((Resnumtmp==(atomA_array[l]->Resnum))&&(iCodetmp==(atomA_array[l]->iCode))&&(strncmp(resnametmp,(atomA_array[l]->resname),3)==0)&&(atomA_array[l]->PDBstring[21]==chaintmp[0])) {
      for (i=0;i<2;i++) atomname[i]=(atomA_array[l]->PDBstring[i+76]);                         
        atomname[i]='\0';
      
        if (strcmp(atomname," H")!=0) {     /*counting of non-H atoms in a nonstandard residue*/
          Atomcount++;
        }
        if ((strcmp(atomname," C")!=0)&&     /*is it a nonstandard atom?*/
            (strcmp(atomname," H")!=0)&&
            (strcmp(atomname," N")!=0)&&
            (strcmp(atomname," O")!=0)&&
            (strcmp(atomname," S")!=0)) {

           nonconsh_nu++;
            if((nonconsh_array = (NONCONSH **)realloc(nonconsh_array, (nonconsh_nu+1)*sizeof(NONCONSH *)))==NULL) {
              fprintf(stderr,"Could not allocate memory\n" );
              return;
            }
            if((nonconsh_array[nonconsh_nu]=(NONCONSH *)malloc(sizeof(NONCONSH)))==NULL) {
              fprintf(stderr,"Could not allocate memory\n" );
               return;
            }

            strncpy((nonconsh_array[nonconsh_nu]->resname),resnametmp,strlen(resnametmp)+1);   /*store nonstandard atoms*/
            strncpy(nonconsh_array[nonconsh_nu]->atom,atomname,strlen(atomname)+1);
            nonconsh_array[nonconsh_nu]->atom[2]='\0';
            nonconsh_array[nonconsh_nu]->resname[3]='\0';
            nonconsh_array[nonconsh_nu]->checked=0;
        }
    }
    nonstd_array[nonstd_nu]->Atomcount=Atomcount; 
 
    }
  }
  }

  };
};

void analyseU(void)
{

int print=1;
double e1=100;
double e2=100;
double e3=100;

double m, c1, c0;
double U11,U22,U33,U12,U23,U13;
int k;

for (k=0; k<adpAcount; k++) {
  U11=adpA_array[k]->U11;
  U22=adpA_array[k]->U22;
  U33=adpA_array[k]->U33;
  U12=adpA_array[k]->U12;
  U13=adpA_array[k]->U13;
  U23=adpA_array[k]->U23;
  
  // Determine coefficients of characteristic poynomial. We write
  //       | a   d   f  |  U11  U12  U13
  //  A =  | d*  b   e  |  U12  U22  U23
  //       | f*  e*  c  |  U13  U23  U33
  double de = U12*U23;                                 // d * e
  double dd = SQR(U12);                                         // d^2
  double ee = SQR(U23);                                         // e^2
  double ff = SQR(U13);                                         // f^2
  m  = U11 + U22 + U33;
  c1 = (U11*U22 + U11*U33 + U22*U33)        // a*b + a*c + b*c - d^2 - e^2 - f^2
          - (dd + ee + ff);
  c0 = U33*dd + U11*ee + U22*ff - U11*U22*U33
            - 2.0 * U13*de;                                     // c*d^2 + a*e^2 + b*f^2 - a*b*c - 2*f*d*e)

  double p, sqrt_p, q, c, s, phi;
  p = SQR(m) - 3.0*c1;
  q = m*(p - (3.0/2.0)*c1) - (27.0/2.0)*c0;
  sqrt_p = sqrt(fabs(p));

  phi = 27.0 * ( 0.25*SQR(c1)*(p - c1) + c0*(q + 27.0/4.0*c0));
  phi = (1.0/3.0) * atan2(sqrt(fabs(phi)), q);
  
  c = sqrt_p*cos(phi);
  s = (1.0/M_SQRT3)*sqrt_p*sin(phi);

  e2  = (1.0/3.0)*(m - c);
  e3  = e2 + s;
  e1 = e2 + c;
  e2 -= s;

 if(((e1<=0)||(e2<=0)||(e3<=0))&&print==1) {
  printf("\n Warning: ADP non positive definite!                                                 e1     e2     e3");
  print=0;
 };
 if((e1<=0)||(e2<=0)||(e3<=0)) {
  printf("\n %s %6.3f %6.3f %6.3f",adpA_array[k]->ADPstring, e1, e2, e3);
  };
  };
}; 

void analysechanges(void)  /*analyse changes during the refinement*/
{

  double x1,x2,y1,y2,z1,z2,dx,dy,dz,dd,dB,B2,O1;
  char Coorstring[9], Bstring[7], Ostring[7], Chstring[1];
  int Resnum2;       
  char Resstring[5];

  int i, k, j, differ=1;

  if (((record2[0]=='A')&&(record2[1]=='T')&&(record2[2]=='O')&&(record2[3]=='M'))||((record2[0]=='H')&&(record2[1]=='E')&&(record2[2]=='T')&&(record2[3]=='A'))) {
    
    for (i=0;i<4;++i) Resstring[i]=record2[i+22];               
    Resstring[i]='\0';
    Resnum2=atof(Resstring);     

    for (k=0;k<atAcount;k++) {  
      if ((atomA_array[k]->Resnum)==Resnum2) {

        for (i=0;i<10;++i) atomA[i]=atomA_array[k]->PDBstring[i+12];                         
        atomA[i]='\0';    

        for (i=0;i<10;++i) atomB[i]=record2[i+12];                         
        atomB[i]='\0';

        differ=strncmp(atomB,atomA,10);

        if (differ==0) {

          for (i=0;i<6;++i) Bstring[i]=record2[i+60];  /* The B factor value is read... */
          Bstring[i]='\0';
          B2=atof(Bstring); 

          for (i=0;i<8;++i) Coorstring[i]=atomA_array[k]->PDBstring[i+30];           /*x*/              
          Coorstring[i]='\0';
          x1=atof(Coorstring); 
 
          for (i=0;i<8;++i) Coorstring[i]=record2[i+30];                            
          Coorstring[i]='\0';
          x2=atof(Coorstring);

          for (i=0;i<8;++i) Coorstring[i]=atomA_array[k]->PDBstring[i+38];   /*y*/
          Coorstring[i]='\0';
          y1=atof(Coorstring); 

          for (i=0;i<8;++i) Coorstring[i]=record2[i+38];
          Coorstring[i]='\0';
          y2=atof(Coorstring); 

          for (i=0;i<8;++i) Coorstring[i]=atomA_array[k]->PDBstring[i+46];    /*z*/
          Coorstring[i]='\0';
          z1=atof(Coorstring); 

          for (i=0;i<8;++i) Coorstring[i]=record2[i+46];
          Coorstring[i]='\0';
          z2=atof(Coorstring);         
    
          for (i=0;i<6;++i) Ostring[i]=atomA_array[k]->PDBstring[i+54];
                                 /* The occupancy value is read... */
          Ostring[i]='\0';
          O1=atof(Ostring);             /* ...converted from a string 
                                    into a number... */
          dx=x1-x2;
          dy=y1-y2;
          dz=z1-z2;
          dd=sqrt(dx*dx+dy*dy+dz*dz);

          dB=(atomA_array[k]->B)-B2;

        break;
     }; /*differ*/
     } /*resnum1=resnum2*/
  } /* for pres pole */  

  
  if (dd>told) {                   /*is the shift larger than the set cutoff?*/
    shiftat_nu++;
    if((shiftat_array = (CHANGEAT **)realloc(shiftat_array, (shiftat_nu+1)*sizeof(CHANGEAT *)))==NULL) {
      fprintf(stderr,"Could not allocate memory\n" );
      return;
    }
    if((shiftat_array[shiftat_nu]=(CHANGEAT *)malloc(sizeof(CHANGEAT)))==NULL) {
      fprintf(stderr,"Could not allocate memory\n" );
      return;
    }
    strncpy(shiftat_array[shiftat_nu]->Atiden,atomA_array[k]->PDBstring,26);     /*store shifted atoms for futher use*/
    shiftat_array[shiftat_nu]->Atiden[26]='\0';
    shiftat_array[shiftat_nu]->dd=dd;
    shiftat_array[shiftat_nu]->dx=dx;
    shiftat_array[shiftat_nu]->dy=dy;
    shiftat_array[shiftat_nu]->dz=dz;
    shiftat_array[shiftat_nu]->dB=dB;
    shiftat_array[shiftat_nu]->B=atomA_array[k]->B;
    shiftat_array[shiftat_nu]->O=O1;
    Chstring[0]=atomA_array[k]->PDBstring[21];
    if (dd>(5*told)) {
     strcpy(shiftat_array[shiftat_nu]->warning,WARN5);
    }
    else if (dd>(4*told)) {
      strcpy(shiftat_array[shiftat_nu]->warning,WARN4);
    }
    else if (dd>(3*told)) {
      strcpy(shiftat_array[shiftat_nu]->warning,WARN3);
    }
      else if (dd>(2*told)) {
      strcpy(shiftat_array[shiftat_nu]->warning,WARN2);
    }
    else {
      strcpy(shiftat_array[shiftat_nu]->warning,WARN1);
    }
    if(ResShiftprev!=Resnum2) {  /* printing a PDB containing shifted atoms */
     for (j=0;j<atAcount;j++) {  
        if((Resnum2==(atomA_array[j]->Resnum))&&(Chstring[0]==(atomA_array[j]->PDBstring[21]))) {
          fprintf(f_shiftoutput ,"%s \n",(atomA_array[j]->PDBstring));
        }
      }
      ResShiftprev=Resnum2; 
    }
  }

  if ((fabs(dB))>tolb) {    /*  B change larger than the set cutoff */
    difB_nu++;
    if((difB_array = (CHANGEAT**)realloc(difB_array, (difB_nu+1)*sizeof(CHANGEAT*)))==NULL) {
      fprintf(stderr,"Could not allocate memory\n" );
      return;
    }
    if((difB_array[difB_nu]=(CHANGEAT*)malloc(sizeof(CHANGEAT)))==NULL){
      fprintf(stderr,"Could not allocate memory\n");
      return;
    }
    strncpy(difB_array[difB_nu]->Atiden,atomA_array[k]->PDBstring,26);   /*store B changes*/
    difB_array[difB_nu]->Atiden[26]='\0';
    difB_array[difB_nu]->dB=dB;
    difB_array[difB_nu]->dd=dd;
    difB_array[difB_nu]->dx=dx;
    difB_array[difB_nu]->dy=dy;
    difB_array[difB_nu]->dz=dz;
    difB_array[difB_nu]->B=atomA_array[k]->B;
    difB_array[difB_nu]->O=O1;
    Chstring[0]=atomA_array[k]->PDBstring[21];
    if (fabs(dB)>(5*tolb)) {
     strcpy(difB_array[difB_nu]->warning,WARN5);
    }
    else if ((fabs(dB))>(4*tolb)) {
      strcpy(difB_array[difB_nu]->warning,WARN4);
    }
    else if ((fabs(dB))>(3*tolb)) {
      strcpy(difB_array[difB_nu]->warning,WARN3);
    }
      else if ((fabs(dB))>(2*tolb)) {
      strcpy(difB_array[difB_nu]->warning,WARN2);
    }
    else {
      strcpy(difB_array[difB_nu]->warning,WARN1);
    }
    if(ResBprev!=Resnum2) {   /* printing a PDB containing atoms with dB>tolb */
     for (j=0;j<atAcount;j++) {  
        if((Resnum2==(atomA_array[j]->Resnum))&&(Chstring[0]==(atomA_array[j]->PDBstring[21]))) {
          fprintf(f_Boutput ,"%s \n",(atomA_array[j]->PDBstring));
        }
      }
      ResBprev=Resnum2; 
    }
  }
    if ((dd>told)&&((fabs(dB))>tolb)) {   /* both shift & B */
    difAl_nu++;
    if((difAl_array = (CHANGEAT **)realloc(difAl_array, (difAl_nu+1)*sizeof(CHANGEAT *)))==NULL) {
      fprintf(stderr,"Could not allocate memory\n");
      return;
    }
    if((difAl_array[difAl_nu]=(CHANGEAT *)malloc(sizeof(CHANGEAT)))==NULL) {
      fprintf(stderr,"Could not allocate memory\n");
      return;
    }
    strncpy(difAl_array[difAl_nu]->Atiden,atomA_array[k]->PDBstring,26);
    difAl_array[difAl_nu]->Atiden[26]='\0';
    difAl_array[difAl_nu]->dd=dd;
    difAl_array[difAl_nu]->dx=dx;
    difAl_array[difAl_nu]->dy=dy;
    difAl_array[difAl_nu]->dz=dz;
    difAl_array[difAl_nu]->dB=dB;
    difAl_array[difAl_nu]->B=atomA_array[k]->B;
    difAl_array[difAl_nu]->O=O1;
    strcpy(difAl_array[difAl_nu]->warning,WARN1);
    } 

  }  /*tady konci ATOM if nebo for k*/

}; /* tady konci funkce */

void sequenceanalysis(void)

{
  char resname[4];
  char iCodetmp='\n';
  int k,i;
  int Resnumtmp=0;
  char Chaintmp='\0';

  for(k=0;k<atAcount;k++) {
    if(atomA_array[k]->PDBstring[0]=='A') {
    if ((Resnumtmp!=atomA_array[k]->Resnum)||(iCodetmp!=atomA_array[k]->iCode)) {
      strncpy(resname,(atomA_array[k]->resname),4);
      Resnumtmp=(atomA_array[k]->Resnum);

      iCodetmp=(atomA_array[k]->iCode);

      if(Chaintmp!=atomA_array[k]->PDBstring[21]) {  /*new chain*/
        i=1;
        sequence_nu++;
        if((sequence_array = (SEQ **)realloc(sequence_array, (sequence_nu+1)*sizeof(SEQ *)))==NULL) {
          fprintf(stderr,"Could not allocate memory\n" );
          return;
        }
        if((sequence_array[sequence_nu]=(SEQ *)malloc(sizeof(SEQ)))==NULL) {
          fprintf(stderr,"Could not allocate memory\n" );
          return;
        }
        if((sequence_array[sequence_nu]->seqstr=(char*)malloc(sizeof(char)))==NULL) {
          fprintf(stderr,"Could not allocate memory\n" );
          sequence_array[sequence_nu]->seqstr[0]='\0';
          return;
        }
    
        (sequence_array[sequence_nu]->chain)=(atomA_array[k]->PDBstring[21]);
        Chaintmp=(atomA_array[k]->PDBstring[21]);
      } /*tady konci if pres chain */
      

      /*convert to one-letter code */
      if (!strncmp(resname, "GLY",3)) {
        if(((sequence_array[sequence_nu]->seqstr) = (char *)realloc((sequence_array[sequence_nu]->seqstr), (i+1)*sizeof(char)))==NULL) {
          fprintf(stderr,"Could not allocate memory\n" );
          return;
        }

        sequence_array[sequence_nu]->chain_size=i;
        sequence_array[sequence_nu]->seqstr[i-1]='G';
        sequence_array[sequence_nu]->seqstr[i]='\0';
        i++;
      }
      else if (!strncmp(resname, "ALA",3)) {
        if(((sequence_array[sequence_nu]->seqstr) = (char *)realloc((sequence_array[sequence_nu]->seqstr), (i+1)*sizeof(char)))==NULL) {
          fprintf(stderr,"Could not allocate memory\n" );
          return;
        }
        
        sequence_array[sequence_nu]->chain_size=i;
        sequence_array[sequence_nu]->seqstr[i-1]='A';
        sequence_array[sequence_nu]->seqstr[i]='\0';
        i++;
      }
      else if (!strncmp(resname, "VAL",3)) {
        if(((sequence_array[sequence_nu]->seqstr) = (char *)realloc((sequence_array[sequence_nu]->seqstr), (i+1)*sizeof(char)))==NULL) {
          fprintf(stderr,"Could not allocate memory\n" );
          return;
        }
       
        sequence_array[sequence_nu]->chain_size=i;
        sequence_array[sequence_nu]->seqstr[i-1]='V';
        sequence_array[sequence_nu]->seqstr[i]='\0';
        i++;
      }
      else if (!strncmp(resname, "LEU",3)) {
        if(((sequence_array[sequence_nu]->seqstr) = (char *)realloc((sequence_array[sequence_nu]->seqstr), (i+1)*sizeof(char)))==NULL) {
          fprintf(stderr,"Could not allocate memory\n" );
          return;
        }
        
        sequence_array[sequence_nu]->chain_size=i;
        sequence_array[sequence_nu]->seqstr[i-1]='L';
        sequence_array[sequence_nu]->seqstr[i]='\0';
        i++;
      }
      else if (!strncmp(resname, "ILE",3)) {
        if(((sequence_array[sequence_nu]->seqstr) = (char *)realloc((sequence_array[sequence_nu]->seqstr), (i+1)*sizeof(char)))==NULL) {
          fprintf(stderr,"Could not allocate memory\n" );
          return;
        }
        
        sequence_array[sequence_nu]->chain_size=i;
        sequence_array[sequence_nu]->seqstr[i-1]='I';
        sequence_array[sequence_nu]->seqstr[i]='\0';
        i++;
      }
      else if (!strncmp(resname, "SER",3)) {
        if(((sequence_array[sequence_nu]->seqstr) = (char *)realloc((sequence_array[sequence_nu]->seqstr), (i+1)*sizeof(char)))==NULL) {
          fprintf(stderr,"Could not allocate memory\n" );
          return;
        }
        
        sequence_array[sequence_nu]->chain_size=i;
        sequence_array[sequence_nu]->seqstr[i-1]='S';
        sequence_array[sequence_nu]->seqstr[i]='\0';
        i++;
      }
      else if (!strncmp(resname, "THR",3)) {
        if(((sequence_array[sequence_nu]->seqstr) = (char *)realloc((sequence_array[sequence_nu]->seqstr), (i+1)*sizeof(char)))==NULL) {
          fprintf(stderr,"Could not allocate memory\n" );
          return;
        }
        
        sequence_array[sequence_nu]->chain_size=i;
        sequence_array[sequence_nu]->seqstr[i-1]='T';
        sequence_array[sequence_nu]->seqstr[i]='\0';
        i++;
      }
      else if (!strncmp(resname, "PHE",3)) {
        if(((sequence_array[sequence_nu]->seqstr) = (char *)realloc((sequence_array[sequence_nu]->seqstr), (i+1)*sizeof(char)))==NULL) {
          fprintf(stderr,"Could not allocate memory\n" );
          return;
        }
       
        sequence_array[sequence_nu]->chain_size=i;
        sequence_array[sequence_nu]->seqstr[i-1]='F';
        sequence_array[sequence_nu]->seqstr[i]='\0';
        i++;
      }
      else if (!strncmp(resname, "TYR",3)) {
        if(((sequence_array[sequence_nu]->seqstr) = (char *)realloc((sequence_array[sequence_nu]->seqstr), (i+1)*sizeof(char)))==NULL) {
          fprintf(stderr,"Could not allocate memory\n" );
          return;
        }
       
        sequence_array[sequence_nu]->chain_size=i;
        sequence_array[sequence_nu]->seqstr[i-1]='Y';
        sequence_array[sequence_nu]->seqstr[i]='\0';
        i++;
      }
      else if (!strncmp(resname, "TRP",3)) {
        if(((sequence_array[sequence_nu]->seqstr) = (char *)realloc((sequence_array[sequence_nu]->seqstr), (i+1)*sizeof(char)))==NULL) {
          fprintf(stderr,"Could not allocate memory\n" );
          return;
        }
        
        sequence_array[sequence_nu]->chain_size=i;
        sequence_array[sequence_nu]->seqstr[i-1]='W';
        sequence_array[sequence_nu]->seqstr[i]='\0';
        i++;
      }
      else if (!strncmp(resname, "CYS",3)) { 
        if(((sequence_array[sequence_nu]->seqstr) = (char *)realloc((sequence_array[sequence_nu]->seqstr), (i+1)*sizeof(char)))==NULL) {
          fprintf(stderr,"Could not allocate memory\n" );
          return;
        }
       
        sequence_array[sequence_nu]->chain_size=i;
        sequence_array[sequence_nu]->seqstr[i-1]='C';
        sequence_array[sequence_nu]->seqstr[i]='\0';
        i++;
      }
      else if (!strncmp(resname, "MET",3)) {
        if(((sequence_array[sequence_nu]->seqstr) = (char *)realloc((sequence_array[sequence_nu]->seqstr), (i+1)*sizeof(char)))==NULL) {
          fprintf(stderr,"Could not allocate memory\n" );
          return;
        }
        
        sequence_array[sequence_nu]->chain_size=i;
        sequence_array[sequence_nu]->seqstr[i-1]='M';
        sequence_array[sequence_nu]->seqstr[i]='\0';
        i++;
      }
      else if (!strncmp(resname, "PRO",3)) {
        if(((sequence_array[sequence_nu]->seqstr) = (char *)realloc((sequence_array[sequence_nu]->seqstr), (i+1)*sizeof(char)))==NULL) {
          fprintf(stderr,"Could not allocate memory\n" );
          return;
        }
      
        sequence_array[sequence_nu]->chain_size=i;
        sequence_array[sequence_nu]->seqstr[i-1]='P';
        sequence_array[sequence_nu]->seqstr[i]='\0';
        i++;
      }
      else if (!strncmp(resname, "HIS",3)) {
        if(((sequence_array[sequence_nu]->seqstr) = (char *)realloc((sequence_array[sequence_nu]->seqstr), (i+1)*sizeof(char)))==NULL) {
          fprintf(stderr,"Could not allocate memory\n" );
          return;
        }
       
        sequence_array[sequence_nu]->chain_size=i;
        sequence_array[sequence_nu]->seqstr[i-1]='H';
        sequence_array[sequence_nu]->seqstr[i]='\0';
        i++;
      }
      else if (!strncmp(resname, "LYS",3)) {
        if(((sequence_array[sequence_nu]->seqstr) = (char *)realloc((sequence_array[sequence_nu]->seqstr), (i+1)*sizeof(char)))==NULL) {
          fprintf(stderr,"Could not allocate memory\n" );
          return;
        }
      
        sequence_array[sequence_nu]->chain_size=i;
        sequence_array[sequence_nu]->seqstr[i-1]='K';
        sequence_array[sequence_nu]->seqstr[i]='\0';
        i++;
      }
      else if (!strncmp(resname, "ARG",3)) {
        if(((sequence_array[sequence_nu]->seqstr) = (char *)realloc((sequence_array[sequence_nu]->seqstr), (i+1)*sizeof(char)))==NULL) {
          fprintf(stderr,"Could not allocate memory\n" );
          return;
        }
       
        sequence_array[sequence_nu]->chain_size=i;
        sequence_array[sequence_nu]->seqstr[i-1]='R';
        sequence_array[sequence_nu]->seqstr[i]='\0';
        i++;
      }
      else if (!strncmp(resname, "GLU",3)) {
        if(((sequence_array[sequence_nu]->seqstr) = (char *)realloc((sequence_array[sequence_nu]->seqstr), (i+1)*sizeof(char)))==NULL) {
          fprintf(stderr,"Could not allocate memory\n" );
          return;
        }
        
        sequence_array[sequence_nu]->chain_size=i;
        sequence_array[sequence_nu]->seqstr[i-1]='E';
        sequence_array[sequence_nu]->seqstr[i]='\0';
        i++;
      }
      else if (!strncmp(resname, "ASP",3)) {
        if(((sequence_array[sequence_nu]->seqstr) = (char *)realloc((sequence_array[sequence_nu]->seqstr), (i+1)*sizeof(char)))==NULL) {
          fprintf(stderr,"Could not allocate memory\n" );
          return;
        }
       
        sequence_array[sequence_nu]->chain_size=i;
        sequence_array[sequence_nu]->seqstr[i-1]='D';
        sequence_array[sequence_nu]->seqstr[i]='\0';
        i++;
      }
      else if (!strncmp(resname, "GLN",3)) {
        if(((sequence_array[sequence_nu]->seqstr) = (char *)realloc((sequence_array[sequence_nu]->seqstr), (i+1)*sizeof(char)))==NULL) {
          fprintf(stderr,"Could not allocate memory\n" );
          return;
        }
       
        sequence_array[sequence_nu]->chain_size=i;
        sequence_array[sequence_nu]->seqstr[i-1]='Q';
        sequence_array[sequence_nu]->seqstr[i]='\0';
        i++;
      }
      else if (!strncmp(resname, "ASN",3)) {
        if(((sequence_array[sequence_nu]->seqstr) = (char *)realloc((sequence_array[sequence_nu]->seqstr), (i+1)*sizeof(char)))==NULL) {
          fprintf(stderr,"Could not allocate memory\n" );
          return;
        }
       
        sequence_array[sequence_nu]->chain_size=i;
        sequence_array[sequence_nu]->seqstr[i-1]='N';
        sequence_array[sequence_nu]->seqstr[i]='\0';
        i++;
      }
      else if (!strncmp(resname, "PYL",3)) {
        if(((sequence_array[sequence_nu]->seqstr) = (char *)realloc((sequence_array[sequence_nu]->seqstr), (i+1)*sizeof(char)))==NULL) {
          fprintf(stderr,"Could not allocate memory\n" );
          return;
        }
       
        sequence_array[sequence_nu]->chain_size=i;
        sequence_array[sequence_nu]->seqstr[i-1]='O';
        sequence_array[sequence_nu]->seqstr[i]='\0';
        i++;
      }
      else if (!strncmp(resname, "SEC",3)) {
        if(((sequence_array[sequence_nu]->seqstr) = (char *)realloc((sequence_array[sequence_nu]->seqstr), (i+1)*sizeof(char)))==NULL) {
          fprintf(stderr,"Could not allocate memory\n" );
          return;
        }
       
        sequence_array[sequence_nu]->chain_size=i;
        sequence_array[sequence_nu]->seqstr[i-1]='U';
        sequence_array[sequence_nu]->seqstr[i]='\0';
        i++;
      }/*else if... other residues? SEM HYP?*/
      else if (!strncmp(resname, "HOH",3)) {
        sequence_array[sequence_nu]->seqstr[i-1]='\0';
        sequence_array[sequence_nu]->seqstr[i]='\0';
        sequence_array[sequence_nu]->chain_size=i;
      }                             
      else {
         if(((sequence_array[sequence_nu]->seqstr) = (char *)realloc((sequence_array[sequence_nu]->seqstr), (i+1)*sizeof(char)))==NULL) {
          fprintf(stderr,"Could not allocate memory\n" );
          return;
        }
       
        sequence_array[sequence_nu]->chain_size=i;
        sequence_array[sequence_nu]->seqstr[i-1]='Z';
        sequence_array[sequence_nu]->seqstr[i]='\0';
        i++;
      }
    } /*konec if - new residue */
    } /* konec if ATOM */
  } /*tady konci for pres pole*/
 
} /*tady konci funkce */

void printresult(void)            /* Self-explanatory */

{
  int i,j,m,l;
  int Resnum, Resnumprev;        /* Current and previous residue 
                                     number indicator */
  char Resstring[4],             /* To get residue number... */
    Confprev;                    /* Conformation of the previous atom 
                                    of the same residue */
  int Nres, Nat;   /*numer of atoms in a nonstandard residue, number of nonstandard atoms*/

  Bav=Bav/atcount;

  printf("\n Results:\n");

  printf("\n Number of residues:\n");
  for(i=0; i <=chain_nu; i++) {
    printf("\n Chain %c ", chain_array[i]->chain[0]);
    printf("  %d ", chain_array[i]->rescount);
  };

  printf("\n Number of waters: %d \n", watercount );

  printf("\n\n Nonstandard residues:");

  for (i=0; i<=nonstd_nu; i++) {
 
  if(nonstd_array[i]->checked!=1) {   /*new nonstandard residue?*/
    
    Nres=1;

    for(j=(i+1);j<=nonstd_nu; j++) {
      
      if(strcmp((nonstd_array[i]->resname),(nonstd_array[j]->resname))==0) {  /* are there more residues with the same name*/
         nonstd_array[j]->checked=1;
         Nres++;
      }
    } 
  
     printf("\n\n%3d x %s:  ",Nres,(nonstd_array[i]->resname) );
   
    for (j=i; j<=nonstd_nu;j++) {
      if(strcmp((nonstd_array[i]->resname),(nonstd_array[j]->resname))==0) {

         printf("%c%d ",(nonstd_array[j]->chain),(nonstd_array[j]->Resnum) );
         
      } 
    }
   
    printf("\n  %d non-H atom(s) |", (nonstd_array[i]->Atomcount) );
  
    for (l=0;l<=nonconsh_nu;++l) {  /*nonstandard atoms in the current residue*/
      
     if ((strcmp((nonstd_array[i]->resname),(nonconsh_array[l]->resname))==0)&&((nonconsh_array[l]->checked)!=1)) {
        Nat=0;
        for(m=l; m<=nonconsh_nu;m++) {
          if((strcmp((nonconsh_array[l]->atom),(nonconsh_array[m]->atom))==0)&&(strcmp((nonconsh_array[l]->resname),(nonconsh_array[m]->resname))==0)) {
            Nat++;
            nonconsh_array[m]->checked=1;
          } 
        } 
        Nat=Nat/Nres;
        printf(" %d x %s |",Nat, nonconsh_array[l]->atom ); 
      }
      
    } 
   
  }
  }
  printf("\n"); 

  printf("\n\n Number of atoms with alternate location: %4ld\n", Aparam);
  char Aistring[7];
  double Ai=0;
  double Aisum=0;
  char Atomid[5];
  char Atomidprev[5];
  strcpy(Atomidprev, "1234");
  Resnumprev=0; Confprev=' ';
  locator=alternate_first;
  while (locator!=NULL){
    for (i=0;i<4;i++) Resstring[i]=locator->PDBstring[i+22];
    Resnum=atoi(Resstring);
    if (Resnum!=Resnumprev) {     /* Is it new residue? */
      if(Aisum>1)printf(" (** sum %4.2f! **) ",Aisum);
      printf("\n ");
      for (i=0;i<9;i++) printf("%c",locator->PDBstring[i+17]);
      printf(": ");
      Confprev=' ';
      Aisum=0;
    };
    for (i=0;i<4;i++) Atomid[i]=locator->PDBstring[i+12];  /* Is it a new atom? */
    if (strcmp(Atomid,Atomidprev)!=0) {
      if(Aisum>1)printf(" (** sum %4.2f **) ",Aisum);
      Confprev=' ';
      Aisum=0;
      strcpy(Atomidprev, Atomid);
      
    };


    if (Confprev!=locator->PDBstring[16]) { 
                                  /* Is it a new alternate location? */
      if (Resnum==Resnumprev) printf(", ");
      printf("%c/",locator->PDBstring[16]);
      for(i=0;i<6;i++) printf("%c", locator->PDBstring[i+54]);
      for (i=0;i<6;i++) Aistring[i]=locator->PDBstring[i+54];
      Ai=atof(Aistring);
      Aisum=Aisum+Ai;
    };

 
    Confprev=locator->PDBstring[16];
    Resnumprev=Resnum;
    locator=locator->nextr;
   
  };
  printf("\n\n All ATOMs:\n");
  printf(" Minimum B: %6.2f\n", Bmin);
  printf(" Maximum B: %6.2f\n", Bmax);
  printf(" Average B: %6.2f\n",Bav);

   for(i=0; i <=chain_nu; i++) {
    chain_array[i]->Bav=(chain_array[i]->Bav)/(chain_array[i]->atcount);
    printf("\n Chain %c ", chain_array[i]->chain[0]);
    printf(" Bmin: %6.2f ", chain_array[i]->Bmin);
    printf(" Bmax: %6.2f ", chain_array[i]->Bmax);
    printf(" Bave: %6.2f ", chain_array[i]->Bav);
  };


  printf("\n\n Minimum occupancy: %6.2f\n", Omin);
  printf(" Maximum occupancy: %6.2f\n", Omax);
  printf("\n Number of atoms with occupancy");
  printf(" differing from 0.00 and 1.00: %4ld\n\n", Oparam); 

  printf("\n Atoms with occupancy close to 0:\n\n"); 
  locator=Zccup_first;
  while (locator!=NULL){
    printf("%s\n",locator->PDBstring);
    locator=locator->nextr;
  }

  printf("\n Atoms with the lowest B factors: \n\n"); 
  for (i=0;i<10;i++) {
        printf(" %s %6.2f \n",lowB_array[i].Atiden,lowB_array[i].B );
      };
  
  printf("\n Atoms with the highest B factors: \n\n"); 
  for (i=0;i<10;i++) {
        printf(" %s %6.2f \n",highB_array[i].Atiden,highB_array[i].B );
      };

  if (LONGLIST) {
    int NEGATIVE=0;
    printf("\n Atoms with occupancy differing from 0.00 and 1.00 \n\n");
    locator=Occup_first;
    while (locator!=NULL){
      printf(" %s \n",locator->PDBstring);
      locator=locator->nextr;
    };
    locator=Occup_first;
    while (locator!=NULL){
      for (i=0;i<6;i++) Aistring[i]=locator->PDBstring[i+54];
      Ai=atof(Aistring);
      if(Ai<0){
        if(!NEGATIVE){
          printf("\n Atoms with negative occupancy \n\n");
          NEGATIVE=1;
        }
        printf(" %s \n",locator->PDBstring);
      }
      locator=locator->nextr;
    };
  };

  if(DIFFER) {
    printf("\n\n Atoms shifted by more than %4.2f A\n",told);
    printf("\n                                   shift   delta_x  delta_y  delta_z  delta_B  final B  occupancy");
    printf("\n                                  =======");
    for(i=0; i <=shiftat_nu; i++) {
      printf("\n %s ", shiftat_array[i]->Atiden);
      printf("%s", shiftat_array[i]->warning);
      printf(" %7.3f ", shiftat_array[i]->dd);
      printf(" %7.3f ", shiftat_array[i]->dx);
      printf(" %7.3f ", shiftat_array[i]->dy);
      printf(" %7.3f ", shiftat_array[i]->dz);
      printf(" %7.2f ", shiftat_array[i]->dB);
      printf(" %7.2f ", shiftat_array[i]->B);
      printf(" %7.2f ", shiftat_array[i]->O);
    };
    printf("\n\n\n Atoms differing in B by more than %4.2f A^2\n",tolb);
    printf("\n                                  delta_B  final B   shift   occupancy");
    printf("\n                                  =======");
    for(i=0; i <=difB_nu; i++) {
      printf("\n %s ",difB_array[i]->Atiden);
      printf("%s", difB_array[i]->warning);
      printf("  %6.2f ", difB_array[i]->dB);
      printf("  %6.2f ", difB_array[i]->B);
      printf(" %7.3f ", difB_array[i]->dd);
      printf("  %6.2f ", difB_array[i]->O);
    };
    printf("\n\n\n Atoms appearing both in shift and B tables");
    printf("\n                               shift   delta_x  delta_y  delta_z  delta_B  final B occupancy");
    for(i=0; i <=difAl_nu; i++) {
      printf("\n %s ",difAl_array[i]->Atiden);
    /* printf("%s", difAl_array[i]->warning); */
      printf(" %7.3f ", difAl_array[i]->dd);
      printf(" %7.3f ", difAl_array[i]->dx);
      printf(" %7.3f ", difAl_array[i]->dy);
      printf(" %7.3f ", difAl_array[i]->dz);
      printf("  %6.2f ", difAl_array[i]->dB);
      printf("  %6.2f ", difAl_array[i]->B);
      printf("  %6.2f ", difAl_array[i]->O);
    };
     printf("\n\n Atoms with position or B change higher than tolerance value were printed into following files:\n shifted_atoms.pdb\n B_atoms.pdb\n");
  };

  if (SEQUENCE) {
    printf("\n Sequence:");
    for(i = 0; i<=sequence_nu; i++) {
    if ((sequence_array[i]->chain_size>0)&&(sequence_array[i]->seqstr[0]!='\0'))
    {
      printf("\n\n>%s chain %c ",fname1,sequence_array[i]->chain);
      fprintf(fasta,"\n>%s chain %c ",fname1,sequence_array[i]->chain);
      for(j=0; j<sequence_array[i]->chain_size;j++) {
        if((j+1)%60==1) { 
          printf("\n ");
          fprintf(fasta,"\n ");
        }
        else if((j+1)%10==1) {
          printf(" ");
          fprintf(fasta," ");
        }
        printf("%c", sequence_array[i]->seqstr[j] );
        fprintf(fasta,"%c", sequence_array[i]->seqstr[j] );
        
      }
    }
    }
    printf("\n Sequence was printed to the pdb.fasta file \n"); 
  }
};

void releasememory(void) /*Selfexplanatory*/

{
  int i=0;
  while (alternate_first!=NULL){           /*Releasing an "alternate" list */
    locator=alternate_first->nextr;
    free(alternate_first->PDBstring);
    free(alternate_first);
    alternate_first=locator;
  };
  while (Occup_first!=NULL){                 /*Releasing a "B=50.00" list */
    locator=Occup_first->nextr;
    /*free(Occup_first->PDBstring);*/
    free(Occup_first);
    Occup_first=locator;
  };
   while (Zccup_first!=NULL){                 /*Releasing a "zero occupancy" list */
    locator=Zccup_first->nextr;
   /* free(Zccup_first->PDBstring);*/
    free(Zccup_first);
    Zccup_first=locator;
  };
  for(i=0; i<=chain_nu; i++) free(chain_array[i]);
  free(chain_array);

  for(i = 0; i <atAcount; i++) {
    /*free(atomA_array[i]->PDBstring);*/
    free(atomA_array[i]);
  }
    free(atomA_array);

  for(i = 0; i <adpAcount; i++) {
    free(adpA_array[i]);
  }
    free(adpA_array);

 for(i=0; i<=nonstd_nu; i++) free(nonstd_array[i]);
 free(nonstd_array);

 for(i=0; i<=nonconsh_nu; i++) free(nonconsh_array[i]);
 free(nonconsh_array);
  

  if(DIFFER) {
    for(i=0; i<=shiftat_nu; i++) free(shiftat_array[i]);
    free(shiftat_array);
    for(i=0; i<=difB_nu;i++) free(difB_array[i]);
    free(difB_array);

    for(i=0; i<=difAl_nu;i++) free(difAl_array[i]);
    free(difAl_array);
  }

  if(SEQUENCE){
    for(i = 0; i <=sequence_nu; i++) {
     free(sequence_array[i]->seqstr);
     free(sequence_array[i]);
  }
    free(sequence_array);
  }
};

int main(int argc, char *argv[])      /*Main routine*/
{
  
  char opt; /* options z prikazove radky */
  int argi;

  argi=argc;
  char *comm;
  comm=argv[0];
  /* prepinace -i -v -h -e*/

  while ((--argi > 0 ) && ((*++argv)[0] == '-')) {

    while ((opt = *++argv[0])) {
      switch (opt) {
        case 'l':
          LONGLIST=1;
          break;
          case 's':
          SEQUENCE=1;
          break; 
        case 'd':
          DIFFER=1;
          break;
  /*      case 'w':
          WATER_CONTACTS=1;
          break; */
        default :
          fprintf(stderr,"Illegal option %c\n",opt);
          printhelp(comm);

          exit(1);
      };
    };
  };

  
  if (argi == 0) {
    fprintf(stderr,"Please specify a pdb file.\n");
    printhelp(comm);
    exit(1);
  };
  
  fname1=*argv++;
  argi--;
  if ((finput1 = fopen(fname1,"r")) == NULL) {
    fprintf(stderr,"Can\'t open %s\n",fname1);
    exit(1);
  };

  if (DIFFER) {
    if (argi == 0) {
      fprintf(stderr,"Please specify a second pdb file.\n");
      printhelp(comm);
      fclose(finput1);
      exit(1);
    };
    
    fname2=*argv++;
    argi--;
    if ((finput2 = fopen(fname2,"r")) == NULL) {
      fprintf(stderr,"Can\'t open %s\n",fname2);
      exit(1);
    };

    if (argi>0){
      argi--;
      told=atof(*argv++);
      }
    else told=TOL_D;
    
    if(argi>0){
        argi--;
        tolb=atof(*argv);
        argv++;
      }
    else tolb=TOL_B;
    }   
  
  initparam();
  printf("\n *** %s - 2015/08/17 very unstable version *** \n",comm);
  printf("\n File name:\t*** %s ***",fname1);
  if(DIFFER) {
    printf("\t(%s)",fname2);
  }
  printf("\n\n");
  getatomsA();
  analyserecord();
  
  if (DIFFER) {
    
    f_shiftoutput = fopen("shifted_atoms.pdb","w");
    f_Boutput = fopen("B_atoms.pdb","w");
    
    while (!feof(finput2)) {
      getrecordB();
      analysechanges();
    }
  }
  if (SEQUENCE) {
    fasta = fopen("pdb.fasta","w");
    sequenceanalysis();
    
  }
  printresult(); 

  if (ANISO==1) analyseU();
  
  
  releasememory();
  fclose(finput1);
  if (DIFFER) {
    fclose(finput2);
    fclose(f_shiftoutput);
    fclose(f_Boutput);
   

  };
 if(SEQUENCE) {
  fclose(fasta);
 }

  printf("\n"); 
  
  exit(0);
};
