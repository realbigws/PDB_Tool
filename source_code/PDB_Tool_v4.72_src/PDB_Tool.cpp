#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include "PDB_Chain_Fold.h"
#include "Confo_Back.h"
#include "Confo_Beta.h"
#include "Confo_Lett.h"
#include "Acc_Surface.h"
#include "Mol_File.h"
#include "Mol_Out.h"
#include "getopt.h"
using namespace std;


//-----------------------------------------------------------------------------------------------------------//
//---- print_help_msg => print help message
void print_help_msg(void) 
{
	cout << "========================================================|" << endl;
	cout << "PDB_Tool  (version 4.81) [2020.07.14]                   |" << endl;
	cout << "          a versatile tool to process PDB file          |" << endl;
	cout << "Usage:   ./PDB_Tool <-i input> <-r range> <-o output>   |" << endl;
	cout << "Or,      ./PDB_Tool <-i inroot> <-L list> <-o outroot>  |" << endl;
	cout << "--------------------------------------------------------|" << endl;
	cout << "For <-i input> input type,                              |" << endl;
	cout << "-i input  : input original PDB file, e.g., 1col.pdb     |" << endl;
	cout << "-r range  : input residue range, e.g., A:1-51           |" << endl;
	cout << "-o output : output PDB file, e.g., 1colA.pdb            |" << endl;
	cout << "--------------------------------------------------------|" << endl;
	cout << "For <-L list> input type,                               |" << endl;
	cout << "-i inroot  : input directory of original PDBs           |" << endl;
	cout << "-o outroot : output directory for generated PDBs        |" << endl;
	cout << "The rows in <list> indicate [input range output]        |" << endl;
	cout << "          e.g.,  1col.pdb A:1-51 1colA.pdb              |" << endl;
	cout << "========================================================|" << endl;
	cout << "                Primary options                         |" << endl;
	cout << "--------------------------------------------------------|" << endl;
	cout << "-M mode : Specify the atoms for output.                 |" << endl;
	cout << "          -2,  output CA+CB                             |" << endl;
	cout << "          -1,  output CA                                |" << endl;
	cout << "           0,  output backbone+CB                       |" << endl;
	cout << "          [1], output full atoms (Set as default)       |" << endl;
	cout << "--------------------------------------------------------|" << endl;
	cout << "-N num  : Specify the numbering type for output.        |" << endl;
	cout << "          -1,  full sequential, for atom and residue    |" << endl;
	cout << "           0,  part sequential, for atom only           |" << endl;
	cout << "          [1], PDB numbering (Set as default)           |" << endl;
	cout << "--------------------------------------------------------|" << endl;
	cout << "                Additional options                      |" << endl;
	cout << "--------------------------------------------------------|" << endl;
	cout << "The following arguments for both input types            |" << endl;
	cout << "        -C 1 to consider non-CA atoms                   |" << endl;
	cout << "        -R 1 to reconstruct missing CB and backbone     |" << endl;
	cout << "        -F 1 for AMI,CLE,SSE; 2 for ACC; 4 for FEAT     |" << endl;
	cout << "           these output files could be combined         |" << endl;
	cout << "           8 for output phi/psi/omega and theta/thor    |" << endl;
	cout << "           9 for output Ca-Ca, Cb-Cb, and hb_value      |" << endl;
	cout << "        -m <mapping_seqres> for missing residues        |" << endl;
	cout << "The following arguments only for <-L list> input type   |" << endl;
	cout << "        -G 1 to output three log files                  |" << endl;
	cout << "========================================================|" << endl;
	exit(-1);
}
//---- WebServer's default input ----//
string INPUT_NAM="";
string INPUT_RAN="_";
string INPUT_OUT="";
string INPUT_LIST="";
string INPUT_MISS="";
int LIST_OR_SINGLE=0; //default: single
int INPUT_MODE=1; //main (default:1)
int INPUT_TYPE=1; //main (default:1)
int INPUT_GLYS=1; //vice (default:1)
int INPUT_NOCA=0; //vice (default:0)
int INPUT_RECO=0; //vice (default:0)
int INPUT_FIFI=0; //vice (default:0)
int INPUT_LOGF=0; //vice (default:0)
int WARN_OUT=1;   //vice (default:1)

//-----------------------------------------------------------------------------------------------------------//
//---- parameter editor ----//
static option long_options[] =
{
	{"input",   no_argument,       NULL, 'i'},
	{"range",   no_argument,       NULL, 'r'},
	{"output",  no_argument,       NULL, 'o'},
	{"List",    no_argument,       NULL, 'L'},
	{"Mode",    no_argument,       NULL, 'M'},
	{"Type",    no_argument,       NULL, 'N'},
	{"Glys",    no_argument,       NULL, 'T'},
	{"Noca",    no_argument,       NULL, 'C'},
	{"Reco",    no_argument,       NULL, 'R'},
	{"Fifi",    no_argument,       NULL, 'F'},
	{"Miss",    no_argument,       NULL, 'm'},
	{"Logf",    no_argument,       NULL, 'G'},
	{"Warn",    no_argument,       NULL, 'W'},
	{0, 0, 0, 0}
};
//-----------------------------------------------------------------------------------------------------------//
//---- process_args => process input parameter args
void process_args(int argc,char** argv) 
{
	string buf;
	int opt;
	if(1==argc)print_help_msg();    
	while(true) 
	{
		int option_index=0;
		opt=getopt_long(argc,argv,"i:r:o:L:M:N:T:C:R:F:m:G:W:",
			   long_options,&option_index);
		if (opt==-1)break;	
		switch(opt) 
		{
			case 'i':
				INPUT_NAM=optarg;
				break;
			case 'r':
				INPUT_RAN=optarg;
				break;
			case 'o':
				INPUT_OUT=optarg;
				break;
			case 'L':
				INPUT_LIST=optarg;
				LIST_OR_SINGLE++;
				break;
			case 'M':
				INPUT_MODE=atoi(optarg);
				break;
			case 'N':
				INPUT_TYPE=atoi(optarg);
				break;
			case 'T':
				INPUT_GLYS=atoi(optarg);
				break;
			case 'C':
				INPUT_NOCA=atoi(optarg);
				break;
			case 'R':
				INPUT_RECO=atoi(optarg);
				break;
			case 'F':
				INPUT_FIFI=atoi(optarg);
				break;
			case 'm':
				INPUT_MISS=optarg;
				break;
			case 'G':
				INPUT_LOGF=atoi(optarg);
				break;
			case 'W':
				WARN_OUT=atoi(optarg);
				break;
			default:
				exit(-1);
		}
	}
}

//------------- check file existance -------//
int Check_File_Existance(string &pdbfile)
{
	ifstream fin;
	fin.open(pdbfile.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"pdbfile %s not found!!\n",pdbfile.c_str());
		return -1;
	}
	return 0;
}

//---------- dynamic programming ----------//
int WWW_Advance_Align_Dyna_Prog_Double(int n1,int n2,const vector<double> &score,
								   double GAP_OPEN1,double GAP_EXT1,double GAP_OPEN2,double GAP_EXT2,
								   double GAP_HEAD1,double GAP_TAIL1,double GAP_HEAD2,double GAP_TAIL2,
								   vector<pair<int,int> > & alignment,double &ali_sco)
{
	int i,j;
	//input
	int m = n1 + 1;  // +1 to account for the extra row,col in
	int n = n2 + 1;  // the DP matrices corresponding to gaps
	int DP_maximal=n;
	int IN_maximal=n2;
	//const value
	const int _H_  = 0;
	const int _S_  = 1;
	const int _V_  = 2;

	//create D and M
	vector <int> D[3];      // the path (directions) matrix
	vector <double> M[3];   // the current scores (values) matrix
	//resize(m,n)
	for (i = 0; i < 3; ++i) 
	{
		D[i].resize(m*n);
		M[i].resize(m*n);
	}
	//init()
	double IN_MIN=-1000000;
	D[_S_][0*DP_maximal+ 0] = -1;
	D[_H_][0*DP_maximal+ 0] = -1;
	D[_V_][0*DP_maximal+ 0] = -1;
	M[_S_][0*DP_maximal+ 0] = 0;
	M[_H_][0*DP_maximal+ 0] = IN_MIN;
	M[_V_][0*DP_maximal+ 0] = IN_MIN;
	for (i = 1; i < m; i++) 
	{
		D[_S_][i*DP_maximal+ 0] = _V_;
		D[_H_][i*DP_maximal+ 0] = _V_;
		D[_V_][i*DP_maximal+ 0] = _V_;
		M[_S_][i*DP_maximal+ 0] = IN_MIN;
		M[_H_][i*DP_maximal+ 0] = IN_MIN;
		M[_V_][i*DP_maximal+ 0] = i*GAP_HEAD1; //-(Params::GAP_OPEN + (i-1)*Params::GAP_EXT);
	}
	for (j = 1; j < n; j++) 
	{
		D[_S_][0*DP_maximal+ j] = _H_;
		D[_H_][0*DP_maximal+ j] = _H_;
		D[_V_][0*DP_maximal+ j] = _H_;
		M[_S_][0*DP_maximal+ j] = IN_MIN;
		M[_H_][0*DP_maximal+ j] = j*GAP_HEAD2; //-(Params::GAP_OPEN + (j-1)*Params::GAP_EXT);
		M[_V_][0*DP_maximal+ j] = IN_MIN;
	}
	//fill(firstSeq, secondSeq, distFunc);
	double gap_open;
	double gap_ext;
	double v1,v2,v3;
	double dist;
	for (i = 1; i < m; i++) 
	{
		for (j = 1; j < n; j++) 
		{
			//condition upper
			if(j==n-1)
			{
				gap_open=GAP_TAIL1;
				gap_ext=GAP_TAIL1;
			}
			else
			{
				gap_open=GAP_OPEN1;
				gap_ext=GAP_EXT1;
			}
			v1 = M[_V_][(i-1)*DP_maximal+ j] + gap_ext;
			v2 = M[_S_][(i-1)*DP_maximal+ j] + gap_open;
			v3 = M[_H_][(i-1)*DP_maximal+ j] + gap_open;
			M[_V_][i*DP_maximal+ j] = std::max(v1, std::max(v2, v3));
			if (M[_V_][i*DP_maximal+ j] == v1) D[_V_][i*DP_maximal+ j] = _V_;
			else if(M[_V_][i*DP_maximal+ j] == v2) D[_V_][i*DP_maximal+ j] = _S_;
			else D[_V_][i*DP_maximal+ j] = _H_;
			//condition left
			if(i==m-1)
			{
				gap_open=GAP_TAIL2;
				gap_ext=GAP_TAIL2;
			}
			else
			{
				gap_open=GAP_OPEN2;
				gap_ext=GAP_EXT2;
			}
			v1 = M[_H_][i*DP_maximal+ j-1] + gap_ext;
			v2 = M[_S_][i*DP_maximal+ j-1] + gap_open;
			v3 = M[_V_][i*DP_maximal+ j-1] + gap_open;
			M[_H_][i*DP_maximal+ j] = std::max(v1, std::max(v2, v3));
			if (M[_H_][i*DP_maximal+ j] == v1) D[_H_][i*DP_maximal+ j] = _H_;
			else if(M[_H_][i*DP_maximal+ j] == v2) D[_H_][i*DP_maximal+ j] = _S_;
			else D[_H_][i*DP_maximal+ j] = _V_;
			//condition diag
			dist = score.at((i-1)*IN_maximal+ j-1);  //Params::K - distFunc(firstSeq[i-1], secondSeq[j-1]);
			v1 = M[_V_][(i-1)*DP_maximal+ j-1] + dist;
			v2 = M[_H_][(i-1)*DP_maximal+ j-1] + dist;
			v3 = M[_S_][(i-1)*DP_maximal+ j-1] + dist;
			M[_S_][i*DP_maximal+ j] = std::max(v1, std::max(v2, v3));
			if (M[_S_][i*DP_maximal+ j] == v3) D[_S_][i*DP_maximal+ j] = _S_;
			else if (M[_S_][i*DP_maximal+ j] == v1) D[_S_][i*DP_maximal+ j] = _V_;
			else D[_S_][i*DP_maximal+ j] = _H_;
		}
	}
	//build(ali, firstSeq, secondSeq, distFunc);
	i = m-1;
	j = n-1;
	v1=M[_V_][i*DP_maximal+ j];
	v2=M[_H_][i*DP_maximal+ j];
	v3=M[_S_][i*DP_maximal+ j];
	double maximal = std::max(v1, std::max(v2, v3));
	int k = -1;
	if(v3==maximal)k = _S_;
	else if(v2==maximal)k = _H_;
	else k = _V_;
	//trace_back
	alignment.clear();
	int count = 0;
	int matches = 0;
	int cur_case=k;
	int pre_case;
	for(;;)
	{
		if(i==0||j==0)break;
		pre_case=D[cur_case][i*DP_maximal+ j];
		switch (cur_case)
		{
			case _S_:
				alignment.push_back(pair<int,int>(i,j)); 
				i--;
				j--;
				++matches;
				break;
			case _V_:
				alignment.push_back(pair<int,int>(i,-j)); 
				i--;
				break;
			case _H_:
				alignment.push_back(pair<int,int>(-i,j)); 
				j--;
				break;
			default:
				cout << "ERROR!! -> advance_global: invalid direction D[" << k << "](" << i << ", " << j << ") = " 
				<< D[k][i*DP_maximal+ j] << endl;
				exit(-1);
		}
		cur_case=pre_case;
		count++;
	}
	while (j> 0) alignment.push_back(pair<int,int>(-i,j)),j--;
	while (i> 0) alignment.push_back(pair<int,int>(i,0)), i--;
	reverse(alignment.begin(), alignment.end());
	ali_sco=maximal;
	return matches;
}
int Seqres_DynaProg(string &seqres,string &ami_,vector <int> &mapping)
{
	//--[0]check
	int len=(int)seqres.length();
	int totnum=(int)ami_.length();

	//--[1]dynamic_programming
	int i,j;
	int head=0;
	int n1=len;    //SEQRES
	int n2=totnum; //ATOM
	vector <double> score;
	score.resize(len*totnum);
	for(i=0;i<n1;i++)
	{
		for(j=0;j<n2;j++)
		{
			if(seqres[i]==ami_[j+head])score.at(i*n2+j)=10;
			else
			{
				if(seqres[i]=='X'||seqres[i]=='Z'||seqres[i]=='.')score.at(i*n2+j)=0;
				else if(ami_[j+head]=='X'||ami_[j+head]=='Z'||ami_[j+head]=='.')score.at(i*n2+j)=0;
				else score.at(i*n2+j)=-15;
			}
		}
	}
	double sco;
	int matchs;
	vector<pair<int,int> > WWW_alignment;
	matchs=WWW_Advance_Align_Dyna_Prog_Double(n1,n2,score,-11,-1,-110,-10,0,0,0,0,
		WWW_alignment,sco);
	int lcmp=(int)WWW_alignment.size();
	
	//extract
	mapping.resize(len);
	for(i=0;i<len;i++)mapping[i]=0; //default: NO
	int first,second;
	int retv=1;
	for(i=0;i<lcmp;i++)
	{
		first=WWW_alignment[i].first;
		second=WWW_alignment[i].second;
		if(first<=0)
		{
			if(second>0)
			{
				retv=-1;
				continue;
			}
		}
		if(first>0 && second>0)
		{
			mapping[first-1]=1;
		}
	}
	return retv;
}

//---------- read FASTA seqres ------------//
int Read_FASTA_SEQRES(string &seqfile,string &seqres,int skip=1) //->from .fasta file
{
	ifstream fin;
	string buf,temp;
	//read
	fin.open(seqfile.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"no such file! %s \n",seqfile.c_str());
		exit(-1);
	}
	//skip
	int i;
	for(i=0;i<skip;i++)
	{
		if(!getline(fin,buf,'\n'))
		{
			fprintf(stderr,"file bad! %s \n",seqfile.c_str());
			exit(-1);
		}
	}
	//process
	temp="";
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		temp+=buf;
	}
	seqres=temp;
	//return
	return (int)seqres.length();
}


//------------- Get_PDB_File_Len -----------//
int Get_PDB_File_Len(string &pdbfile) //-> only suitable for pdb_BC100 pdb_file
{
	//--- list for mapping ---//
	map<string, int > ws_mapping;
	map<string, int>::iterator iter;
	//read
	ifstream fin;
	string buf,temp;
	fin.open(pdbfile.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"pdbfile %s not found!!\n",pdbfile.c_str());
		exit(-1);
	}
	//process
	int len;
	int count=0;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		len=(int)buf.length();
		if(len<3)continue;
		//check ATOM
		if(len<4)continue;
		temp=buf.substr(0,4);
		if(temp!="ATOM" && temp!="HETA")continue;
		//check CA
		temp=buf.substr(13,2);
		if(temp!="CA")continue;
		//record name
		temp=buf.substr(21,6);
		iter = ws_mapping.find(temp);
		if(iter != ws_mapping.end())continue;
		count++;
		ws_mapping.insert(map < string, int >::value_type(temp, count));
	}
	//return
	return count;
}

//------ output protein feature files -------//2015_02_20//
//--> output AMI_SSE_CLE
void Output_Protein_Features_AMI_SSE_CLE(
	string &outroot,string &outname,int moln,PDB_Residue *pdb,
	vector <int> &mapping,string &miss_seq)
{
	//init
	int i;
	FILE *fpp;
	string file;
	PDB_Chain_Fold chain_fold;
	chain_fold.initialize_simple(moln,' ');
	for(i=0;i<moln;i++)chain_fold.set_residue(i,pdb[i]);
	int retv;
	//calc
	retv=chain_fold.calculate_CLE();
	if(retv!=0)
	{
		fprintf(stderr,"[%s]CLE_BAD!!!\n",outname.c_str());
		exit(-1);
	}
	retv=chain_fold.calculate_SSE();
	if(retv!=0)
	{
		fprintf(stderr,"[%s]SSE_BAD!!!\n",outname.c_str());
		exit(-1);
	}
	//------ output AMI,CLE,SSE ------//
	string AMI,CLE,SSE;
	AMI=chain_fold.get_sequence();
	CLE=chain_fold.get_CLE();
	SSE=chain_fold.get_SSE();
	//output AMI
	file="";
	file=file+outroot+"/"+outname+".ami";
	fpp=fopen(file.c_str(),"wb");
	if(fpp==0)
	{
		fprintf(stderr,"ERROR: file %s can't be opened. \n",file.c_str());
	}
	else
	{
		fprintf(fpp,">%s\n",outname.c_str());
		int cur=0;
		for(i=0;i<(int)mapping.size();i++)
		{
			if(mapping[i]==1)
			{
				fprintf(fpp,"%c",AMI[cur]);
				cur++;
			}
			else
			{
				fprintf(fpp,"-");
			}
		}
		fprintf(fpp,"\n");
		fclose(fpp);
	}
	//output CLE
	file="";
	file=file+outroot+"/"+outname+".cle";
	fpp=fopen(file.c_str(),"wb");
	if(fpp==0)
	{
		fprintf(stderr,"ERROR: file %s can't be opened. \n",file.c_str());
	}
	else
	{
		fprintf(fpp,">%s\n",outname.c_str());
		int cur=0;
		for(i=0;i<(int)mapping.size();i++)
		{
			if(mapping[i]==1)
			{
				fprintf(fpp,"%c",CLE[cur]);
				cur++;
			}
			else
			{
				fprintf(fpp,"-");
			}
		}
		fprintf(fpp,"\n");
		fclose(fpp);
	}
	//output SSE
	file="";
	file=file+outroot+"/"+outname+".sse";
	fpp=fopen(file.c_str(),"wb");
	if(fpp==0)
	{
		fprintf(stderr,"ERROR: file %s can't be opened. \n",file.c_str());
	}
	else
	{
		fprintf(fpp,">%s\n",outname.c_str());
		int cur=0;
		for(i=0;i<(int)mapping.size();i++)
		{
			if(mapping[i]==1)
			{
				fprintf(fpp,"%c",SSE[cur]);
				cur++;
			}
			else
			{
				fprintf(fpp,"-");
			}
		}
		fprintf(fpp,"\n");
		fclose(fpp);
	}
}

//--> output ACC and ACC_Value
void Output_Protein_Features_ACC(
	string &outroot,string &outname,Acc_Surface *acc_surface,
	int moln,PDB_Residue *pdb,XYZ **mcc,int *mcc_side,char *ami,char *acc,
	vector <int> &mapping,string &miss_seq)
{
	//init
	int i;
	FILE *fpp;
	string file;
	PDB_Chain_Fold chain_fold;
	chain_fold.initialize_simple(moln,' ');
	for(i=0;i<moln;i++)chain_fold.set_residue(i,pdb[i]);
	//------ output ACC_Code and ACC_Value ------//
	for(i=0;i<moln;i++)pdb[i].get_XYZ_array(mcc[i],mcc_side[i]);
	acc_surface->AC_Calc_SolvAcc(mcc,ami,moln,acc,mcc_side);
	//output ACC_Code
	file="";
	file=file+outroot+"/"+outname+".acc";
	fpp=fopen(file.c_str(),"wb");
	if(fpp==0)
	{
		fprintf(stderr,"ERROR: file %s can't be opened. \n",file.c_str());
	}
	else
	{
		fprintf(fpp,">%s\n",outname.c_str());
		int cur=0;
		for(i=0;i<(int)mapping.size();i++)
		{
			if(mapping[i]==1)
			{
				fprintf(fpp,"%c",acc[cur]);
				cur++;
			}
			else
			{
				fprintf(fpp,"-");
			}
		}
		fprintf(fpp,"\n");
		fclose(fpp);
	}
	//output ACC_Value
	file="";
	file=file+outroot+"/"+outname+".acc_value";
	fpp=fopen(file.c_str(),"wb");
	if(fpp==0)
	{
		fprintf(stderr,"ERROR: file %s can't be opened. \n",file.c_str());
	}
	else
	{
		int cur=0;
		for(i=0;i<(int)mapping.size();i++)
		{
			if(mapping[i]==1)
			{
				char c=ami[cur];
				fprintf(fpp,"%4d %c %3d %3d\n",i+1,c,acc_surface->AC_normal[cur],acc_surface->AC_output[cur]);
				cur++;
			}
			else
			{
				char c='-';
				if(miss_seq!="")c=miss_seq[i];
				fprintf(fpp,"%4d %c   -   -\n",i+1,c);
			}
		}
		fclose(fpp);
	}
}

//------ output protein feature in one file -------//2015_02_20//
// file format should be consistent with TPL file
/*
//////////// Features
  Num Res  Missing   SSE    CLE   ACC   pACC  CNa CNb   Xca       Yca       Zca       Xcb       Ycb       Zcb
   1   E      0       L      1     2     87   2   1     19.400     4.600    31.600    17.889     4.542    31.872
   2   N      0       E      5     2     48   4   3     20.500     7.600    33.700    21.326     8.630    32.887
   3   I      0       E      5     1     18   3   5     18.900     9.100    36.800    18.496     8.208    37.931
   4   E      0       E      5     2     53   3   2     19.900    12.700    37.700    19.630    13.781    36.734
   5   V      0       E      5     0      0   6   9     20.200    13.600    41.400    20.814    12.502    42.276
   6   H      0       E      5     1     32   6   3     20.700    17.200    42.800    19.590    18.263    42.565
   7   M      0       E      5     0      0   9   8     22.700    17.900    45.900    24.187    17.551    45.967
   8   L      0       E      5     1     12   5   4     21.100    20.900    47.700    19.558    20.886    47.424
   9   N      0       E      5     1     33   4   3     21.400    23.000    50.900    22.004    24.410    50.855
*/
//-------- ACC<->Int -------//
int ACC_To_Int(char c)
{
	switch(c)
	{
		case 'B': return 0;
		case 'M': return 1;
		case 'E': return 2;
		default: return 1;
	}
}
//----- output protein features -----//
void Output_Protein_Features(
	string &outroot,string &outname,Acc_Surface *acc_surface,XYZ *mol,XYZ *mcb,
	int moln,PDB_Residue *pdb,XYZ **mcc,int *mcc_side,char *ami,char *acc,
	vector <int> &mapping,string &miss_seq)
{
	//init
	int i,j;
	int retv;
	PDB_Chain_Fold chain_fold;
	chain_fold.initialize_simple(moln,' ');
	for(i=0;i<moln;i++)chain_fold.set_residue(i,pdb[i]);

	//------ calculate AMI,CLE,SSE ------//
	retv=chain_fold.calculate_CLE();
	if(retv!=0)
	{
		fprintf(stderr,"[%s]CLE_BAD!!!\n",outname.c_str());
		exit(-1);
	}
	retv=chain_fold.calculate_SSE();
	if(retv!=0)
	{
		fprintf(stderr,"[%s]SSE_BAD!!!\n",outname.c_str());
		exit(-1);
	}
	string AMI,CLE,SSE;
	AMI=chain_fold.get_sequence();
	CLE=chain_fold.get_CLE();
	SSE=chain_fold.get_SSE();

	//------ calculate ACC_Code and ACC_Value ------//
	for(i=0;i<moln;i++)pdb[i].get_XYZ_array(mcc[i],mcc_side[i]);
	acc_surface->AC_Calc_SolvAcc(mcc,ami,moln,acc,mcc_side);

	//------ calculate contact number for CA and CB ---//
	double DIST_CUTOFF = 64; //-> 8.0A (note that in older version of TPL, the CA/CB contact cutoff is 7.0A
	vector <int> cn_ca(moln,0);
	vector <int> cn_cb(moln,0);
	for(i=0;i<moln;i++)
	{
		for(j=0;j<moln;j++)
		{
			double dist_ca=mol[i].distance_square(mol[j]);
			double dist_cb=mcb[i].distance_square(mcb[j]);
			//-> normal condition
			if( abs(i-j)>=4 )
			{
				if(dist_ca <= DIST_CUTOFF)cn_ca[i]++;
				if(dist_cb <= DIST_CUTOFF)cn_cb[i]++;
			}
			//-> chain broken
			else if( abs(i-j)>=1 )
			{
				//check PDB_residue
				string str1,str2;
				pdb[i].get_PDB_residue_number(str1);
				pdb[j].get_PDB_residue_number(str2);
				str1=str1.substr(1,4);
				str2=str2.substr(1,4);
				int pos1=atoi(str1.c_str());
				int pos2=atoi(str2.c_str());
				if( abs(pos1-pos2)>=4 )
				{
					if(dist_ca <= DIST_CUTOFF)cn_ca[i]++;
					if(dist_cb <= DIST_CUTOFF)cn_cb[i]++;
				}
			}
		}
	}

/*
	//------ calculate core region ------//
	int MinHelixLen = 4;
	int MinBetaLen = 3;
	int MinContact = 1;
	vector <int> Core(moln,1);
	for(i=0;i<moln;i++)
	{
		//calculate core
		Core[i] = 1;
		if(SSE[i]!='H' && SSE[i]!='E') continue;
		Core[i] = 2;
		int sslen = 1;
		for(j=i-1;j>0;j--)
		{
			if(SSE[j]!=SSE[i]) break;
			sslen++;
		}
		for(j=i+1;j<moln;j++)
		{
			if(SSE[j]!=SSE[i]) break;
			sslen++;
		}
		if(SSE[i]=='H' && sslen >= MinHelixLen && cn_ca[i]>MinContact)Core[i] = 5;
		if(SSE[i]=='E' && sslen >= MinBetaLen  && cn_ca[i]>MinContact)Core[i] = 5;
	}
*/

	//------ output feature files -------//
	string file=outroot+"/"+outname+".feature";
	FILE *fpp=fopen(file.c_str(),"wb");
	if(fpp==0)
	{
		fprintf(stderr,"ERROR: file %s can't be opened. \n",file.c_str());
	}
	else
	{
		fprintf(fpp,"  Num Res  Missing   SSE    CLE   ACC   pACC  CNa CNb   Xca       Yca       Zca       Xcb       Ycb       Zcb\n");
		int cur=0;
		for(i=0;i<(int)mapping.size();i++)
		{
			if(mapping[i]==1)
			{
				char c=ACC_To_Int(acc[cur]);
				fprintf(fpp,"%4d   %c      0       %c      %c     %1d    %3d  %2d  %2d   %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f\n",
					i+1,AMI[cur],SSE[cur],CLE[cur],c,acc_surface->AC_normal[cur],cn_ca[cur],cn_cb[cur],mol[cur].X,mol[cur].Y,mol[cur].Z,mcb[cur].X,mcb[cur].Y,mcb[cur].Z);
				cur++;
			}
			else
			{
				char c='-';
				if(miss_seq!="")c=miss_seq[i];
				fprintf(fpp,"%4d   %c      1       L      R     2     45   0   0  \n",i+1,c);
			}
		}
	}
	fclose(fpp);
}



//----- output protein angles -----//
//-> output format
/*
>POS  PHI  PSI  OMEGA  THETA  THOR
....

*/
void Output_Protein_Angles(
	string &outroot,string &outname,
	int moln,PDB_Residue *pdb,char *ami,
	vector <int> &mapping,string &miss_seq)
{
	//init
	int i;
	int retv;
	PDB_Chain_Fold chain_fold;
	chain_fold.initialize_simple(moln,' ');
	for(i=0;i<moln;i++)chain_fold.set_residue(i,pdb[i]);
	//calculate phi_psi_omega
	if(chain_fold.calculate_phi_psi_omega()!=0)
	{
		fprintf(stderr,"phi_psi_omega calculation cailed \n");
		exit(-1);
	}
	//calculate theta_tor
	if(chain_fold.calculate_theta_tau()!=0)
	{
		fprintf(stderr,"theta_tau calculation cailed \n");
		exit(-1);
	}
	//output phi_psi_omega
	vector <vector <double> > phi_psi_omega_out;
	chain_fold.print_phi_psi_omega(phi_psi_omega_out);
	//output theta_tor
	vector <vector <double> > theta_tor_out;
	chain_fold.print_theta_tau(theta_tor_out);

	//------ output feature files -------//
	string file=outroot+"/"+outname+".angles";
	FILE *fpp=fopen(file.c_str(),"wb");
	if(fpp==0)
	{
		fprintf(stderr,"ERROR: file %s can't be opened. \n",file.c_str());
	}
	else
	{
		fprintf(fpp,">POS Res  PHI        PSI       OMEGA     THETA       THOR  \n");
		int cur=0;
		for(i=0;i<(int)mapping.size();i++)
		{
			if(mapping[i]==1)
			{
				//print
				stringstream os;
				os<<setw(4)<<i + 1<<" ";
				os<<setw(1)<<ami[cur]<<" ";
				os<<setw(10)<<phi_psi_omega_out.at(cur).at(0) / M_PI * 180<<" ";
				os<<setw(10)<<phi_psi_omega_out.at(cur).at(1) / M_PI * 180<<" ";
				os<<setw(10)<<phi_psi_omega_out.at(cur).at(2) / M_PI * 180<<" ";
				os<<setw(10)<<theta_tor_out.at(cur).at(0) / M_PI * 180<<" ";
				os<<setw(10)<<theta_tor_out.at(cur).at(1) / M_PI * 180<<" ";
				//out
				string buf=os.str();
				fprintf(fpp,"%s\n",buf.c_str());
				cur++;
			}
			else
			{
				char c='-';
				if(miss_seq!="")c=miss_seq[i];
				//print
				stringstream os;
				os<<setw(4)<<i + 1<<" ";
				os<<setw(1)<<c<<" ";
				os<<setw(10)<<" - "<<" ";
				os<<setw(10)<<" - "<<" ";
				os<<setw(10)<<" - "<<" ";
				os<<setw(10)<<" - "<<" ";
				os<<setw(10)<<" - "<<" ";
				//out
				string buf=os.str();
				fprintf(fpp,"%s\n",buf.c_str());
			}
		}
	}
	fclose(fpp);
}


//---------- output Jinbo'style distance matrix ------------//__180520__//
//-> example (in Hai-cang format) [starting from 1] [no gap]
// CaCaMatrix = np.ones( (fullSeqLen, fullSeqLen), np.float16) * (-1)
// CbCbMatrix = np.ones( (fullSeqLen, fullSeqLen), np.float16) * (-1)
// CgCgMatrix = np.ones( (fullSeqLen, fullSeqLen), np.float16) * (-1)  -> not consider now
// CaCgMatrix = np.ones( (fullSeqLen, fullSeqLen), np.float16) * (-1)  -> not consider now
// NOMatrix   = np.ones( (fullSeqLen, fullSeqLen), np.float16) * (-1)
/*
i j <Ca-Ca> <Cb-Cb> <hydro_val>
......

[legend]:

the first and second column is the two positions in 1-base (i.e., starting from 1).
the third column is the Ca-Ca distance. (-1 means not available)
the fourth column is the Cb-Cb distance. (-1 means not available)
the fifth column is the hydro_bond value. (-1 means not available)

#----------- equation of hydro_bond value ----------------#
#-> in the order of N,Ca,C,O,H
double Hydro_Bond::HB_Calc_Single(XYZ **HB_mol,int i,int j)
{
        double dho,dhc,dnc,dno;
        dho=HB_mol[i][4].distance(HB_mol[j][3]);
        dhc=HB_mol[i][4].distance(HB_mol[j][2]);
        dno=HB_mol[i][0].distance(HB_mol[j][3]);
        dnc=HB_mol[i][0].distance(HB_mol[j][2]);
        return (332.0*0.42*0.2*(1.0/dno+1.0/dhc-1.0/dho-1.0/dnc));
}

*/

void Output_Protein_Distances(
	string &outroot,string &outname,Hydro_Bond *hydro_bond,
	int moln,PDB_Residue *pdb,XYZ **mcc, int *mcc_side,char *ami,
	vector <int> &mapping,string &miss_seq)
{
	//init
	int i,j;
	//calc hydro_bond
	for(i=0;i<moln;i++)pdb[i].get_XYZ_array(mcc[i],mcc_side[i]);
	vector <vector <double> > hb_mat;
	hydro_bond->HB_Input_Mol(mcc,ami,moln);
	hydro_bond->HB_Calc_Hydro_Bond(hb_mat);

	//------ output feature files -------//
	string file=outroot+"/"+outname+".distances";
	FILE *fpp=fopen(file.c_str(),"wb");
	if(fpp==0)
	{
		fprintf(stderr,"ERROR: file %s can't be opened. \n",file.c_str());
	}
	else
	{
		fprintf(fpp,">PS1  PS2 A A   CA_CA      CB_CB      HB_val      \n");
		int cur1=0;
		for(i=0;i<(int)mapping.size();i++)
		{
			if(mapping[i]==0)
			{
				for(j=0;j<(int)mapping.size();j++)
				{
					char c1='-';
					if(miss_seq!="")c1=miss_seq[i];
					char c2='-';
					if(miss_seq!="")c2=miss_seq[j];
					//print
					stringstream os;
					os<<setw(4)<<i + 1<<" ";
					os<<setw(4)<<j + 1<<" ";
					os<<setw(1)<<c1<<" ";
					os<<setw(1)<<c2<<" ";
					os<<setw(10)<<" - "<<" ";
					os<<setw(10)<<" - "<<" ";
					os<<setw(10)<<" - "<<" ";
					//out
					string buf=os.str();
					fprintf(fpp,"%s\n",buf.c_str());
				}
				continue;
			}
			int cur2=0;
			for(j=0;j<(int)mapping.size();j++)
			{
				if(mapping[j]==1)
				{
					//get atoms
					XYZ ca_i,ca_j;
					pdb[cur1].get_backbone_atom( "CA ",ca_i );
					pdb[cur2].get_backbone_atom( "CA ",ca_j );
					XYZ cb_i,cb_j;
					pdb[cur1].get_sidechain_atom( "CB ",cb_i );
					pdb[cur2].get_sidechain_atom( "CB ",cb_j );
					//get distance
					double distance;
					vector <double> distances;
					//-> Ca-Ca distance (symmetric)
					distance=ca_i.distance(ca_j);
					distances.push_back(distance);
					//-> Cb-Cb distance (symmetric)
					distance=cb_i.distance(cb_j);
					distances.push_back(distance);
					//-> HB value
					double hb_val=hb_mat[cur1][cur2];
					//print
					stringstream os;
					os<<setw(4)<<i + 1<<" ";
					os<<setw(4)<<j + 1<<" ";
					os<<setw(1)<<ami[cur1]<<" ";
					os<<setw(1)<<ami[cur2]<<" ";
					os<<setw(10)<<distances.at(0)<<" ";
					os<<setw(10)<<distances.at(1)<<" ";
					os<<setw(10)<<hb_val<<" ";
					//out
					string buf=os.str();
					fprintf(fpp,"%s\n",buf.c_str());
					cur2++;
				}
				else
				{
					char c1='-';
					if(miss_seq!="")c1=miss_seq[i];
					char c2='-';
					if(miss_seq!="")c2=miss_seq[j];
					//print
					stringstream os;
					os<<setw(4)<<i + 1<<" ";
					os<<setw(4)<<j + 1<<" ";
					os<<setw(1)<<c1<<" ";
					os<<setw(1)<<c2<<" ";
					os<<setw(10)<<" - "<<" ";
					os<<setw(10)<<" - "<<" ";
					os<<setw(10)<<" - "<<" ";
					//out
					string buf=os.str();
					fprintf(fpp,"%s\n",buf.c_str());
				}
			}
			cur1++;
		}
	}
	fclose(fpp);
}





//==================== PDB_Back_Process ===============// (process list)
//[list_style]
//first_line: input_dir
//second_line: output_dir
//others:
//input range output-> e.g., 1col.pdb A:1-51 1colA.res
//[note] -> should only contain two  ' '
int PDB_Back_Process(string &input_dir,string &list,string &output_dir,
	int OutType,int OutMode,int OutGlys,int OutNoca,int OutReco,int OutFifi,int OutLogf,int OutWarn)
{
	int totlen=30000;
	//class
	Mol_File mol_input;
	Mol_Out mol_output;
	Confo_Lett confo_lett;
	Confo_Beta *confo_beta=new Confo_Beta(totlen);
	Confo_Back *confo_back=new Confo_Back(totlen);
	Acc_Surface *acc_surface=new Acc_Surface(totlen);
	Hydro_Bond *hydro_bond=new Hydro_Bond(totlen);
	//init
	ifstream fin;
	string buf;
	string input,chain,output;
	string file,temp;
	string path,outa;
	char c;
	int wwscount=0;
	int i;
	int len;
	int pos[10];
	int wscount;
	int ret_val;
	int moln;
	PDB_Residue *pdb=new PDB_Residue[totlen];
	FILE *fp=0;
	FILE *fq=0;
	if(OutLogf==1)
	{
		fp=fopen("ws_real_list","wb");
		fq=fopen("ws_valid_list","wb");
		//judge
		if(fp==0 || fq==0)
		{
			fprintf(stderr,"ERROR: logfile can't be opened. \n");
			OutLogf=0;
		}
	}
	string TER="TER                                                                             ";


	//recon_related
	XYZ *mol=new XYZ[totlen];   //CA
	XYZ *mcb=new XYZ[totlen];   //CB
	XYZ **mbb;                  //BackBone (N,CA,C,O,CB)
	XYZ **mcc;                  //BackBone (N,CA,C,O,CB,...,)
	NewArray2D(&mbb,totlen,5);
	NewArray2D(&mcc,totlen,15);
	int *mcc_side=new int[totlen];
	char *ami=new char[totlen+1];
	char *cle=new char[totlen+1];
	char *acc=new char[totlen+1];


	//open
	fin.open(list.c_str(),ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"Process_List Not Found!![%s]\n",list.c_str());
		fin.close();
		fin.clear();
		return -1;
	}
	if(OutLogf==1)
	{
		mol_input.flog=fopen("ws_pdb_log","wb");
		mol_input.LOGOUT=1;
	}
	mol_input.CbBACK=1;
	if(OutNoca==1)mol_input.CaONLY=0; //consider Non-CA atoms !!//__110408__//
	//macro
	mol_input.OUTP_MODE=OutType;
	mol_input.PROC_MODE=OutMode;
	mol_input.GLYC_MODE=OutGlys;
	//memory limit
	mol_input.WARNING_out=OutWarn;
	mol_input.MEMORY_LIMIT=totlen;


	//--- get input and output
	path=input_dir;
	outa=output_dir;

/*
	//---just_temp---//[get path]
	{
		if(!getline(fin,buf,'\n'))
		{
			fprintf(stderr,"NO_INPUT_DIRECTORY[%s]!!!\n",buf.c_str());
			exit(-1);
		}
		path=buf;
		if(!getline(fin,buf,'\n'))
		{
			fprintf(stderr,"NO_OUTPUT_DIRECTORY[%s]!!!\n",buf.c_str());
			exit(-1);
		}
		outa=buf;
	}
*/
	//process
	int zero=-1;
	int count=0;
	for(;;)
	{
		//read
		if(!getline(fin,buf,'\n'))break;
		getline_end(buf,0x0D);
		//check
		len=(int)buf.length();
		wscount=0;
		for(i=0;i<len;i++)
		{
			if(buf[i]==' ')
			{
				pos[wscount]=i;
				wscount++;
			}
		}
		if(wscount!=2)
		{
			fprintf(stderr,"LIST_Format ERROR!!![%s]\n",buf.c_str());
			exit(-2);
		}
		//calc
		input=buf.substr(0,pos[0]);
		chain=buf.substr(pos[0]+1,pos[1]-pos[0]-1);
		c=chain[0];
		output=buf.substr(pos[1]+1,len-pos[1]-1);
		file=path+'/'+input;
		temp="0 "+file+' '+chain;

		//check the existance of file
		ret_val=Check_File_Existance(file);
		if(ret_val==-1)
		{
			if(OutLogf==1)fprintf(fq,"%s X -1\n",buf.c_str());
			continue;
		}

		//pre_process
		ret_val=mol_input.XYZ_Tranform(temp,moln,0,0,0,0,0,0);
		if(ret_val<=0)
		{
			if(ret_val!=-12345)
			{
				if(OutLogf==1)fprintf(fq,"%s @ %5d\n",buf.c_str(),zero);
				continue;
			}
		}
		//check memory
		if(ret_val==-12345)
		{
			//add memory
			if(moln>totlen)
			{
				totlen=moln;
				//delete
				delete [] pdb;
				delete [] mol;
				delete [] mcb;
				delete [] ami;
				delete [] cle;
				delete [] acc;
				DeleteArray2D(&mbb,totlen);
				DeleteArray2D(&mcc,totlen);
				delete [] mcc_side;
				delete confo_beta;
				delete confo_back;
				delete acc_surface;
				delete hydro_bond;
				//create
				pdb=new PDB_Residue[totlen];
				mol=new XYZ[totlen];
				mcb=new XYZ[totlen];
				ami=new char[totlen+1];
				cle=new char[totlen+1];
				acc=new char[totlen+1];
				NewArray2D(&mbb,totlen,5);
				NewArray2D(&mcc,totlen,15);
				mcc_side=new int[totlen];
				confo_beta=new Confo_Beta(totlen);
				confo_back=new Confo_Back(totlen);
				acc_surface=new Acc_Surface(totlen);
				hydro_bond=new Hydro_Bond(totlen);
				//memory limit
				mol_input.MEMORY_LIMIT=totlen;
			}
		}

		//reload
		int PRE_LOAD_=mol_input.PRE_LOAD;
		int WARNING_out_=mol_input.WARNING_out;
		mol_input.PRE_LOAD=1;
		mol_input.WARNING_out=0;
		mol_input.XYZ_Tranform(temp,moln,0,mol,ami,0,0,pdb);
		mol_input.PRE_LOAD=PRE_LOAD_;
		mol_input.WARNING_out=WARNING_out_;


		//check
		{
			int i,j;
			int correct;
			int iret;
			correct=1; //default:OK
			for(i=0;i<moln;i++)
			{
				iret=pdb[i].PDB_residue_backbone_check(4);
				if(iret!=1)
				{
					correct=0;
					break;
				}
				iret=pdb[i].PDB_residue_CB_check();
				if(iret!=1)
				{
					correct=0;
					break;
				}
			}

			//judge
			if(correct==1)
			{
				wwscount++;
				if(OutLogf==1)
				{
					fprintf(fp,"%s\n",output.c_str());
					fprintf(fq,"%s 1 %5d\n",buf.c_str(),moln);
				}
			}
			else
			{
				wwscount++;
				if(OutLogf==1)
				{
					fprintf(fp,"%s\n",output.c_str());
					fprintf(fq,"%s 0 %5d\n",buf.c_str(),moln);
				}
				if(OutReco==1)
				{
					//[1]recon
					confo_lett.btb_ori(0,0,0,moln,mol,cle);
					confo_back->Recon_Back_Main(mol,cle,moln,mbb);      //given CA, recon BackBone (N,CA,C,O,CB)
					confo_beta->Recon_Beta_21(mol,mcb,moln,ami,cle);    //given CA, recon CB
					//[2]assign
					for(i=0;i<moln;i++)
					{
						//CA missing!!
						if(pdb[i].get_backbone_part_index(1)==0)
						{
							//backbone (N,CA,C,O)
							for(j=0;j<4;j++)pdb[i].set_backbone_atom(j,mbb[i][j]);
							//sidechain (CB)
							if(i==0||i==moln-1)pdb[i].set_sidechain_atom(0,mbb[i][4]);
							else pdb[i].set_sidechain_atom(0,mcb[i]);
						}
						else
						{
							//backbone (N,CA,C,O)
							for(j=0;j<4;j++)
							{
								if(pdb[i].get_backbone_part_index(j)==0)
								{
									pdb[i].set_backbone_atom(j,mbb[i][j]);
								}
							}
							//sidechain (CB)
							if(pdb[i].get_sidechain_part_index(0)==0)
							{
								if(i==0||i==moln-1)pdb[i].set_sidechain_atom(0,mbb[i][4]);
								else pdb[i].set_sidechain_atom(0,mcb[i]);
							}
						}
					}
				}
			}

			//get CB
			for(i=0;i<moln;i++)pdb[i].get_sidechain_atom( "CB ",mcb[i] );

			//get missing map
			vector <int> mapping;
			mapping.resize(moln);
			for(i=0;i<moln;i++)mapping[i]=1;
			string miss_seq="";

			//output
			FILE *fpdb;
			file=outa+"/"+output+".pdb";
			fpdb=fopen(file.c_str(),"wb");
			if(fpdb==0)
			{
				fprintf(stderr,"ERROR: file %s can't be opened. \n",file.c_str());
			}
			else
			{
				mol_output.Output_PDB_III(fpdb,moln,pdb,c,OutType,OutMode,OutGlys);
				fclose(fpdb);
			}

			//output others
			if(OutFifi!=0)
			{
				//-> basic output
				if(OutFifi==1 || OutFifi==3 || OutFifi==5 || OutFifi==7)
					Output_Protein_Features_AMI_SSE_CLE(outa,output,moln,pdb,mapping,miss_seq);
				if(OutFifi==2 || OutFifi==3 || OutFifi==6 || OutFifi==7)
					Output_Protein_Features_ACC(outa,output,acc_surface,moln,pdb,mcc,mcc_side,ami,acc,mapping,miss_seq);
				if(OutFifi==4 || OutFifi==5 || OutFifi==6 || OutFifi==7)
					Output_Protein_Features(outa,output,acc_surface,mol,mcb,moln,pdb,mcc,mcc_side,ami,acc,mapping,miss_seq);
				//-> additionals
				//--| angles
				if(OutFifi==8)
					Output_Protein_Angles(outa,output,moln,pdb,ami,mapping,miss_seq);
				//--| distances
				if(OutFifi==9)
					Output_Protein_Distances(outa,output,hydro_bond,moln,pdb,mcc,mcc_side,ami,mapping,miss_seq);
			}
		}

		//printf
		count++;
	}
	//terminal
	delete [] pdb;
	delete [] ami;
	delete [] cle;
	delete [] acc;
	delete [] mol;
	delete [] mcb;
	DeleteArray2D(&mbb,totlen);
	DeleteArray2D(&mcc,totlen);
	delete [] mcc_side;
	delete confo_beta;
	delete confo_back;
	delete acc_surface;
	delete hydro_bond;
	fin.close();
	fin.clear();
	return wwscount;
}

//===================== single_process ===================//__110710__//
//[single_style]
void PDB_Back_Process_Single(string &input,string &range,string &output,
	int OutType,int OutMode,int OutGlys,int OutNoca,int OutReco,int OutFifi,int OutLogf,int OutWarn,
	string &input_miss)
{
	//get length
	int ret_val;
	int totlen=Get_PDB_File_Len(input);
	if(totlen<=0)
	{
		fprintf(stderr,"pdbfile %s length error!!\n",input.c_str());
		exit(-1);
	}
	//get seqres_miss
	string miss_seq="";
	if(input_miss!="")
	{
		int miss_len=Read_FASTA_SEQRES(input_miss,miss_seq);
		if(miss_len<=0)
		{
			fprintf(stderr,"input_miss %s length error!!\n",input_miss.c_str());
			exit(-1);
		}
	}
	//class
	Mol_File mol_input;
	Mol_Out mol_output;
	Confo_Lett confo_lett;
	Confo_Beta *confo_beta=new Confo_Beta(totlen);
	Confo_Back *confo_back=new Confo_Back(totlen);
	Acc_Surface *acc_surface=new Acc_Surface(totlen);
	Hydro_Bond *hydro_bond=new Hydro_Bond(totlen);
	mol_input.MODRES=1;
	//init
	int moln;
	PDB_Residue *pdb=new PDB_Residue[totlen];
	string TER="TER                                                                             ";

	//recon_related
	XYZ *mol=new XYZ[totlen];   //CA
	XYZ *mcb=new XYZ[totlen];   //CB
	XYZ **mbb;                  //BackBone (N,CA,C,O,CB)
	XYZ **mcc;                  //BackBone (N,CA,C,O,CB,...,)
	NewArray2D(&mbb,totlen,5);
	NewArray2D(&mcc,totlen,15);
	int *mcc_side=new int[totlen];
	char *ami=new char[totlen+1];
	char *cle=new char[totlen+1];
	char *acc=new char[totlen+1];

	//open
	mol_input.CbBACK=1;
	if(OutNoca==1)mol_input.CaONLY=0; //consider Non-CA atoms !!//__110408__//
	//macro
	mol_input.OUTP_MODE=OutType;
	mol_input.PROC_MODE=OutMode;
	mol_input.GLYC_MODE=OutGlys;
	//memory limit
	mol_input.WARNING_out=OutWarn;
	mol_input.MEMORY_LIMIT=totlen;

	//process
	{
		//pre_process
		ret_val=mol_input.XYZ_Input(input,range,0,moln,0,0,0,0,0);
		if(ret_val<=0)
		{
			if(ret_val!=-12345)goto end;
		}
		if(ret_val!=1)goto end;
		//check memory
		if(ret_val==-12345)
		{
			//add memory
			if(moln>totlen)
			{
				totlen=moln;
				//delete
				delete [] pdb;
				delete [] mol;
				delete [] mcb;
				delete [] ami;
				delete [] cle;
				delete [] acc;
				DeleteArray2D(&mbb,totlen);
				DeleteArray2D(&mcc,totlen);
				delete [] mcc_side;
				delete confo_beta;
				delete confo_back;
				delete acc_surface;
				delete hydro_bond;
				//create
				pdb=new PDB_Residue[totlen];
				mol=new XYZ[totlen];
				mcb=new XYZ[totlen];
				ami=new char[totlen+1];
				cle=new char[totlen+1];
				acc=new char[totlen+1];
				NewArray2D(&mbb,totlen,5);
				NewArray2D(&mcc,totlen,15);
				mcc_side=new int[totlen];
				confo_beta=new Confo_Beta(totlen);
				confo_back=new Confo_Back(totlen);
				acc_surface=new Acc_Surface(totlen);
				hydro_bond=new Hydro_Bond(totlen);
				//memory limit
				mol_input.MEMORY_LIMIT=totlen;
			}
		}
		//reload
		int PRE_LOAD_=mol_input.PRE_LOAD;
		int WARNING_out_=mol_input.WARNING_out;
		mol_input.PRE_LOAD=1;
		mol_input.WARNING_out=0;
		mol_input.XYZ_Input(input,range,0,moln,mol,ami,0,0,pdb);
		mol_input.PRE_LOAD=PRE_LOAD_;
		mol_input.WARNING_out=WARNING_out_;

		//check
		{
			int i,j;
			int correct;
			int iret;
			correct=1; //default:OK
			for(i=0;i<moln;i++)
			{
				iret=pdb[i].PDB_residue_backbone_check(4);
				if(iret!=1)
				{
					correct=0;
					break;
				}
				iret=pdb[i].PDB_residue_CB_check();
				if(iret!=1)
				{
					correct=0;
					break;
				}
			}

			//judge
			if(correct!=1)
			{
				if(OutReco==1)
				{
					//[1]recon
					confo_lett.btb_ori(0,0,0,moln,mol,cle);
					confo_back->Recon_Back_Main(mol,cle,moln,mbb);      //given CA, recon BackBone (N,CA,C,O,CB)
					confo_beta->Recon_Beta_21(mol,mcb,moln,ami,cle);    //given CA, recon CB
					//[2]assign
					for(i=0;i<moln;i++)
					{
						//CA missing!!
						if(pdb[i].get_backbone_part_index(1)==0)
						{
							//backbone (N,CA,C,O)
							for(j=0;j<4;j++)pdb[i].set_backbone_atom(j,mbb[i][j]);
							//sidechain (CB)
							if(i==0||i==moln-1)pdb[i].set_sidechain_atom(0,mbb[i][4]);
							else pdb[i].set_sidechain_atom(0,mcb[i]);
						}
						else
						{
							//backbone (N,CA,C,O)
							for(j=0;j<4;j++)
							{
								if(pdb[i].get_backbone_part_index(j)==0)
								{
									pdb[i].set_backbone_atom(j,mbb[i][j]);
								}
							}
							//sidechain (CB)
							if(pdb[i].get_sidechain_part_index(0)==0)
							{
								if(i==0||i==moln-1)pdb[i].set_sidechain_atom(0,mbb[i][4]);
								else pdb[i].set_sidechain_atom(0,mcb[i]);
							}
						}
					}
				}
			}

			//get CB
			for(i=0;i<moln;i++) pdb[i].get_sidechain_atom( "CB ",mcb[i] );

			//get missing map
			vector <int> mapping;
			if(miss_seq=="")
			{
				mapping.resize(moln);
				for(i=0;i<moln;i++)mapping[i]=1;
			}
			else
			{
				string AMI(ami);
				Seqres_DynaProg(miss_seq,AMI,mapping);
			}

			//output
			FILE *fpdb;
			fpdb=fopen(output.c_str(),"wb");
			if(fpdb==0)
			{
				fprintf(stderr,"ERROR: file %s can't be opened. \n",output.c_str());
			}
			else
			{
				mol_output.Output_PDB_III(fpdb,moln,pdb,'_',OutType,OutMode,OutGlys);
				fclose(fpdb);
			}

			//output others
			if(OutFifi!=0)
			{
				//-> get name and root
				string outroot,outname;
				getBaseName(output,outname,'/','.');
				getRootName(output,outroot,'/');
				//-> basic output
				if(OutFifi==1 || OutFifi==3 || OutFifi==5 || OutFifi==7)
					Output_Protein_Features_AMI_SSE_CLE(outroot,outname,moln,pdb,mapping,miss_seq);
				if(OutFifi==2 || OutFifi==3 || OutFifi==6 || OutFifi==7)
					Output_Protein_Features_ACC(outroot,outname,acc_surface,moln,pdb,mcc,mcc_side,ami,acc,mapping,miss_seq);
				if(OutFifi==4 || OutFifi==5 || OutFifi==6 || OutFifi==7)
					Output_Protein_Features(outroot,outname,acc_surface,mol,mcb,moln,pdb,mcc,mcc_side,ami,acc,mapping,miss_seq);
				//-> additionals
				//--| angles
				if(OutFifi==8)
					Output_Protein_Angles(outroot,outname,moln,pdb,ami,mapping,miss_seq);
				//--| distances
				if(OutFifi==9)
					Output_Protein_Distances(outroot,outname,hydro_bond,moln,pdb,mcc,mcc_side,ami,mapping,miss_seq);
			}
		}
	}
end:
	//terminal
	delete [] pdb;
	delete [] ami;
	delete [] cle;
	delete [] acc;
	delete [] mol;
	delete [] mcb;
	DeleteArray2D(&mbb,totlen);
	DeleteArray2D(&mcc,totlen);
	delete [] mcc_side;
	delete confo_beta;
	delete confo_back;
	delete acc_surface;
	delete hydro_bond;
}

//============== main ===============//
int main(int argc, char** argv)
{
	//---- BackBone ----//(process list or single) PDB_Backbone//
	{
		process_args(argc,argv);
		//-> if set with -F, then automatically set -R
		if(INPUT_FIFI!=0)INPUT_RECO=1;
		//-> list or single
		if(LIST_OR_SINGLE>0)PDB_Back_Process(INPUT_NAM,INPUT_LIST,INPUT_OUT,INPUT_TYPE,INPUT_MODE,INPUT_GLYS,INPUT_NOCA,INPUT_RECO,INPUT_FIFI,INPUT_LOGF,WARN_OUT);
		else PDB_Back_Process_Single(INPUT_NAM,INPUT_RAN,INPUT_OUT,INPUT_TYPE,INPUT_MODE,INPUT_GLYS,INPUT_NOCA,INPUT_RECO,INPUT_FIFI,INPUT_LOGF,WARN_OUT,INPUT_MISS);
		exit(0);
	}
}
