#include <iostream>
#include <stdio.h>
#include <cmath>
#include <vector>
#include "mgmres.hpp"

using namespace std;
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//							Definição da classe célula
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class Cell{
public:
  	double x,y,rN,r,rA,uptake;
  	int time,state,prev_state;
  	void set(double X,double Y,double RN,double R,double RA,double Uptake, int Time,int State){x=X;y=Y;rN=RN;r=R;rA=RA;uptake = Uptake;time=Time;state=State;}
  	void print();
};

void diference(vector<Cell>* Celulas, double X2[], double D, double h, double h_f, int N, int M, int iteracao, double EGF_backg, double conc_EGF_front);
void reaction(vector<Cell>* Celulas, int N,int M,double h,double h_f,double EGF[],double EGF_backg);

int main ( int argc, char *argv[] ){  	
  	//Parâmetros do domínio
  	double D=1e-12	, h=10.0, h_f=1.0, N=340.0, M=340.0, EGF_backg = 0.000, conc_EGF_front = 0.0/*0.263717*/;
  	int n = (N/h)+1, m = (M/h)+1;
  	double X2[n*m];
  	
  	//Parâmetros dos agentes
  	double rN = 5.295,r = 9.953, rA = 1.214*r,EGF_cell = 2.65;
  	vector<Cell> Celulas;
  	int time=0;

  	//Condição inicial
  	Cell a; 
  	a.set(34.0,34.0,rN,r,rA,EGF_cell,time,2);
  	Celulas.push_back(a);  
  	a.set(34.0,306.0,rN,r,rA,EGF_cell,time,1);
  	Celulas.push_back(a);  
  	a.set(306.0,306.0,rN,r,rA,EGF_cell,time,2);
  	Celulas.push_back(a);  
  	a.set(306.0,34.0,rN,r,rA,EGF_cell,time,1);
  	Celulas.push_back(a);  
  	/*
  	a.set(220.0,260.0,rN,r,rA,EGF_cell,time,2);
  	Celulas.push_back(a);  
  	a.set(190.0,190.0,rN,r,rA,EGF_cell,time,1);
  	Celulas.push_back(a);  
  	a.set(200.0,200.0,rN,r,rA,EGF_cell,time,2);
  	Celulas.push_back(a);  
  	a.set(210.0,210.0,rN,r,rA,EGF_cell,time,1);
  	Celulas.push_back(a); */ 
  	  
  	diference(&Celulas, X2, D, h, h_f, N, M,time,EGF_backg,conc_EGF_front);

  	return 0;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//					Cálculo da Concentração de EGF no microambiente
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void diferenceEGF(const vector<Cell>& Celulas, double X2[], double D, double h, double h_f, double N, double M,int iteracao,double EGF_backg,double conc_EGF_front){
	double tol = 1.0e-6;
	double delta_t = 0.1;
	double temp = 0;
	int n = (N/h)+1;
	int m = (M/h)+1;
	double sigma = delta_t/pow(h,2);
	double * B = new double[n*m];
	double * X1 = new double[n*m]; 
	int nz_num = 4*3 + 2*(n-2)*4 + 2*(m-2)*4 + (n-2)*(m-2)*5;
	int row_ptr[n*m+1],col_ind[nz_num];
	double val[nz_num];
	row_ptr[0]=0;
	
	double * EGF = new double[n*m]; 
	
	reaction(Celulas, N,M, h, h_f, EGF,EGF_backg);	
	
	//Define as coordenadas
	double XX[n],YY[m];
	for(int i=0;i<n;i++) XX[i]= (i*h);
	for(int i=0;i<m;i++) YY[i]= (i*h);
	
	//Monta a matriz A em CRS
	int k=0;
  	for (int i = 0; i < n*m; i++){ 
       		for (int j = 0; j < n*m; j++){ 
            		if (i==j){
            			if (row_ptr[i] != k) row_ptr[i+1] = k+1;
            			col_ind[k] = j;
				val[k]= 1+2*D*sigma;
				k++;
				continue;
			}
	  		if ( (i%n!=0 && i-j==1) || (i%n!=n-1 && j-i==1) ){
	  			if (row_ptr[i] != k) row_ptr[i+1] = k+1;
            			col_ind[k] = j;
				val[k]= -D*sigma*0.5;
				k++;
				continue;
	  		}
	  		if (i-j==n || j-i==n){
	  			if (row_ptr[i] != k) row_ptr[i+1] = k+1;
            			col_ind[k] = j;
				val[k]= -D*sigma*0.5;
				k++;
	  		}
       		}
  	}
        
    //Condição inicial
	for(int i=0;i<n*m;i++){
		X1[i] = conc_EGF_front;		
		X2[i] = conc_EGF_front;
	}
	//Condição de Contorno
	for(int i=0;i<n*m;i++){
		if (i%n == 0 || i%n == n-1 || i/n == 0 || i/n==m-1){
			for(int j=row_ptr[i];j<row_ptr[i+1];j++){
				if (col_ind[j] == i) val[j]=1.0;
				else val[j]=0.0;
			}
		}
		
	}
	
	//Parametros para mgmres
	int itr_max=1000, mr=5;
	double tol_abs=1.0e-5, tol_rel=1.0e-5;
	
	//Loop no tempo
	while (temp< 10.0){        	
        //Atualiza B[ ]
 	  	for(int i=0;i<n*m;i++){
			B[i] = (1-2*D*sigma)*X2[i]+EGF[i];
			if (i%n!=0)
				B[i]+= (D*sigma/2)*X2[i-1];
			if (i%n!=n-1)
				B[i]+= (D*sigma/2)*X2[i+1];
			if ((i-n)>=0)
				B[i]+= (D*sigma/2)*X2[i-n];
			if ((i+n)<n*m)
				B[i]+= (D*sigma/2)*X2[i+n];
	  	}
		//Condição de Contorno
	  	for(int i=0;i<n*m;i++){
			if (i%n == 0 || i%n == n-1 || i/n == 0 || i/n==m-1){
				B[i]=conc_EGF_front;
			}
	  	}
		//Calcula X[]
		pmgmres_ilu_cr ( n*m, nz_num, row_ptr, col_ind, val, X2, B, itr_max, mr, tol_abs, tol_rel );
		
		for(int i=0;i<n*m;i++) X1[i]= X2[i];
		printf("TEMPO = %f \n",temp);
		
		//Atualiza o tempo
 		temp += delta_t ;
	}
	
	//Cria arquivo egf*.dat
	FILE *arq;
	char name[100];
  	sprintf(name,"egf-%05d.dat",iteracao);
	arq = fopen(name,"w");
	fprintf(arq,"# %d %d \n",N,M);
	fprintf(arq,"# %d %d \n",n,m);
	fprintf(arq,"# %ld %d \n",(*Celulas).size(),iteracao);
	for(int i=0;i<n*m;i++){
		fprintf(arq,"%lf %lf %lf \n",XX[i%n],YY[i/n],X2[i]);
	}
	fclose(arq);
	
	//Libera memória
	delete [] B;
	delete [] X1;
	delete [] EGF;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//				Cálculo do termo reativo para equação de concentração de nutriente
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void reaction(vector<Cell>* Celulas, int N,int M,double h,double h_f,double EGF[], double EGF_backg){
	//Domínio N X M com malha refinada de espaçamento h_f
	int m = M/h_f+1;
	int n = N/h_f+1;
	double uptake_malhaf[m][n];
	
	//Número de células
	int num_cell = (*Celulas).size();
	
	//Define malha fina
	double XX_f[n],YY_f[m];
	for(int i=0;i< n;i++) XX_f[i]= i*h_f;
	for(int i=0;i< m;i++) YY_f[i]= i*h_f;
	
	
	//Inicializa termo de captação
	for(int i=0;i< m;i++){
		for(int j=0;j< n;j++){
			uptake_malhaf[i][j] = 0.0;
		}
	}
	
	//Calcula termo de captação 
	for(int k=0;k<num_cell;k++){
		for(int i=0;i< m;i++){
			for(int j=0;j< n;j++){
				if ( pow(pow(XX_f[j]-(*Celulas)[k].x,2) + pow(YY_f[i]-(*Celulas)[k].y,2),0.5) <= (*Celulas)[k].r){
					uptake_malhaf[i][j] += (*Celulas)[k].uptake;
				}
			}
		}
	}
	
	//Termo de captação background
	for(int i=0;i< m;i++){
		for(int j=0;j< n;j++){
			if (uptake_malhaf[i][j] == 0.0)
				uptake_malhaf[i][j] = EGF_backg;
		}
	}
	
	//Domínio N X M com malha refinada de espaçamento h
	m = M/h+1;
	n = N/h+1;
	
	//Calcula o termo reação na malha do nutriente	
	for(int i=0;i<n*m;i++){
		double EGF_temp=0.0;
		//Nós fora da fronteira
		if (i%n != 0 && i%n != n-1 && i/n !=0 && i/n != m-1){
			for(int j=(i%n)*h - h/2 ;j<= (i%n)*h + h/2;j++){
				for(int l=(i/n)*h - h/2;l<= (i/n)*h + h/2;l++){
					EGF_temp += uptake_malhaf[l][j];
				}
			}
			EGF[i] = EGF_temp/pow(h+1,2);
			continue;
		}
		//Nó (0,0)
		if (i%n==0 && i/n ==0){
			for(int j= 0;j<= h/2;j++){
				for(int l= 0;l<= h/2;l++){
					EGF_temp += uptake_malhaf[l][j];
				}
			}
			EGF[i] = EGF_temp/pow(h/2+1,2);
			continue;
		}
		//Nó (N,0)
		if (i%n==n-1 && i/n ==0){
			for(int j= (i%n)*h - h/2;j<= (i%n)*h;j++){
				for(int l= 0;l<= h/2;l++){
					EGF_temp += uptake_malhaf[l][j];
				}
			}
			EGF[i] = EGF_temp/pow(h/2+1,2);
			continue;
		}
		//Nó (0,M)
		if (i%n==0 && i/n ==m-1){
			for(int j= 0;j<= h/2;j++){
				for(int l= (i/n)*h - h/2;l<= (i/n)*h;l++){
					EGF_temp += uptake_malhaf[l][j];
				}
			}
			EGF[i] = EGF_temp/pow(h/2+1,2);
			continue;
		}
		//Nó (N,M)
		if (i%n==n-1 && i/n ==m-1){
			for(int j= (i%n)*h - h/2;j<= (i%n)*h;j++){
				for(int l= (i/n)*h - h/2;l<= (i/n)*h;l++){
					EGF_temp += uptake_malhaf[l][j];
				}
			}
			EGF[i] = EGF_temp/pow(h/2+1,2);
			continue;
		}
		//Nós (0,X)
		if (i%n==0 && i/n != m-1 && i/n != 0){
			for(int j= 0;j<= h/2;j++){
				for(int l= (i/n)*h - h/2;l<= (i/n)*h + h/2;l++){
					EGF_temp += uptake_malhaf[l][j];
				}
			}
			EGF[i] = EGF_temp/((h/2+1)*(h+1));
			continue;
		}
		//Nós (X,0)
		if (i/n==0 && i%n != 0 && i%n != n-1){
			for(int j= (i%n)*h - h/2;j<= (i%n)*h + h/2;j++){
				for(int l= 0;l<= h/2;l++){
					EGF_temp += uptake_malhaf[l][j];
				}
			}
			EGF[i] = EGF_temp/((h/2+1)*(h+1));
			continue;
		}
		//Nós (N,X)
		if (i%n==n-1 && i/n != m-1 && i/n != 0){
			for(int j= (i%n)*h - h/2;j<= (i%n)*h;j++){
				for(int l= (i/n)*h - h/2;l<= (i/n)*h + h/2;l++){
					EGF_temp += uptake_malhaf[l][j];
				}
			}
			EGF[i] = EGF_temp/((h/2+1)*(h+1));
			continue;
		}
		//Nós (X,M)
		if (i/n==m-1 && i%n != n-1 && i%n != 0){
			for(int j= (i%n)*h - h/2;j<= (i%n)*h + h/2;j++){
				for(int l= (i/n)*h - h/2;l<= (i/n)*h;l++){
					EGF_temp += uptake_malhaf[l][j];
				}
			}
			EGF[i] = EGF_temp/((h/2+1)*(h+1));
			continue;
		}

	}
}
