//Program ranwalk1.cpp

//Made:5/12/09;Modified:

//This program simulates the random walk in 3D in a geometry consisting of parallel cylinders

  #include <iostream>
  #include <fstream>
  #include <stdlib.h>
  #include <math.h>
  #include <iomanip>
  #include <time.h>
  #include <cstdlib> 
  #include <algorithm>
  #include <string>	
  #include <complex>
  #include <string>



using namespace std;

#define Pi 3.14159265

/***********************RANDOM NUMBER GENERATOR (KISS) BEGIN****************************/
//RNG JKISS from http://www0.cs.ucl.ac.uk/staff/d.jones/GoodPracticeRNG.pdf

//Define structure used in RNG
struct RNGvalues
{
	unsigned int x;	   //seed
	unsigned int y;    //seed
	unsigned int z;    //seed
	unsigned int c;    //seed
	double r;   //random number output from RNG
};


//Seed initialization using entropy pool in /dev/urandom.  This should only be called once.
unsigned int devrand()
{
	unsigned int seed;
	FILE *urandom;

	urandom = fopen ("/dev/random", "r");
	fread (&seed, sizeof (seed), 1, urandom);
	fclose(urandom);
	cout<<"seed is " << seed<<endl;

	return seed;
}

//Initialize RNG
RNGvalues init_KISS()
{
	
	RNGvalues v;

	unsigned int x,y,z,c;
	unsigned long long t;

	x=devrand();
	while (!(y=devrand())); /* y must not be zero! */ z = devrand();
	/* We don’t really need to set c as well but let's anyway... */
	/* NOTE: offset c by 1 to avoid z=c=0 */
	c = devrand() % 698769068 + 1; /* Should be less than 698769069 */
	
	v.x = x;
	v.y = y;
	v.z = z;
	v.c = c;
	return v;
}

//RNG script

RNGvalues JKISS(RNGvalues v){

unsigned long long t;

v.x = 314527869 * v.x + 1234567;
v.y ^= v.y << 5; v.y ^= v.y >> 7; v.y ^= v.y << 22;
t = 4294584393ULL * v.z + v.c; v.c = t >> 32; v.z = t;

v.r = (v.x+v.y+v.z)/ 4294967295.0;
//return x + y + z;
return v;
}

/***********************RANDOM NUMBER GENERATOR (KISS) END ******************************/



double diffclock(clock_t clock1,clock_t clock2)
{
	double diffticks=clock1-clock2;
	double diffms=(diffticks*10)/CLOCKS_PER_SEC;
	return diffms;
}

double max ( double a, double b ) {
  return (b<a)?a:b;     // or: return comp(b,a)?a:b; for the comp version
}
double min ( double a, double b ) {
  return (b>a)?a:b;     // or: return comp(b,a)?a:b; for the comp version
}
/*double abs ( double a) {
  return sqrt(a*a)  ;   
}*/
// Simulation constants
// defining resolution
#define dt 7.5e5  // time step in ps
const int TN=16e5;  //number of time steps
const double TT=TN*dt;  //Total time
//const int TN=100*D*1e12/(4*P*P*dt);			// total time steps
const long NP =100;		// number of points
const int timepoints=10000;  //number of time points to record
const int Tstep=TN/timepoints;


#define Din 0.5e-9					// Diffusion coefficient inside the cylinders in m^2/s
#define Dout 2e-9					// Diffusion coefficient outside the cylinders in m^2/
#define P 0					// permeability in m/s
const double stepIN=sqrt(6*dt*Din);	//in micrometers
const double stepOUT=sqrt(6*dt*Dout);//in micrometers
double pIE=stepIN*1e-6*P/Din; // may be wrong
double pEI=stepOUT*1e-6*P/Dout; // may be wrong

const int sizetheta=8;

	  

  int main(int argc, char *argv[] ){
 	// Declaration of the indices
	int i,j,ii,jj,ti,tj,k,jk,kl,t1,t2;
	unsigned long long t;
	clock_t begin=clock();
	clock_t end=clock();
	
	RNGvalues v;
	v = init_KISS();  //Initialize the RNG

	cout<<"START"<<endl;
	// generate or load phantom geometry
	cout<<"Generating phantom:"<<endl;
	begin=clock();
	double radius = 5.0;		//in micrometers
	double density = 0.5;
	double res = 0.1;			// in micrometers
	double FOV = 100;
    double *xm,*ym,*r, Xmax, Ymax;
	int Nmax;
	unsigned int **A;
	char geometry='I';
	if (geometry=='I'){//import geometry from file
		//double density;
		ifstream myfile1 ("phantom_density.out", ios::in);
		myfile1>>density;
		myfile1.close();
		//double radius;
		ifstream myfile2 ("phantom_radius.out", ios::in);
		myfile2>>radius;
		myfile2.close();
		//double res;
		ifstream myfile3 ("phantom_res.out", ios::in);
		myfile3>>res;
		myfile3.close();
		//double Xmax;
		ifstream myfile4 ("phantom_Xmax.out", ios::in);
		myfile4>>Xmax;
		myfile4.close();
		//double Ymax;
		ifstream myfile5 ("phantom_Ymax.out", ios::in);
		myfile5>>Ymax;
		myfile5.close();
		int Ncyl;
		ifstream myfile6 ("phantom_N.out", ios::in);
		myfile6>>Ncyl;
		myfile6.close();
		cout<<Ncyl<<" cylinders (radius "<<radius<<" micron) in a geometry with size "<< Xmax<<" micron ..."<<endl;
		ifstream myfile11 ("phantom_Nmax.out", ios::in);
		myfile11>>Nmax;
		myfile11.close();
		xm=new double [Ncyl+1];
		ym=new double [Ncyl+1];
		r=new double [Ncyl+1];
		ifstream myfile7 ("phantom_r.out", ios::in);
		i=1;
		while (!myfile7.eof()&i<=Ncyl){
			myfile7>>r[i];
			//cout<<" "<<i<<" "<<r[i];
			i++;
			}
		ifstream myfile8 ("phantom_xm.out", ios::in);
		i=1;
		while (!myfile8.eof()&i<=Ncyl){
			myfile8>>xm[i];
			//cout<<" "<<i<<" "<<xm[i]<<endl;
			i++;
			}
		ifstream myfile9 ("phantom_ym.out", ios::in);
		i=1;
		while (!myfile9.eof()&i<=Ncyl){
			myfile9>>ym[i];
			//cout<<" "<<i<<" "<<ym[i]<<endl;
			i++;
			}
		int MSX=ceil(Xmax/res)+1;
		int MSY=ceil(Ymax/res)+1;
		A = new unsigned int *[MSX];
		for (i = 0; i < MSX; ++i){
		A[i]=new unsigned int[MSY];
			for (j = 0;j<MSY;j++)
				A[i][j]=0;
			}
		ifstream myfile10 ("phantom_A.out", ios::in);
		//while (!myfile10.eof()){
			for (i=0;i<MSY;i++){
				for (j=0;j<MSX;j++){
				myfile10>>A[j][i];
				//if (A[i][j]!=0){
				//	cout<<" "<<A[i][j]<<" ";
				//	}
				}
				}
		/*ofstream Amatrix ("phantom_Aoutput.out");
		//while (!myfile10.eof()){
			for (i=0;i<MSY;i++){
				for (j=0;j<MSX;j++){
				Amatrix<<A[j][i]<<" ";
				//if (A[i][j]!=0){
				//	cout<<" "<<A[i][j]<<" ";
				//	}
				//}
				}
		}*/
	}
//	geometry=='S';
	/*if (geometry=='S'){
		//const int Ncyl=2000;
		//double FOV=sqrt(Ncyl*pi*radius*radius/density);
		int Ncyl=density/Pi*(FOV/radius)*(FOV/radius);
		double Pmax=Pi/4;
		double Rmax=sqrt(Pmax/Pi/Ncyl)*FOV;
		int Nx=(FOV/(2*Rmax));
		int Ny=Nx;
		Xmax=ceil(Nx*2*Rmax/res)*res;
		Ymax=ceil(Ny*2*Rmax/res)*res;
		Ncyl=Nx*Ny;
		cout<<Ncyl<<" cylinders (radius "<<radius<<" micron) in a square geometry with size "<< Xmax<<" micron ..."<<endl;
		FOV=Xmax;
		xm=new double [Ncyl];
		ym=new double [Ncyl];
		r=new double [Ncyl];
		int teller=0;
		for (int k=0;k<(Nx-1);k++){
			for (int l=0;l<(Ny-1);l++){
				r[teller]=radius;
				xm[teller]=(1+2*k)*Rmax;
				ym[teller]=(1+2*l)*Rmax;
				teller++;
				//cout<<"xm = "<<xm[teller-1]<<"  ";
				}
			}
		int MSX=ceil(Xmax/res)+1;
		int MSY=ceil(Ymax/res)+1;
		A = new int *[MSX];
		for (i = 0; i < MSX; ++i){
		A[i]=new int[MSY];
			for (j = 0;j<MSY;j++)
				A[i][j]=0;
			}
		for (i=0;i<teller-2;i++){
			for (ii=floor((xm[i]-r[i])/res);ii<=ceil((xm[i]+r[i])/res);ii++){
				if(ii>MSX-1)
					ti=ii-MSX;
				else if (ii<0)
					ti=ii+MSX;
				else
					ti=ii;
				for (jj=floor((ym[i]-r[i])/res);jj<=ceil((ym[i]+r[i])/res);jj++){
					if(jj>MSY-1)
						tj=jj-MSY;
					else
						tj=jj;
					if(((ii*res-xm[i])*(ii*res-xm[i])+(jj*res-ym[i])*(jj*res-ym[i]))<=(r[i]+res)*(r[i]+res))		
						A[ti][tj]=i;
					else 
						A[ti][tj]=0;
			//cout<<A[ti][tj]<<endl;
					}
				}
			}
		}*/
	end=clock();
	cout << "Done! Elpased time "<<double(diffclock(end,begin))/10 << " s"<< endl;
	
	// start simulation
	cout << "Start of simulation " << endl;
	cout << "Number of tracer points " <<  NP << endl;
	cout << "Time points " <<  TN << endl;
//	cout << "pIE " <<  pIE << endl;
//	cout << "pEI " <<  pEI << endl;
//	cout << "stepIN [mum] " <<  stepIN << endl;
	cout << "stepOUT [mum]" <<  stepOUT << endl;

	// generating seed points
	cout<<"Tracer points...";
	cout << "Nmax:  "<<Nmax<<endl;
	begin=clock();
	double xip[NP],yip[NP];
	double xp,yp,xpr,ypr,tru,stap,vs,vt,sinvs,cosvs,sinvt,xpp,ypp,zpp,xnp,ynp,znp,xB,yB,zB,xC,yC,zC,X,Y,R,a,b,c,d,s,sn,srest,amin,aplus,xmin,ymin,pex;
	double CminX,XminC,CminY,YminC,cosphi,sinphi,cos2phi,sin2phi,xD,yD,zD,cosvss,sinvss;
	double vp[3];
	int xfl,yfl,xce,yce;
	bool intrusion, intrusion1, intrusion2;
	//int t;
	k=0;
		while(k<NP){
		//xp=Xmax/3+random_uniform_0_1()*Xmax/3;
		//yp=Ymax/3+random_uniform_0_1()*Ymax/3;
		//cout << "random seed 2"<< random_uniform_0_1()<< endl;
		//cout<<"k: "<<k<<endl;
		v=JKISS(v);
		xp=v.r*Xmax;
		//cout<<"randx is "<< v.r <<endl;
		//cout<<"idum is "<< idum << endl;
		//cout<<"k is " << k << endl;
		//cout<<"xp is " << xp << endl;
		v = JKISS(v);
		yp=v.r*Ymax;
		xpr=xp/res;
		//cout<<"randy is "<<v.r<< endl;
		ypr=yp/res;
		// choose to take only outside points SHOULD ALWAYS BE THE CASE
		xfl=floor(xpr);yfl=floor(ypr);xce=ceil(xpr);yce=ceil(ypr);
		/*cout<< "xfl is " << xfl << endl;
		cout<< "yfl is " << yfl << endl;
		cout << "xce is " << xce << endl;
		cout << "yce is " << yce <<endl;
		cout << "A[xfl][yfl]:  " << A[xfl][yfl] <<endl;
		cout << "A[xce][yce]:  " << A[xce][yce] << endl;
		cout << "max A floor/ceil:  " << max(A[xfl][yfl],A[xce][yce]) << endl;
		cout << "A[xfl][yce]:  " << A[xfl][yce] << endl;
		cout << "A[xce][yfl]:  " << A[xce][yfl] << endl;
		cout << "max A floor and ceil:  " << max(A[xfl][yce],A[xce][yfl])<< endl;
		*/
		//cout << "max all " <<max(max(A[xfl][yfl],A[xce][yce]),max(A[xfl][yce],A[xce][yfl]))<<endl;
		t=max(max(A[xfl][yfl],A[xce][yce]),max(A[xfl][yce],A[xce][yfl]));
		//cout << "t is "<< t << endl;

		if (t<Nmax)
			intrusion=(t!=0)&&((min(abs(xp-xm[t]),min(abs(xp-xm[t]-Xmax),\
			(xp-xm[t]+Xmax))))*(min(abs(xp-xm[t]),min(abs(xp-xm[t]-Xmax),\
			(xp-xm[t]+Xmax))))+(min(abs(yp-ym[t]),min(abs(yp-ym[t]-Ymax),\
			(yp-ym[t]+Ymax))))*(min(abs(yp-ym[t]),min(abs(yp-ym[t]-Ymax),\
			(yp-ym[t]+Ymax))))<r[t]*r[t]);
		else{
			cout<<"270: t>Nmax: "<<t<<endl;
			t1=floor(t/Nmax);
			//cout<<"t1 is " << t1<<endl;
			t2=t-t1*Nmax;
			
			intrusion1=(t1!=0)&&((min(abs(xp-xm[t1]),min(abs(xp-xm[t1]-Xmax),\
			(xp-xm[t1]+Xmax))))*(min(abs(xp-xm[t1]),min(abs(xp-xm[t1]-Xmax),\
			(xp-xm[t1]+Xmax))))+(min(abs(yp-ym[t1]),min(abs(yp-ym[t1]-Ymax),\
			(yp-ym[t1]+Ymax))))*(min(abs(yp-ym[t1]),min(abs(yp-ym[t1]-Ymax),\
			(yp-ym[t1]+Ymax))))<r[t1]*r[t1]);
			//cout << "intrusion1 is " << intrusion1 << endl;
			//cout << "t2 is " << t2 << endl;
			//cout << "xm[t2] is "<< xm[t2] << endl;
			//cout << "ym [t2] is " << ym[t2] << endl;
			intrusion2=(t2!=0)&&((min(abs(xp-xm[t2]),min(abs(xp-xm[t2]-Xmax),\
			(xp-xm[t2]+Xmax))))*(min(abs(xp-xm[t2]),min(abs(xp-xm[t2]-Xmax),\
			(xp-xm[t2]+Xmax))))+(min(abs(yp-ym[t2]),min(abs(yp-ym[t2]-Ymax),\
			(yp-ym[t2]+Ymax))))*(min(abs(yp-ym[t2]),min(abs(yp-ym[t2]-Ymax),\
			(yp-ym[t2]+Ymax))))<r[t2]*r[t2]);
			//cout << "intrustion2 is "<< intrusion2 << endl;
			intrusion=intrusion1||intrusion2;
			t=intrusion1*t1+intrusion2*t2;
			t=max(intrusion1*t1,intrusion2*t2);
//trusion=(t!=0)&&((min(abs(xp-xm[t]),min(abs(xp-xm[t]-Xmax),(xp-xm[t]+Xmax))))*(min(abs(xp-xm[t]),min(abs(xp-xm[t]-Xmax),(xp-xm[t]+Xmax))))+(min(abs(yp-ym[t]),min(abs(yp-ym[t]-Ymax),(yp-ym[t]+Ymax))))*(min(abs(yp-ym[t]),min(abs(yp-ym[t]-Ymax),(yp-ym[t]+Ymax))))<r[t]*r[t]);
		}
		if (!intrusion){
			xip[k]=xp;
			yip[k]=yp;
			k++;
			}
		}
	end=clock();
	cout << "done! Elpased time "<<double(diffclock(end,begin))/10 << " s"<< endl;
 
	// simulating pseudo random walk
	cout<<"simulating pseudo random walk...";
//	cout << "random seed 3"<< random_uniform_0_1()<< endl;
	begin=clock();
	int ErrorPoints=0;
	//double xstore[NP][timepoints], ystore[NP][timepoints], zstore[NP][timepoints];
	double MOM2[timepoints][sizetheta],MOM4[timepoints][sizetheta];
	for(jk=0;jk<timepoints;jk++){
		for(kl=0;kl<sizetheta;kl++){
			MOM2[jk][kl]=0;
			MOM4[jk][kl]=0;
			}
		}
	double theta [sizetheta];
	for(i=0;i<sizetheta;i++)
		theta[i]=i*Pi/sizetheta; //do not understand

	double timesim[timepoints];
	for (i=0;i<timepoints;i++){
		timesim[i]=TT/timepoints*i; //do not understand
		}
	for (k=0;k<NP;k++){		//loop over the different points
		cout<<"k: "<<k<<endl;
		double dXtemp[timepoints],dYtemp[timepoints],dZtemp[timepoints];
		int N=1;
		bool errorcodeee=1;
		bool position,position1,position2;
		intrusion=0;
		xpp=xip[k];
		ypp=yip[k];
		zpp=0;
		int NT=0;
		int Tt=0;
		int xjm,yjm,xjmtemp,yjmtemp,xjm1,yjm1;
		int xjmp = 0;
		int yjmp = 0;
		for (tru=dt;tru<=TT;tru+=dt){//loop over the different timepoints
			xpr=xpp/res;
			ypr=ypp/res;
			if (tru==dt){
				xfl=floor(xpr);yfl=floor(ypr);xce=ceil(xpr);yce=ceil(ypr);
				t=max(max(A[xfl][yfl],A[xce][yce]),max(A[xfl][yce],A[xce][yce]));
				}
			position=0; //define if point is INSIDE or OUTSIDE cylinders -> can be omitted in the stripped version
		/*	if (t<Nmax)
				position=(t!=0)&&((min(abs(xpp-xm[t]),min(abs(xpp-xm[t]-Xmax),(xpp-xm[t]+Xmax))))*(min(abs(xpp-xm[t]),min(abs(xpp-xm[t]-Xmax),(xpp-xm[t]+Xmax))))+(min(abs(ypp-ym[t]),min(abs(ypp-ym[t]-Ymax),(ypp-ym[t]+Ymax))))*(min(abs(ypp-ym[t]),min(abs(ypp-ym[t]-Ymax),(ypp-ym[t]+Ymax))))<r[t]*r[t]);
			else {
				t1=floor(t/1000);
				t2=t-t1*1000;
				position1=(t1!=0)&&((min(abs(xpp-xm[t1]),min(abs(xpp-xm[t1]-Xmax),(xpp-xm[t1]+Xmax))))*(min(abs(xpp-xm[t1]),min(abs(xpp-xm[t1]-Xmax),(xpp-xm[t1]+Xmax))))+(min(abs(ypp-ym[t1]),min(abs(ypp-ym[t1]-Ymax),(ypp-ym[t1]+Ymax))))*(min(abs(ypp-ym[t1]),min(abs(ypp-ym[t1]-Ymax),(ypp-ym[t1]+Ymax))))<r[t1]*r[t1]);
				position2=(t2!=0)&&((min(abs(xpp-xm[t2]),min(abs(xpp-xm[t2]-Xmax),(xpp-xm[t2]+Xmax))))*(min(abs(xpp-xm[t2]),min(abs(xpp-xm[t2]-Xmax),(xpp-xm[t2]+Xmax))))+(min(abs(ypp-ym[t2]),min(abs(ypp-ym[t2]-Ymax),(ypp-ym[t2]+Ymax))))*(min(abs(ypp-ym[t2]),min(abs(ypp-ym[t2]-Ymax),(ypp-ym[t2]+Ymax))))<r[t2]*r[t2]);
				position=position1||position2;
				t=position1*t1+position2*t2; // assuming one of the two is zero
//				position=(t!=0)&&((min(abs(xpp-xm[t]),min(abs(xpp-xm[t]-Xmax),(xpp-xm[t]+Xmax))))*(min(abs(xpp-xm[t]),min(abs(xpp-xm[t]-Xmax),(xpp-xm[t]+Xmax))))+(min(abs(ypp-ym[t]),min(abs(ypp-ym[t]-Ymax),(ypp-ym[t]+Ymax))))*(min(abs(ypp-ym[t]),min(abs(ypp-ym[t]-Ymax),(ypp-ym[t]+Ymax))))<r[t]*r[t]);
			}*/

		// Calculation of point B, new point in a random direction
			/*if (position){
				cout<<"error t:"<<max(max(A[xfl][yfl],A[xce][yce]),max(A[xfl][yce],A[xce][yce]))<<endl;
				stap=stepIN;}
			else*/
			stap=stepOUT;
			v = JKISS(v);// The way to assign random direction is wrong, Feb 15th, 2017
			vs=v.r*2*Pi;
			v = JKISS(v);
			vt=acos(1-v.r*2); // modified Feb 15th, 2017
			cosvs=cos(vs);sinvs=sin(vs);sinvt=sin(vt);
			vp[0]=sinvt*cosvs;
			vp[1]=sinvt*sinvs;
			vp[2]=cos(vt);
			xB=xpp+stap*vp[0];
			yB=ypp+stap*vp[1];
			zB=zpp+stap*vp[2];//mistake ends here
			// check for jumps at the edges
			xjm=0;// It can be wrong if do not accumlate jumps in the following steps
			yjm=0;
			if (xB>Xmax){
				xB=xB-Xmax;xjm=1;}
			if(xB<0){
				xB=xB+Xmax;xjm=-1;}
			if (yB>Ymax){
				yB=yB-Ymax;yjm=1;}
			if (yB<0){
				yB=yB+Ymax;yjm=-1;}

			xpr=xB/res;ypr=yB/res;
			if (position){
				intrusion=(xB-xm[t])*(xB-xm[t])+(yB-ym[t])*(yB-ym[t])>=r[t]*r[t];
				}
			else{
				xfl=floor(xpr);yfl=floor(ypr);xce=ceil(xpr);yce=ceil(ypr);
				t=max(max(A[xfl][yfl],A[xce][yce]),max(A[xfl][yce],A[xce][yfl]));
				if (t<Nmax)
					intrusion=(t!=0)&&((min(abs(xB-xm[t]),min(abs(xB-xm[t]-Xmax),(xB-xm[t]+Xmax))))*(min(abs(xB-xm[t]),min(abs(xB-xm[t]-Xmax),(xB-xm[t]+Xmax))))+(min(abs(yB-ym[t]),min(abs(yB-ym[t]-Ymax),(yB-ym[t]+Ymax))))*(min(abs(yB-ym[t]),min(abs(yB-ym[t]-Ymax),(yB-ym[t]+Ymax))))<r[t]*r[t]);
				else {
	//				cout<<"378: t>Nmax: "<<t<<endl;
					t1=floor(t/Nmax);
					t2=t-t1*Nmax;
					intrusion1=(t1!=0)&&((min(abs(xB-xm[t1]),min(abs(xB-xm[t1]-Xmax),(xB-xm[t1]+Xmax))))*(min(abs(xB-xm[t1]),min(abs(xB-xm[t1]-Xmax),(xB-xm[t1]+Xmax))))+(min(abs(yB-ym[t1]),min(abs(yB-ym[t1]-Ymax),(yB-ym[t1]+Ymax))))*(min(abs(yB-ym[t1]),min(abs(yB-ym[t1]-Ymax),(yB-ym[t1]+Ymax))))<r[t1]*r[t1]);
					intrusion2=(t2!=0)&&((min(abs(xB-xm[t2]),min(abs(xB-xm[t2]-Xmax),(xB-xm[t2]+Xmax))))*(min(abs(xB-xm[t2]),min(abs(xB-xm[t2]-Xmax),(xB-xm[t2]+Xmax))))+(min(abs(yB-ym[t2]),min(abs(yB-ym[t2]-Ymax),(yB-ym[t2]+Ymax))))*(min(abs(yB-ym[t2]),min(abs(yB-ym[t2]-Ymax),(yB-ym[t2]+Ymax))))<r[t2]*r[t2]);
					//t=intrusion1*t1+intrusion2*t2;
					t=max(intrusion1*t1,intrusion2*t2);					
					intrusion=intrusion1||intrusion2;
	//				cout<<"387: t chosen: "<<t<<endl;
					//intrusion=(t!=0)&&((min(abs(xB-xm[t]),min(abs(xB-xm[t]-Xmax),(xB-xm[t]+Xmax))))*(min(abs(xB-xm[t]),min(abs(xB-xm[t]-Xmax),(xB-xm[t]+Xmax))))+(min(abs(yB-ym[t]),min(abs(yB-ym[t]-Ymax),(yB-ym[t]+Ymax))))*(min(abs(yB-ym[t]),min(abs(yB-ym[t]-Ymax),(yB-ym[t]+Ymax))))<r[t]*r[t]);
					}
				}
			if(!intrusion){
				xnp=xB;ynp=yB;znp=zB;
				}
			else {//calculation of point C, where the particles crosses the boundary
				xjmtemp=0;
				yjmtemp=0;
				while(intrusion && errorcodeee){
					xmin=xpp-xm[t];
					X=xmin;
					if (abs(xmin+Xmax)<min(abs(xmin),abs(xmin-Xmax)))
						X=xmin+Xmax;
					else if (abs(xmin-Xmax)<min(abs(xmin),abs(xmin+Xmax)))
						X=xmin-Xmax;
					ymin=ypp-ym[t];
					Y=ymin;
					if (abs(ymin+Ymax)<min(abs(ymin),abs(ymin-Ymax)))
						Y=ymin+Ymax;
					else if (abs(ymin-Ymax)<min(abs(ymin),abs(ymin+Ymax)))
						Y=ymin-Ymax;
					R=r[t];
					a=vp[0]*X+vp[1]*Y;
					b=X*X+Y*Y-R*R;
					c=vp[0]*vp[0]+vp[1]*vp[1];
					d=a*a-b*c;
					amin=-a-sqrt(d);
					aplus=-a+sqrt(d);
					/*if (position){
						s=max(aplus,amin);
					}
					else{*/
						if (abs(aplus)<abs(amin))
							s=aplus/c;
						else
							s=amin/c;
					//	}
					xC=xpp+s*vp[0];
					yC=ypp+s*vp[1];
					zC=zpp+s*vp[2];
					xjm=0;yjm=0;
					// check for jumps at the edges
					if (xC>Xmax){
						xC=xC-Xmax;xjm=1;}
					if (xC<0){
						xC=xC+Xmax;xjm=-1;}
					if (yC>Ymax){
						yC=yC-Ymax;yjm=1;}
					if (yC<0){
						yC=yC+Ymax;yjm=-1;}
					srest=stap-s;
					/*if (abs(s)<1e-13){
						errorcodeee=0;
						//cout<<"problem ";
						break;
						}*/
					// point C at the boundary
					// calculation of new point on travel path point D: intrusion or reflection
					//pex=random_uniform_0_1();
					sn=abs(R-sqrt(X*X+Y*Y));
					//pIE=2*sn*1e-6*P/Din;// the normal distance to the surface should be taken into account here
					//pEI=2*sn*1e-6*P/Dout;
					/*if(position && pex<=pIE){//particle crosses membrane from inside to outside
						srest=srest*sqrt(Dout/Din);
						position=0;
						}
					else if(!position && pex<=pEI){//particle crosses membrane from outside to inside
						srest=srest*sqrt(Din/Dout);
						position=1;
						}
					else{*/
					//elastic collision
						CminX=xC-xm[t];
                        XminC = CminX;
						if (abs(CminX+Xmax)<min(abs(CminX),abs(CminX-Xmax)))
							XminC=CminX+Xmax;
						else if(abs(CminX-Xmax)<min(abs(CminX),abs(CminX+Xmax)))
							XminC=CminX-Xmax;
						CminY=yC-ym[t];
                        YminC = CminY;
						if (abs(CminY+Ymax)<min(abs(CminY),abs(CminY-Ymax)))
							YminC=CminY+Ymax;
						else if(abs(CminY-Ymax)<min(abs(CminY),abs(CminY+Ymax)))
							YminC=CminY-Ymax;
						cosphi=XminC/r[t];sinphi=YminC/r[t];
						cos2phi=cosphi*cosphi-sinphi*sinphi;
						sin2phi=2*cosphi*sinphi;
						sinvt=-sinvt;
						cosvss=cos2phi*cosvs+sin2phi*sinvs;
						sinvss=sin2phi*cosvs-cos2phi*sinvs;
						cosvs=cosvss;sinvs=sinvss;
						vp[0]=cosvs*sinvt;
						vp[1]=sinvs*sinvt;
					//	}
					xD=xC+srest*vp[0];
					yD=yC+srest*vp[1];
					zD=zC+srest*vp[2];
					xjm1=0;yjm1=0; //check again for jumps at the edges
					if (xD>Xmax){
						xD=xD-Xmax;xjm1=1;}
					if (xD<0){
						xD=xD+Xmax;xjm1=-1;}
					if (yD>Ymax){
						yD=yD-Ymax;yjm1=1;}
					if (yD<0){
						yD=yD+Ymax;yjm1=-1;}
					xpr=xD/res;
					ypr=yD/res; //modified Feb. 15th, 2017
				/*	if (position){
						intrusion=(xD-xm[t])*(xD-xm[t])+(yD-ym[t])*(yD-ym[t])>=r[t]*r[t];
						}
					else{*/
						xfl=floor(xpr);yfl=floor(ypr);xce=ceil(xpr);yce=ceil(ypr);
						t=max(max(A[xfl][yfl],A[xce][yce]),max(A[xfl][yce],A[xce][yfl]));
						if (t<Nmax)
						intrusion=(t!=0)&&((min(abs(xD-xm[t]),min(abs(xD-xm[t]-Xmax),(xD-xm[t]+Xmax))))*(min(abs(xD-xm[t]),min(abs(xD-xm[t]-Xmax),(xD-xm[t]+Xmax))))+(min(abs(yD-ym[t]),min(abs(yD-ym[t]-Ymax),(yD-ym[t]+Ymax))))*(min(abs(yD-ym[t]),min(abs(yD-ym[t]-Ymax),(yD-ym[t]+Ymax))))<r[t]*r[t]);
						else {
							t1=floor(t/Nmax);
							t2=t-t1*Nmax;
							intrusion1=(t1!=0)&&((min(abs(xD-xm[t1]),min(abs(xD-xm[t1]-Xmax),(xD-xm[t1]+Xmax))))*(min(abs(xD-xm[t1]),min(abs(xD-xm[t1]-Xmax),(xD-xm[t1]+Xmax))))+(min(abs(yD-ym[t1]),min(abs(yD-ym[t1]-Ymax),(yD-ym[t1]+Ymax))))*(min(abs(yD-ym[t1]),min(abs(yD-ym[t1]-Ymax),(yD-ym[t1]+Ymax))))<r[t1]*r[t1]);
							intrusion2=(t2!=0)&&((min(abs(xD-xm[t2]),min(abs(xD-xm[t2]-Xmax),(xD-xm[t2]+Xmax))))*(min(abs(xD-xm[t2]),min(abs(xD-xm[t2]-Xmax),(xD-xm[t2]+Xmax))))+(min(abs(yD-ym[t2]),min(abs(yD-ym[t2]-Ymax),(yD-ym[t2]+Ymax))))*(min(abs(yD-ym[t2]),min(abs(yD-ym[t2]-Ymax),(yD-ym[t2]+Ymax))))<r[t2]*r[t2]);
							intrusion=intrusion1||intrusion2;
							//t=intrusion1*t1+intrusion2*t2;
							t=max(intrusion1*t1,intrusion2*t2);
						}
				//	}
					if(!intrusion){//particle does not cross a membrane
						xnp=xD;ynp=yD;znp=zD;
						xjm=xjm+xjm1;yjm=yjm+yjm1;
						}
					else{//particle crosses a membrane
						stap=srest;
						xpp=xC;
						ypp=yC;
						xjmtemp=xjmtemp+xjm;
						yjmtemp=yjmtemp+yjm;
						}
					xjm=xjmtemp+xjm;
					yjm=yjmtemp+yjm;
					}
				}			
				xjmp=xjmp+xjm;
				yjmp=yjmp+yjm;
				xpp=xnp;
				ypp=ynp;
				zpp=znp;
				/*if (NT<0){
					cout<<"problem"<<endl;
					}*/
				if (NT/Tstep==Tt){
					dXtemp[Tt]=xpp-xip[k]+Xmax*xjmp;
					dYtemp[Tt]=ypp-yip[k]+Ymax*yjmp;
					dZtemp[Tt]=zpp;
					//xstore[k][Tt]=xnp;
					//ystore[k][Tt]=ynp;
					//zstore[k][Tt]=znp;
					Tt++;
					}
				NT++;
				/*if (!errorcodeee){
					cout<<"problem"<<endl;
				}*/
				}//end of time travel loop
			
				if(errorcodeee){
			for(jk=0;jk<timepoints;jk++){
				for(kl=0;kl<sizetheta;kl++){
					MOM2[jk][kl]=MOM2[jk][kl]+(dZtemp[jk]*cos(theta[kl])+dXtemp[jk]*sin(theta[kl]))*(dZtemp[jk]*cos(theta[kl])+dXtemp[jk]*sin(theta[kl]));
					MOM4[jk][kl]=MOM4[jk][kl]+(dZtemp[jk]*cos(theta[kl])+dXtemp[jk]*sin(theta[kl]))*(dZtemp[jk]*cos(theta[kl])+dXtemp[jk]*sin(theta[kl]))*(dZtemp[jk]*cos(theta[kl])+dXtemp[jk]*sin(theta[kl]))*(dZtemp[jk]*cos(theta[kl])+dXtemp[jk]*sin(theta[kl]));
					}
				}
			}
		else 
			ErrorPoints++;	
		}// end of for-loop over the different seed points
		int NPfinal=NP-ErrorPoints;
		end=clock();
	cout << "Done! Elpased time "<<double(diffclock(end,begin))/10 << " s"<< endl;
	cout << "Number of error points: "<<ErrorPoints<<endl;
	cout << "Writing data to file...";

	begin=clock();
	char MOM2filename[100];
	char MOM4filename[100];
	char Xfilename[100];
	char Yfilename[100];
	char Zfilename[100];


	int Pfile=P*1e8;
	int TTfile=TT*1e-9;
	int Rfile = radius;		//in micrometers
	int Densfile = density*100;;
	v = JKISS(v);
	//cout<< "v.x is " << v.x << endl;
	double timeSIM1 = JKISS(v).r*1000000;
	unsigned int timeSIM2 = v.x;

	sprintf(MOM2filename,"MOM2simTT%dTN%dNP%ldP%dSQUAREradius%dDensity%dtime%f_%u.txt",TTfile,TN,NP,Pfile,Rfile,Densfile,timeSIM1,timeSIM2);
	sprintf(MOM4filename,"MOM4simTT%dTN%dNP%ldP%dSQUAREradius%dDensity%dtime%f_%u.txt",TTfile,TN,NP,Pfile,Rfile,Densfile,timeSIM1,timeSIM2);
	//sprintf(Xfilename,"XsimTT%dTN%dNP%dP%dSQUAREradius%dDensity%d.txt",TTfile,TN,NP,Pfile,Rfile,Densfile);
	//sprintf(Yfilename,"YsimTT%dTN%dNP%dP%dSQUAREradius%dDensity%d.txt",TTfile,TN,NP,Pfile,Rfile,Densfile);
	//sprintf(Zfilename,"ZsimTT%dTN%dNP%dP%dSQUAREradius%dDensity%d.txt",TTfile,TN,NP,Pfile,Rfile,Densfile);

	cout.setf(ios_base::internal,ios_base::adjustfield);
	ofstream MOM2file(MOM2filename);
	ofstream MOM4file(MOM4filename);
	for(jk=0;jk<timepoints;jk++){
			MOM2file<<setw(4)<<setfill(' ')<<jk*Tstep*dt/1e9<<"ms ";
			MOM4file<<setw(4)<<setfill(' ')<<jk*Tstep*dt/1e9<<"ms ";
			for(kl=0;kl<sizetheta;kl++){
				MOM2[jk][kl]=MOM2[jk][kl]/NPfinal;
				MOM4[jk][kl]=MOM4[jk][kl]/NPfinal;
				MOM2file<<setw(32)<<setfill(' ')<<MOM2[jk][kl];
				MOM4file<<setw(30)<<setfill(' ')<<MOM4[jk][kl];
				}
			MOM2file<<endl;
			MOM4file<<endl;
			}
	MOM2file.close();
	MOM4file.close();




	/*ofstream Xfile(Xfilename);
	ofstream Yfile(Yfilename);
	ofstream Zfile(Zfilename);
	for(jk=0;jk<timepoints;jk++){
		Xfile<<setw(4)<<setfill(' ')<<jk*Tstep*dt/1e9<<"ms ";
		Yfile<<setw(4)<<setfill(' ')<<jk*Tstep*dt/1e9<<"ms ";
		Zfile<<setw(4)<<setfill(' ')<<jk*Tstep*dt/1e9<<"ms ";
		for(kl=0;kl<NPfinal;kl++){
			Xfile<<setw(12)<<setfill(' ')<<xstore[kl][jk]<<" ";
			Yfile<<setw(12)<<setfill(' ')<<ystore[kl][jk]<<" ";
			Zfile<<setw(12)<<setfill(' ')<<zstore[kl][jk]<<" ";
			}
		Xfile<<endl;
		Yfile<<endl;
		Zfile<<endl;
		}
	Xfile.close();
	Yfile.close();
	Zfile.close();*/
	end=clock();
	cout << "Done! Elpased time "<<double(diffclock(end,begin))/10 << " s"<< endl;
    cout << "\n Press exit sto terminate program\n";
	//delete pointers
	delete A,xm,ym,r;


    return 0;
  }



