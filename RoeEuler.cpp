#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>
using namespace std;

#define npts 400
#define maxt 5500
#define gamma 1.4
#define G1 (gamma-1)/(2*gamma)
#define G2 (gamma+1)/(2*gamma)
#define G3 (2*gamma)/(gamma-1)
#define G4 2.0/(gamma-1)
#define G5 2.0/(gamma+1)
#define G6 (gamma-1)/(gamma+1)
#define G7 (gamma-1)/2.0
#define G8 (gamma-1.0)
#define len 100
#define MIN(a,b) (((a)<(b)) ? (a):(b))
#define MAX(a,b) (((a)>(b)) ? (a):(b))
#define K -1.0

int InitArrays(double (&d)[npts],double (&u)[npts],double (&p)[npts],double (&c)[npts],double (&unknd)[npts],double (&unkn_dv)[npts],double (&unkn_de)[npts],double dx);
double CalcCFL(double d[npts],double u[npts],double p[npts],double c[npts],double dx);
int RoeRiemann(double (&d)[npts],double (&u)[npts],double (&p)[npts],double (&c)[npts],double (&unkn_d)[npts],double (&unkn_dv)[npts],double (&unkn_de)[npts],double (&fluxDensity)[npts],double (&fluxDensityVelocity)[npts],double (&fluxDensityEnergy)[npts]);
int pressurefunc(double &f,double &fd,double ptemp,double dk,double pk,double ck);
int update(double (&d)[npts],double (&u)[npts],double (&p)[npts],double (&c)[npts],double (&unkn_d)[npts],double (&unkn_dv)[npts],double (&unkn_de)[npts],double (&fluxDensity)[npts],double (&fluxDensityVelocity)[npts],double (&fluxDensityEnergy)[npts],double dt,double dx);
int RoeFlux(double (&d)[npts],double (&u)[npts],double (&p)[npts],double (&c)[npts],double (&unkn_d)[npts],double (&unkn_dv)[npts],double (&unkn_de)[npts],double (&fluxDensity)[npts],double (&fluxDensityVelocity)[npts],double (&fluxDensityEnergy)[npts]);
double superbee(double r);
double calcdLcontribution(double d[npts],int i);
double calcdRcontribution(double d[npts],int i);
double calcuLcontribution(double u[npts],int i);
double calcuRcontribution(double u[npts],int i);
double calcpLcontribution(double p[npts],int i);
double calcpRcontribution(double p[npts],int i);
double sgn(double x){
	// this is the sign function needed for superbee //
	if(x>=0){
		return 1;	
	}	
	if(x<0){
		return -1;	
	}
	
}


int Boundary(double (&d)[npts],double (&u)[npts],double (&p)[npts],double (&unkn_d)[npts],double (&unkn_dv)[npts],double (&unkn_de)[npts]){

	// this function assigns the boundary conditions of the domain //
	d[0]=d[1];
	u[0]=u[1];
	p[0]=p[1];
	unkn_d[0]=d[0];
	unkn_dv[0]=d[0]*u[0];
	unkn_de[0]=(p[0]/G8)+(0.5*d[0]*u[0]*u[0]);
	
	d[npts-1]=d[npts-2];
	u[npts-1]=u[npts-2];
	p[npts-1]=p[npts-2];
	unkn_d[npts-1]=d[npts-2];
	unkn_dv[npts-1]=d[npts-2]*u[npts-2];
	unkn_de[npts-1]=(p[npts-2]/G8)+(0.5*d[npts-2]*u[npts-2]*u[npts-2]);

	return 0;
}

int main(void){

	// this is the main entry point of the code //
	double d[npts],u[npts],p[npts],c[npts],e[npts],fluxDensity[npts],fluxDensityVelocity[npts],fluxDensityEnergy[npts];
	double unkn_d[npts],unkn_dv[npts],unkn_de[npts];
	double t,dt,dx;
	ofstream myfile;
	myfile.open("dump");
	t=0.0;
	dx=(double)len/npts;
	
	/*put initial conditions */
	InitArrays(d,u,p,c,unkn_d,unkn_dv,unkn_de,dx);
	
	for(int i=0;i<maxt;i++){
		// correct boundary values //
		Boundary(d,u,p,unkn_d,unkn_dv,unkn_de);
		/*calculate CFL condition */
		dt=CalcCFL(d,u,p,c,dx);
		/*solve riemann problem and obtain fluxes for euler equation*/
		RoeRiemann(d,u,p,c,unkn_d,unkn_dv,unkn_de,fluxDensity,fluxDensityVelocity,fluxDensityEnergy);
		/*update unknowns */	
		update(d,u,p,c,unkn_d,unkn_dv,unkn_de,fluxDensity,fluxDensityVelocity,fluxDensityEnergy,dt,dx);
		
	}//end time loop 
	
	/*print results */
	double x;
	x=0;
	for(int i=0;i<npts;i++){
		x+=dx;
		myfile<<x<<"    "<<u[i]<<endl;	
	}
	myfile.close();

	return 0;
}

double superbee(double r){
	double fi,min1,min2,max1;
	min1=MIN(r,2.0);
	min2=MIN(2*r,1.0);
	max1=MAX(min1,min2);
	fi=MAX(0,max1);	
	

	return fi;
}

double calcdLcontribution(double d[npts],int i){
		// this function calculates the left density condition //

		double duplus,duminus,k,fiL,fiL2,fiR,fiR2,rL,rR,duplus2,value;	
	
		
		duplus=(1e-9+fabs(d[i+1]-d[i]))*sgn((d[i+1]-d[i]));
		
		duminus=(1e-9+fabs(d[i]-d[i-1]))*sgn((d[i]-d[i-1]));
		if(i==npts-1){
			duplus2=duplus;
		}
		else {
		
			duplus2=(1e-9+fabs(d[i+2]-d[i+1]))*sgn((d[i+2]-d[i+1]));			
		}
		rL=duplus/(duminus);
		rR=duplus/(duplus2);
		fiL=superbee(rL);
		fiL2=superbee(1/rL);
		value=0.5*duminus*(0.5*((1-K)*fiL+(1+K)*rR*fiL2));
		
	return value;
}

double calcuLcontribution(double u[npts],int i){
		// this function calculates the left velocity u  condition //
		double duplus,duminus,k,fiL,fiL2,fiR,fiR2,rL,rR,duplus2,value;	
	
	
		duplus=(1e-9+fabs(u[i+1]-u[i]))*sgn((u[i+1]-u[i]));		
		duminus=(1e-9+fabs(u[i]-u[i-1]))*sgn((u[i]-u[i-1]));
		if(i==npts-1){
			duplus2=duplus;
		}
		else {
			
			duplus2=(1e-9+fabs(u[i+2]-u[i+1]))*sgn((u[i+2]-u[i+1]));			
		}
		rL=duplus/(duminus);
		rR=duplus/(duplus2);
		fiL=superbee(rL);
		fiL2=superbee(1/rL);
		value=0.5*duminus*(0.5*((1-K)*fiL+(1+K)*rR*fiL2));
		
	return value;
}
double calcpLcontribution(double p[npts],int i){
		// this function calculates the left pressure condition //
		double duplus,duminus,k,fiL,fiL2,fiR,fiR2,rL,rR,duplus2,value;	
	
		
		duplus=(1e-9+fabs(p[i+1]-p[i]))*sgn((p[i+1]-p[i]));		
		duminus=(1e-9+fabs(p[i]-p[i-1]))*sgn((p[i]-p[i-1]));
		if(i==npts-1){
			duplus2=duplus;
		}
		else {			
			duplus2=(1e-9+fabs(p[i+2]-p[i+1]))*sgn((p[i+2]-p[i+1]));			
		}
		rL=duplus/(duminus);
		rR=duplus/(duplus2);
		fiL=superbee(rL);
		fiL2=superbee(1/rL);
		value=0.5*duminus*(0.5*((1-K)*fiL+(1+K)*rR*fiL2));
		
	return value;
}

double calcdRcontribution(double d[npts],int i){
	// this function calculates the right density condition //
		double duplus,duminus,k,fiL,fiL2,fiR,fiR2,rL,rR,duplus2,value;	
		duplus=(1e-9+fabs(d[i+1]-d[i]))*sgn((d[i+1]-d[i]));
		duminus=(1e-9+fabs(d[i]-d[i-1]))*sgn((d[i]-d[i-1]));
		if(i==npts-1){
			duplus2=duplus;
		}
		else {

			duplus2=(1e-9+fabs(d[i+2]-d[i+1]))*sgn((d[i+2]-d[i+1]));			
		}
		rL=duplus/(duminus);
		rR=duplus/(duplus2);
		fiR=superbee(rR);
		fiR2=superbee(1/rR);
		value=0.5*duplus*(0.5*((1-K)*(1/(rR))*fiR+(1+K)*fiR2));

	return value;
}

double calcuRcontribution(double u[npts],int i){
		// this function calculates the right velocity u condition //
		double duplus,duminus,k,fiL,fiL2,fiR,fiR2,rL,rR,duplus2,value;	

		duplus=(1e-9+fabs(u[i+1]-u[i]))*sgn((u[i+1]-u[i]));
		duminus=(1e-9+fabs(u[i]-u[i-1]))*sgn((u[i]-u[i-1]));
		if(i==npts-1){
			duplus2=duplus;
		}
		else {
			//duplus2=u[i+2]-u[i+1];
			duplus2=(1e-9+fabs(u[i+2]-u[i+1]))*sgn((u[i+2]-u[i+1]));			
		}
		rL=duplus/(duminus);
		rR=duplus/(duplus2);
		fiR=superbee(rR);
		fiR2=superbee(1/rR);
		value=0.5*duplus*(0.5*((1-K)*(1/(rR))*fiR+(1+K)*fiR2));
		
	return value;
}

double calcpRcontribution(double p[npts],int i){
		// this function calculates the right pressure condition //
		double duplus,duminus,k,fiL,fiL2,fiR,fiR2,rL,rR,duplus2,value;	
		
		duplus=(1e-9+fabs(p[i+1]-p[i]))*sgn((p[i+1]-p[i]));
		duminus=(1e-9+fabs(p[i]-p[i-1]))*sgn((p[i]-p[i-1]));
		if(i==npts-1){
			duplus2=duplus;
		}
		else {
			
			duplus2=(1e-9+fabs(p[i+2]-p[i+1]))*sgn((p[i+2]-p[i+1]));		
		}
		rL=duplus/(duminus);
		rR=duplus/(duplus2);
		fiR=superbee(rR);
		fiR2=superbee(1/rR);
		value=0.5*duplus*(0.5*((1-K)*(1/(rR))*fiR+(1+K)*fiR2));
		
	return value;
}

int RoeRiemann(double (&d)[npts],double (&u)[npts],double (&p)[npts],double (&c)[npts],double (&unkn_d)[npts],double (&unkn_dv)[npts],double (&unkn_de)[npts],double (&fluxDensity)[npts],double (&fluxDensityVelocity)[npts],double (&fluxDensityEnergy)[npts]){
	double dl,dr,ul,ur,pl,pr,pguess,pold,udiff,ptemp,fl,fld,fr,frd,diff,cl,cr,utemp,xovert,hl,hr,ratio;
	double davg,uavg,havg,cavg,drho,du,dp,dv1,dv2,dv3,ws1,ws2,ws3,da,favg1,favg2,favg3,flux1,flux2,flux3;
	double Rmatrix[3][3],lamda[3];
	
	double duplus,duminus,k,fiL,fiL2,fiR,fiR2,rL,rR,duplus2;
	xovert=0.0;	
	for(int i=0;i<npts-1;i++){
		dl=d[i]+calcdLcontribution(d,i);
		dr=d[i+1]-calcdRcontribution(d,i);
		ul=u[i]+calcuLcontribution(u,i);
		ur=u[i+1]-calcuRcontribution(u,i);
		pl=p[i]+calcpLcontribution(p,i);
		pr=p[i+1]-calcpRcontribution(p,i);
/*		dl=d[i];
		dr=d[i+1];
		ul=u[i];
		ur=u[i+1];
		pl=p[i];
		pr=p[i+1];
	*/	//cl=c[i];
		//cr=c[i+1];
		cl=sqrt(gamma*pl/dl);
		cr=sqrt(gamma*pr/dr);
		hl=(ul+pl)/dl;
		hr=(ur+pr)/dr;		
		

		/*compute rhoe average values */
		ratio=sqrt(dr/dl);
		davg=ratio*dl;
		uavg=(ul+ratio*ur)/(1+ratio);
		havg=(hl+ratio*hr)/(1+ratio);
		cavg=(gamma-1)*(havg-0.5*uavg*uavg);
		

		/*differences in primitive variables */
		drho=dr-dl;
		du=ur-ul;
		dp=pr-pl;
		
		


		/* eigen values at interface i+1/2*/
		lamda[0]=fabs(uavg);
		lamda[1]=fabs(uavg+cavg);
		lamda[2]=fabs(uavg-cavg);

		


		double lamdal[3],lamdar[3];
		double eta1,eta2,eta3;
		
		lamdal[0]=fabs(u[i]);
		lamdal[1]=fabs(u[i]+c[i]);
		lamdal[2]=fabs(u[i]-c[i]);
		lamdar[0]=fabs(u[i+1]);
		lamdar[1]=fabs(u[i+1]+c[i+1]);
		lamdar[2]=fabs(u[i+1]-c[i+1]);
			
		
		eta1=MAX((lamda[0]-lamdal[0]),(lamdar[0]-lamda[0]));
		eta1=MAX(eta1,0);
		eta2=MAX((lamda[1]-lamdal[1]),(lamdar[1]-lamda[1]));
		eta2=MAX(eta2,0);
		eta3=MAX((lamda[2]-lamdal[2]),(lamdar[2]-lamda[2]));
		eta3=MAX(eta3,0);
		if(lamda[0]<eta1){
			lamda[0]=eta1;		
		}	
		if(lamda[1]<eta2){
			lamda[1]=eta2;		
		}
		if(lamda[2]<eta3){
			lamda[2]=eta3;		
		}
	
		double delw[3],r1[3],r2[3],r3[3];
		r1[0]=1.0;
		r1[1]=uavg;
		r1[2]=uavg*uavg/2.0;
		r2[0]=0.5*davg/cavg*1.0;
		r2[1]=0.5*(davg/cavg)*(uavg+cavg);
		r2[2]=0.5*(davg/cavg)*(havg+(uavg*cavg));
		r3[0]=-0.5*(davg/cavg)*1.0;
		r3[1]=-0.5*(davg/cavg)*(uavg-cavg);
		r3[2]=-0.5*(davg/cavg)*(havg-(uavg-cavg));
		
		delw[0]=drho-(dp/(cavg*cavg));
		delw[1]=du+(dp/(davg*cavg));
		delw[2]=du-(dp/(davg*cavg));
	
		/*compute average flux */
		favg1=0.5*(dl*ul+dr*ur);
		favg2=0.5*((dl*ul*ul+pl)+(dr*ur*ur+pr));
		favg3=0.5*((dl*ul*(cl/(gamma-1)+0.5*ul*ul))+(dr*ur*(cr/(gamma-1)+0.5*ur*ur)));
	
		/*add matrix dissipation term to complete roe flux */
		
			
		favg1=favg1-0.5*((lamda[0]*delw[0]*r1[0])+(lamda[1]*delw[1]*r2[0])+(lamda[2]*delw[2]*r3[0]));
		favg2=favg2-0.5*((lamda[0]*delw[0]*r1[1])+(lamda[1]*delw[1]*r2[1])+(lamda[2]*delw[2]*r3[1]));
		favg3=favg3-0.5*((lamda[0]*delw[0]*r1[2])+(lamda[1]*delw[1]*r2[2])+(lamda[2]*delw[2]*r3[2]));

		//update flux values for gridpoints //
		fluxDensity[i]=favg1;
		fluxDensityVelocity[i]=favg2;
		fluxDensityEnergy[i]=favg3;

			
	}//end for loop for points

		
	return 0;
}


double CalcCFL(double d[npts],double u[npts],double p[npts],double c[npts],double dx){
	
	// this function calculates the CFL values //
	double smax,smaxtemp,cflcoeff,dt;
	cflcoeff=0.05;
	smax=0.0;

	for(int i=0;i<npts;i++){
		c[i]=sqrt(gamma*p[i]/d[i]);	
		smaxtemp=fabs(u[i]+c[i]);
		if(smaxtemp>smax){
			smax=smaxtemp;	
			
		}
		
	}
	dt=cflcoeff*dx/smax;	

	return dt;
}

int update(double (&d)[npts],double (&u)[npts],double (&p)[npts],double (&c)[npts],double (&unkn_d)[npts],double (&unkn_dv)[npts],double (&unkn_de)[npts],double (&fluxDensity)[npts],double (&fluxDensityVelocity)[npts],double (&fluxDensityEnergy)[npts],double dt,double dx){
	//this function updates the internal node points //

	// similar to initialize update the interior nodes only. //
	/*first update unknowns */
	for(int i=1;i<npts-1;i++){
		unkn_d[i]=unkn_d[i]+dt/dx*(fluxDensity[i-1]-fluxDensity[i]);
		unkn_dv[i]=unkn_dv[i]+dt/dx*(fluxDensityVelocity[i-1]-fluxDensityVelocity[i]);
		unkn_de[i]=unkn_de[i]+dt/dx*(fluxDensityEnergy[i-1]-fluxDensityEnergy[i]);	
	}	
	/*update godunov variables */
	for(int i=1;i<npts-1;i++){
		d[i]=unkn_d[i];
		u[i]=unkn_dv[i]/d[i];
		p[i]=G8*(unkn_de[i]-0.5*unkn_dv[i]*u[i]);	
	}
	return 0;
}

int InitArrays(double (&d)[npts],double (&u)[npts],double (&p)[npts],double (&c)[npts],double (&unkn_d)[npts], double (&unkn_dv)[npts],double (&unkn_de)[npts],double dx){

	//this function initilizes the domain with initial condition //

	// loop through the interior nodes 1.....(npts-1). the corner nodes treated in boundary conditions function // 

	for(int i=1;i<npts-1;i++){		
		if((i*dx)<=(len/2)){ //initialize left side of membrane
			unkn_d[i]=1.0;
			unkn_dv[i]=0.0;
			unkn_de[i]=2.5;		
			d[i]=unkn_d[i];
			u[i]=unkn_dv[i]/d[i];
			p[i]=G8*(unkn_de[i]-0.5*unkn_dv[i]*u[i]);		
			c[i]=gamma*p[i]/d[i];	
		}
		else { //initialize right side of membrane 
			unkn_d[i]=0.1;
			unkn_dv[i]=0.0;
			unkn_de[i]=0.25;		
			d[i]=unkn_d[i];
			u[i]=unkn_dv[i]/d[i];
			p[i]=G8*(unkn_de[i]-0.5*unkn_dv[i]*u[i]);
			c[i]=gamma*p[i]/d[i];
		}
		
		
	}

	return 0;
}

/*
fl=load('results1d');
fl=fl';
scatter(fl(1,:),fl(2,:));
plot(fl(1,:),fl(2,:));
*/
