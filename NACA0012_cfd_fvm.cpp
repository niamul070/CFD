/* A 2D Euler solver for unstructured grid flow over NACA 0012 airfoil at transonic flows. Two cases are ran Ma 0.8 at 8 deg angle
/  and Ma 0.4 at 4 deg angle of attack. Structure and Vector data structures of C++ is used to store unknowns. Cell Centered Finite Volume scheme is used.          */

#include<iostream>
#include <math.h>
#include <fstream>
#include<string.h>
#include <cstdlib>
#include<sstream>
#include <vector>
using namespace std;

#define ndim 2
#define gamma 1.4
#define pinf 1.0
#define rhoinf 1.4
//#define alpha .1396263401595463661 //value in radians 8
#define alpha 0.0698131701
#define Ma 0.4



/* This global variable stores the element/cell wise geometric and flow properties */
struct Cell{
	/* geometric properties */
	int id;
	int elemid;
	double area;
	double centroidx;
	double centroidy;
	double normal12x;
	double normal12y;
	double normal23x;
	double normal23y;
	double normal31x;
	double normal31y;
	double normal12unit;
	double normal23unit;
	double normal31unit;
	int bctype;
	int bface[2];
	int wface[2];
	double inwardnormalx;
	double inwardnormaly;
	double elemTimeStep;
	vector <int> neigh;
	/* flow properties */
	double pressure;
	double velocityu;
	double velocityv;
	double speed;
	double speedofsound;
	double density;
	double densityVelocityu;
	double densityVelocityv;
	double densityVelocity;
	double densityEnergy;

	double fluxDensityxContribution;
	double fluxDensityyContribution;
	double fluxDensityVelocityuxContribution;
	double fluxDensityVelocityuyContribution;
	double fluxDensityVelocityvxContribution;
	double fluxDensityVelocityvyContribution;
	double fluxDensityEnergyxContribution;
	double fluxDensityEnergyyContribution;

	/*fluxes */
	double fluxDensity;
	double fluxDensityVelocityu;
	double fluxDensityVelocityv;
	double fluxDensityEnergy;
	
	/*old values */
	double oldDensity;
	double oldDensityVelocityu;
	double oldDensityVelocityv;
	double oldDensityVelocity;
	double oldDensityEnergy;
	
	
};

/* The following global variables stores the mesh properties */
int connectivity[11143*3];
double cordinate[5696*2];
int bnode[249*2];

vector <Cell> triangle;
double ff_densityvelocity,ff_densityvelocityx,ff_densityvelocityy,ff_densityEnergy,ff_densityEnergyxContribution,ff_densityEnergyyContribution,ff_velocityx,ff_velocityy;


/*Name and Parameterlist of the functions used in the code */
double Areacalc(int a,int b, int c);
double Centroidcalcx(int a,int b,int c);
double Centroidcalcy(int a,int b,int c);
int FindNeighbor(int ielemNode,int nelem);
int ComputeNormal(int ielemNode,int nelem);
double ComputeDensityFlux(int ielem,int nelem);
double ComputeDensityVelocityuFlux(int ielem,int nelem);
double ComputeDensityVelocityvFlux(int ielem,int nelem);
double ComputeDensityEnergyFlux(int ielem,int nelem);
double ArtificialViscosityDensity(int ielem,int nelem);
double ArtificialViscosityDensityVelocityu(int ielem,int nelem);
double ArtificialViscosityDensityVelocityv(int ielem,int nelem);
double ArtificialViscosityDensityEnergy(int ielem,int nelem);
int CorrectBoundaryValues(int i, int nelem,int node1,int node2,int node3);
int calcInwardNormal(int ielem,int a1,int b1);
int ReadMeshFile(int npoin,int nelem,int bcnodes);
double convertToDouble(std::string const& s);
int CalcBoundaryFaces(int npoin,int nelem,int bcnodes);
int WriteDatatoFile(int npoin,int nelem);



/* This is the main method. Entry point of the code */

int main(void){

	int npoin,nelem,i,j,a,b,c,RK,bcnodes;
	double x1,y1,x2,y2,x3,y3;
	npoin=5696;
	nelem=11143;
	bcnodes=249;
	RK=3;
	double n1,n2,n3,p;
	int count = 0;
	i =0;
	std::string buffer;
	
	int p1=0;
	int p2=0;
	int p3=0;
	int q1;
	double vinf,speedOfSound,vx,vy,pressure,ff_energy;
	cout<<"starting "<<endl;
	
	/* Read meshfile and fomulate the global variables : coordinate array, connectivity array, boundary nodes etc...  */
	ReadMeshFile(npoin,nelem,bcnodes);
	Cell tempvec; 
	/* Formulate the cells geometric parameters by looping through the connectivity array */
	for(i=0;i<nelem*3;i+=3){
		tempvec.id=i/3;
		a=connectivity[i];
		b=connectivity[i+1];
		c=connectivity[i+2];
		tempvec.area=Areacalc(a,b,c);
		tempvec.centroidx=Centroidcalcx(a,b,c);
		tempvec.centroidy=Centroidcalcy(a,b,c);
		tempvec.elemid=(i/3)+1;	
		tempvec.bctype=3;		
		triangle.push_back(tempvec);	
	}
	
	
	
	/*starting to formulate neighbor elements */
	for(i=0;i<nelem*3;i+=3){
		a=connectivity[i];
		b=connectivity[i+1];
		c=connectivity[i+2];
		FindNeighbor(i,nelem);
	}


	/*compute normals for each elements */
	for(i=0;i<nelem*3;i+=3){
		a=connectivity[i];
		b=connectivity[i+1];
		c=connectivity[i+2];
		ComputeNormal(i,nelem);
	}
	
	/* compute free stream/ farfield parameters */
	speedOfSound=sqrt(gamma*pinf/rhoinf);
	vx=Ma*speedOfSound*cos(alpha);
	vy=Ma*speedOfSound*sin(alpha);
	vinf=sqrt(vx*vx+vy*vy);
	ff_energy = (0.5*(vinf*vinf)) + (pinf/(rhoinf*(gamma-1)));
	ff_densityvelocity=rhoinf*vinf;
	ff_densityvelocityx=vx*rhoinf;
	ff_densityvelocityy=rhoinf*vy;
	ff_velocityx=vx;
	ff_velocityy=vy;
	ff_densityEnergy=rhoinf*ff_energy;
	ff_densityEnergyxContribution=(ff_densityEnergy+pinf)*vx;
	ff_densityEnergyyContribution=(ff_densityEnergy+pinf)*vy;
	
	
	/*put initial conditions to all cells */
	for(i=0;i<nelem*3;i+=3){
		a=connectivity[i];
		b=connectivity[i+1];
		c=connectivity[i+2];		
		triangle[i/3].pressure=pinf;
		triangle[i/3].velocityu=vx;
		triangle[i/3].velocityv=vy;
		triangle[i/3].speed=sqrt(vx*vx+vy*vy);
		triangle[i/3].speedofsound=speedOfSound;
		triangle[i/3].density=rhoinf;
		triangle[i/3].densityVelocityu=rhoinf*vx;
		triangle[i/3].densityVelocityv=rhoinf*vy;
		triangle[i/3].densityVelocity=rhoinf*vinf;
		triangle[i/3].densityEnergy=rhoinf*ff_energy;
		/* initially put all fluxes to zero */
		triangle[i/3].fluxDensity=0.0;
		triangle[i/3].fluxDensityVelocityu=0.0;
		triangle[i/3].fluxDensityVelocityv=0.0;
		triangle[i/3].fluxDensityEnergy=0.0;
		triangle[i/3].fluxDensityxContribution=0.0;
		triangle[i/3].fluxDensityyContribution=0.0;
		triangle[i/3].fluxDensityVelocityuxContribution=0.0;
		triangle[i/3].fluxDensityVelocityuyContribution=0.0;
		triangle[i/3].fluxDensityVelocityvxContribution=0.0;
		triangle[i/3].fluxDensityVelocityvyContribution=0.0;
		triangle[i/3].fluxDensityEnergyxContribution=0.0;
		triangle[i/3].fluxDensityEnergyyContribution=0.0;
		/*set old timestep values to zero */
		triangle[i/3].oldDensity=0.0;
		triangle[i/3].oldDensityVelocityu=0.0;
		triangle[i/3].oldDensityVelocityv=0.0;
		triangle[i/3].oldDensityEnergy=0.0;
	}

	/* obtain the boundary faces from the boundary elements */
	CalcBoundaryFaces(npoin,nelem,bcnodes);
	
	
	cout<<"running it will take a while"<<endl;
	/* now start timeloop  */	
	for(i=0;i<1000;i++){
		//get old values 
		for(int m=0;m<nelem*3;m+=3){
			triangle[m/3].oldDensity=triangle[m/3].density;
			triangle[m/3].oldDensityVelocityu=triangle[m/3].densityVelocityu;
			triangle[m/3].oldDensityVelocityv=triangle[m/3].densityVelocityv;
			triangle[m/3].oldDensityEnergy=triangle[m/3].densityEnergy;	
			triangle[m/3].elemTimeStep=(0.5 * sqrt(triangle[m/3].area)) / (triangle[m/3].speed + triangle[m/3].speedofsound);			
		}
		
		//do RK stages set value to desired one 
		for(j=0;j<RK;j++){

			/*compute flow properties */
			for(int k=0;k<nelem*3;k+=3){			
				triangle[k/3].velocityu=triangle[k/3].densityVelocityu/triangle[k/3].density;
				triangle[k/3].velocityv=triangle[k/3].densityVelocityv/triangle[k/3].density;
				triangle[k/3].speed=sqrt(triangle[k/3].velocityu*triangle[k/3].velocityu+triangle[k/3].velocityv*triangle[k/3].velocityv);
				triangle[k/3].pressure=(gamma-1)*(triangle[k/3].densityEnergy-(0.5*triangle[k/3].density*triangle[k/3].speed*triangle[k/3].speed));
				triangle[k/3].speedofsound=sqrt(gamma*triangle[k/3].pressure/triangle[k/3].density);				
			}
			/*compute flux contributions */
			for(int k=0;k<nelem*3;k+=3){
				triangle[k/3].fluxDensityxContribution=triangle[k/3].densityVelocityu;
				triangle[k/3].fluxDensityyContribution=triangle[k/3].densityVelocityv;
			}
			
			for(int k=0;k<nelem*3;k+=3){
				triangle[k/3].fluxDensityVelocityuxContribution=triangle[k/3].densityVelocityu*triangle[k/3].velocityu+triangle[k/3].pressure;
				triangle[k/3].fluxDensityVelocityuyContribution=triangle[k/3].densityVelocityu*triangle[k/3].velocityv;
			}
			
			for(int k=0;k<nelem*3;k+=3){
				triangle[k/3].fluxDensityVelocityvxContribution=triangle[k/3].densityVelocityv*triangle[k/3].velocityv+triangle[k/3].pressure;
				triangle[k/3].fluxDensityVelocityvyContribution=triangle[k/3].densityVelocityv*triangle[k/3].velocityu;				
			}
			

			for(int k=0;k<nelem*3;k+=3){
				triangle[k/3].fluxDensityEnergyxContribution=(triangle[k/3].densityEnergy+triangle[k/3].pressure)*triangle[k/3].velocityu;
				triangle[k/3].fluxDensityEnergyyContribution=(triangle[k/3].densityEnergy+triangle[k/3].pressure)*triangle[k/3].velocityv;	
					
			}
			
			/* compute fluxes */
			for(int k=0;k<nelem*3;k+=3){				
				/*now compute fluxes in each element */
				triangle[k/3].fluxDensity=ComputeDensityFlux(k,nelem);
				//cout<<i<<"  itr rk: "<<j<<"  elem "<<k/3<<"\tflux density\t"<<triangle[k/3].fluxDensity<<endl;
		
			}
			for(int k=0;k<nelem*3;k+=3){
				triangle[k/3].fluxDensityVelocityu=ComputeDensityVelocityuFlux(k,nelem);				
			}
			for(int k=0;k<nelem*3;k+=3){
				triangle[k/3].fluxDensityVelocityv=ComputeDensityVelocityvFlux(k,nelem);
				
			}
			for(int k=0;k<nelem*3;k+=3){
				triangle[k/3].fluxDensityEnergy=ComputeDensityEnergyFlux(k,nelem);				
			}
			for(int k=0;k<nelem*3;k+=3){
				/* correct the flux by putting artificial viscosity */
				triangle[k/3].fluxDensity=triangle[k/3].fluxDensity-ArtificialViscosityDensity(k,nelem);
				
				triangle[k/3].fluxDensityVelocityu=triangle[k/3].fluxDensityVelocityu-ArtificialViscosityDensityVelocityu(k,nelem);
				triangle[k/3].fluxDensityVelocityv=triangle[k/3].fluxDensityVelocityv-ArtificialViscosityDensityVelocityv(k,nelem);
				triangle[k/3].fluxDensityEnergy=triangle[k/3].fluxDensityEnergy-ArtificialViscosityDensityEnergy(k,nelem);

			}//end artificial viscosity loop
			
			
			for(int k=0;k<nelem*3;k+=3){
				/*add old values to obtain new values */
				triangle[k/3].density = triangle[k/3].oldDensity + (triangle[k/3].elemTimeStep/(RK+1-j))*(triangle[k/3].fluxDensity/triangle[k/3].area);				
			}
			for(int k=0;k<nelem*3;k+=3){	
				triangle[k/3].densityVelocityu = triangle[k/3].oldDensityVelocityu + (triangle[k/3].elemTimeStep/(RK+1-j))*(triangle[k/3].fluxDensityVelocityu/triangle[k/3].area);				
			}
			for(int k=0;k<nelem*3;k+=3){
				triangle[k/3].densityVelocityv = triangle[k/3].oldDensityVelocityv + (triangle[k/3].elemTimeStep/(RK+1-j))*(triangle[k/3].fluxDensityVelocityv/triangle[k/3].area);				
			}
			for(int k=0;k<nelem*3;k+=3){
				triangle[k/3].densityEnergy = triangle[k/3].oldDensityEnergy + (triangle[k/3].elemTimeStep/(RK+1-j))*(triangle[k/3].fluxDensityEnergy/triangle[k/3].area);

			}	//end all element level loops

		//CorrectBoundaryValues(int i, int nelem,int node1,int node2,int node3);

		}//end RK stage loop		
	
	}//end steady state iteration for loop

	/* write output values to file for paraview */
	WriteDatatoFile(npoin,nelem);
	
	cout<<"finised"<<endl;
	return 0;
}
//end main 



double ComputeDensityFlux(int ielem,int nelem){
	double flux,flux1,flux2,flux3,flux4,flux5,flux6;
	flux=flux1=flux2=flux3=flux4=flux5=flux6=0.0;
	int node1,node2,node3,neibor;
	node1=connectivity[ielem];
	node2=connectivity[ielem+1];
	node3=connectivity[ielem+2];
	/* this function calculates the density flux. First it loops through the neighbors obtains the common edges/faces between cells and then calculates the fluxes. Special consideration is given for fluxes of the elements/cells that are at the farfield boundary */
	for (int it=0;it<triangle[ielem/3].neigh.size();it++){
		neibor=triangle[ielem/3].neigh[it];
		if((node1==connectivity[(neibor-1)*3] or node1==connectivity[((neibor-1)*3)+1] or node1==connectivity[((neibor-1)*3)+2]) and (node2==connectivity[(neibor-1)*3] or node2==connectivity[((neibor-1)*3)+1] or node2==connectivity[((neibor-1)*3)+2])){

				flux1=0.5*(triangle[ielem/3].normal12x*(triangle[ielem/3].fluxDensityxContribution+triangle[neibor-1].fluxDensityxContribution))+0.5*(triangle[ielem/3].normal12y*(triangle[ielem/3].fluxDensityyContribution+triangle[neibor-1].fluxDensityyContribution));				
		}
		else if((node2==connectivity[(neibor-1)*3] or node2==connectivity[((neibor-1)*3)+1] or node2==connectivity[((neibor-1)*3)+2]) and (node3==connectivity[(neibor-1)*3] or node3==connectivity[((neibor-1)*3)+1] or node3==connectivity[((neibor-1)*3)+2])){
				flux2=0.5*(triangle[ielem/3].normal23x*(triangle[ielem/3].fluxDensityxContribution+triangle[neibor-1].fluxDensityxContribution))+0.5*(triangle[ielem/3].normal23y*(triangle[ielem/3].fluxDensityyContribution+triangle[neibor-1].fluxDensityyContribution));				
			}
			else if((node1==connectivity[(neibor-1)*3] or node1==connectivity[((neibor-1)*3)+1] or node1==connectivity[((neibor-1)*3)+2]) and (node3==connectivity[(neibor-1)*3] or node3==connectivity[((neibor-1)*3)+1] or node3==connectivity[((neibor-1)*3)+2])){
			flux3=0.5*(triangle[ielem/3].normal31x*(triangle[ielem/3].fluxDensityxContribution+triangle[neibor-1].fluxDensityxContribution))+0.5*(triangle[ielem/3].normal31y*(triangle[ielem/3].fluxDensityyContribution+triangle[neibor-1].fluxDensityyContribution));			
			} 
	}//end neighbor for loop
 
	if (triangle[ielem/3].bctype==4){ /*for a farfield boundary elem */
			
			if((node1==triangle[ielem/3].bface[0] or node1==triangle[ielem/3].bface[1]) and (node2==triangle[ielem/3].bface[0] or node2==triangle[ielem/3].bface[1])){
				
				flux4=0.5*(triangle[ielem/3].normal12x*(triangle[ielem/3].fluxDensityxContribution+ff_densityvelocityx))+0.5*(triangle[ielem/3].normal12y*(triangle[ielem/3].fluxDensityyContribution+ff_densityvelocityy));				
			}
			else if((node2==triangle[ielem/3].bface[0] or node2==triangle[ielem/3].bface[1]) and (node3==triangle[ielem/3].bface[0] or node3==triangle[ielem/3].bface[1])){
				
				flux5=0.5*(triangle[ielem/3].normal23x*(triangle[ielem/3].fluxDensityxContribution+ff_densityvelocityx))+0.5*(triangle[ielem/3].normal23y*(triangle[ielem/3].fluxDensityyContribution+ff_densityvelocityy));				
			}
			else if((node1==triangle[ielem/3].bface[0]or node1==triangle[ielem/3].bface[1]) and (node3==triangle[ielem/3].bface[0] or node3==triangle[ielem/3].bface[1])){
				flux6=0.5*(triangle[ielem/3].normal31x*(triangle[ielem/3].fluxDensityxContribution+ff_densityvelocityx))+0.5*(triangle[ielem/3].normal31y*(triangle[ielem/3].fluxDensityyContribution+ff_densityvelocityy));				
			}
	
	}// end if for farfield boundary

	flux=flux1+flux2+flux3+flux4+flux5+flux6; 	
	return flux;
}

double ComputeDensityVelocityuFlux(int ielem,int nelem){
	double flux;
	flux=0.0;
	int node1,node2,node3,neibor,nei1,nei2,nei3;
	node1=connectivity[ielem];
	node2=connectivity[ielem+1];
	node3=connectivity[ielem+2];
	/* this function calculates the density velocity u flux. First it loops through the neighbors obtains the common edges/faces between cells and then calculates the fluxes. Special consideration is given for fluxes of the elements/cells that are at the farfield boundary and at the wing boundary. */

	for (int it=0;it<triangle[ielem/3].neigh.size();it++){
		neibor=triangle[ielem/3].neigh[it];
		/* for general element */
		if((node1==connectivity[(neibor-1)*3] or node1==connectivity[((neibor-1)*3)+1] or node1==connectivity[((neibor-1)*3)+2]) and (node2==connectivity[(neibor-1)*3] or node2==connectivity[((neibor-1)*3)+1] or node2==connectivity[((neibor-1)*3)+2])){
				flux=flux+0.5*(triangle[ielem/3].normal12x*(triangle[ielem/3].fluxDensityVelocityuxContribution+triangle[neibor-1].fluxDensityVelocityuxContribution))+0.5*(triangle[ielem/3].normal12y*(triangle[ielem/3].fluxDensityVelocityuyContribution+triangle[neibor-1].fluxDensityVelocityuyContribution));
				
		}
		else if((node2==connectivity[(neibor-1)*3] or node2==connectivity[((neibor-1)*3)+1] or node2==connectivity[((neibor-1)*3)+2]) and (node3==connectivity[(neibor-1)*3] or node3==connectivity[((neibor-1)*3)+1] or node3==connectivity[((neibor-1)*3)+2])){
			flux=flux+0.5*(triangle[ielem/3].normal23x*(triangle[ielem/3].fluxDensityVelocityuxContribution+triangle[neibor-1].fluxDensityVelocityuxContribution))+0.5*(triangle[ielem/3].normal23y*(triangle[ielem/3].fluxDensityVelocityuyContribution+triangle[neibor-1].fluxDensityVelocityuyContribution));
			
		}
		else if((node1==connectivity[(neibor-1)*3] or node1==connectivity[((neibor-1)*3)+1] or node1==connectivity[((neibor-1)*3)+2]) and (node3==connectivity[(neibor-1)*3] or node3==connectivity[((neibor-1)*3)+1] or node3==connectivity[((neibor-1)*3)+2])){
			flux=flux+0.5*(triangle[ielem/3].normal31x*(triangle[ielem/3].fluxDensityVelocityuxContribution+triangle[neibor-1].fluxDensityVelocityuxContribution))+0.5*(triangle[ielem/3].normal31y*(triangle[ielem/3].fluxDensityVelocityuyContribution+triangle[neibor-1].fluxDensityVelocityuyContribution));
			
		} 
	}//end neighbor for loop 

	
	if (triangle[ielem/3].bctype==4){ /*for a farfield boundary element */
			//cout<<"FARFIELD "<<endl;
			if((node1==triangle[ielem/3].bface[0] or node1==triangle[ielem/3].bface[1]) and (node2==triangle[ielem/3].bface[0] or node2==triangle[ielem/3].bface[1])){
				flux=flux+0.5*(triangle[ielem/3].normal12x*(triangle[ielem/3].fluxDensityVelocityuxContribution+ff_densityvelocityx*ff_velocityx+pinf))+0.5*(triangle[ielem/3].normal12y*(triangle[ielem/3].fluxDensityVelocityuyContribution+ff_densityvelocityx*ff_velocityy));
				
			}
			else if((node2==triangle[ielem/3].bface[0] or node2==triangle[ielem/3].bface[1]) and (node3==triangle[ielem/3].bface[0] or node3==triangle[ielem/3].bface[1])){
				flux=flux+0.5*(triangle[ielem/3].normal23x*(triangle[ielem/3].fluxDensityVelocityuxContribution+(ff_densityvelocityx*ff_velocityx+pinf)))+0.5*(triangle[ielem/3].normal23y*(triangle[ielem/3].fluxDensityVelocityuyContribution+(ff_densityvelocityx*ff_velocityy)));
			
			}
			else if((node1==triangle[ielem/3].bface[0]or node1==triangle[ielem/3].bface[1]) and (node3==triangle[ielem/3].bface[0] or node3==triangle[ielem/3].bface[1])){
				flux=flux+0.5*(triangle[ielem/3].normal31x*(triangle[ielem/3].fluxDensityVelocityuxContribution+ff_densityvelocityx*ff_velocityx+pinf))+0.5*(triangle[ielem/3].normal31y*(triangle[ielem/3].fluxDensityVelocityuyContribution+ff_densityvelocityx*ff_velocityy));
				
			}

	}// end if for farfield boundary  
	
	if(triangle[ielem/3].bctype==0){ /* for a wing boundary element */
			//cout<<"WING "<<ielem/3<<endl;
			if((node1==triangle[ielem/3].wface[0] or node1==triangle[ielem/3].wface[1]) and (node2==triangle[ielem/3].wface[0] or node2==triangle[ielem/3].wface[1])){

				flux=flux+(triangle[ielem/3].normal12x*(triangle[ielem/3].pressure));
			
			}
			else if((node2==triangle[ielem/3].wface[0] or node2==triangle[ielem/3].wface[1]) and (node3==triangle[ielem/3].wface[0] or node3==triangle[ielem/3].wface[1])){

				flux=flux+(triangle[ielem/3].normal23x*(triangle[ielem/3].pressure));
			
			}
			else if((node1==triangle[ielem/3].wface[0]or node1==triangle[ielem/3].wface[1]) and (node3==triangle[ielem/3].wface[0] or node3==triangle[ielem/3].wface[1])){

				flux=flux+(triangle[ielem/3].normal31x*(triangle[ielem/3].pressure));
				
			}			
			
	}// end if for wing boundary
	
	
	return flux;
}

double ComputeDensityVelocityvFlux(int ielem,int nelem){
	double flux,flux1,flux2,flux3,flux4,flux5,flux6,flux7,flux8,flux9;
	flux=flux1=flux2=flux3=flux4=flux5=flux6=flux7=flux8=flux9=0.0;
	int node1,node2,node3,neibor,nei1,nei2,nei3;
	node1=connectivity[ielem];
	node2=connectivity[ielem+1];
	node3=connectivity[ielem+2];
	/* this function calculates the density velocity v flux. First it loops through the neighbors obtains the common edges/faces between cells and then calculates the fluxes. Special consideration is given for fluxes of the elements/cells that are at the farfield boundary */	

	for (int it=0;it<triangle[ielem/3].neigh.size();it++){
		neibor=triangle[ielem/3].neigh[it];
		if((node1==connectivity[(neibor-1)*3] or node1==connectivity[((neibor-1)*3)+1] or node1==connectivity[((neibor-1)*3)+2]) and (node2==connectivity[(neibor-1)*3] or node2==connectivity[((neibor-1)*3)+1] or node2==connectivity[((neibor-1)*3)+2])){
				flux1=0.5*(triangle[ielem/3].normal12y*(triangle[ielem/3].fluxDensityVelocityvxContribution+triangle[neibor-1].fluxDensityVelocityvxContribution))+0.5*(triangle[ielem/3].normal12x*(triangle[ielem/3].fluxDensityVelocityvyContribution+triangle[neibor-1].fluxDensityVelocityvyContribution));
				
		}
		else if((node2==connectivity[(neibor-1)*3] or node2==connectivity[((neibor-1)*3)+1] or node2==connectivity[((neibor-1)*3)+2]) and (node3==connectivity[(neibor-1)*3] or node3==connectivity[((neibor-1)*3)+1] or node3==connectivity[((neibor-1)*3)+2])){
				flux2=0.5*(triangle[ielem/3].normal23y*(triangle[ielem/3].fluxDensityVelocityvxContribution+triangle[neibor-1].fluxDensityVelocityvxContribution))+0.5*(triangle[ielem/3].normal23x*(triangle[ielem/3].fluxDensityVelocityvyContribution+triangle[neibor-1].fluxDensityVelocityvyContribution));
			
		}
		else if((node1==connectivity[(neibor-1)*3] or node1==connectivity[((neibor-1)*3)+1] or node1==connectivity[((neibor-1)*3)+2]) and (node3==connectivity[(neibor-1)*3] or node3==connectivity[((neibor-1)*3)+1] or node3==connectivity[((neibor-1)*3)+2])){
				flux3=0.5*(triangle[ielem/3].normal31y*(triangle[ielem/3].fluxDensityVelocityvxContribution+triangle[neibor-1].fluxDensityVelocityvxContribution))+0.5*(triangle[ielem/3].normal31x*(triangle[ielem/3].fluxDensityVelocityvyContribution+triangle[neibor-1].fluxDensityVelocityvyContribution));
			
		}
	}//end neighbor for loop 
	
	if(triangle[ielem/3].bctype==0){ /* for a wing boundary element */
			if((node1==triangle[ielem/3].wface[0] or node1==triangle[ielem/3].wface[1]) and (node2==triangle[ielem/3].wface[0] or node2==triangle[ielem/3].wface[1])){

				flux4=(triangle[ielem/3].normal12y*(triangle[ielem/3].pressure));
				
			}
			else if((node2==triangle[ielem/3].wface[0] or node2==triangle[ielem/3].wface[1]) and (node3==triangle[ielem/3].wface[0] or node3==triangle[ielem/3].wface[1])){

				flux5=(triangle[ielem/3].normal23y*(triangle[ielem/3].pressure));
				
			}
			else if((node1==triangle[ielem/3].wface[0]or node1==triangle[ielem/3].wface[1]) and (node3==triangle[ielem/3].wface[0] or node3==triangle[ielem/3].wface[1])){

				flux6=(triangle[ielem/3].normal31y*(triangle[ielem/3].pressure));
				
			}			
			
	}// end if for wing boundary

	if (triangle[ielem/3].bctype==4){ /*for a farfield boundary */
			if((node1==triangle[ielem/3].bface[0] or node1==triangle[ielem/3].bface[1]) and (node2==triangle[ielem/3].bface[0] or node2==triangle[ielem/3].bface[1])){
				flux7=0.5*(triangle[ielem/3].normal12y*(triangle[ielem/3].fluxDensityVelocityvxContribution+ff_densityvelocityy*ff_velocityy+pinf))+0.5*(triangle[ielem/3].normal12x*(triangle[ielem/3].fluxDensityVelocityvyContribution+ff_densityvelocityx*ff_velocityy));
				
			}
			else if((node2==triangle[ielem/3].bface[0] or node2==triangle[ielem/3].bface[1]) and (node3==triangle[ielem/3].bface[0] or node3==triangle[ielem/3].bface[1])){
				flux8=0.5*(triangle[ielem/3].normal23y*(triangle[ielem/3].fluxDensityVelocityvxContribution+ff_densityvelocityy*ff_velocityy+pinf))+0.5*(triangle[ielem/3].normal23x*(triangle[ielem/3].fluxDensityVelocityvyContribution+ff_densityvelocityx*ff_velocityy));
			
			}
			else if((node1==triangle[ielem/3].bface[0]or node1==triangle[ielem/3].bface[1]) and (node3==triangle[ielem/3].bface[0] or node3==triangle[ielem/3].bface[1])){
				flux9=0.5*(triangle[ielem/3].normal31y*(triangle[ielem/3].fluxDensityVelocityvxContribution+ff_densityvelocityy*ff_velocityy+pinf))+0.5*(triangle[ielem/3].normal31x*(triangle[ielem/3].fluxDensityVelocityvyContribution+ff_densityvelocityx*ff_velocityy));
				
			}

	}// end if for farfield boundary 
	flux=flux1+flux2+flux3+flux4+flux5+flux6+flux7+flux8+flux9;
	
	return flux;
}


double ComputeDensityEnergyFlux(int ielem,int nelem){

	double flux,flux1,flux2,flux3,flux4,flux5,flux6,flux7,flux8,flux9;
	flux=flux1=flux2=flux3=flux4=flux5=flux6=flux7=flux8=flux9=0.0;
	int node1,node2,node3,neibor;
	node1=connectivity[ielem];
	node2=connectivity[ielem+1];
	node3=connectivity[ielem+2];
	
	/* this function calculates the density energy flux. First it loops through the neighbors obtains the common edges/faces between cells and then calculates the fluxes. Special consideration is given for fluxes of the elements/cells that are at the farfield boundary */
	for (int it=0;it<triangle[ielem/3].neigh.size();it++){
		neibor=triangle[ielem/3].neigh[it];
		if((node1==connectivity[(neibor-1)*3] or node1==connectivity[((neibor-1)*3)+1] or node1==connectivity[((neibor-1)*3)+2]) and (node2==connectivity[(neibor-1)*3] or node2==connectivity[((neibor-1)*3)+1] or node2==connectivity[((neibor-1)*3)+2])){

				flux1=0.5*(triangle[ielem/3].normal12x*(triangle[ielem/3].fluxDensityEnergyxContribution+triangle[neibor-1].fluxDensityEnergyxContribution))+0.5*(triangle[ielem/3].normal12y*(triangle[ielem/3].fluxDensityEnergyyContribution+triangle[neibor-1].fluxDensityEnergyyContribution));
				
		}
		else if((node2==connectivity[(neibor-1)*3] or node2==connectivity[((neibor-1)*3)+1] or node2==connectivity[((neibor-1)*3)+2]) and (node3==connectivity[(neibor-1)*3] or node3==connectivity[((neibor-1)*3)+1] or node3==connectivity[((neibor-1)*3)+2])){
				flux2=0.5*(triangle[ielem/3].normal23x*(triangle[ielem/3].fluxDensityEnergyxContribution+triangle[neibor-1].fluxDensityEnergyxContribution))+0.5*(triangle[ielem/3].normal23y*(triangle[ielem/3].fluxDensityEnergyyContribution+triangle[neibor-1].fluxDensityEnergyyContribution));
				
			}
			else if((node1==connectivity[(neibor-1)*3] or node1==connectivity[((neibor-1)*3)+1] or node1==connectivity[((neibor-1)*3)+2]) and (node3==connectivity[(neibor-1)*3] or node3==connectivity[((neibor-1)*3)+1] or node3==connectivity[((neibor-1)*3)+2])){
			flux3=0.5*(triangle[ielem/3].normal31x*(triangle[ielem/3].fluxDensityEnergyxContribution+triangle[neibor-1].fluxDensityEnergyxContribution))+0.5*(triangle[ielem/3].normal31y*(triangle[ielem/3].fluxDensityEnergyyContribution+triangle[neibor-1].fluxDensityEnergyyContribution));
		
			} 
	}//end neighbor for loop
 
	if (triangle[ielem/3].bctype==4){ /*for a farfield boundary elem */
			if((node1==triangle[ielem/3].bface[0] or node1==triangle[ielem/3].bface[1]) and (node2==triangle[ielem/3].bface[0] or node2==triangle[ielem/3].bface[1])){
				
				flux4=0.5*(triangle[ielem/3].normal12x*(triangle[ielem/3].fluxDensityEnergyxContribution+ff_densityEnergyxContribution))+0.5*(triangle[ielem/3].normal12y*(triangle[ielem/3].fluxDensityEnergyyContribution+ff_densityEnergyyContribution));
				
			}
			else if((node2==triangle[ielem/3].bface[0] or node2==triangle[ielem/3].bface[1]) and (node3==triangle[ielem/3].bface[0] or node3==triangle[ielem/3].bface[1])){
				
				flux5=0.5*(triangle[ielem/3].normal23x*(triangle[ielem/3].fluxDensityEnergyxContribution+ff_densityEnergyxContribution))+0.5*(triangle[ielem/3].normal23y*(triangle[ielem/3].fluxDensityEnergyyContribution+ff_densityEnergyyContribution));
				
			}
			else if((node1==triangle[ielem/3].bface[0]or node1==triangle[ielem/3].bface[1]) and (node3==triangle[ielem/3].bface[0] or node3==triangle[ielem/3].bface[1])){
				flux6=0.5*(triangle[ielem/3].normal31x*(triangle[ielem/3].fluxDensityEnergyxContribution+ff_densityEnergyxContribution))+0.5*(triangle[ielem/3].normal31y*(triangle[ielem/3].fluxDensityEnergyyContribution+ff_densityEnergyyContribution));
				
			}
	
	}// end if for farfield boundary 
	flux=flux1+flux2+flux3+flux4+flux5+flux6;
	

	return flux;
	
}




double Areacalc(int a,int b,int c){
	/* this function calculates the area of the triangles in unstructured grid */
	double cross,x1,y1,z1,x2,y2,z2,x3,y3,z3=0.0;
	double ab,ac,bc,s=0.0;
	x1=cordinate[(a-1)*ndim];
	y1=cordinate[(a-1)*ndim+1];
	x2=cordinate[(b-1)*ndim];
	y2=cordinate[(b-1)*ndim+1];
	x3=cordinate[(c-1)*ndim];
	y3=cordinate[(c-1)*ndim+1];
	ab=(x2-x1)*(y3-y1);
	ac=(x3-x1)*(y2-y1);
	cross=(fabs(ab-ac)/2);

return cross;
}

double Centroidcalcx(int a,int b,int c){
	/*this sub calculates the x coordinate of the centroid of the triangles. (not used ) */
	double cross,x1,y1,z1,x2,y2,z2,x3,y3,z3;
	double xn,yn,zn;
	x1=cordinate[(a-1)*ndim];
	y1=cordinate[(a-1)*ndim+1];
	x2=cordinate[(b-1)*ndim];
	y2=cordinate[(b-1)*ndim+1];
	x3=cordinate[(c-1)*ndim];
	y3=cordinate[(c-1)*ndim+1];
	
	
	xn=(x1+x2+x3)/3.0;

	
	
	return xn;
}

double Centroidcalcy(int a,int b,int c){
	/*this sub calculates the y coordinate of the centroid of the triangles. (not used ) */
	double cross,x1,y1,z1,x2,y2,z2,x3,y3,z3;
	double xn,yn,zn;
	x1=cordinate[(a-1)*ndim];
	y1=cordinate[(a-1)*ndim+1];
	x2=cordinate[(b-1)*ndim];
	y2=cordinate[(b-1)*ndim+1];
	x3=cordinate[(c-1)*ndim];
	y3=cordinate[(c-1)*ndim+1];
	
	
	yn=(y1+y2+y3)/3.0;
	
	
	return yn;
}


int FindNeighbor(int ielemNode,int nelem){
/* this function calculates the neighbors of an elements (element surrounding elements) through brute force approach. stores them
in structure */	
	
	int count1=0;
	int node1,node2,node3=0;
	ofstream myfile1;

	node1=connectivity[ielemNode];
	node2=connectivity[ielemNode+1];
	node3=connectivity[ielemNode+2];
	
	count1=ielemNode;
	
	for(int ik=0;ik<nelem*3;ik+=3){
		
		if(ik==ielemNode) continue;
		if((node1==connectivity[ik] and node2==connectivity[ik+1]) or (node1==connectivity[ik]  and node2==connectivity[ik+2]) or (node1==connectivity[ik+1] and node2==connectivity[ik]) or (node1==connectivity[ik+1] and node2==connectivity[ik+2]) or (node1==connectivity[ik+2] and node2==connectivity[ik]) or (node1==connectivity[ik+2] and node2==connectivity[ik+1])){
			triangle[ielemNode/3].neigh.push_back(ik/3+1);
		}
		
		if((node1==connectivity[ik] and node3==connectivity[ik+1]) or (node1==connectivity[ik]  and node3==connectivity[ik+2]) or (node1==connectivity[ik+1] and node3==connectivity[ik]) or (node1==connectivity[ik+1] and node3==connectivity[ik+2]) or (node1==connectivity[ik+2] and node3==connectivity[ik]) or (node1==connectivity[ik+2] and node3==connectivity[ik+1])){
			triangle[ielemNode/3].neigh.push_back(ik/3+1);
		}
		
		 if((node3==connectivity[ik] and node2==connectivity[ik+1]) or (node3==connectivity[ik]  and node2==connectivity[ik+2]) or (node3==connectivity[ik+1] and node2==connectivity[ik]) or (node3==connectivity[ik+1] and node2==connectivity[ik+2]) or (node3==connectivity[ik+2] and node2==connectivity[ik]) or (node3==connectivity[ik+2] and node2==connectivity[ik+1])){
			triangle[ielemNode/3].neigh.push_back(ik/3+1);
		}
		
		count1++;
		
		
	}
	return 0;
}


int ComputeNormal(int ielemNode,int nelem){

	double dx1,dx2,dx3,dy1,dy2,dy3,drl1,drl2,drl3;
	int node1,node2,node3;
	node1=node2=node3=0;
	dx1=dx2=dx3=dy1=dy2=dy3=drl1=drl2=drl3=0.0;
	node1=connectivity[ielemNode];
	node2=connectivity[ielemNode+1];
	node3=connectivity[ielemNode+2];
	/* this function calculates the Normals for each faces/edges of the triangles */
	
	dy1=cordinate[(node1-1)*ndim]-cordinate[(node2-1)*ndim];
	dy2=cordinate[(node2-1)*ndim]-cordinate[(node3-1)*ndim];
	dy3=cordinate[(node3-1)*ndim]-cordinate[(node1-1)*ndim];
	dx1=cordinate[(node2-1)*ndim+1]-cordinate[(node1-1)*ndim+1];
	dx2=cordinate[(node3-1)*ndim+1]-cordinate[(node2-1)*ndim+1];
	dx3=cordinate[(node1-1)*ndim+1]-cordinate[(node3-1)*ndim+1];
	drl1=sqrt(dx1*dx1+dy1*dy1);
	drl2=sqrt(dx2*dx2+dy2*dy2);
	drl3=sqrt(dx3*dx3+dy3*dy3);
	triangle[ielemNode/3].normal12x=-1*dx1;
	triangle[ielemNode/3].normal12y=-1*dy1;
	triangle[ielemNode/3].normal23x=-1*dx2;
	triangle[ielemNode/3].normal23y=-1*dy2;
	triangle[ielemNode/3].normal31x=-1*dx3;
	triangle[ielemNode/3].normal31y=-1*dy3;
	triangle[ielemNode/3].normal12unit=drl1;
	triangle[ielemNode/3].normal23unit=drl2;
	triangle[ielemNode/3].normal31unit=drl3;

	return 0;
}

int calcInwardNormal(int ielem,int a1, int b1){
	/* this function calculates the inwards Normals for each faces/edges of the triangles (not used) */
	double normal,x1,x2,y1,y2,dx1,dy1,drl1;
	x1=cordinate[(a1-1)*ndim];
	y1=cordinate[(a1-1)*ndim+1];
	x2=cordinate[(b1-1)*ndim];
	y2=cordinate[(b1-1)*ndim+1];
	dx1=cordinate[(a1-1)*ndim+1]-cordinate[(b1-1)*ndim+1];
	dy1=cordinate[(b1-1)*ndim]-cordinate[(a1-1)*ndim];
	drl1=sqrt(dx1*dx1+dy1*dy1);
	triangle[ielem/3].inwardnormalx=dx1/drl1;
	triangle[ielem/3].inwardnormaly=dy1/drl1;

	return 0;
}


int CorrectBoundaryValues(int i, int nelem,int node1,int node2,int node3){

	double dot_product,my_speed,my_mach_number,normalx,normaly,ff_velocity_n,my_velocity_n,ff_speedofsound,boundary_energy;
	double avgSpeedofSound,avgSoundDensity,boundary_pressure,boundary_density,boundary_vn,boundary_velocity_x,boundary_velocity_y;
	double boundaryCorrectedV,boundary_vtx,boundary_vty,boundary_vt,boundary_speed,tangentx,tangenty,unit_n;
	
	dot_product = 0.0;
	my_speed = 0.0;
	normalx=0.0;
	normaly=0.0;
	ff_speedofsound=sqrt(gamma*pinf/rhoinf);

	if((node1==triangle[i/3].bface[0] or node1==triangle[i/3].bface[1]) and (node2==triangle[i/3].bface[0] or node2==triangle[i/3].bface[1])){
			unit_n=sqrt(triangle[i/3].normal12x*triangle[i/3].normal12x	+triangle[i/3].normal12y*triangle[i/3].normal12y);	
			dot_product=(-1*triangle[i/3].normal12x)*triangle[i/3].velocityu*unit_n+(-1*triangle[i/3].normal12y)*triangle[i/3].velocityv*unit_n;
			normalx=-1*triangle[i/3].normal12x;
			normaly=-1*triangle[i/3].normal12y;
	}
	else if((node2==triangle[i/3].bface[0] or node2==triangle[i/3].bface[1]) and (node3==triangle[i/3].bface[0] or node3==triangle[i/3].bface[1])){
			unit_n=sqrt(triangle[i/3].normal23x*triangle[i/3].normal23x	+triangle[i/3].normal23y*triangle[i/3].normal23y);
			dot_product=(-1*triangle[i/3].normal23x)*triangle[i/3].velocityu*unit_n+(-1*triangle[i/3].normal23y)*triangle[i/3].velocityv*unit_n;
			normalx=-1*triangle[i/3].normal23x;
			normaly=-1*triangle[i/3].normal23y;
	}
	else if((node1==triangle[i/3].bface[0]or node1==triangle[i/3].bface[1]) and (node3==triangle[i/3].bface[0] or node3==triangle[i/3].bface[1])){
			unit_n=sqrt(triangle[i/3].normal31x*triangle[i/3].normal31x	+triangle[i/3].normal31y*triangle[i/3].normal31y);
			dot_product=(-1*triangle[i/3].normal31x)*triangle[i/3].velocityu*unit_n+(-1*triangle[i/3].normal31y)*triangle[i/3].velocityv*unit_n;
			normalx=-1*triangle[i/3].normal31x;
			normaly=-1*triangle[i/3].normal31y;	
	}
	//cout<<"  elem "<<i/3<<"\tbefore flux density\t"<<triangle[i/3].velocityu<<endl;
	tangentx=-1*normaly/unit_n;
	tangenty=-1*normalx/unit_n;			
	/* determine whether this i sub or supersonic flow */
	my_mach_number = triangle[i/3].speed / triangle[i/3].speedofsound;
	if(my_mach_number < 1.0){	 /* subsonic */
		ff_velocity_n = ff_velocityx*normalx + ff_velocityy*normaly;
		my_velocity_n = dot_product;
		avgSpeedofSound = 0.5*(triangle[i/3].speedofsound + ff_speedofsound);
		avgSoundDensity = 0.5*(rhoinf*ff_speedofsound + triangle[i/3].density*triangle[i/3].speedofsound);
		boundary_vtx=(triangle[i/3].velocityu-(my_velocity_n*normalx/unit_n))/(tangentx);
		boundary_vty=(triangle[i/3].velocityv-(my_velocity_n*normaly/unit_n))/tangenty;
		boundary_vt=sqrt(boundary_vtx*boundary_vtx+boundary_vty*boundary_vty);
		//cout<<"at vt calc "<<triangle[i/3].speed<<"\t"<<my_velocity_n<<endl;		
		if(dot_product > 0){	/* subsonic inflow */
			//cout<<"inside subsonic inflow "<<endl;		
			boundary_pressure = 0.5*(pinf+ triangle[i/3].pressure+ (ff_velocity_n-my_velocity_n)*avgSoundDensity);
			boundary_density = rhoinf + ((boundary_pressure-pinf)/(avgSpeedofSound*avgSpeedofSound));
			boundary_vn = ff_velocity_n + (pinf - boundary_pressure)/avgSoundDensity;
			//cout<<i/3<<"  "<<boundary_pressure<<"\t"<<boundary_density<<"\t"<<boundary_vn<<"\t"<<boundary_vt<<endl;
		}
		else{ /* subsonic outflow */
			//cout<<"inside subsonic outflow "<<endl;		
			boundary_pressure = pinf;
			boundary_density = triangle[i/3].density + (pinf - triangle[i/3].pressure)/(avgSpeedofSound*avgSpeedofSound);
			boundary_vn = ff_velocity_n + (pinf - triangle[i/3].pressure)/avgSoundDensity;
			//cout<<i/3<<"  "<<boundary_pressure<<"\t"<<boundary_density<<"\t"<<boundary_vn<<"\t"<<boundary_vt<<endl;
		}
		boundaryCorrectedV=sqrt(boundary_vt*boundary_vt+boundary_vn*boundary_vn);
		boundary_velocity_x = boundaryCorrectedV*cos(alpha);
		boundary_velocity_y = boundaryCorrectedV*sin(alpha);
		//cout<<"corrected v "<<boundaryCorrectedV<<endl;
	}
	else { /* supersonic */
		if(dot_product > 0){	/* inflow */
			/* all quantities come from infinity */
			//cout<<"inside supersonic inflow "<<endl;
			boundary_pressure = pinf;
			boundary_density = rhoinf;
			boundary_velocity_x = ff_velocityx;
			boundary_velocity_y = ff_velocityy;
			//cout<<i/3<<"  "<<boundary_pressure<<"\t"<<boundary_density<<"\t"<<boundary_velocity_x<<"   "<<boundary_velocity_y<<endl;
		}
		else {/* outflow */
			//cout<<"inside supersonic outflow "<<endl;
			boundary_pressure = triangle[i/3].pressure;
			boundary_density = triangle[i/3].density;
			boundary_velocity_x = triangle[i/3].velocityu;
			boundary_velocity_y = triangle[i/3].velocityv;
			//cout<<i/3<<"  "<<boundary_velocity_x<<"   "<<boundary_velocity_y<<endl;	
		}	
	}//end supersonic else
				
	boundary_energy = (boundary_pressure / (boundary_density*gamma-1)) + 0.5*(boundary_velocity_x*boundary_velocity_x + boundary_velocity_y*boundary_velocity_y);
	//cout<<"  elem "<<i/3<<"\tbefore flux density\t"<<triangle[i/3].densityVelocityv<<endl;			
	/* set the variables to their corrected values */
	triangle[i/3].pressure=boundary_pressure;
	triangle[i/3].density= boundary_density;
	triangle[i/3].densityVelocityu = boundary_density*boundary_velocity_x;
	triangle[i/3].densityVelocityv = boundary_density*boundary_velocity_y;
	triangle[i/3].densityEnergy = boundary_density*boundary_energy;
	//cout<<"  elem "<<i/3<<"\tafter flux density\t"<<triangle[i/3].densityVelocityv<<endl;		

	return 0;
}
	



double ArtificialViscosityDensity(int ielem,int nelem){
	/* this sub calculates the artificial viscosity due to density flux */
	double artViscosity,lamdamax,smoothparam;
	artViscosity=0.0;
	smoothparam=0.08;
	int node1,node2,node3,neibor;
	node1=connectivity[ielem];
	node2=connectivity[ielem+1];
	node3=connectivity[ielem+2];	
	for (int it=0;it<triangle[ielem/3].neigh.size();it++){
		neibor=triangle[ielem/3].neigh[it];
		//if(triangle[ielem/3].bctype !=0 and triangle[ielem/3].bctype !=4){
			lamdamax=0.5*(triangle[ielem/3].speed+triangle[ielem/3].speedofsound+triangle[neibor-1].speed+triangle[neibor-1].speedofsound);
	
			
			if((node1==connectivity[(neibor-1)*3] or node1==connectivity[((neibor-1)*3)+1] or node1==connectivity[((neibor-1)*3)+2]) and (node2==connectivity[(neibor-1)*3] or node2==connectivity[((neibor-1)*3)+1] or node2==connectivity[((neibor-1)*3)+2])){
				artViscosity+=triangle[ielem/3].normal12unit*smoothparam*lamdamax*(triangle[ielem/3].density-triangle[neibor-1].density);							
				
			}
			else if((node2==connectivity[(neibor-1)*3] or node2==connectivity[((neibor-1)*3)+1] or node2==connectivity[((neibor-1)*3)+2]) and (node3==connectivity[(neibor-1)*3] or node3==connectivity[((neibor-1)*3)+1] or node3==connectivity[((neibor-1)*3)+2])){
				artViscosity+=triangle[ielem/3].normal23unit*smoothparam*lamdamax*(triangle[ielem/3].density-triangle[neibor-1].density);
				
			}
			else if((node1==connectivity[(neibor-1)*3] or node1==connectivity[((neibor-1)*3)+1] or node1==connectivity[((neibor-1)*3)+2]) and (node3==connectivity[(neibor-1)*3] or node3==connectivity[((neibor-1)*3)+1] or node3==connectivity[((neibor-1)*3)+2])){

				artViscosity+=triangle[ielem/3].normal31unit*smoothparam*lamdamax*(triangle[ielem/3].density-triangle[neibor-1].density);
				
			}
		//}
		
	}
	
	return artViscosity;
}

double ArtificialViscosityDensityVelocityu(int ielem,int nelem){
	/* this sub calculates the artificial viscosity due to density velocity u flux */
	double artViscosity,lamdamax,smoothparam;
	artViscosity=0.0;
	smoothparam=0.08;
	int node1,node2,node3,neibor;
	node1=connectivity[ielem];
	node2=connectivity[ielem+1];
	node3=connectivity[ielem+2];	
	for (int it=0;it<triangle[ielem/3].neigh.size();it++){
		neibor=triangle[ielem/3].neigh[it];
		//if(triangle[ielem/3].bctype !=0 and triangle[ielem/3].bctype !=4){

			lamdamax=0.5*(triangle[ielem/3].speed+triangle[ielem/3].speedofsound+triangle[neibor-1].speed+triangle[neibor-1].speedofsound);

			
			if((node1==connectivity[(neibor-1)*3] or node1==connectivity[((neibor-1)*3)+1] or node1==connectivity[((neibor-1)*3)+2]) and (node2==connectivity[(neibor-1)*3] or node2==connectivity[((neibor-1)*3)+1] or node2==connectivity[((neibor-1)*3)+2])){

				artViscosity+=triangle[ielem/3].normal12unit*smoothparam*lamdamax*(triangle[ielem/3].densityVelocityu-triangle[neibor-1].densityVelocityu);
			}
			else if((node2==connectivity[(neibor-1)*3] or node2==connectivity[((neibor-1)*3)+1] or node2==connectivity[((neibor-1)*3)+2]) and (node3==connectivity[(neibor-1)*3] or node3==connectivity[((neibor-1)*3)+1] or node3==connectivity[((neibor-1)*3)+2])){

				artViscosity+=triangle[ielem/3].normal23unit*smoothparam*lamdamax*(triangle[ielem/3].densityVelocityu-triangle[neibor-1].densityVelocityu);
			}
			else if((node1==connectivity[(neibor-1)*3] or node1==connectivity[((neibor-1)*3)+1] or node1==connectivity[((neibor-1)*3)+2]) and (node3==connectivity[(neibor-1)*3] or node3==connectivity[((neibor-1)*3)+1] or node3==connectivity[((neibor-1)*3)+2])){

				artViscosity+=triangle[ielem/3].normal31unit*smoothparam*lamdamax*(triangle[ielem/3].densityVelocityu-triangle[neibor-1].densityVelocityu);
			}
		//}
		
	}
		
	return artViscosity;
}

double ArtificialViscosityDensityVelocityv(int ielem,int nelem){
	/* this sub calculates the artificial viscosity due to density velocity v flux */
	double artViscosity,lamdamax,smoothparam;
	artViscosity=0.0;
	smoothparam=0.08;
	int node1,node2,node3,neibor;
	node1=connectivity[ielem];
	node2=connectivity[ielem+1];
	node3=connectivity[ielem+2];	
	for (int it=0;it<triangle[ielem/3].neigh.size();it++){
		neibor=triangle[ielem/3].neigh[it];
		//if(triangle[ielem/3].bctype !=0 and triangle[ielem/3].bctype !=4){

			lamdamax=0.5*(triangle[ielem/3].speed+triangle[ielem/3].speedofsound+triangle[neibor-1].speed+triangle[neibor-1].speedofsound);
	
			
			if((node1==connectivity[(neibor-1)*3] or node1==connectivity[((neibor-1)*3)+1] or node1==connectivity[((neibor-1)*3)+2]) and (node2==connectivity[(neibor-1)*3] or node2==connectivity[((neibor-1)*3)+1] or node2==connectivity[((neibor-1)*3)+2])){

				artViscosity+=triangle[ielem/3].normal12unit*smoothparam*lamdamax*(triangle[ielem/3].densityVelocityv-triangle[neibor-1].densityVelocityv);
			}
			else if((node2==connectivity[(neibor-1)*3] or node2==connectivity[((neibor-1)*3)+1] or node2==connectivity[((neibor-1)*3)+2]) and (node3==connectivity[(neibor-1)*3] or node3==connectivity[((neibor-1)*3)+1] or node3==connectivity[((neibor-1)*3)+2])){

				artViscosity+=triangle[ielem/3].normal23unit*smoothparam*lamdamax*(triangle[ielem/3].densityVelocityv-triangle[neibor-1].densityVelocityv);
			}
			else if((node1==connectivity[(neibor-1)*3] or node1==connectivity[((neibor-1)*3)+1] or node1==connectivity[((neibor-1)*3)+2]) and (node3==connectivity[(neibor-1)*3] or node3==connectivity[((neibor-1)*3)+1] or node3==connectivity[((neibor-1)*3)+2])){

				artViscosity+=triangle[ielem/3].normal31unit*smoothparam*lamdamax*(triangle[ielem/3].densityVelocityv-triangle[neibor-1].densityVelocityv);
			}
		//}
		
	}
		
	return artViscosity;
}

double ArtificialViscosityDensityEnergy(int ielem,int nelem){
	/* this sub calculates the artificial viscosity due to density energy flux */

	double artViscosity,lamdamax,smoothparam;
	artViscosity=0.0;
	smoothparam=0.2;
	int node1,node2,node3,neibor;
	node1=connectivity[ielem];
	node2=connectivity[ielem+1];
	node3=connectivity[ielem+2];	
	for (int it=0;it<triangle[ielem/3].neigh.size();it++){
		neibor=triangle[ielem/3].neigh[it];
		//if(triangle[ielem/3].bctype !=0 and triangle[ielem/3].bctype !=4){

			lamdamax=0.5*(triangle[ielem/3].speed+triangle[ielem/3].speedofsound+triangle[neibor-1].speed+triangle[neibor-1].speedofsound);
		
			if((node1==connectivity[(neibor-1)*3] or node1==connectivity[((neibor-1)*3)+1] or node1==connectivity[((neibor-1)*3)+2]) and (node2==connectivity[(neibor-1)*3] or node2==connectivity[((neibor-1)*3)+1] or node2==connectivity[((neibor-1)*3)+2])){

				artViscosity+=triangle[ielem/3].normal12unit*smoothparam*lamdamax*(triangle[ielem/3].densityEnergy-triangle[neibor-1].densityEnergy);
			}
			else if((node2==connectivity[(neibor-1)*3] or node2==connectivity[((neibor-1)*3)+1] or node2==connectivity[((neibor-1)*3)+2]) and (node3==connectivity[(neibor-1)*3] or node3==connectivity[((neibor-1)*3)+1] or node3==connectivity[((neibor-1)*3)+2])){

				artViscosity+=triangle[ielem/3].normal23unit*smoothparam*lamdamax*(triangle[ielem/3].densityEnergy-triangle[neibor-1].densityEnergy);
			}
			else if((node1==connectivity[(neibor-1)*3] or node1==connectivity[((neibor-1)*3)+1] or node1==connectivity[((neibor-1)*3)+2]) and (node3==connectivity[(neibor-1)*3] or node3==connectivity[((neibor-1)*3)+1] or node3==connectivity[((neibor-1)*3)+2])){

				artViscosity+=triangle[ielem/3].normal31unit*smoothparam*lamdamax*(triangle[ielem/3].densityEnergy-triangle[neibor-1].densityEnergy);
			}
		//}
		
	}
	
	
	return artViscosity;
}

int ReadMeshFile(int npoin, int nelem,int bcnodes){
	/* this sub reads the mesh and formulates coordinate,connectivity and bnode arrays */
	int i,j,a,b,c;
	double x1,y1,x2,y2,x3,y3;
	
	double n1,n2,n3,p;
	int count = 0;
	i =0;
	std::string buffer;
	
	int p1=0;
	int p2=0;
	int p3=0;
	int q1;
	
	ifstream fin("naca0012.mesh.fine");
	
	while (!fin.eof()) {
		getline(fin, buffer, '\n');
		count++;  
		
		if(count>=10 and count <=(nelem+9)){
			i=buffer.find_first_not_of(' ');
			n1 = convertToDouble (buffer.substr(8,14));
			n2=convertToDouble(buffer.substr(16,21));
			n3=convertToDouble(buffer.substr(24,30));
		    connectivity[p1]=n1;
			connectivity[p1+1]=n2;
			connectivity[p1+2]=n3;
			p1+=3;
		}
		
		else if(count>=11154 and count<16850){
			n1 = convertToDouble (buffer.substr(8,29));
			i=buffer.find_first_not_of(' ',35);
			n2=convertToDouble(buffer.substr(34,60));
			//cout<<n1<<"   "<<n2<<endl;
			cordinate[p2]=n1;
			cordinate[p2+1]=n2;
			//cout<<cordinate[p2]<<"    "<<cordinate[p2+1]<<endl;
			p2+=2;				
		}
		else if(count>22547 and count<=(22547+bcnodes)){
			n1 = convertToDouble (buffer.substr(0,8));
			i=buffer.find_first_not_of(' ',35);
			n2=convertToDouble(buffer.substr(8,9));
			bnode[p3]=n1;
			bnode[p3+1]=n2;	
			//cout<<n1<<"   "<<n2<<endl;
			p3+=2;			
		}
	  
	} //end while eof
	
	return 0;
}


int CalcBoundaryFaces(int npoin, int nelem,int bcnodes){
	int j1,a1,b1,k,l,i,j;
	/*this sub obtaines which faces are at the boundary of an element */
	for(i=0;i<nelem*3;i+=3){		
		//if(triangle[i/3].bctype==4){						
			for(j=0;j<3;j++){
				j1=(j+1)%3;
				a1=connectivity[3*i/3+j];
				b1=connectivity[3*i/3+j1];
				for(k=0;k<bcnodes*2;k+=2){
					if(a1==bnode[k] and bnode[k+1]==4){
						for(l=0;l<bcnodes*2;l+=2){
							if(b1==bnode[l] and bnode[l+1]==4){
								triangle[i/3].bface[0]=a1;
								triangle[i/3].bface[1]=b1;
								triangle[i/3].bctype=4;
								calcInwardNormal(i,a1,b1);
								break;						
							}						
						}					
					}					
				}			
			}
		//}
	}
	
	for(i=0;i<nelem*3;i+=3){		
		//if(triangle[i/3].bctype==0){						
			for(j=0;j<3;j++){
				j1=(j+1)%3;
				a1=connectivity[3*i/3+j];
				b1=connectivity[3*i/3+j1];
				for(k=0;k<bcnodes*2;k+=2){
					if(a1==bnode[k] and bnode[k+1]==0){
						for(l=0;l<bcnodes*2;l+=2){
							if(b1==bnode[l]and bnode[l+1]==0){
								triangle[i/3].wface[0]=a1;
								triangle[i/3].wface[1]=b1;
								triangle[i/3].bctype=0;
								//calcInwardNormal(i,a1,b1);
								break;						
							}						
						}					
					}					
				}			
			}
		//}
	}
	


	return 0;
}

int WriteDatatoFile(int npoin,int nelem){
	/* this sub writes the output values to file */
	int i,j;	
	ofstream myfile;
	myfile.open("outputs4deg.vtk");
	//int esup1[npoin*nelem];
	int *esup1=(int*)malloc((npoin*nelem)*sizeof(int));
	//int esup2[npoin+1];
	int *esup2=(int*)malloc((npoin+1)*sizeof(int));
	double avgpressure;
	avgpressure=0.0;
	int ipoi1;
	int start,end;
	double point_values[npoin];	
	int ipoin, istor;
	
	for(j = 0; j < nelem; j++){
		for(i = 0; i < 3; i++){
			ipoi1=connectivity[3*j+i];
			esup2[ipoi1]=esup2[ipoi1]+1;
		}		
	}
	
	for(i = 1; i < npoin+1; i++){
		esup2[i] += esup2[i-1];		
	}

	
	
	for(j = 0; j < nelem; j++){
		for(i = 0; i < 3; i++)		{
			ipoin=(connectivity[3*j+i])-1;
			istor=esup2[ipoin]+1;
			esup2[ipoin]=istor;
			esup1[istor-1]=j+1;			
		}
	}
	
	for(i = npoin+1; i > 0; i--)
	{
		esup2[i] = esup2[i-1];
	}
	esup2[0] = 0;

	
	myfile<<"# vtk DataFile Version 2.0"<<endl<<endl;
	myfile<<"ASCII"<<endl<<"DATASET UNSTRUCTURED_GRID"<<endl;
	myfile<<endl<<"POINTS"<<"    "<<npoin<<"    "<<"float"<<endl;
	for(i=0;i<npoin*2;i+=2){
		myfile<<cordinate[i]<<"    "<<cordinate[i+1]<<"    "<<0.0<<endl;	
	}
	myfile<<endl<<endl;
	myfile<<"CELLS"<<"    "<<nelem<<"    "<<4*nelem<<endl;
	for(i=0;i<nelem*3;i+=3){
		myfile<<3.0<<"    "<<connectivity[i]-1<<"    "<<connectivity[i+1]-1<<"    "<<connectivity[i+2]-1<<endl;	
	}
	myfile<<endl;
	myfile<<"CELL_TYPES"<<"    "<<nelem<<endl;
	
	for(i=0;i<nelem;i++){
		myfile<<5<<endl;
		
	}


	myfile<<endl<<"POINT_DATA"<<"    "<<npoin<<endl<<"VECTORS"<<"    "<<"Velocity float  "<<endl;
	double vxvalues[npoin],vyvalues[npoin];	
	for(j = 0; j < npoin; j++){				
		start = esup2[j], end = esup2[j+1];
		vxvalues[j] = 0.0;
		vyvalues[j]=0.0;		
		for(i = start; i < end; i++){						
			vxvalues[j] += triangle[esup1[i]-1].densityVelocityu/triangle[esup1[i]-1].density;
			vyvalues[j]+= triangle[esup1[i]-1].densityVelocityv/triangle[esup1[i]-1].density; 
			//cout<<"point "<<j<<"\tesup"<<esup1[i]-1<<"\tpressure\t"<<triangle[esup1[i]-1].pressure<<endl;
			
		}
		vxvalues[j] /= (double)(end-start);
		vyvalues[j]/=(double)(end-start);
		cout<<"average value "<<point_values[j]<<endl;
		myfile<<vxvalues[j]<<"    "<<vyvalues[j]<<"    "<<0.0<<endl;		
		
	}	

	myfile<<endl<<"SCALARS"<<"  "<<"Pressure  float "<<1<<endl<<"LOOKUP_TABLE default"<<endl;
	
	for(j = 0; j < npoin; j++){				
		start = esup2[j], end = esup2[j+1];
		point_values[j] = 0.0;		
		for(i = start; i < end; i++){						
			point_values[j] += triangle[esup1[i]-1].pressure;  
			cout<<"point "<<j<<"\tesup"<<esup1[i]-1<<"\tpressure\t"<<triangle[esup1[i]-1].pressure<<endl;
			
		}
		point_values[j] /= (double)(end-start);
		cout<<"average value "<<point_values[j]<<endl;
		myfile<<point_values[j]<<endl;		
		
	}

	myfile<<endl<<"SCALARS"<<"  "<<"Mach_no  float "<<1<<endl<<"LOOKUP_TABLE default"<<endl;
	double mach[npoin];
	for(j = 0; j < npoin; j++){				
		start = esup2[j], end = esup2[j+1];
		mach[j] = 0.0;		
		for(i = start; i < end; i++){						
			mach[j] += triangle[esup1[i]-1].speed/triangle[esup1[i]-1].speedofsound;  
			
		}
		mach[j] /= (double)(end-start);
		//cout<<"average value "<<point_values[j]<<endl;
		myfile<<mach[j]<<endl;		
		
	}

	myfile<<endl<<"SCALARS"<<"  "<<"Energy  float "<<1<<endl<<"LOOKUP_TABLE default"<<endl;
	double energy[npoin];
	for(j = 0; j < npoin; j++){				
		start = esup2[j], end = esup2[j+1];
		energy[j] = 0.0;		
		for(i = start; i < end; i++){						
			energy[j] += triangle[esup1[i]-1].densityEnergy/triangle[esup1[i]-1].density;  
			
		}
		energy[j] /= (double)(end-start);
		//cout<<"average value "<<point_values[j]<<endl;
		myfile<<energy[j]<<endl;		
		
	}

	myfile<<endl<<"SCALARS"<<"  "<<"Density  float "<<1<<endl<<"LOOKUP_TABLE default"<<endl;
	double density[npoin];
	for(j = 0; j < npoin; j++){				
		start = esup2[j], end = esup2[j+1];
		density[j] = 0.0;		
		for(i = start; i < end; i++){						
			density[j] += triangle[esup1[i]-1].density;  
			
		}
		density[j] /= (double)(end-start);
		//cout<<"average value "<<point_values[j]<<endl;
		myfile<<density[j]<<endl;		
		
	}
	
	
	myfile.close();
	
	
	return 0;
}

double convertToDouble(std::string const& s)
{
	std::istringstream i(s);
	double x;
	if (!(i >> x))
	  return -1;
	return x;
}




