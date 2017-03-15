#include<iostream>
#include<vector>
#include<math.h>
#include<fstream>
#define npoint 36
#define tolerance 0.09 
#define MIN(a,b) (((a)<(b)) ? (a):(b))

using namespace std;

vector <double> coordinates;
vector <double> triangles;
vector <double> faces;
vector <double> activepoints;

int checktri(double face1,double face2){
	// this function checks if a particular face is present more than two time in the triangle list. //
	// it was needed so same faces cannot enter the face list more than two times //
	// one face can only be shared by two triangles //
	int count;
	count=0;
	for(int i=0;i<triangles.size();i+=3){
		if((triangles[i]==face1 and triangles[i+1]==face2) or (triangles[i]==face2 and triangles[i+1]==face1)){	
			count++;					
		}
		if((triangles[i+1]==face1 and triangles[i+2]==face2) or (triangles[i+1]==face2 and triangles[i+2]==face1)){			
			count++;		
		}
		if((triangles[i]==face1 and triangles[i+2]==face2) or (triangles[i]==face2 and triangles[i+2]==face1)){			
			count++;		
		}	
		if(count>=2){
			count=0;
			return 1;	
		}
				
	}	
	return 0;	
}


int checkface(double face1, double replace){

	// this function checks the lenght of any face //

	double x,y,dist;
	x=coordinates[(face1-1)*2-1]-coordinates[(replace-1)*2-1];
	y=coordinates[(face1-1)*2]-coordinates[(replace-1)*2];
	dist=sqrt(x*x+y*y);
	for(int i=0;i<faces.size();i+=2){
		if((faces[i]==face1 and faces[i+1]==replace) or (faces[i]==replace and faces[i+1]==face1)){			
			return 1;		
		}			
	}
	if(dist>=0.5){
		return 1;	
	}
	return 2;
}

int initboundaries(){

	// this function initializes the coordinates of initial front manually along the boundary //

	for(int i=0;i<(npoint*2);i+=2){
		if(i/2<8){
			coordinates.push_back(i/2*0.1);
			coordinates.push_back(0.0);	
		}
		if(i/2>=8 and i/2<14){
			coordinates.push_back(0.7);
			coordinates.push_back((i/2-7)*0.1);				
			
		}
	}		
	/* manually insert rest of the points */
	coordinates.push_back(0.6);
	coordinates.push_back(0.6);
	coordinates.push_back(0.5);
	coordinates.push_back(0.6);
	coordinates.push_back(0.4);
	coordinates.push_back(0.6);
	coordinates.push_back(0.3);
	coordinates.push_back(0.6);
	coordinates.push_back(0.2);
	coordinates.push_back(0.6);
	coordinates.push_back(0.1);
	coordinates.push_back(0.6);
	coordinates.push_back(0.0);
	coordinates.push_back(0.6);			
		
	coordinates.push_back(0.0);
	coordinates.push_back(0.5);
	coordinates.push_back(0.0);
	coordinates.push_back(0.4);
	coordinates.push_back(0.0);
	coordinates.push_back(0.3);
	coordinates.push_back(0.0);
	coordinates.push_back(0.2);
	coordinates.push_back(0.0);
	coordinates.push_back(0.1);
			
	
	
	coordinates.push_back(0.2);
	coordinates.push_back(0.40);
	coordinates.push_back(0.30);
	coordinates.push_back(0.40);
	coordinates.push_back(0.40);
	coordinates.push_back(0.40);
	coordinates.push_back(0.50);
	coordinates.push_back(0.40);
			
	coordinates.push_back(0.5);
	coordinates.push_back(0.30);

	coordinates.push_back(0.50);
	coordinates.push_back(0.20);
	coordinates.push_back(0.40);
	coordinates.push_back(0.20);
	coordinates.push_back(0.30);
	coordinates.push_back(0.20);
	coordinates.push_back(0.20);
	coordinates.push_back(0.20);
			
	coordinates.push_back(0.2);
	coordinates.push_back(0.30);		
			
	
	return 0;
}

int initfaces(){
	// this function constructs the initial face list from the initial front //

	for(int i=1;i<(npoint+1);i++){
		if(i==26 and i+1==27){
			faces.push_back(i);
			faces.push_back(1);		
		}
		else if(i==36){
			faces.push_back(36);
			faces.push_back(27);		
		}
		else{
		faces.push_back(i);
		faces.push_back(i+1);
		}
	}
	return 0;
}

int writetofile(){

	// this function writes the coordinates and connectivities to the file //
	ofstream myfile;
	int i;
	myfile.open("triangles.vtk");
	myfile<<"# vtk DataFile Version 2.0"<<endl<<endl;
	myfile<<"ASCII"<<endl<<"DATASET UNSTRUCTURED_GRID"<<endl;
	myfile<<endl<<"POINTS"<<"    "<<coordinates.size()/2<<"    "<<"float"<<endl;


	for(i=0;i<coordinates.size();i+=2){
		myfile<<coordinates[i]<<"    "<<coordinates[i+1]<<"    "<<0.0<<endl;	
	}
	myfile<<endl<<endl;

	myfile<<"CELLS"<<"    "<<triangles.size()/3<<"    "<<4*triangles.size()/3<<endl;
	for(i=0;i<triangles.size();i+=3){
		myfile<<3.0<<"    "<<triangles[i]-1<<"    "<<triangles[i+1]-1<<"    "<<triangles[i+2]-1<<endl;	
	}
	myfile<<endl;
	myfile<<"CELL_TYPES"<<"    "<<triangles.size()/3<<endl;
	
	for(i=0;i<triangles.size()/3;i++){
		myfile<<5<<endl;
		
	}
	myfile.close();

	return 0;
}
int CheckConflict(double xc,double yc,double face1,double face2){
	
	// this function checks if a proposed face is creating conflict with any other existing face //
 	
	double ax,ay,bx,by,cx,cy,dx,dy,v1x,v1y,v2x,v2y,v3x,v3y,v4x,v4y,v5x,v5y,v6x,v6y,vax,vay,vaz,vbx,vby,vbz,vcx,vcy,vcz;
	double va,vb,vc,vd,vdx,vdy,vdz;
	double trifacep1,trifacep2,trifacep3;
	
	
	if(xc<0.0 or xc>0.7 or yc<0.0 or yc>0.6){
		
		//return 1;	
	}
	if(xc>0.2 and xc<0.5 and yc>0.2 and yc<0.4){
		
		//return 1;	
	}


	for(int i=0;i<triangles.size();i+=3){
		trifacep1=triangles[i];
		trifacep2=triangles[i+1];
		trifacep3=triangles[i+2];
		
		//check between 1-ipnew and triface 12//
		ax=coordinates[(face1-1)*2];
		ay=coordinates[((face1-1)*2)+1];
		bx=xc;
		by=yc;
		cx=coordinates[(trifacep1-1)*2];
		cy=coordinates[(trifacep1-1)*2+1];
		dx=coordinates[(trifacep2-1)*2];
		dy=coordinates[(trifacep2-1)*2+1];
		v1x=bx-ax;
		v1y=by-ay;
		v2x=cx-ax;
		v2y=cy-ay;
		v3x=dx-ax;
		v3y=dy-ay;
		v4x=dx-cx;
		v4y=dy-cy;
		v5x=-v2x;
		v5y=-v2y;
		v6x=bx-cx;
		v6y=by-cy;
		vaz=v1x*v2y-v2x*v1y;
		vbz=v1x*v3y-v1y*v3x;
		vcz=v4x*v5y-v4y*v5x;
		vdz=v4x*v6y-v4y*v6x;			
		if(vaz*vbz<0.0 and vcz*vdz<0.0){
			return 1;		
		}
		//check between 1-ipnew and triface 23//
		ax=coordinates[(face1-1)*2];
		ay=coordinates[((face1-1)*2)+1];
		bx=xc;
		by=yc;
		cx=coordinates[(trifacep2-1)*2];
		cy=coordinates[(trifacep2-1)*2+1];
		dx=coordinates[(trifacep3-1)*2];
		dy=coordinates[(trifacep3-1)*2+1];
		v1x=bx-ax;
		v1y=by-ay;
		v2x=cx-ax;
		v2y=cy-ay;
		v3x=dx-ax;
		v3y=dy-ay;
		v4x=dx-cx;
		v4y=dy-cy;
		v5x=-v2x;
		v5y=-v2y;
		v6x=bx-cx;
		v6y=by-cy;
		vaz=v1x*v2y-v2x*v1y;
		vbz=v1x*v3y-v1y*v3x;
		vcz=v4x*v5y-v4y*v5x;
		vdz=v4x*v6y-v4y*v6x;			
		if(vaz*vbz<0.0 and vcz*vdz<0.0){
			return 1;		
		}
		//check between 1-ipnew and triface 13	
		ax=coordinates[(face1-1)*2];
		ay=coordinates[((face1-1)*2)+1];
		bx=xc;
		by=yc;		
		cx=coordinates[(trifacep1-1)*2];
		cy=coordinates[(trifacep1-1)*2+1];
		dx=coordinates[(trifacep3-1)*2];
		dy=coordinates[(trifacep3-1)*2+1];
		v1x=bx-ax;
		v1y=by-ay;
		v2x=cx-ax;
		v2y=cy-ay;
		v3x=dx-ax;
		v3y=dy-ay;
		v4x=dx-cx;
		v4y=dy-cy;
		v5x=-v2x;
		v5y=-v2y;
		v6x=bx-cx;
		v6y=by-cy;
		vaz=v1x*v2y-v2x*v1y;
		vbz=v1x*v3y-v1y*v3x;
		vcz=v4x*v5y-v4y*v5x;
		vdz=v4x*v6y-v4y*v6x;			
		if(vaz*vbz<0.0 and vcz*vdz<0.0){
			return 1;		
		}
			
			
		//check between ipnew-2 and triface 12
		bx=coordinates[(face2-1)*2];
		by=coordinates[((face2-1)*2)+1];
		ax=xc;
		ay=yc;
		cx=coordinates[(trifacep1-1)*2];
		cy=coordinates[(trifacep1-1)*2+1];
		dx=coordinates[(trifacep2-1)*2];
		dy=coordinates[(trifacep2-1)*2+1];
		v1x=bx-ax;
		v1y=by-ay;
		v2x=cx-ax;
		v2y=cy-ay;
		v3x=dx-ax;
		v3y=dy-ay;
		v4x=dx-cx;
		v4y=dy-cy;
		v5x=-v2x;
		v5y=-v2y;
		v6x=bx-cx;
		v6y=by-cy;
		vaz=v1x*v2y-v2x*v1y;
		vbz=v1x*v3y-v1y*v3x;
		vcz=v4x*v5y-v4y*v5x;
		vdz=v4x*v6y-v4y*v6x;			
		if(vaz*vbz<0.0 and vcz*vdz<0.0){
			return 1;		
		}
		//check between ipnew-2 and triface 23//
		bx=coordinates[(face2-1)*2];
		by=coordinates[((face2-1)*2)+1];
		ax=xc;
		ay=yc;
		cx=coordinates[(trifacep2-1)*2];
		cy=coordinates[(trifacep2-1)*2+1];
		dx=coordinates[(trifacep3-1)*2];
		dy=coordinates[(trifacep3-1)*2+1];
		v1x=bx-ax;
		v1y=by-ay;
		v2x=cx-ax;
		v2y=cy-ay;
		v3x=dx-ax;
		v3y=dy-ay;
		v4x=dx-cx;
		v4y=dy-cy;
		v5x=-v2x;
		v5y=-v2y;
		v6x=bx-cx;
		v6y=by-cy;
		vaz=v1x*v2y-v2x*v1y;
		vbz=v1x*v3y-v1y*v3x;
		vcz=v4x*v5y-v4y*v5x;
		vdz=v4x*v6y-v4y*v6x;			
		if(vaz*vbz<0.0 and vcz*vdz<0.0){
			return 1;		
		}
		//check between ipnew-2 and triface 13
		bx=coordinates[(face2-1)*2];
		by=coordinates[((face2-1)*2)+1];
		ax=xc;
		ay=yc;
		cx=coordinates[(trifacep1-1)*2];
		cy=coordinates[(trifacep1-1)*2+1];
		dx=coordinates[(trifacep3-1)*2];
		dy=coordinates[(trifacep3-1)*2+1];
		v1x=bx-ax;
		v1y=by-ay;
		v2x=cx-ax;
		v2y=cy-ay;
		v3x=dx-ax;
		v3y=dy-ay;
		v4x=dx-cx;
		v4y=dy-cy;
		v5x=-v2x;
		v5y=-v2y;
		v6x=bx-cx;
		v6y=by-cy;
		vaz=v1x*v2y-v2x*v1y;
		vbz=v1x*v3y-v1y*v3x;
		vcz=v4x*v5y-v4y*v5x;
		vdz=v4x*v6y-v4y*v6x;			
		if(vaz*vbz<0.0 and vcz*vdz<0.0){
			return 1;		
		}
		if(face1==27 and face2==23){
			return 0;
		}		
	}	

	return 0;
}


int addnewElems(int points);

int main(){	
	// this is the main entry point of the code //
	int j;
	int points;
	points=npoint;
	initboundaries();
	initfaces();
	addnewElems(points);
	writetofile();

	return 0;

}

int addnewElems(int points){
	
	// this function creates the new	
	
	int j;
	j=0;
	int k,tvalue;
	double dely, delx,face1,face2,normal,slope,nx,ny,dotproduct,a,area,h;
	double point1x,point1y,point2x,point2y,tx,ty;
	double x0,y0,xc,yc,xd,yd,dcd,ipnew;	
	tvalue=0;
	while(!faces.empty()){
		//loop until the face list is empty //		
		//pop face i.e. extract two vertices of the face and delete them from the list // 
		
		//an if condition set just for debugging purposes. it forcefully terminates the loop //
		if(j>664){
			//break;		
		}
		//extract and delete the top face from the list //		
		face1=faces[0];
		face2=faces[1];
		faces.erase(faces.begin());		
		faces.erase(faces.begin());
		//check if face1 face2 is already part of any triangle //
		tvalue=checktri(face1,face2);
		if(tvalue==1){
			continue;	
		}
		//calculate ipnew check if it creates any conflict //
		k=j/2;
		point1x=coordinates[(face1-1)*2];
		point2x=coordinates[(face2-1)*2];
		point1y=coordinates[((face1-1)*2)+1];
		point2y=coordinates[((face2-1)*2)+1];
		x0=0.5*(point2x+point1x);
		y0=0.5*(point2y+point1y);
		tx=point2x-point1x;
		ty=point2y-point1y;
		a=sqrt(tx*tx+ty*ty);
		h=0.866*a;
		tx=tx/a;
		ty=ty/a;
		nx=-ty;
		ny=tx;
		xc=x0+h*nx;
		yc=y0+h*ny;
		
		//see if ipnew can be replaced with any existing point //
		double xx,yy,replace,side12x,side12y,side23x,side23y,side31x,side31y,s12,s23,s31,s3,warea,sarea1,temp,dist;
		int flag,retValue,fvalue;
		fvalue=10;
		flag=10;
		for(int m=0;m<coordinates.size();m+=2){				
			xx=coordinates[m];
			yy=coordinates[m+1];
			temp=sqrt((xx-xc)*(xx-xc)+(yy-yc)*(yy-yc));			
			if(temp<=tolerance){
				xc=xx;
				yc=yy;
				flag=1;
				replace=m/2+1;
				if(face1==27 and face2==23){
					replace--;								
				}					
			}
		}
		//special case for four corner nodes //
		if((face1==26 and face2==1) or (face1==8 and face2==9) or(face1==14 and face2==15) or (face1==21 and face2==22)){
			continue;
		}
		
		if(flag ==10){ //ipnew is selected as no conflict reported. add new faces in facelist and triangles in triangle list 	
			retValue=CheckConflict(xc,yc,face1,face2);			
			if(retValue==0){
				coordinates.push_back(xc);
				coordinates.push_back(yc);
				points+=1;
				faces.push_back(face1);
				faces.push_back(points);
				fvalue=10;
				faces.push_back(points);
				faces.push_back(face2);
				triangles.push_back(face1);
				triangles.push_back(points);
				triangles.push_back(face2);
			}
										
		}//end if condition of flag 
		if (flag==1){  //ipnew is replaced by existing point  //			
			retValue=CheckConflict(xc,yc,face1,face2);
			if(retValue==0){								
				if((face1==7 and face2==8)or(face1==13 and face2==14) or(face1==20 and face2==21)or(face1==26 and face2==1)){
					//check if any of the faces is already in the face list //
					fvalue=checkface(face1,replace);
					if(fvalue==2){
						faces.push_back(face1);
						faces.push_back(replace);
						fvalue=10;						
					}
					triangles.push_back(face1);
					triangles.push_back(replace);
					triangles.push_back(face2);									
				}
				else if(face1==1 and face2==2){
					fvalue=checkface(face1,replace);
					if(fvalue==2){
						faces.push_back(face1);
						faces.push_back(replace);
						fvalue=10;						
					}
					triangles.push_back(face1);
					triangles.push_back(replace);
					triangles.push_back(face2);				
				}
				else{
					fvalue=checkface(face1,replace);
					if(fvalue==2){
						faces.push_back(face1);
						faces.push_back(replace);
						fvalue=10;						
					}
					fvalue=checkface(replace,face2);
					if(fvalue==2){
						faces.push_back(replace);
						faces.push_back(face2);	
						fvalue=10;	
					}					
					triangles.push_back(face1);
					triangles.push_back(replace);
					triangles.push_back(face2);
				}//end else of face check			
			}//end else of existing point check									
		}//end else flag condition		
		j+=2;	
	} // end while faces list is not empty
	return 0;
}


