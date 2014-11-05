#include <stdio.h>
#include <stdlib.h>


int main(){

  FILE *fp,*fq;

  int n[3];
  int types,anz=0;
  double lat1[3],lat2[3],lat3[3];
  double p[3],p_new[3];
  double a;
  int i,j1,j2,j3,cells;
  char line1[200],line2[200],line3[100];



  if((fp = fopen("cell.dat","r")) == NULL){
    printf("\n    File cell.dat not found\n\n");
    printf("    Sample file cell.dat has been written:\n\n");
    fq = fopen("cell.dat","w");
    fprintf(fq,"%d %d %d\n%d\n",2,2,1,4);
    fclose(fq);
    printf("%d %d %d \n%d\n\n",2,2,1,4);
    printf("    Program stops!\n\n");
    exit(1);
  }
  else{
    printf("\n    file cell.dat found\n\n");
  }
  
  fscanf(fp,"%d %d %d",&n[0],&n[1],&n[2]);
  fscanf(fp,"%d",&types);
  fclose(fp);
  cells = n[0]*n[1]*n[2];

  int num[types];

  fp = fopen("POSCAR","r");
  fq = fopen("POSCAR_SC","w");

  fgets(line1,sizeof line1,fp);
  sprintf(line2,"%dx%dx%d_supercell %s",n[0],n[1],n[2],line1);
  fprintf(fq,"%s",line2);

  fscanf(fp,"%le",&a);
  fprintf(fq,"%le\n",a);

  fscanf(fp,"%le %le %le",&lat1[0],&lat1[1],&lat1[2]);
  fscanf(fp,"%le %le %le",&lat2[0],&lat2[1],&lat2[2]);
  fscanf(fp,"%le %le %le",&lat3[0],&lat3[1],&lat3[2]);
  fprintf(fq,"%le %le %le\n",(double)n[0]*lat1[0],(double)n[0]*lat1[1],(double)n[0]*lat1[2]);
  fprintf(fq,"%le %le %le\n",(double)n[1]*lat2[0],(double)n[1]*lat2[1],(double)n[1]*lat2[2]);
  fprintf(fq,"%le %le %le\n",(double)n[2]*lat3[0],(double)n[2]*lat3[1],(double)n[2]*lat3[2]);
  
  for(i=0;i<types;i++){
    fscanf(fp,"%d",&num[i]);
    anz+=num[i];
    if(i==0){
      fprintf(fq,"%d",num[i]*cells);
    }
    else{
      fprintf(fq," %d",num[i]*cells);
    }
  }
  fprintf(fq,"\n");

  fscanf(fp,"%s",line3);
  fprintf(fq,"%s\n",line3);
    
  for(i=0;i<anz;i++){
    
    fscanf(fp,"%le %le %le",&p[0],&p[1],&p[2]);

    for(j1=0;j1<n[0];j1++){
      for(j2=0;j2<n[1];j2++){
	for(j3=0;j3<n[2];j3++){

	  p_new[0] = (p[0] + (double)j1)/(double)n[0];
	  p_new[1] = (p[1] + (double)j2)/(double)n[1];
	  p_new[2] = (p[2] + (double)j3)/(double)n[2];

	  fprintf(fq,"%le %le %le\n",p_new[0],p_new[1],p_new[2]);
	  
	}
      }
    }

  }

  fclose(fp);
  fclose(fq);


  return 0;

}
