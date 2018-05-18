#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <omp.h>
#include <memory.h>

#define LEADER 0

/*
 * run with
 * mpicc -fopenmp -Wall main.c -o elhaouari.exe
 * mpirun -np 4  elhaouari.exe matrice vecteur
 */
struct tablo {
    int * tab;
    int size;
};
struct matrice{
    int **tab;
    int size;
};


int multiplyvectors(int *v,int *b,int size){
    int result=0;
    #pragma omp parallel for reduction (+:result)
    for (int i = 0; i < size ; ++i) {
        result += v[i] * b[i];
    }
    return result;
}


int mod(int a, int b)
{
    int r = a % b;
    return r < 0 ? r + b : r;
}

struct matrice * allocateMatrice(int size){

    struct matrice *mat = malloc(sizeof(struct matrice));
    mat->size=size;
    mat->tab = (int **) malloc(size * sizeof(int *));
    for(int i = 0; i < size; i++)
    {
        mat->tab[i] = calloc(size,sizeof(int));
    }
    return mat;
}

struct tablo * allocateTablo(int size) {
    struct tablo * tmp = malloc(sizeof(struct tablo));
    tmp->size = size;
    tmp->tab = calloc(size,sizeof(int));
    //printArray(tmp);
    return tmp;
}


int getSizeN(char *fileName) {

   FILE * fichier =NULL;
    fichier = fopen(fileName, "r+");

    if(fichier == NULL){
        printf("Error couldn't read file null ");
        MPI_Finalize();
    }
    int token = 0;
    int size = 0;
    while(1){
        int ret = fscanf(fichier,"%d",&token);
        if(ret==1) {
            size++;
        } else if(ret == EOF)
            break;
    }
    return  size;
}



struct tablo *getVector(char *fileName, int size){
    FILE * fichier =NULL;
    fichier = fopen(fileName, "r+");

    if(fichier == NULL){
        printf("Error couldn't read file null ");
        MPI_Finalize();
    }

    struct tablo *v = allocateTablo(size);
    int token = 0;
    int i = 0;
    while(1){
        int ret = fscanf(fichier,"%d",&token);
        if(ret==1) {
            v->tab[i]= token;
            i++;
        } else if(ret == EOF)
            break;
    }

    return v;

}

struct matrice *getMatrice(char *fileName, int n){

    struct matrice *matrix = allocateMatrice(n);


    FILE * fichier =NULL;
    fichier = fopen(fileName, "r+");

    if(fichier == NULL){
        printf("Error couldn't read file null ");
        MPI_Finalize();
    }

    int token = 0;
    int i = 0;
    while(1){
        int ret = fscanf(fichier,"%d",&token);
        if(ret==1) {
            matrix->tab[i/n][i%n]= token;
            i++;
        } else if(ret == EOF)
            break;
    }

    return matrix;

}



struct tablo *filltablo(char *filename){
    int n = getSizeN(filename);

    struct tablo *t = getVector(filename,n);
    return t;
}


struct matrice *fillmatrice(char *filename, int size){

    struct matrice *m = getMatrice(filename,size );
    return m;
}



struct tablo *broadcast(int currentrank,int world_size,struct tablo *vecteur, char *vecteurfile){


    if(currentrank==LEADER){

        vecteur = filltablo(vecteurfile);

        MPI_Send(&(vecteur->size),1,MPI_INT,
                 (currentrank+1)%world_size, //send to next
                 0,MPI_COMM_WORLD);
        MPI_Send(vecteur->tab,vecteur->size,MPI_INT,
                 (currentrank+1)%world_size, //send to
                 0,MPI_COMM_WORLD);

        printf("process %d sending to process %d \n", currentrank,(currentrank+1)%world_size );

    } else{
        int vectsize;

        MPI_Recv(&vectsize,1,MPI_INT,currentrank -1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        printf("process %d received token %d from process %d \n", currentrank,vectsize,currentrank-1 );

        vecteur = allocateTablo(vectsize);


        MPI_Recv(vecteur->tab,vecteur->size,MPI_INT,currentrank -1, 0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        printf("process %d received vecteur %d %d  %d %d   from process %d \n", currentrank,
               vecteur->tab[0],
               vecteur->tab[1],
               vecteur->tab[2],
               vecteur->tab[3]
                ,currentrank-1 );

        if(currentrank!= world_size-1){
            MPI_Send(&(vecteur->size),1,MPI_INT,(currentrank+1)%world_size,0,MPI_COMM_WORLD);
            MPI_Send(vecteur->tab,vecteur->size,MPI_INT,(currentrank+1)%world_size,0,MPI_COMM_WORLD);

        }

    }
    return vecteur;

}


struct matrice *scatter(struct matrice *A, struct tablo *vecteur, int currentrank, int world_size, char *matricefilename){
    if(currentrank==LEADER) {
        A   = fillmatrice(matricefilename, vecteur->size);

        //   MPI_Send(&(A->size), 1, MPI_INT, (world_rank + 1) % world_size, 0, MPI_COMM_WORLD);
        for(int i =1*(vecteur->size/world_size); i<vecteur->size ;i++){
            MPI_Send(A->tab[i], A->size, MPI_INT, (currentrank + 1) % world_size, 0, MPI_COMM_WORLD);
        }
    } else{
        int start= currentrank*(vecteur->size / world_size);
        int end= vecteur->size;//(world_rank+1)*(vecteur->size / world_size);
        A = allocateMatrice(vecteur->size);

        for(int i =start;i<end;i++){

            MPI_Recv(A->tab[i],vecteur->size,MPI_INT,currentrank -1, 0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            printf("process %d received matrice %d %d  %d %d  from process %d \n", currentrank,
                   A->tab[i][0],
                   A->tab[i][1],
                   A->tab[i][2],
                   A->tab[i][3]
                    ,currentrank-1 );

        }


        int next = (currentrank+1)*(vecteur->size / world_size);
        for (int i = next; i <end ; i++) {
            MPI_Send(A->tab[i], A->size, MPI_INT, (currentrank + 1) % world_size, 0, MPI_COMM_WORLD);


        }
    }

    return A;


}



struct tablo *computes(struct tablo *vecteur,int  currentrank, int world_size, struct matrice *A){
    struct tablo *computed = allocateTablo(vecteur->size);



    int lastindex =(currentrank+1)*(vecteur->size / world_size);
    int firstindex = currentrank*(vecteur->size / world_size);
    //pour gérer le cas N%P != 0
    if((vecteur->size%world_size != 0)&&(currentrank == world_size -1))lastindex = vecteur->size;

    for(int i = firstindex ;
        i<lastindex;
        i++){
        computed->tab[i] = multiplyvectors(vecteur->tab,A->tab[i],vecteur->size);
        printf("in processor %d calculated %d\n", currentrank,computed->tab[i]);
    }
    return  computed;
}





struct tablo *gather(struct tablo *computed, int currentrank, int world_size ){

    if(currentrank == LEADER) {

        for (int i = 1 ; i < world_size ; i++) {

         int index_in_result = mod(currentrank-i,world_size )*(computed->size / world_size);
         int size_of_data = (i==1) ?
                            computed->size/world_size +  computed->size%world_size
                            : computed->size/world_size;


   printf("size : %d, index %d  \n", size_of_data, index_in_result);

            MPI_Recv(computed->tab+index_in_result,
                     size_of_data,
                     MPI_INT, mod(currentrank - 1 + world_size, world_size), //precedent
                     0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            printf("received %d in index %d \n ", computed->tab[index_in_result],
                   index_in_result
            );
        }


    } else{
        int size_of_data = (currentrank==world_size-1) ?
                           computed->size/world_size +  computed->size%world_size
                                  : computed->size/world_size;


        int index = currentrank*(computed->size / world_size);


        MPI_Send(computed->tab+index,
                 size_of_data,MPI_INT,(currentrank+1)%world_size,0,MPI_COMM_WORLD
        );


        printf("\n Size of data %d for rank %d  with index %d +++ current rank +1 %d \n",size_of_data,currentrank,index , (currentrank+1)%world_size);

        for (int i = 0; i < ( currentrank -1 + world_size) % world_size; i++) {

            MPI_Recv(computed->tab+currentrank*(computed->size / world_size), (computed->size/world_size), MPI_INT, (currentrank - 1 + world_size) % world_size,
                     0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            MPI_Send(computed->tab+currentrank*(computed->size / world_size), (computed->size/world_size), MPI_INT, (currentrank + 1) % world_size,
                     0, MPI_COMM_WORLD);
        }

    }

    return computed;

}



int main(int argc, char** argv) {

    //Init mpi environment
    MPI_Init(NULL,NULL);

    //
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD,&world_size);


    //info transmise



    struct tablo *vecteur =NULL;


    vecteur = broadcast(world_rank,world_size,vecteur,argv[2]);

    printf("vecteur %d %d %d %d ", vecteur->tab[0],vecteur->tab[1],vecteur->tab[2],vecteur->tab[3]);

    struct matrice * A =NULL;

    A = scatter(A,vecteur,world_rank,world_size,argv[1]);


    printf("vecteur %d %d %d %d \n", A->tab[3][0],A->tab[3][1],A->tab[3][2],A->tab[3][3]);


    struct tablo *result = computes(vecteur,world_rank,world_size,A);


    result=gather(result,world_rank,world_size);


    if (world_rank == LEADER) {
        for (int i = 0; i < result->size; i++) {
            printf("%d\n", result->tab[i]);
        }
    }


    MPI_Finalize();

}