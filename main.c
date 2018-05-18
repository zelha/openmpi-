#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <omp.h>

#define LEADER 0

/*
 * run with
 * mpicc -fopenmp -Wall main.c -o main.exe
 * mpirun -np 4 a.out mat2_1.txt
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
struct tablo *filltablo(int size){
    struct tablo *t = allocateTablo(4);
    t->tab[0]=-10;
    t->tab[1]=46;
    t->tab[2]=-25;
    t->tab[3]=-23;
    return t;
}


struct matrice *fillmatrice(){
    struct matrice *m = allocateMatrice(4);
    m->tab[0][0]=-10;
    m->tab[0][1]=46;
    m->tab[0][2]=-25;
    m->tab[0][3]=-23;

    m->tab[1][0]=-17;
    m->tab[1][1]=33;
    m->tab[1][2]=1;
    m->tab[1][3]=-9;

    m->tab[2][0]=47;
    m->tab[2][1]=0;
    m->tab[2][2]=33;
    m->tab[2][3]=-23;

    m->tab[3][0]=-20;
    m->tab[3][1]=30;
    m->tab[3][2]=-32;
    m->tab[3][3]=-18;
    return m;
}



struct tablo *broadcast(int currentrank,int world_size, struct tablo *vecteur){


    if(currentrank==LEADER){
        vecteur = filltablo(4); //changer ça

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


struct matrice *scatter(struct matrice *A, struct tablo *vecteur, int currentrank, int world_size){
    if(currentrank==LEADER) {
        A   = fillmatrice();

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

        //on gere le dernier tout seul.

  /*      MPI_Recv(computed->tab+(mod(currentrank-1,world_size )*(computed->size / world_size)),
                 (computed->size/world_size + computed->size % world_size),
                 MPI_INT, mod( currentrank - 1 + world_size, world_size),
                 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//TODO clean code


        printf("received %d  in index %d \n ",
               computed->tab[( mod(currentrank-1, world_size ))*(computed->size / world_size)]
               ,(mod(currentrank-1, world_size ))*(computed->size / world_size)
        );
*/

        for (int i = 1 ; i < world_size ; i++) {

       //     printf(" %d %d  \n ",(currentrank-i),world_size );

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

/*    } else if (currentrank == world_size - 1 ){
        MPI_Send(computed->tab+currentrank*(computed->size / world_size),
                 (computed->size/world_size + computed->size % world_size),MPI_INT,(currentrank+1)%world_size,0,MPI_COMM_WORLD
        );

        for (int i = 0; i < ( currentrank -1 + world_size) % world_size; i++) {

            MPI_Recv(computed->tab+currentrank*(computed->size / world_size), (computed->size/world_size), MPI_INT, (currentrank - 1 + world_size) % world_size,
                     0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            MPI_Send(computed->tab+currentrank*(computed->size / world_size), (computed->size/world_size), MPI_INT, (currentrank + 1) % world_size,
                     0, MPI_COMM_WORLD);
        }*/
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
    int token;

    struct tablo *vecteur;

    struct tablo *computed;

    vecteur = broadcast(world_rank,world_size,vecteur);

    printf("vecteur %d %d %d %d ", vecteur->tab[0],vecteur->tab[1],vecteur->tab[2],vecteur->tab[3]);

    struct matrice * A;

    A = scatter(A,vecteur,world_rank,world_size);


    printf("vecteur %d %d %d %d \n", A->tab[3][0],A->tab[3][1],A->tab[3][2],A->tab[3][3]);


    struct tablo *result = computes(vecteur,world_rank,world_size,A);


    result=gather(result,world_rank,world_size);


    if (world_rank == LEADER) {
        printf("final vector");
        for (int i = 0; i < result->size; i++) {
            printf("%d\n", result->tab[i]);
        }
    }


    MPI_Finalize();

}