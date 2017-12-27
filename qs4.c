#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <mpi.h>
#include <string.h>


void swap(int* a, int*b){
	int t = *a;
	*a = *b;
	*b = t;
}

int partition(int arr[], int low, int high){

	int pivot = arr[high];
	int i = (low-1);
	int j=0;
	for(j=low;j<=high-1;j++){
		if(arr[j]<=pivot){
			i++;
			swap(&arr[i],&arr[j]);
		}
	}
	swap(&arr[i+1], &arr[high]);
	return(i+1);
}

void quickSort(int arr[], int low, int high){
	if(low<high){
		int pi = partition(arr,low,high);
		quickSort(arr,low,pi-1);
		quickSort(arr,pi+1,high);
	} 
}


int main(int argc, char* argv[]){

	// Fixed array - can be changed for dynamic memory allocation
	int GArray[64] = {3,14,15,12,9,7,5,10,20,30,22,16,90,110,150,67,120,87,5,4,26,15,2156,12315623,15,64134,154512,516,1315,21513,512,63341,126,315,2121,622,132,13,646,421,311,54,7894,4654,5415,41,46545,11,64654,465,546,88313,654,87489,456,54643,654641,88913,443,6461,74646,798,78975};
	int globalSize = sizeof(GArray)/sizeof(GArray[0]);
	int rank,size,chunkSize,pivot,localPivot,breakPoint,flag=0;

	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);

	// Initialize local array and compute the chunkSize
	pivot = ceil(globalSize/2);
	pivot = GArray[pivot];
	chunkSize = globalSize / size ;

	int LArray[chunkSize],i=0,j=0;
	int SampleArr[size];
	int FullSampleArr[size*size];



	/* ------------------------------- BUILD LOCAL ARRAYS ------------------------------------------------------*/

	MPI_Scatter(GArray,chunkSize,MPI_INT,LArray,chunkSize,MPI_INT,0,MPI_COMM_WORLD);

	// Sort the individual arrays
	quickSort(LArray,0,chunkSize-1);

	MPI_Barrier(MPI_COMM_WORLD);

	for(i=0;i<chunkSize;i++)
		printf("%d ",LArray[i]);
	printf("\n");


	/* -------------------------- BUILD SAMPLE SET ---------------------------------------------------------------*/

	// Collect Sample set and sort it
	for(i=0;i<size;i++)
		SampleArr[i]=LArray[i*globalSize/(size*size)];
	

	/*------------------------- AGGREGATE SAMPLE SET AND COLLECT INTO A SINGLE ARRAY ------------------------------*/
	MPI_Gather(SampleArr,size,MPI_INT,FullSampleArr,size,MPI_INT,0,MPI_COMM_WORLD);
	quickSort(FullSampleArr,0,size*size-1);


	MPI_Bcast(FullSampleArr,size*size,MPI_INT,0,MPI_COMM_WORLD);

	localPivot = FullSampleArr[(int)ceil(sizeof(FullSampleArr)/sizeof(FullSampleArr[0]))/2];

	

	/*---------------------------- BREAK POINT AND SUBARRAY GENERATION -------------------------------------------*/

	// Calculate the break point of each local array
	for(breakPoint=0;breakPoint<chunkSize;breakPoint++){ 
		if(localPivot<=LArray[breakPoint]){
			if(localPivot==LArray[breakPoint])
				flag=1;
			break;
		}
	}


	int leftSub[chunkSize];
	int rightSub[chunkSize];


	// SET VALUES FOR THE LEFT SUBARRAY
	for(i=0;i<breakPoint;i++)
		leftSub[i]=LArray[i];	// Set the values upto breakpoint for the initial values

	for(i=breakPoint ; i<chunkSize ; i++)
		leftSub[i]=-32767;	// Set remaining values to -32767

	

	// SET VALUES FOR THE RIGHT SUBARRAY

	if(flag==0){
		for(i=breakPoint;i<chunkSize;i++)
			rightSub[i-breakPoint]=LArray[i];
		for(i=chunkSize-breakPoint;i<chunkSize;i++)
			rightSub[i]=-32767;
	}

	else{
		breakPoint++; // Skip the broadcasted pivot
		for(i=breakPoint ; i<chunkSize ; i++){
			rightSub[i-breakPoint]=LArray[i];
		}
		for(i=chunkSize-breakPoint;i<breakPoint+1;i++)
			rightSub[i]=-32767;
	}


	MPI_Barrier(MPI_COMM_WORLD);


	int leftSubSize = sizeof(leftSub)/sizeof(int);
	int rightSubSize = sizeof(rightSub)/sizeof(int);


	MPI_Barrier(MPI_COMM_WORLD);

	printf("\n");
	for(i=0;i<leftSubSize;i++)
		printf("%d ",leftSub[i]);

	MPI_Barrier(MPI_COMM_WORLD);

	rightSubSize = sizeof(rightSub)/sizeof(int);

	printf("\n");
		for(i=0;i<rightSubSize;i++)
		printf("%d ",rightSub[i]);
	printf("\n");

	
/*------------------------------ ALL to ALL SECTION ---------------------*/


	int Sorted1[globalSize];
	int Sorted2[globalSize];
	
//	MPI_Alltoall(leftSub,chunkSize,MPI_INT,Sorted1,chunkSize,MPI_INT,MPI_COMM_WORLD);
//	MPI_Alltoall(rightSub,chunkSize,MPI_INT,Sorted2,chunkSize,MPI_INT,MPI_COMM_WORLD);
	MPI_Gather(leftSub,chunkSize,MPI_INT,Sorted1,chunkSize,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Gather(rightSub,chunkSize,MPI_INT,Sorted1,chunkSize,MPI_INT,1,MPI_COMM_WORLD);

	//if(rank!=0)
	//	MPI_Send(Sorted2,chunkSize,MPI_INT,0,0,MPI_COMM_WORLD);
	//if(rank==0)
	//	MPI_Recv(Sorted2,chunkSize,MPI_INT,1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

	if(rank==0){
		printf("\n");
		quickSort(Sorted1,0, globalSize-1);
		for(i=0;i<globalSize;i++)
			printf("%d ",Sorted1[i]);
		printf("%d\n", localPivot);
		printf("\n");
	}
	if(rank==1){
		printf("\n");
		quickSort(Sorted1, 0, globalSize-1);
		for(i=0;i<globalSize;i++)           
			printf("%d ",Sorted1[i]);
	printf("\n");
	}
	int finalArray[globalSize];

	if(rank == 0){
		MPI_Status status;
		MPI_Recv(Sorted2, globalSize, MPI_INT, 1, 0, MPI_COMM_WORLD, &status);
	}
	if(rank == 1){
		MPI_Send(Sorted1, globalSize, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}
	
	if(rank == 0){
		int j = 0;
		for(int i = 0; i < globalSize; i++){
			if(Sorted1[i]!=-32767){
				finalArray[j] = Sorted1[i];
				j++;
			}
		}

		//printf("%d ", localPivot);
		finalArray[j] = localPivot;
		j++;

		for(int i = 0; i < globalSize; i++){
			if(Sorted2[i]!=-32767){
				finalArray[j] = Sorted2[i];
				j++;
			}
		}

		for(int i = 0; i < globalSize; i++){
			printf("%d ",finalArray[i]);
		}
	}
	MPI_Finalize();
}
