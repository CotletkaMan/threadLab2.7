#include "paralel.h"
#include <string.h>
#include <iostream>

int** createBlockOrder(int* order , int y , int x , int countBlocks){
	int** intervals = new int*[countBlocks + countBlocks - 1];
	for(int i = 0 ; i < countBlocks + countBlocks - 1 ; i++)
		intervals[i] = new int[2];
	int span = 0 , step = (y + 1) / countBlocks;

	intervals[0][0] = 1;
	for(int i = 0 ; i < y ; i++){
		if((i + 1) % step){
			for(int j = 0 ; j < x ; j++)
				order[i * x + j] -= span * x;
		}
		else{
			intervals[2 * span + 1][0] = intervals[2 * span][1] = intervals[2 * span][0] + (step - 1) * x;
			intervals[2  * span + 2][0] = intervals[2 * span + 1][1] = intervals[2 * span + 1][0] + x;
			for(int j = 0 ; j < x ; j++)
				order[i * x + j] = y * x + 1 + (span - countBlocks + 1) * x + j;
			span++;
		}
	}
	intervals[2 * span][1] = x * y + 1;
	return intervals;
}

FILE* openFile(char* path){
	return fopen(path , "w+");
}

void closeFile(FILE* file){
	fclose(file);
}

void writeToFile(FILE* file , double** array , int* order,  int y , int x){
	for(int i = 0 ; i < y ; i++ ){
		fprintf(file, "%d %f\n", i , *array[order[i + 1]]);	
	}
	fprintf(file , "\n\n");
}

void createScript(int X , int Y , int T , double max){
	FILE *file = fopen("plotScript.sh\0" , "w");
	fprintf(file, "i = 0\nn = %d\nset terminal gif animate\nset output 'report.gif'\nset xrange [0:%d]\nset yrange [0:%f]\nload 'animateGnu'", T , Y, max);
	fclose(file);
}