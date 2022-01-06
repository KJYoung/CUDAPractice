/*
 * First version of Striped Arrayed Ungapped Alignment ( Not optimized )
 * 2021/01/06/13:17
 */

#include <stdio.h>
#include <iostream>
#include <fstream>

using namespace std;

int getCodonTableIndex(char aa){
    switch (aa) {
        case 'A':
            return 0;
        case 'R':
            return 1;
        case 'N':
            return 2;
        case 'D':
            return 3;
        case 'C':
            return 4;
        case 'Q':
            return 5;
        case 'E':
            return 6;
        case 'G':
            return 7;
        case 'H':
            return 8;
        case 'I':
            return 9;
        case 'L':
            return 10;
        case 'K':
            return 11;
        case 'M':
            return 12;
        case 'F':
            return 13;
        case 'P':
            return 14;
        case 'S':
            return 15;
        case 'T':
            return 16;
        case 'W':
            return 17;
        case 'Y':
            return 18;
        case 'V':
            return 19;
        case 'B':
            return 20;
        case 'Z':
            return 21;
        case 'X':
            return 22;
        case '*': default:
            return 23;
    }
}
int getSubMatrix(char targetAA, char queryAA){
    int targetIdx = getCodonTableIndex(targetAA);
    int queryIdx = getCodonTableIndex(queryAA);
    int largerIdx = (targetIdx > queryIdx) ? targetIdx : queryIdx;
    int smallerIdx = (targetIdx > queryIdx) ? queryIdx : targetIdx;

    static int codon_table[50][50]={
            { 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,-2,-1,-1,-1,-4},
            {-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,-1,-2, 0,-1,-4},
            {-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3, 4,-3, 0,-1,-4},
            {-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3, 4,-3, 1,-1,-4},
            { 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-1,-3,-1,-4},
            {-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2, 0,-2, 4,-1,-4},
            {-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2, 1,-3, 4,-1,-4},
            { 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,-1,-4,-2,-1,-4},
            {-2, 0, 1,-1,-3,-0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3, 0,-3, 0,-1,-4},
            {-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,-3, 3,-3,-1,-4},
            {-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,-4, 3,-3,-1,-4},
            {-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2, 0,-3, 1,-1,-4},
            {-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,-3, 2,-1,-1,-4},
            {-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,-3, 0,-3,-1,-4},
            {-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,-2,-3,-1,-1,-4},
            { 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2, 0,-2, 0,-1,-4},
            { 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0,-1,-1,-1,-1,-4},
            {-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,-4,-2,-2,-1,-4},
            {-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,-3,-1,-2,-1,-4},
            { 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,-3, 2,-2,-1,-4},
            {-2,-1, 4, 4,-3, 0, 1,-1, 0,-3,-4, 0,-3,-3,-2, 0,-1,-4,-3,-3, 4,-3, 0,-1,-4},
            {-1,-2,-3,-3,-1,-2,-3,-4,-3, 3, 3,-3, 2, 0,-3,-2,-1,-2,-1, 2,-3, 3,-3,-1,-4},
            {-1, 0, 0, 1,-3, 4, 4,-2, 0,-3,-3, 1,-1,-3,-1, 0,-1,-2,-2,-2, 0,-3, 4,-1,-4},
            {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-4},
            {-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4, 1}
    };

    return codon_table[smallerIdx][largerIdx];
}
__device__
int getCodonTableIndexD(char aa){
    switch (aa) {
        case 'A':
            return 0;
        case 'R':
            return 1;
        case 'N':
            return 2;
        case 'D':
            return 3;
        case 'C':
            return 4;
        case 'Q':
            return 5;
        case 'E':
            return 6;
        case 'G':
            return 7;
        case 'H':
            return 8;
        case 'I':
            return 9;
        case 'L':
            return 10;
        case 'K':
            return 11;
        case 'M':
            return 12;
        case 'F':
            return 13;
        case 'P':
            return 14;
        case 'S':
            return 15;
        case 'T':
            return 16;
        case 'W':
            return 17;
        case 'Y':
            return 18;
        case 'V':
            return 19;
        case 'B':
            return 20;
        case 'Z':
            return 21;
        case 'X':
            return 22;
        case '*': default:
            return 23;
    }
}
__device__
int getSubMatrixD(char targetAA, char queryAA){
    int targetIdx = getCodonTableIndexD(targetAA);
    int queryIdx = getCodonTableIndexD(queryAA);
    int largerIdx = (targetIdx > queryIdx) ? targetIdx : queryIdx;
    int smallerIdx = (targetIdx > queryIdx) ? queryIdx : targetIdx;

    static int codon_table[50][50]={
            { 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,-2,-1,-1,-1,-4},
            {-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,-1,-2, 0,-1,-4},
            {-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3, 4,-3, 0,-1,-4},
            {-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3, 4,-3, 1,-1,-4},
            { 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-1,-3,-1,-4},
            {-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2, 0,-2, 4,-1,-4},
            {-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2, 1,-3, 4,-1,-4},
            { 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,-1,-4,-2,-1,-4},
            {-2, 0, 1,-1,-3,-0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3, 0,-3, 0,-1,-4},
            {-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,-3, 3,-3,-1,-4},
            {-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,-4, 3,-3,-1,-4},
            {-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2, 0,-3, 1,-1,-4},
            {-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,-3, 2,-1,-1,-4},
            {-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,-3, 0,-3,-1,-4},
            {-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,-2,-3,-1,-1,-4},
            { 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2, 0,-2, 0,-1,-4},
            { 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0,-1,-1,-1,-1,-4},
            {-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,-4,-2,-2,-1,-4},
            {-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,-3,-1,-2,-1,-4},
            { 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,-3, 2,-2,-1,-4},
            {-2,-1, 4, 4,-3, 0, 1,-1, 0,-3,-4, 0,-3,-3,-2, 0,-1,-4,-3,-3, 4,-3, 0,-1,-4},
            {-1,-2,-3,-3,-1,-2,-3,-4,-3, 3, 3,-3, 2, 0,-3,-2,-1,-2,-1, 2,-3, 3,-3,-1,-4},
            {-1, 0, 0, 1,-3, 4, 4,-2, 0,-3,-3, 1,-1,-3,-1, 0,-1,-2,-2,-2, 0,-3, 4,-1,-4},
            {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-4},
            {-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4, 1}
    };

    return codon_table[smallerIdx][largerIdx];
}

void fastaParser(string FilePath, string& parsedSeq){
    ifstream fastaStream(FilePath);
    string line;
    getline(fastaStream, line); // pass the first line.
    while(getline(fastaStream, line)){
        parsedSeq.append(line);
    }
    cout << "[" << parsedSeq << "]" << endl;
}

__global__ void unGappedAlignGPU_stripedArr(char* vectorCurr, char* vectorPrev, char queryi, char* target, char* max, int targetLeng)
{
    int i = threadIdx.x;
    int jump = blockDim.x;


    if(i == 0){
        char subVal = getSubMatrixD(queryi, target[0]);
        *vectorCurr = subVal > 0 ? subVal : 0;
        max[i] = (*vectorCurr > max[i]) ? *vectorCurr : max[i];
    }


    for (int j = i+1; j < targetLeng; j+= jump) {
        // < PARALLELIZE > //
        char tempA = *(vectorPrev + j - 1) + getSubMatrixD(queryi, target[j]);
        *(vectorCurr + j) = (tempA > 0) ? tempA : 0;
        max[j] = (tempA > max[j]) ? tempA : max[j];
    }
}

void unGappedAlignCPU(char* vectorCurr, char* vectorPrev, char queryi, char* target, char* max, int targetLeng)
{
    *vectorCurr = std::max(getSubMatrix(queryi, target[0]), 0);
    *max = std::max(*vectorCurr, *max);
    for (int j = 1; j < targetLeng; j++) {
        *(vectorCurr + j) = std::max(*(vectorPrev + j - 1) + getSubMatrix(queryi, target[j]), 0);
        *max = std::max(*(vectorCurr + j), *max);
    }
}

#define HUMANHEMALPHA5 "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYRMVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYRMVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYRMVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYRMVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR"
#define HUMANHEMBETA5  "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYHMVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYHMVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYHMVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYHMVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"
#define HUMANHEMALPHA "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR"
#define HUMANHEMBETA "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"

#define HUMANHEMALPHA0 "MVLSPADKTNVKA"
#define HUMANHEMBETA0  "MVHLTPEEKSA"

#define QUERYSEQUENCE HUMANHEMALPHA5
#define TARGETSEQUENCE HUMANHEMBETA5
int main() {
    using namespace std;
    string hemA, hemB;

    int iter = 100;
    // CPU task.
    double CPUtime, GPUtime;
    {
        char humanHemAlpha[] = QUERYSEQUENCE;
        char humanHemBeta[] = TARGETSEQUENCE;

        char *query = humanHemAlpha;
        char *target = humanHemBeta;

        int targetSize = strlen(target);
        printf("Start : with target size %d \n", targetSize);

        const clock_t begin_time = clock();
        char max = 0;

        char *vectorCurr = (char *) malloc(targetSize * sizeof(char));
        char *vectorPrev = (char *) malloc(targetSize * sizeof(char));
        memset(vectorPrev, 0, targetSize * sizeof(char));


        for (int tx = 0; tx < iter; tx++) {
            for (int i = 0; query[i] != '\0'; i++) {
                unGappedAlignCPU(vectorCurr, vectorPrev, query[i], target, &max, targetSize);

                char *vectorTmp = vectorPrev;
                vectorPrev = vectorCurr;
                vectorCurr = vectorTmp;
            }
        }
        CPUtime = clock() - begin_time;
        cout << "CPU result : " << (int) max << endl;
        cout << "Elapsed time : " << CPUtime / CLOCKS_PER_SEC << endl;
    }
    //GPU task with arrayed striped
    {
        int NUM_CUDA_THREAD = 256;

        char humanHemAlpha[] = QUERYSEQUENCE;
        char humanHemBeta[] = TARGETSEQUENCE;

        int querySize = strlen(humanHemAlpha);
        int targetSize = strlen(humanHemBeta);

        char *query, *target;
        query = humanHemAlpha;

        const clock_t begin_time = clock();

        cudaMalloc((void **) &target, targetSize * sizeof(char));
        cudaMemcpy(target, humanHemBeta, targetSize, cudaMemcpyHostToDevice);

        char *vectorCurr, *vectorPrev;
        cudaMalloc((void **) &vectorCurr, targetSize * sizeof(char));
        cudaMalloc((void **) &vectorPrev, targetSize * sizeof(char));
        cudaMemset(vectorCurr, 0, targetSize*sizeof(char));
        cudaMemset(vectorPrev, 0, targetSize*sizeof(char));

        char *max;
        cudaMalloc((void **) &max, NUM_CUDA_THREAD * sizeof(char));
        cudaMemset(max, 0, NUM_CUDA_THREAD*sizeof(char));

        char maxRes = 0;
        char *maxHost = (char*)malloc(NUM_CUDA_THREAD * sizeof(char));
        for (int tx = 0; tx < iter; tx++) {
            for(int i=0; i<querySize; i++){
                cudaMemset(vectorCurr, 0, sizeof(char));
                unGappedAlignGPU_stripedArr<<<1,NUM_CUDA_THREAD>>>(vectorCurr, vectorPrev, query[i], target, max, targetSize);
                cudaDeviceSynchronize();

                cudaMemcpy(maxHost, max, NUM_CUDA_THREAD * sizeof(char), cudaMemcpyDeviceToHost);

                for(int ini=0; ini<NUM_CUDA_THREAD; ini++){
                    if(maxRes < maxHost[ini]) {
                        maxRes = maxHost[ini];
                        //printf("Report [%d] at query %d.%d\n", maxRes, i, ini);
                    }
                }

                char* vectorTmp = vectorPrev;
                vectorPrev = vectorCurr;
                vectorCurr = vectorTmp;
            }
        }
        GPUtime = clock() - begin_time;
        cout << "GPU result : " << (int) maxRes << "[striped with Array overhead]" << endl;
        cout << "Elapsed time : " << GPUtime / CLOCKS_PER_SEC << endl;

        cudaFree(target); cudaFree(vectorCurr); cudaFree(vectorPrev); cudaFree(max);
        delete[] maxHost;
    }
    cout << "Speed Ratio : G/C = " << GPUtime / CPUtime << endl;
    return 0;
}
