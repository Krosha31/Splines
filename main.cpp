#include <iostream>
#include <math.h>
#include <stdio.h>


struct spline {
    int N = 0;
    float *x = nullptr;
    float *y = nullptr;
    float *h = nullptr;
    float *l = nullptr;
    float* delta = nullptr;
    float *lambda = nullptr;
    float *c = nullptr;
    float *b = nullptr;
    float *d = nullptr;
};

void CountNumberOfPoints(FILE* InFile, spline* spl){
    bool flag = false;
    do{
        flag = false;
        while(fgetc(InFile)!='\n' && !feof(InFile))
            flag = true;
        if(flag)
            spl->N++;
    }while(!feof(InFile));
    spl->N--;
}


void ReadPoints(FILE* InFile, spline* spl){
    int i=0;
    for(i=0; i<spl->N+1; i++){
        fscanf(InFile, "%f", &spl->x[i]);
        fscanf(InFile, "%f", &spl->y[i]);
    }
}


void AllocSpline(spline* spl){
    spl->x = new float[spl->N+1];
    spl->y = new float[spl->N+1];
    spl->h = new float[spl->N+1];
    spl->l = new float[spl->N+1];
    spl->delta = new float[spl->N+1];
    spl->lambda = new float[spl->N+1];
    spl->c = new float[spl->N+1];
    spl->d = new float[spl->N+1];
    spl->b = new float[spl->N+1];
}


void SplineDestructor(spline* spl){
    delete [] spl->x;
    delete [] spl->y;
    delete [] spl->h;
    delete [] spl->l;
    delete [] spl->delta;
    delete [] spl->lambda;
    delete [] spl->c;
    delete [] spl->d;
    delete [] spl->b;
}




spline* MakeSpline(FILE* InFile) {
    spline* spl = new spline;
    CountNumberOfPoints(InFile, spl);
    rewind(InFile);
    AllocSpline(spl);
    ReadPoints(InFile, spl);
    int k;
    for(k=1; k<=spl->N; k++){
        spl->h[k] = spl->x[k] - spl->x[k-1];
        if(spl->h[k]==0){
            printf("\nError, x[%d]=x[%d]\n", k, k-1);
            return 0;
        }
        spl->l[k] = (spl->y[k] - spl->y[k-1])/spl->h[k];
    }
    spl->delta[1] = - spl->h[2]/(2*(spl->h[1]+spl->h[2]));
    spl->lambda[1] = 1.5*(spl->l[2] - spl->l[1])/(spl->h[1]+spl->h[2]);
    for(k=3; k<=spl->N; k++){
        spl->delta[k-1] = - spl->h[k]/(2*spl->h[k-1] + 2*spl->h[k] + spl->h[k-1]*spl->delta[k-2]);
        spl->lambda[k-1] = (3*spl->l[k] - 3*spl->l[k-1] - spl->h[k-1]*spl->lambda[k-2]) /
                      (2*spl->h[k-1] + 2*spl->h[k] + spl->h[k-1]*spl->delta[k-2]);
    }
    spl->c[0] = 0;
    spl->c[spl->N] = 0;
    for(k=spl->N; k>=2; k--){
        spl->c[k-1] = spl->delta[k-1]*spl->c[k] + spl->lambda[k-1];
    }
    for(k=1; k<=spl->N; k++){
        spl->d[k] = (spl->c[k] - spl->c[k-1])/(3*spl->h[k]);
        spl->b[k] = spl->l[k] + (2*spl->c[k]*spl->h[k] + spl->h[k]*spl->c[k-1])/3;
    }
}


double GetPoint(spline* spl, float x) {
    int left = 0, right = spl->N;
    while (left + 1 < right) {
        int c = (left + right) / 2;
        if (x > spl->x[c]) {
            left = c;
        }
        else {
            right = c;
        }
    }
    return spl->y[right] + spl->b[right]*(x - spl->x[right]) + spl->c[right]*pow(x - spl->x[right], 2) + spl->d[right]*pow(x - spl->x[right], 3);
}


float FindPoint(spline* spl1, spline* spl2, bool &flag_point) {
    float left = std::max(spl1->x[0], spl2->x[0]), right = std::min(spl1->x[spl1->N], spl2->x[spl2->N]);
    float eps = 0.01;
    if (((GetPoint(spl1, left) - GetPoint(spl2, left)) * (GetPoint(spl1, right) - GetPoint(spl2, right))) > eps) {
        flag_point = false;
        return 0;
    }
    while (right - left > eps) {
        float c = (left + right) / 2;
        if (abs(GetPoint(spl1, c) - GetPoint(spl2, c)) < eps)
            return c;
        if ((GetPoint(spl1, c) - GetPoint(spl2, c) * (GetPoint(spl1, right) - GetPoint(spl2, right))) < 0) {
            left = c;
        }
        else {
            right = c;
        }
    }
    return (left + right) / 2;
}


void Test(spline* spl){
    float start = spl->x[0];
    float end = spl->x[spl->N];
    float step = (end - start)/20;
    for(float s = start; s<=end; s+= step){
        //find k, where s in [x_k-1; x_k]
        int k;
        for(k=1; k<=spl->N; k++){
            if(s>=spl->x[k-1] && s<=spl->x[k]){
                break;
            }
        }
        float F = spl->y[k] + spl->b[k]*(s-spl->x[k]) + spl->c[k]*pow(s-spl->x[k], 2) + spl->d[k]*pow(s-spl->x[k], 3);
        printf("%f\t%f\t%f\n", s, GetPoint(spl, s), F);
    }
}


int main(){
    char filename[256];
    FILE* InFile = nullptr;
    while (InFile == nullptr) {
        printf("\nInput filename: ");
        scanf("%s", filename);
        InFile = fopen(filename, "rt");
    }
    spline* spl1 = MakeSpline(InFile);
    InFile = nullptr;
    while (InFile == nullptr) {
        printf("\nInput filename: ");
        scanf("%s", filename);
        InFile = fopen(filename, "rt");
    }
    spline* spl2 = MakeSpline(InFile);
    bool flag_point = true;
    float point = FindPoint(spl1, spl2, flag_point);
    if (!flag_point) {
       printf("There's no intersection point");
    }
    else {
        printf("Such point is x=%f, y=%f", point, GetPoint(spl1, point));
    }
    //Test(spl2);
    SplineDestructor(spl1);
    SplineDestructor(spl2);
    return 0;
}