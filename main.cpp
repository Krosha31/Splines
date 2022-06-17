#include <iostream>
#include <math.h>
#include <stdio.h>

const float eps = 0.01;


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
            return nullptr;
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
    return spl;
}


int FindInterval(spline* spl, float x) {
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
    return right;
}


float GetPoint(spline* spl, float x) {
    int right = FindInterval(spl, x);
    return spl->y[right] + spl->b[right]*(x - spl->x[right]) + spl->c[right]*pow(x - spl->x[right], 2) + spl->d[right]*pow(x - spl->x[right], 3);
}

//Производная
float GetDiff1(spline* spl, float x) {
    int right = FindInterval(spl, x);
    return spl->b[right] + 2 * spl->c[right] * (x - spl->x[right]) + 3 * spl->d[right] * pow(x - spl->x[right], 2);
}


float GetDiff2(spline* spl, float x) {
    int right = FindInterval(spl, x);
    return 2 * spl->c[right] + 6 * spl->d[right] * (x - spl->x[right]);
}


float FindPoint(spline* spl1, spline* spl2, bool &flag_point) {
    float left = std::max(spl1->x[0], spl2->x[0]), right = std::min(spl1->x[spl1->N], spl2->x[spl2->N]);
    // проверяем, пересекаются ли сплайны
    if (((GetPoint(spl1, left) - GetPoint(spl2, left)) * (GetPoint(spl1, right) - GetPoint(spl2, right))) > eps) {
        flag_point = false;
        return 0;
    }
    // ищем точку пересечения бинарным поиском
    while (right - left > eps) {
        float c = (left + right) / 2;
        //если нашли
        if (abs(GetPoint(spl1, c) - GetPoint(spl2, c)) < eps)
            return c;
        //если не нашли
        if (((GetPoint(spl1, c) - GetPoint(spl2, c)) * (GetPoint(spl1, right) - GetPoint(spl2, right))) < 0) {
            left = c;
        }
        else {
            right = c;
        }
    }
    //Возвращаем найденный х
    return (left + right) / 2;
}

//Функция, которую мы исследуем для нахождения расстояния. Это функция расстояний двух точек разных сплайнов. Соответственно,
// мы ищем ее минимум
float MinFunc(spline* spl1, spline* spl2, float x1, float x2) {
    return pow((x1 - x2), 2) + pow((GetPoint(spl1, x1) - GetPoint(spl2, x2)), 2);
}
//Частные производные функции, исследуемой на наимеьшее значение
float GetDiffMinFuncX1(spline* spl1, spline* spl2, float x1, float x2) {
    return 2 * (GetDiff1(spl1, x1) * (GetPoint(spl1, x1) - GetPoint(spl2, x2)) + x1 - x2);
}

float GetDiffMinFuncX2(spline* spl1, spline* spl2, float x1, float x2) {
    return 2 * (-GetDiff1(spl2, x2) * (GetPoint(spl1, x1) - GetPoint(spl2, x2)) - x1 + x2);
}

//Функция, которая исследуется для нахождения длины шага
float FuncForLambda(spline* spl1, spline* spl2, float x[], float grad[], float lambda) {
    float dop[2];
    dop[0] = x[0] - lambda * grad[0];
    dop[1] = x[1] - lambda * grad[1];
    return MinFunc(spl1, spl2, dop[0], dop[1]);
}

//Находим оптимальную длину шага методом наискорейшего спуска
float FindMin(spline* spl1, spline* spl2, float a, float b, float x[]) {
    const float fi=1.6180339887;
    float y_c[2], x_c[2];
    x_c[0] = b - ((b - a)/fi);
    x_c[1] = a + ((b - a)/fi);
    float grad[2];
    grad[0] = GetDiffMinFuncX1(spl1, spl2, x[0], x_c[1]);
    grad[1] = GetDiffMinFuncX2(spl1, spl2, x[0], x_c[1]);
    y_c[0] = FuncForLambda(spl1, spl2, x, grad, x_c[0]);
    y_c[1] = FuncForLambda(spl1, spl2, x, grad, x_c[1]);
    //Метод золотого сечения
    while (abs(b - a) > eps) {
        if (y_c[0] <= y_c[1]) {
            b = x_c[1];
            x_c[1] = x_c[0];
            x_c[0] = b - ((b - a) / fi);
            y_c[1] = y_c[0];
            y_c[0] = FuncForLambda(spl1, spl2, x, grad, x_c[0]);
        }
        else {
            a = x_c[0];
            x_c[0] = x_c[1];
            x_c[1] = a + ((b - a) / fi);
            y_c[0] = y_c[1];
            y_c[1] = FuncForLambda(spl1, spl2, x, grad, x_c[1]);
        }
    }
    return (a + b) / 2;
}

// Находим расстояние при помощи градиентного спуска
float DistanceBetweenSplines(spline* spl1, spline* spl2, float x1, float x2) {
    float pred = 1000;
    float x[2];
    x[0] = x1;
    x[1] = x2;
    while (abs(pred - MinFunc(spl1, spl2, x[0], x[1])) > eps) {
       float predx[2], grad[2];
       predx[0] = x[0];
       predx[1] = x[1];
       pred = MinFunc(spl1, spl2, predx[0], predx[1]);
       grad[0] = GetDiffMinFuncX1(spl1, spl2, predx[0], predx[1]);
       grad[1] = GetDiffMinFuncX2(spl1, spl2, predx[0], predx[1]);
       float lambda = FindMin(spl1, spl2, 0, 0.05, predx);
       x[0] = predx[0] - lambda * grad[0];
       x[1] = predx[0] - lambda * grad[1];
    }
    return pow(MinFunc(spl1, spl2, x[0], x[1]), 0.5);
}

//Выводим значения в файл для построения графиков в питоне
void Test(spline* spl, int number){
    float start = spl->x[0];
    float end = spl->x[spl->N];
    float step = 0.01;
    FILE* OutFile = nullptr;
    if (number == 1)
        OutFile = fopen("spline1.txt", "wt");
    else if (number == 2)
        OutFile = fopen("spline2.txt", "wt");
    for(float s = start; s<=end; s+= step){
        fprintf(OutFile, "%f\t%f\n", s, GetPoint(spl, s));
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
       printf("There's no intersection point\n");
        printf("Distance is %f", DistanceBetweenSplines(spl1, spl2, -1, -1));
    }
    else {
        printf("Such point is x=%f, y=%f\n", point, GetPoint(spl1, point));
    }
    Test(spl1, 1);
    Test(spl2, 2);
    SplineDestructor(spl1);
    SplineDestructor(spl2);
    return 0;
}