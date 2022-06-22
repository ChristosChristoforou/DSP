#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

//Number of samples ex. q = 3 => 2^3 points
#define q	3
//N-point FFT, iFFT ex. 1D = 000000000001B => 1<<3 => 000000001000B = 8D
#define N	(1<<q)		/* N-point FFT, iFFT */
#define PI  atan2(0, -1)

struct complex_polar{
    double magnitude;
    double angle;
};

struct complex_cartisian{
    double real;
    double imaginary;
};

typedef struct complex_number{
    struct complex_cartisian cartisian;
    struct complex_polar polar;
}complex;

typedef struct complex_array{
    size_t len;
    complex* array;
}cplx_buff;

//Complex Numbers

void cartisian_to_polar(complex* num);
void polar_to_cartisian(complex* num);
double randf(double min, double max);
void print_vector(const char* title, cplx_buff x);

//FFT

void _FFT(complex* v, size_t n, complex* tmp);
void FFT(cplx_buff* vector, cplx_buff* fft_vector);

//IFFT

void _IFFT(complex* v, size_t n, complex* tmp);
void IFFT(cplx_buff* vector, cplx_buff* fft_vector);

//Convolution

cplx_buff Convolution(cplx_buff* vector_1, cplx_buff* vector_2);

//Correlation

cplx_buff Correlation(cplx_buff* vector_1, cplx_buff* vector_2);

//Spectral density

cplx_buff Spectral_Density(cplx_buff* vector_1, cplx_buff* vector_2);

//Coherence

cplx_buff Coherence(cplx_buff* vector_1, cplx_buff* vector_2);

//Scaled correlation

cplx_buff Scaled_Corelation(cplx_buff* vector_1, cplx_buff* vector_2, size_t scale);


int main(){

    cplx_buff arr_1, arr_2;
    arr_1.len = N;
    arr_1.array = (complex*)malloc(sizeof(complex)*arr_1.len);
    arr_2.len = N;
    arr_2.array = (complex*)malloc(sizeof(complex)*arr_2.len);

    
    for (size_t i = 0; i < arr_1.len; i++){
        arr_1.array[i].cartisian.real = 0.125 * cos(2 * PI * i/(double)N);
        arr_1.array[i].cartisian.imaginary = 0;
        arr_2.array[i].cartisian.real = -0.125 * sin(2 * PI * i/(double)N);
        arr_2.array[i].cartisian.imaginary = 0;
        cartisian_to_polar(&arr_1.array[i]);
        cartisian_to_polar(&arr_2.array[i]);
    }

    cplx_buff conv = Scaled_Corelation(&arr_1, &arr_2, 2);
    print_vector("Convolution", conv);


    free(arr_1.array);
    free(arr_2.array);
    free(conv.array);
    return 0;
}


void cartisian_to_polar(complex* num){
    num->polar.magnitude = sqrt(pow(num->cartisian.real, 2) + pow(num->cartisian.imaginary, 2));
    num->polar.angle = atan2(num->cartisian.imaginary, num->cartisian.real) * 180/PI;
}

void polar_to_cartisian(complex* num){
    num->cartisian.real = num->polar.magnitude * cos(num->polar.angle);
    num->cartisian.real = num->polar.magnitude * sin(num->polar.angle);
}

double randf(double min, double max){
    double range = (max - min); 
    double div = RAND_MAX / range;
    return min + (rand() / div);
}

void print_vector(const char* title, cplx_buff x){
    printf("%s (dim = %ld):\n", title, x.len);
    for(size_t i = 0; i < x.len; i++ ){
        cartisian_to_polar(&x.array[i]);
        printf("Polar:     %5.2f exp(%5.2f i)\n\n", x.array[i].polar.magnitude, x.array[i].polar.angle);
    }
    printf("\n");
    return;
}

void _FFT(complex* v, size_t n, complex* tmp){
    if (n > 1){			/* otherwise, do nothing and return */
        complex z, w, *vo, *ve;
        ve = tmp;
        vo = tmp + n/2;
        for (size_t k = 0; k < n/2; k++){
            ve[k] = v[2*k];
            vo[k] = v[2*k + 1];
        }
        _FFT( ve, n/2, v );		/* FFT on even-indexed elements of v[] */
        _FFT( vo, n/2, v );		/* FFT on odd-indexed elements of v[] */
        for (size_t m = 0; m < n/2; m++) {
            w.cartisian.real = cos(2 * PI * m/(double)n);
            w.cartisian.imaginary = -sin(2 * PI * m/(double)n);
            z.cartisian.real = w.cartisian.real * vo[m].cartisian.real - w.cartisian.imaginary * vo[m].cartisian.imaginary;	/* Re(w*vo[m]) */
            z.cartisian.imaginary = w.cartisian.real * vo[m].cartisian.imaginary + w.cartisian.imaginary * vo[m].cartisian.real;	/* Im(w*vo[m]) */
            v[m].cartisian.real = ve[m].cartisian.real + z.cartisian.real;
            v[m].cartisian.imaginary = ve[m].cartisian.imaginary + z.cartisian.imaginary;
            v[m + n/2].cartisian.real = ve[m].cartisian.real - z.cartisian.real;
            v[m + n/2].cartisian.imaginary = ve[m].cartisian.imaginary - z.cartisian.imaginary;
        }
    }
}

void FFT(cplx_buff* vector, cplx_buff* fft_vector){
    cplx_buff proxy, tmp;

    proxy.len = vector->len;
    proxy.array = (complex*)malloc(sizeof(complex) * proxy.len);
    for (size_t i = 0; i < proxy.len; i++){
        proxy.array[i] = vector->array[i];
    }
    
    tmp.len = vector->len;
    tmp.array = (complex*)malloc(sizeof(complex) * tmp.len);

    _FFT(proxy.array, proxy.len, tmp.array);

    for (size_t i = 0; i < fft_vector->len; i++){
        fft_vector->array[i] = proxy.array[i];
        cartisian_to_polar(&fft_vector->array[i]);
    }
    free(proxy.array);
    free(tmp.array);

}

void _IFFT(complex* v, size_t n, complex* tmp){
  if(n > 1) {			/* otherwise, do nothing and return */
    complex z, w, *vo, *ve;
    ve = tmp;
    vo = tmp + n/2;
    for(size_t k = 0; k < n/2; k++) {
      ve[k] = v[2*k];
      vo[k] = v[2*k + 1];
    }
    _IFFT( ve, n/2, v );		/* FFT on even-indexed elements of v[] */
    _IFFT( vo, n/2, v );		/* FFT on odd-indexed elements of v[] */
    for(size_t m = 0; m < n/2; m++) {
      w.cartisian.real = cos(2 * PI * m/(double)n);
      w.cartisian.imaginary = sin(2 * PI * m/(double)n);
      z.cartisian.real = w.cartisian.real * vo[m].cartisian.real - w.cartisian.imaginary * vo[m].cartisian.imaginary;	/* Re(w*vo[m]) */
      z.cartisian.imaginary = w.cartisian.real * vo[m].cartisian.imaginary + w.cartisian.imaginary * vo[m].cartisian.real;	/* Im(w*vo[m]) */
      v[m].cartisian.real = ve[m].cartisian.real + z.cartisian.real;
      v[m].cartisian.imaginary = ve[m].cartisian.imaginary + z.cartisian.imaginary;
      v[m + n/2].cartisian.real = ve[m].cartisian.real - z.cartisian.real;
      v[m + n/2].cartisian.imaginary = ve[m].cartisian.imaginary - z.cartisian.imaginary;
    }
  }
}

void IFFT(cplx_buff* vector, cplx_buff* ifft_vector){
    cplx_buff proxy, tmp;

    proxy.len = vector->len;
    proxy.array = (complex*)malloc(sizeof(complex) * proxy.len);
    for (size_t i = 0; i < proxy.len; i++){
        proxy.array[i] = vector->array[i];
    }
    
    tmp.len = vector->len;
    tmp.array = (complex*)malloc(sizeof(complex) * tmp.len);

    _IFFT(proxy.array, proxy.len, tmp.array);

    for (size_t i = 0; i < ifft_vector->len; i++){
        ifft_vector->array[i] = proxy.array[i];
        cartisian_to_polar(&ifft_vector->array[i]);
    }
    free(proxy.array);
    free(tmp.array);

}

cplx_buff Convolution(cplx_buff* vector_1, cplx_buff* vector_2){
    cplx_buff conv;
    conv.len = vector_1->len;
    conv.array = (complex*)malloc(sizeof(complex) * conv.len);

    cplx_buff proxy_1, proxy_2, prod;
    proxy_1.len = vector_1->len;
    proxy_1.array = (complex*)malloc(sizeof(complex) * proxy_1.len);
    proxy_2.len = vector_2->len;
    proxy_2.array = (complex*)malloc(sizeof(complex) * proxy_2.len);
    prod.len = vector_1->len;
    prod.array = (complex*)malloc(sizeof(complex) * prod.len);


    FFT(vector_1, &proxy_1);
    FFT(vector_2, &proxy_2);

    for (size_t i = 0; i < prod.len; i++){
        prod.array[i].cartisian.real = proxy_1.array[i].cartisian.real * proxy_2.array[i].cartisian.real - proxy_1.array[i].cartisian.imaginary * proxy_2.array[i].cartisian.imaginary; /* Re(proxy_1[i]*proxy_2[i]) */
        prod.array[i].cartisian.imaginary = proxy_1.array[i].cartisian.real * proxy_2.array[i].cartisian.imaginary + proxy_1.array[i].cartisian.imaginary * proxy_2.array[i].cartisian.real; /* Im(proxy_1[i]*proxy_2[i]) */
    }
    
    IFFT(&prod, &conv);

    free(proxy_1.array);
    free(proxy_2.array);
    free(prod.array);
    return conv;
}

cplx_buff Correlation(cplx_buff* vector_1, cplx_buff* vector_2){
    cplx_buff corr;
    corr.len = vector_1->len;
    corr.array = (complex*)malloc(sizeof(complex) * corr.len);

    cplx_buff proxy_1, proxy_2, prod, tmp;
    proxy_1.len = vector_1->len;
    proxy_1.array = (complex*)malloc(sizeof(complex) * proxy_1.len);
    proxy_2.len = vector_2->len;
    proxy_2.array = (complex*)malloc(sizeof(complex) * proxy_2.len);
    prod.len = vector_1->len;
    prod.array = (complex*)malloc(sizeof(complex) * prod.len);


    FFT(vector_1, &proxy_1);
    FFT(vector_2, &proxy_2);

    for (size_t i = 0; i < prod.len; i++){
        prod.array[i].cartisian.real = proxy_1.array[i].cartisian.real * proxy_2.array[i].cartisian.real + proxy_1.array[i].cartisian.imaginary * proxy_2.array[i].cartisian.imaginary; /* Re(~{proxy_1[i]}*proxy_2[i]) */
        prod.array[i].cartisian.imaginary = proxy_1.array[i].cartisian.real * proxy_2.array[i].cartisian.imaginary - proxy_1.array[i].cartisian.imaginary * proxy_2.array[i].cartisian.real; /* Im(~{proxy_1[i]}*proxy_2[i]) */
    }
    
    IFFT(&prod, &corr);

    free(proxy_1.array);
    free(proxy_2.array);
    free(prod.array);
    return corr;
}

cplx_buff Spectral_Density(cplx_buff* vector_1, cplx_buff* vector_2){
    cplx_buff Spectr;
    Spectr.len = vector_1->len;
    Spectr.array = (complex*)malloc(sizeof(complex) * Spectr.len);

    cplx_buff proxy_1, proxy_2, prod, tmp;
    proxy_1.len = vector_1->len;
    proxy_1.array = (complex*)malloc(sizeof(complex) * proxy_1.len);
    proxy_2.len = vector_2->len;
    proxy_2.array = (complex*)malloc(sizeof(complex) * proxy_2.len);
    prod.len = vector_1->len;
    prod.array = (complex*)malloc(sizeof(complex) * prod.len);
    tmp.len = vector_1->len;
    tmp.array = (complex*)malloc(sizeof(complex) * tmp.len);

    FFT(vector_1, &proxy_1);
    FFT(vector_2, &proxy_2);

    for (size_t i = 0; i < prod.len; i++){
        prod.array[i].cartisian.real = proxy_1.array[i].cartisian.real * proxy_2.array[i].cartisian.real + proxy_1.array[i].cartisian.imaginary * proxy_2.array[i].cartisian.imaginary; /* Re(~{proxy_1[i]}*proxy_2[i]) */
        prod.array[i].cartisian.imaginary = proxy_1.array[i].cartisian.real * proxy_2.array[i].cartisian.imaginary - proxy_1.array[i].cartisian.imaginary * proxy_2.array[i].cartisian.real; /* Im(~{proxy_1[i]}*proxy_2[i]) */
    }
    
    IFFT(&prod, &Spectr);
    _FFT(Spectr.array, Spectr.len, tmp.array);

    free(proxy_1.array);
    free(proxy_2.array);
    free(prod.array);
    free(tmp.array);
    return Spectr;
}

cplx_buff Coherence(cplx_buff* vector_1, cplx_buff* vector_2){
    cplx_buff cohr;
    cohr.len = vector_1->len;
    cohr.array = (complex*)malloc(sizeof(complex) * cohr.len);

    cplx_buff spctr_xy, spctr_xx, spctr_yy, tmp;
    spctr_xy = Spectral_Density(vector_1, vector_2);
    spctr_xx = Spectral_Density(vector_1, vector_1);
    spctr_yy = Spectral_Density(vector_2, vector_2);
    tmp.len = vector_1->len;
    tmp.array = (complex*)malloc(sizeof(complex) * tmp.len);
    
    for (size_t k = 0; k < spctr_xy.len; k++){
        //Numerator
        spctr_xy.array[k].cartisian.real = pow(sqrt(pow(spctr_xy.array[k].cartisian.real, 2) + pow(spctr_xy.array[k].cartisian.imaginary, 2)), 2);
        spctr_xy.array[k].cartisian.real = 0;

        //Denominator
        tmp.array[k].cartisian.real = spctr_xx.array[k].cartisian.real * spctr_yy.array[k].cartisian.real - spctr_xx.array[k].cartisian.imaginary * spctr_yy.array[k].cartisian.imaginary;	/* Re(spctr_xx * spctr_yy) */
        tmp.array[k].cartisian.imaginary = spctr_xx.array[k].cartisian.real * spctr_yy.array[k].cartisian.imaginary + spctr_xx.array[k].cartisian.imaginary * spctr_yy.array[k].cartisian.real; /* Im(spctr_xx * spctr_yy) */
        cartisian_to_polar(&tmp.array[k]);
        tmp.array[k].polar.magnitude = 1/tmp.array[k].polar.magnitude;
        tmp.array[k].polar.angle = -tmp.array[k].polar.angle;
        polar_to_cartisian(&tmp.array[k]);

    }

  
    for (size_t k = 0; k < cohr.len; k++){
        cohr.array[k].cartisian.real = spctr_xy.array[k].cartisian.real * tmp.array[k].cartisian.real - spctr_xy.array[k].cartisian.imaginary * tmp.array[k].cartisian.imaginary;	/* Re(|spctr_xy|^2 * tmp) */
        cohr.array[k].cartisian.imaginary = spctr_xy.array[k].cartisian.real * tmp.array[k].cartisian.imaginary + spctr_xy.array[k].cartisian.imaginary * tmp.array[k].cartisian.real; /* Im(|spctr_xy|^2 * tmp) */
    }


    free(spctr_xy.array);
    free(spctr_xx.array);
    free(spctr_yy.array);
    free(tmp.array);
    return cohr;
}

cplx_buff Scaled_Corelation(cplx_buff* vector_1, cplx_buff* vector_2, size_t scale){

    size_t K = round(vector_1->len/scale);

    cplx_buff scorel;
    scorel.len = vector_1->len - vector_1->len % K;
    scorel.array = (complex*)malloc(sizeof(complex) * scorel.len);

    // https://en.wikipedia.org/wiki/Pearson_correlation_coefficient#For_a_sample

    complex mean_1, mean_2, tmp1, tmp2, median_prod, ssd_prod;
    cplx_buff sx, sy;
    sx.len = scale;
    sx.array = (complex*)malloc(sizeof(complex) * sx.len);
    sy.len = K;
    sy.array = (complex*)malloc(sizeof(complex) * sy.len);

    for (size_t k = 0; k < K; k++){
        mean_1.cartisian.real = 0;
        mean_1.cartisian.imaginary = 0;
        mean_2.cartisian.real = 0;
        mean_2.cartisian.imaginary = 0;
        for (size_t i = 0; i < scale; i++){
            mean_1.cartisian.real += vector_1->array[i + k].cartisian.real;
            mean_1.cartisian.imaginary += vector_1->array[i + k].cartisian.imaginary;
            mean_2.cartisian.real += vector_2->array[i + k].cartisian.real;
            mean_2.cartisian.imaginary += vector_2->array[i + k].cartisian.imaginary;
        }
        mean_1.cartisian.real = mean_1.cartisian.real / scale;
        mean_1.cartisian.imaginary = mean_1.cartisian.imaginary / scale;
        mean_2.cartisian.real = mean_2.cartisian.real / scale;
        mean_2.cartisian.imaginary = mean_2.cartisian.imaginary / scale;

        sx.array[k].cartisian.real = 0;
        sx.array[k].cartisian.imaginary = 0;
        sy.array[k].cartisian.real = 0;
        sy.array[k].cartisian.imaginary = 0;
        for (size_t i = 0; i < scale; i++){
            tmp1.cartisian.real = vector_1->array[i + k].cartisian.real - mean_1.cartisian.real;
            tmp1.cartisian.imaginary = vector_1->array[i + k].cartisian.imaginary - mean_1.cartisian.imaginary;

            sx.array[k].cartisian.real += tmp1.cartisian.real * tmp1.cartisian.real - tmp1.cartisian.imaginary * tmp1.cartisian.imaginary; /* Re((xi-X)^2) */
            sx.array[k].cartisian.imaginary += tmp1.cartisian.real * tmp1.cartisian.imaginary + tmp1.cartisian.imaginary * tmp1.cartisian.real; /* Im((xi-X)^2) */

            tmp1.cartisian.real = vector_2->array[i + k].cartisian.real - mean_2.cartisian.real;
            tmp1.cartisian.imaginary = vector_2->array[i + k].cartisian.imaginary - mean_2.cartisian.imaginary;
            
            sy.array[k].cartisian.real += tmp1.cartisian.real * tmp1.cartisian.real - tmp1.cartisian.imaginary * tmp1.cartisian.imaginary; /* Re((yi-Y)^2) */
            sy.array[k].cartisian.imaginary += tmp1.cartisian.real * tmp1.cartisian.imaginary + tmp1.cartisian.imaginary * tmp1.cartisian.real; /* Im((yi-Y)^2) */
        }

        // sqrt(sx/n-1)
        sx.array[k].cartisian.real = sx.array[k].cartisian.real/(scale - 1);
        sx.array[k].cartisian.imaginary = sx.array[k].cartisian.imaginary/(scale - 1);
        cartisian_to_polar(&sx.array[k]);
        sx.array[k].polar.magnitude = sqrt(sx.array[k].polar.magnitude);
        sx.array[k].polar.angle = sx.array[k].polar.angle / 2;
        polar_to_cartisian(&sx.array[k]);
        
        // sqrt(sy/n-1)
        sy.array[k].cartisian.real = sy.array[k].cartisian.real/(scale - 1);
        sy.array[k].cartisian.imaginary = sy.array[k].cartisian.imaginary/(scale - 1);
        cartisian_to_polar(&sy.array[k]);
        sy.array[k].polar.magnitude = sqrt(sy.array[k].polar.magnitude);
        sy.array[k].polar.angle = sy.array[k].polar.angle / 2;
        polar_to_cartisian(&sy.array[k]);

        // (n-1)sx*sy
        ssd_prod.cartisian.real = (scale - 1) * (sx.array[K].cartisian.real * sy.array[K].cartisian.real - sx.array[K].cartisian.imaginary * sy.array[K].cartisian.imaginary); /* (scale - 1) * Re(sx*sy) */
        ssd_prod.cartisian.imaginary = (scale - 1) * (sx.array[K].cartisian.real * sy.array[K].cartisian.imaginary + sx.array[K].cartisian.imaginary * sy.array[K].cartisian.real); /* scale *Im((yi-Y)^2) */

        //Sum(xiyi)
        tmp2.cartisian.real = 0;
        tmp2.cartisian.imaginary = 0;
        for (size_t i = 0; i < scale; i++){
            tmp2.cartisian.real += vector_1->array[i + k].cartisian.real * vector_2->array[i + k].cartisian.real - vector_1->array[i + k].cartisian.imaginary * vector_2->array[i + k].cartisian.imaginary; /* Re((yi-Y)^2) */
            tmp2.cartisian.imaginary += vector_1->array[i + k].cartisian.real * vector_2->array[i + k].cartisian.imaginary + vector_1->array[i + k].cartisian.imaginary * vector_2->array[i + k].cartisian.real; /* Im((yi-Y)^2) */
        }
        median_prod.cartisian.real = scale * (mean_1.cartisian.real * mean_2.cartisian.real - mean_1.cartisian.imaginary * mean_2.cartisian.imaginary); /* scale * Re(XY) */
        median_prod.cartisian.imaginary = scale * (mean_1.cartisian.real * mean_2.cartisian.imaginary + mean_1.cartisian.imaginary * mean_2.cartisian.real); /* scale *Im((yi-Y)^2) */

        scorel.array[k].cartisian.real = tmp2.cartisian.real - median_prod.cartisian.real;
        scorel.array[k].cartisian.imaginary = tmp2.cartisian.imaginary - median_prod.cartisian.imaginary;
        cartisian_to_polar(&scorel.array[k]);
        cartisian_to_polar(&ssd_prod);
        scorel.array[k].polar.magnitude = scorel.array[k].polar.magnitude / ssd_prod.polar.magnitude;
        scorel.array[k].polar.angle = scorel.array[k].polar.angle - ssd_prod.polar.angle;
        polar_to_cartisian(&scorel.array[k]);
    }
    
    free(sx.array);
    free(sy.array);
    return scorel;
}