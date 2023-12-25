#include <iostream>
#include <vector>
#include <math.h>
#define PI 3.1415926

using namespace std;

class complex{
	public:
		float real;
		float imag;
		
		complex() : real(0.0), imag(0.0) {}
		complex(float real, float imag) : real(real), imag(imag) {}
		complex operator-(complex &other);		
		complex operator+(complex &other);		
		friend complex mul(complex, complex);
		complex operator*(complex &other);
		float twoNorm();
};

complex mul(complex a, complex b)
{
	complex result;
	result.real = a.real * b.real - a.imag * b.imag;
	result.imag = a.real * b.imag + a.imag * b.real;

	return result;
}

complex complex::operator+(complex &other)
{
	return complex(real + other.real, imag + other.imag);
}

complex complex::operator-(complex &other)
{
	return complex(real - other.real, imag - other.imag);
}

complex complex::operator*(complex &other)
{
	return complex(real*other.real - imag*other.imag, imag*other.real + real*other.imag);
}

float complex::twoNorm()
{
	return (real * real + imag * imag);
}

uint32_t bit_reverse(uint32_t x)
{
	uint32_t n = x;
	n = ((n & 0xffff0000) >> 16) | ((n & 0x0000ffff) << 16);
    n = ((n & 0xff00ff00) >>  8) | ((n & 0x00ff00ff) <<  8);
    n = ((n & 0xf0f0f0f0) >>  4) | ((n & 0x0f0f0f0f) <<  4);
    n = ((n & 0xcccccccc) >>  2) | ((n & 0x33333333) <<  2);
    n = ((n & 0xaaaaaaaa) >>  1) | ((n & 0x55555555) <<  1);
	n = n >> 28;
    return n;
}

void scrambling(vector<complex> &a, vector<complex> &A)
{
	int n = a.size();

	for(uint32_t k = 0; k < n - 1; k++)
	{
		uint32_t f = bit_reverse(k);
		A.at(k) = a.at(f);
		A.at(f) = a.at(k);
	}
}

vector<complex> iter_fft(vector<complex> a)
{
	int n = a.size();
	vector<complex> A(n);
	scrambling(a, A);

	for(int s = 1; s <= 4; s++)
	{
		int m = 1, g = s;
		while(g != 0)
		{
			m *= 2;
			g--;
		}
		complex wm(cos(2*PI/m), sin(2*PI/m));
		for(int k = 0; k <= n - 1; k+=m)
		{
			complex w(1.0, 0.0);
			for(int j = 0; j <= m/2 - 1; j++)
			{
				complex t, u;
				t = mul(w, A.at(k+j+m/2));
				u = A.at(k+j);
				A.at(k+j) = u + t;
				A.at(k+j+m/2) = u - t;
				w = mul(w, wm);
			}
		}
	}

	return A;
}

float float_to_fixed (float fp, int intWord, int fracWord) {
	float upperBound = pow(2.0, intWord - 1) - pow(2.0, -fracWord);
	float lowerBound = -pow(2.0, intWord - 1);
	float scaleFactor = pow(2.0, fracWord);

	float temp = fp * scaleFactor;
	float result = floor(fp * scaleFactor) / scaleFactor;

	if(result > upperBound)
		result = upperBound;
	else if(result < lowerBound)
		result = lowerBound;
	else
		result = result;

	return result;
}

vector<complex> recur_fft(vector<complex> a, bool Isfixed, int fracWord)
{
	int fftSize = a.size();

	if(fftSize == 1)
	{
		return a;
	}
	else
	{
		complex twiddleN(cos(2*PI/fftSize), -sin(2*PI/fftSize));
		if(Isfixed)
		{
			float wn_re = float_to_fixed(cos(2*PI/fftSize), 1, fracWord);
			float wn_im = float_to_fixed(-sin(2*PI/fftSize), 1, fracWord);
			complex temp(wn_re, wn_im);
			twiddleN = temp;
		}
		complex w(1.0, 0.0);
		vector<complex> aEven(fftSize/2);
		vector<complex> aOdd(fftSize/2);
		vector<complex> yEven(fftSize/2);
		vector<complex> yOdd(fftSize/2);
		for(int i = 0; i < fftSize; i+=2)
		{
			aEven.at(i/2) = a.at(i);
			aOdd.at(i/2) = a.at(i+1);
		}
		yEven = recur_fft(aEven, Isfixed, fracWord);
		yOdd = recur_fft(aOdd, Isfixed, fracWord);
		vector<complex> y(fftSize);
		for(int i = 0; i < fftSize/2; i++)
		{
			complex temp;
			temp = mul(w, yOdd.at(i));
			y.at(i) = temp + yEven.at(i);
			y.at(i+fftSize/2) = yEven.at(i) - temp;
			w = mul(w, twiddleN);
		}

		return y;
	}
}

int main(int argc, char* argv[])
{
	int N = 16;
	vector<complex> x(N), y(N), qX(N), qY(N), fY(N);
	vector<float> snr1(N), snr2(N);
	
	for(int j = 1; j < N; j++)
	{
		float fPower = 0, diffPowerQ = 0, diffPowerF = 0;
		for(int i = 0; i < N; i++)
		{
			// world length starts from 2 bits, contained 1 sign bit, 1 integer bit and j fractional bits
			float t = float_to_fixed(cos(2*PI*i/N), 1, j);
			complex temp(cos(2*PI*i/N), 0.0);
			x.at(i) = temp;
			complex temp1(t, 0.0);
			qX.at(i) = temp1;
			//cout << x.at(i).real << ", " << x.at(i).imag << endl;
		}
		y = recur_fft(x, false, 0);
		qY = recur_fft(qX, false, 0);
		fY = recur_fft(qX, true, j);
		for(int i = 0; i < N; i++)
		{
			complex temp;
			fPower += y.at(i).twoNorm();
			temp = y.at(i) - qY.at(i);
			diffPowerQ += temp.twoNorm();
			temp = y.at(i) - fY.at(i);
			diffPowerF += temp.twoNorm();
		}
		//cout << fPower << endl;
		cout << 10 * log10((fPower - diffPowerQ) / fPower) << ", " << 10 * log10((fPower - diffPowerF) / fPower) << endl;
	}
	return 0;
}
