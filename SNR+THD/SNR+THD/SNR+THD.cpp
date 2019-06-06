// SNR+THD.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include <fstream>
#include <cmath>
#include <fftw3.h>

typedef struct wavhead {
	char RIFF_ID[4];
	uint32_t RIFF_Size;
	char RIFF_Type[4];
	char Format_ID[4];
	uint32_t Format_Size;
	uint16_t AudioFormat;
	uint16_t NumChannels;
	uint32_t SampleRate;
	uint32_t ByteRate;
	uint16_t BlockAlign;
	uint16_t BitsPerSample;
	uint16_t BlockAlign2;
	char Data_ID[4];
	uint16_t Data_Size_Low16;
	uint16_t Data_Size_High16;
}WavInfo;

const double pi = 3.141592653589793238462643383279502884197;

// ----------------------------------------------------------------
// the Hilbert Transform
// input: 
// output:
// ----------------------------------------------------------------
void Hilbert_r2c(const double* in, fftw_complex* out, const uint32_t N) {
	uint32_t i = 0;
	for (i = 0; i < N; i++) {
		out[i][0] = in[i];
		out[i][1] = 0;
	}
	fftw_execute(fftw_plan_dft_1d(N, out, out, FFTW_FORWARD, FFTW_ESTIMATE));
	for (i = 0; i < N; i++) {
		if (i == 0 || 2 * i == N)continue;
		if (2 * i < N) {
			out[i][0] *= 2.0;
			out[i][1] *= 2.0;
		}
		else {
			out[i][0] = 0.0;
			out[i][1] = 0.0;
		}
	}
	fftw_execute(fftw_plan_dft_1d(N, out, out, FFTW_BACKWARD, FFTW_ESTIMATE));
}

// ----------------------------------------------------------------
// the 0-order modified Bessel function of the first kind
// ----------------------------------------------------------------
double bessel_i0(double x) {
	double r = 1.0, s = 1.0, i = 1.0;
	x /= 2.0;
	x *= x;
	for (i = 1.0; s > 1e-10; i++) {
		s *= x;
		s /= (i * i);
		r += s;
	}
	return r;
}

// ----------------------------------------------------------------
// ----------------------------------------------------------------
double kaiser_win(uint32_t n, uint32_t N) {
	double a = 2.0 * n / (N - 1) - 1;
	return bessel_i0(38 * sqrt(1.0 - a * a)) / bessel_i0(38);
}

// ----------------------------------------------------------------
// ----------------------------------------------------------------
void Power(double& Ps, double& Pn, double& Ph, double* PSD, const double* in, const uint32_t N) {

	uint32_t Np = N / 2 + 1;
	double* in_windowed = (double*)fftw_malloc(N * sizeof(double));
	fftw_complex* s = (fftw_complex*)fftw_malloc((N / 2 + 1) * sizeof(fftw_complex));

	int i = 0, j = 0, imax = 0, il = 0, ir = 0;
	double si = 0.0, sw = 0.0, sw2 = 0.0, ff = 0.0, w = 0.0, rbw = 0.0;
	for (i = 0; i < N; i++) {
		si += in[i];
	}
	si /= N;
	for (i = 0; i < N; i++) {
		w = kaiser_win(i, N);
		sw += w;
		sw2 += w * w;
		in_windowed[i] = (in[i] - si) * w;
	}
	rbw = sw2 / (sw * sw);
	fftw_execute(fftw_plan_dft_r2c_1d(N, in_windowed, s, FFTW_ESTIMATE));
	for (i = 0; i < N / 2 + 1; i++) {
		PSD[i] = 2.0 * (s[i][0] * s[i][0] + s[i][1] * s[i][1]) / (N * sw2);
		Pn += PSD[i];
		if (PSD[i] > PSD[imax])imax = i;
	}
	
	double P[7] = { 0.0 };
	Ps = Ph = 0.0;
	for (j = 0; j < 7; j++) {
		switch (j) {
		case 0: i = 0; break;
		case 1: i = imax; break;
		default: i = (uint32_t)round(j * ff); break;
		}
		int ic = i;
		il = (i == 0) ? 0 : (i - 1);
		ir = (i == N / 2) ? (N / 2) : (i + 1);
		i = (PSD[i] > PSD[il]) ? ((PSD[i] > PSD[ir]) ? i : ir) : ((PSD[il] > PSD[ir]) ? il : ir);
		il = i - 1;
		ir = i + 1;
		while (il >= 0 && PSD[il] < PSD[il + 1])il--;
		il++;
		while (ir <= N / 2 && PSD[ir] < PSD[ir - 1])ir++;
		ir--;

		for (i = il; i <= ir; i++) {
			P[j] += PSD[i];
			Pn -= PSD[i];
			Np--;
			if (j == 1) {
				ff += i * PSD[i];
			}
		}
		if (j == 1)ff /= P[1];
		if (P[j] < rbw * PSD[idx])P[j] = rbw * Psd[ic];
		if (j > 1)Ph += P[j];
	}
	Ps = P[1];
	Pn *= (N / 2 + 1);
	Pn /= (double)Np;
}

int main()
{
	//16-bit mono PCM wav only
	FILE* wavfile;
	char path[256] = "";
	WavInfo info1;
	std::cout <<  "Wav file path (16-bit mono PCM wav only):" << std::endl;
	std::cin >> path;
	fopen_s(&wavfile, path, "rb");
	fread((void*)& info1, sizeof(WavInfo), 1, wavfile);

	/**/
	std::cout << info1.RIFF_ID[0] << info1.RIFF_ID[1] << info1.RIFF_ID[2] << info1.RIFF_ID[3] << std::endl
		<< info1.RIFF_Size << std::endl
		<< info1.RIFF_Type[0] << info1.RIFF_Type[1] << info1.RIFF_Type[2] << info1.RIFF_Type[3] << std::endl
		<< info1.Format_ID[0] << info1.Format_ID[1] << info1.Format_ID[2] << info1.Format_ID[3] << std::endl
		<< info1.Format_Size << std::endl
		<< info1.AudioFormat << std::endl
		<< info1.NumChannels << std::endl
		<< info1.SampleRate << std::endl
		<< info1.ByteRate << std::endl
		<< info1.BlockAlign << std::endl
		<< info1.BitsPerSample << std::endl
		<< info1.BlockAlign2 << std::endl
		<< info1.Data_ID[0] << info1.Data_ID[1] << info1.Data_ID[2] << info1.Data_ID[3] << std::endl
		<< info1.Data_Size_Low16 << std::endl
		<< info1.Data_Size_High16 << std::endl;
	std::cout << std::endl;
	/**/

	//Number of samples
	uint32_t N = (uint32_t)info1.Data_Size_High16;
	N <<= 16;
	N += (uint32_t)info1.Data_Size_Low16;
	N /= sizeof(int16_t);

	int16_t* rawdata = (int16_t*)fftw_malloc(N * sizeof(int16_t));
	double* formdata = (double*)fftw_malloc(N * sizeof(double));
	//fftw_complex* s = (fftw_complex*)fftw_malloc((N / 2 + 1) * sizeof(fftw_complex));
	double* PSD = (double*)fftw_malloc((N / 2 + 1) * sizeof(double));

	fread((void*)rawdata, sizeof(int16_t), N, wavfile);
	uint32_t i;
	for (i = 0; i < N; i++) {
		formdata[i] = (double)rawdata[i] / (double)(1 << 15);
	}
	fftw_free((void*)rawdata); rawdata = NULL;
	//fftw_execute(fftw_plan_dft_r2c_1d(N, formdata, s, FFTW_ESTIMATE));
	
	double Ps = 0.0, Pn = 0.0, Ph = 0.0;
	Power(Ps, Pn, Ph, PSD, formdata, N);
	double SNR = 10.0 * log10(Ps / Pn), THD = 10.0 * log10(Ph / Ps);

	/*
	double P = 0, Pm = 0, Ps = 0, Pn = 0, Ph = 0, FundFreq = 0;
	uint32_t fi = 1;
	for (i = 1; i <= (N / 2 + 1); i++) {
		P = (s[i][0] * s[i][0] + s[i][1] * s[i][1]);
		Pn += P;
		if (P > Pm) {
			Pm = P;
			fi = i;
		}
	}
	FundFreq = fi * (double)info1.SampleRate / (double)N;
	Ps = (s[fi][0] * s[fi][0] + s[fi][1] * s[fi][1]);
	Pn -= Ps;
	for (i = 2 * fi; i <= 6 * fi; i += fi) {
		Ph += (s[i][0] * s[i][0] + s[i][1] * s[i][1]); std::cout << i / fi << "\t" << (s[i][0] * s[i][0] + s[i][1] * s[i][1])/Ps << std::endl;
	}
	Pn -= Ph;
	double SNR = 10 * log10(Ps / Pn), THD = 10 * log10(Ph / Ps);
	*/

	fftw_free((void*)PSD); PSD = NULL;
	fftw_free((void*)formdata); formdata = NULL;
	fclose(wavfile);

	//std::cout << "wu" << std::endl;
	std::cout << std::endl;
	std::cout << "Signal Power: " << Ps << std::endl;// 10 * log10(Ps) << "dB" << std::endl;
	std::cout << "Noise Power: " << Pn << std::endl;// 10 * log10(Pn) << "dB" << std::endl;
	std::cout << "Harmonic Power: " << Ph << std::endl;// 10 * log10(Ph) << "dB" << std::endl;
	//std::cout << "Fundamental Frequency: " << FundFreq << "Hz" << std::endl;
	std::cout << "SNR: " << SNR << "dB" << std::endl;
	std::cout << "THD: " << THD << "dB" << std::endl << std::endl;
	system("pause");
	return 0;
}

// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门提示: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
