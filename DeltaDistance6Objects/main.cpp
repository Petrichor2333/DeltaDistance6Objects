// DeltaDistance6Objects.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include "pch.h"
#include "config.h"
#include "DDFunc.h"
#include "SolvingMinimum.h"
#include "SolvingMinimumMove.h"
#include "GenerateRandom.h"

int main() {

	srand((unsigned)time(NULL));

	double time_cost;
	double sum_time_cost = 0;
	int num_of_iter = 100000;
	int num_of_failed = 0;
	double interval_boundary = 5;
	double Velocity_boundary = 3;

	int *NumIters = new int[num_of_iter];
	double *TimeCosts = new double[num_of_iter];


	vector<int> NumRange;
	NumRange.push_back(0);
	NumRange.push_back(5);
	NumRange.push_back(10);
	NumRange.push_back(15);
	NumRange.push_back(20);
	NumRange.push_back(30);
	NumRange.push_back(50);
	NumRange.push_back(55);

	vector<int> NumIter(NumRange.size() - 1, 0);


	for (int i = 0;i < num_of_iter;i++) {

		double R1[3][3], R2[3][3];
		double p1[3], p2[3];
		double v1[3], v2[3];

		GenerateRandomRotTra(R1, p1, interval_boundary);
		GenerateRandomRotTra(R2, p2, interval_boundary);
		GenerateRandomVel(v1, Velocity_boundary);
		GenerateRandomVel(v2, Velocity_boundary);

		//Select 2 objects two be tested

		ObjectDDFuncInfo* objA = &Object1DDFuncInfo;
		//ObjectDDFuncInfo* objA = &Object2DDFuncInfo;
		//ObjectDDFuncInfo* objA = &Object3DDFuncInfo;
		//ObjectDDFuncInfo* objA = &Object4DDFuncInfo;
		//ObjectDDFuncInfo* objA = &Object5DDFuncInfo;
		//ObjectDDFuncInfo* objA = &Object6DDFuncInfo;

		//ObjectDDFuncInfo* objB = &Object1DDFuncInfo;
		ObjectDDFuncInfo* objB = &Object2DDFuncInfo;
		//ObjectDDFuncInfo* objB = &Object3DDFuncInfo;
		//ObjectDDFuncInfo* objB = &Object4DDFuncInfo;
		//ObjectDDFuncInfo* objB = &Object5DDFuncInfo;
		//ObjectDDFuncInfo* objB = &Object6DDFuncInfo;

		double x[3];
		x[0] = (p1[0] + p2[0]) / 2;
		x[1] = (p1[1] + p2[1]) / 2;
		x[2] = (p1[2] + p2[2]) / 2;


		ResultsOfSolvingMinimum Result;				//Collision detection for stationary objects
		//ResultsOfSolvingMinimumMove Result;		//Collision detection for moving objects

		auto start = chrono::system_clock::now();
		bool SolSucc = SolveForCollision(objA, objB, x, &Result, R1, p1, R2, p2);					//Collision detection for stationary objects
		//bool SolSucc = SolveForCollisionMove(objA, objB, x, &Result, R1, p1, R2, p2, v1, v2);		//Collision detection for moving objects
		auto end = chrono::system_clock::now();
		auto duration = chrono::duration_cast<chrono::microseconds>(end - start);
		time_cost = double(duration.count());
		sum_time_cost = sum_time_cost + time_cost;

		NumIters[i] = Result.num_iter;
		TimeCosts[i] = time_cost;


		for (int k = 0;k < NumIter.size();k++) {
			if (Result.num_iter <= NumRange[k + 1]) {
				NumIter[k]++;
			}
		}

		if (Result.if_success == 0) {
			num_of_failed++;
		}


		if (i % 10000 == 0 && i > 0) {
			cout << "iteration: " << i << endl;
		}

	}

	double SumNumIters = 0, SumTimeCosts = 0;
	for (int i = 0;i < num_of_iter;i++) {
		SumNumIters += (double)(NumIters[i]);
		SumTimeCosts += TimeCosts[i];
	}
	double AveNumIters = SumNumIters / num_of_iter;
	double AveTimeCosts = SumTimeCosts / num_of_iter;

	double SDNumIters = 0, SDTimeCosts = 0;
	for (int i = 0;i < num_of_iter;i++) {
		SDNumIters += ((double)(NumIters[i]) - AveNumIters)*((double)(NumIters[i]) - AveNumIters);
		SDTimeCosts += (TimeCosts[i] - AveTimeCosts)*(TimeCosts[i] - AveTimeCosts);
	}
	SDNumIters = sqrt(SDNumIters / num_of_iter);
	SDTimeCosts = sqrt(SDTimeCosts / num_of_iter);

	cout << endl;
	cout << "Number of iteration:  Average, standard deviation: " << endl;
	cout << AveNumIters << ", " << SDNumIters << endl;
	cout << "Time cost:  Average, standard deviation: " << endl;
	cout << AveTimeCosts << ", " << SDTimeCosts << endl;


	for (int k = 0;k < NumIter.size();k++) {
		cout << "i" << NumRange[k + 1] << ": \t" << NumIter[k] << " \t" << 100 * (double)(NumIter[k]) / num_of_iter << "%" << endl;
	}

	cout << endl;
	cout << "Time cost: " << sum_time_cost / num_of_iter << endl;
	cout << endl;

	cout << endl;
	cout << "Number of failed: " << num_of_failed << endl;
	cout << endl;


	return 0;
}

// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

