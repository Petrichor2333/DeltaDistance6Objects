// DDistCollDetec.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include "pch.h"
#include "config.h"
#include "DDFunc.h"
#include "SolvingMinimum.h"
#include "SolvingMinimumMove.h"
#include "GenerateRandom.h"

int main() {

	srand((unsigned)time(NULL));

	double interval_boundary = 5;
	double Velocity_boundary = 3;
	int num_of_iter = 100000;	

	//double time_cost;
	//double sum_time_cost = 0;
	//double sum_time_cost_move = 0;
	
	

	int *NumIters = new int[num_of_iter];
	double *TimeCosts = new double[num_of_iter];

	int *NumItersMove = new int[num_of_iter];
	double *TimeCostsMove = new double[num_of_iter];

	/*
	vector<int> NumRange;
	NumRange.push_back(0);
	NumRange.push_back(5);
	NumRange.push_back(10);
	NumRange.push_back(15);
	NumRange.push_back(20);
	NumRange.push_back(30);
	NumRange.push_back(50
	NumRange.push_back(55);*/

	//vector<int> NumIter(NumRange.size() - 1, 0);

	ObjElData ObjElParaData1, ObjElParaData2;
	Obj1MSData Obj1MSParaData1, Obj1MSParaData2;
	Obj2PCData Obj2PCParaData1, Obj2PCParaData2;
	Obj3DMData Obj3DMParaData1, Obj3DMParaData2;
	Obj4CKData Obj4CKParaData1, Obj4CKParaData2;
	Obj5PLData Obj5PLParaData1, Obj5PLParaData2;
	Obj6HNData Obj6HNParaData1, Obj6HNParaData2;

	//Select 1 objects two be tested	
	ObjectDDFuncInfo objA0 = { &ObjElParaData1, DDFuncSqGrad_Obj_EL ,DDFuncSqHess_Obj_EL ,DDFuncSqValu_Obj_EL ,DDFuncSqValuIndep_Obj_EL };
	ObjectDDFuncInfo objA1 = { &Obj1MSParaData1, DDFuncSqGrad_Obj1_MS ,DDFuncSqHess_Obj1_MS ,DDFuncSqValu_Obj1_MS ,DDFuncSqValuIndep_Obj1_MS };
	ObjectDDFuncInfo objA2 = { &Obj2PCParaData1, DDFuncSqGrad_Obj2_PC ,DDFuncSqHess_Obj2_PC ,DDFuncSqValu_Obj2_PC ,DDFuncSqValuIndep_Obj2_PC };	
	ObjectDDFuncInfo objA3 = { &Obj3DMParaData1, DDFuncSqGrad_Obj3_DM ,DDFuncSqHess_Obj3_DM ,DDFuncSqValu_Obj3_DM ,DDFuncSqValuIndep_Obj3_DM };	
	ObjectDDFuncInfo objA4 = { &Obj4CKParaData1, DDFuncSqGrad_Obj4_CK ,DDFuncSqHess_Obj4_CK ,DDFuncSqValu_Obj4_CK ,DDFuncSqValuIndep_Obj4_CK };	
	ObjectDDFuncInfo objA5 = { &Obj5PLParaData1, DDFuncSqGrad_Obj5_PL ,DDFuncSqHess_Obj5_PL ,DDFuncSqValu_Obj5_PL ,DDFuncSqValuIndep_Obj5_PL };	
	ObjectDDFuncInfo objA6 = { &Obj6HNParaData1, DDFuncSqGrad_Obj6_HN ,DDFuncSqHess_Obj6_HN ,DDFuncSqValu_Obj6_HN ,DDFuncSqValuIndep_Obj6_HN };

	//Select 2 objects two be tested
	ObjectDDFuncInfo objB0 = { &ObjElParaData2, DDFuncSqGrad_Obj_EL ,DDFuncSqHess_Obj_EL ,DDFuncSqValu_Obj_EL ,DDFuncSqValuIndep_Obj_EL };
	ObjectDDFuncInfo objB1 = { &Obj1MSParaData2, DDFuncSqGrad_Obj1_MS ,DDFuncSqHess_Obj1_MS ,DDFuncSqValu_Obj1_MS ,DDFuncSqValuIndep_Obj1_MS };
	ObjectDDFuncInfo objB2 = { &Obj2PCParaData2, DDFuncSqGrad_Obj2_PC ,DDFuncSqHess_Obj2_PC ,DDFuncSqValu_Obj2_PC ,DDFuncSqValuIndep_Obj2_PC };
	ObjectDDFuncInfo objB3 = { &Obj3DMParaData2, DDFuncSqGrad_Obj3_DM ,DDFuncSqHess_Obj3_DM ,DDFuncSqValu_Obj3_DM ,DDFuncSqValuIndep_Obj3_DM };
	ObjectDDFuncInfo objB4 = { &Obj4CKParaData2, DDFuncSqGrad_Obj4_CK ,DDFuncSqHess_Obj4_CK ,DDFuncSqValu_Obj4_CK ,DDFuncSqValuIndep_Obj4_CK };
	ObjectDDFuncInfo objB5 = { &Obj5PLParaData2, DDFuncSqGrad_Obj5_PL ,DDFuncSqHess_Obj5_PL ,DDFuncSqValu_Obj5_PL ,DDFuncSqValuIndep_Obj5_PL };
	ObjectDDFuncInfo objB6 = { &Obj6HNParaData2, DDFuncSqGrad_Obj6_HN ,DDFuncSqHess_Obj6_HN ,DDFuncSqValu_Obj6_HN ,DDFuncSqValuIndep_Obj6_HN };

	ObjectDDFuncInfo* pobjA[7];
	ObjectDDFuncInfo* pobjB[7];
	pobjA[0] = &objA0, pobjA[1] = &objA1, pobjA[2] = &objA2, pobjA[3] = &objA3, pobjA[4] = &objA4, pobjA[5] = &objA5, pobjA[6] = &objA6;
	pobjB[0] = &objB0, pobjB[1] = &objB1, pobjB[2] = &objB2, pobjB[3] = &objB3,	pobjB[4] = &objB4, pobjB[5] = &objB5, pobjB[6] = &objB6;

	for (int obj_i = 0;obj_i < 7;obj_i++) {
		for (int obj_j = obj_i;obj_j < 7;obj_j++) {

			if (obj_i == 0 && obj_j != 0) {
				continue;
			}

			int num_of_failed = 0;
			int num_of_failed_move = 0;

			for (int i = 0;i < num_of_iter;i++) {

				double R1[3][3], R2[3][3];
				double p1[3], p2[3];
				double R1m[3][3], R2m[3][3];
				double p1m[3], p2m[3];
				double v1[3], v2[3];

				GenerateRandomRotTra(R1, p1, interval_boundary);
				GenerateRandomRotTra(R2, p2, interval_boundary);
				GenerateRandomRotTra(R1m, p1m, interval_boundary);
				GenerateRandomRotTra(R2m, p2m, interval_boundary);
				GenerateRandomVel(v1, Velocity_boundary);
				GenerateRandomVel(v2, Velocity_boundary);

				ResultsOfSolvingMinimum Result;				//Collision detection for stationary objects
				ResultsOfSolvingMinimumMove ResultMove;		//Collision detection for moving objects

				auto start = chrono::system_clock::now();
				bool SolSucc = SolveForCollision(pobjA[obj_i], pobjB[obj_j], &Result, R1, p1, R2, p2);					//Collision detection for stationary objects
				auto end = chrono::system_clock::now();
				auto duration = chrono::duration_cast<chrono::microseconds>(end - start);
				double time_cost = double(duration.count());
				NumIters[i] = Result.num_iter;
				TimeCosts[i] = time_cost;
				//sum_time_cost = sum_time_cost + time_cost;

				auto startMove = chrono::system_clock::now();
				bool SolSuccMove = SolveForCollisionMove(pobjA[obj_i], pobjB[obj_j], &ResultMove, R1m, p1m, R2m, p2m, v1, v2);		//Collision detection for moving objects
				auto endMove = chrono::system_clock::now();
				auto durationMove = chrono::duration_cast<chrono::microseconds>(endMove - startMove);
				double time_cost_move = double(durationMove.count());
				//sum_time_cost_move = sum_time_cost_move + time_cost;
				NumItersMove[i] = ResultMove.num_iter;
				TimeCostsMove[i] = time_cost_move;

				/*
				for (int k = 0;k < NumIter.size();k++) {
					if (Result.num_iter <= NumRange[k + 1]) {
						NumIter[k]++;
					}
				}*/

				if (Result.if_success == 0) {
					num_of_failed++;
				}
				if (ResultMove.if_success == 0) {
					num_of_failed_move++;
				}

				/*
				if (i % 10000 == 0 && i > 0) {
					cout << "iteration: " << i << endl;
				}
				*/
			}

			double SumNumIters = 0, SumTimeCosts = 0;
			double SumNumItersMove = 0, SumTimeCostsMove = 0;
			for (int i = 0;i < num_of_iter;i++) {
				SumNumIters += (double)(NumIters[i]);
				SumTimeCosts += TimeCosts[i];
				SumNumItersMove += (double)(NumItersMove[i]);
				SumTimeCostsMove += TimeCostsMove[i];
			}
			double AveNumIters = SumNumIters / num_of_iter;
			double AveTimeCosts = SumTimeCosts / num_of_iter;
			double AveNumItersMove = SumNumItersMove / num_of_iter;
			double AveTimeCostsMove = SumTimeCostsMove / num_of_iter;

			double SDNumIters = 0, SDTimeCosts = 0;
			double SDNumItersMove = 0, SDTimeCostsMove = 0;
			for (int i = 0;i < num_of_iter;i++) {
				SDNumIters += ((double)(NumIters[i]) - AveNumIters)*((double)(NumIters[i]) - AveNumIters);
				SDTimeCosts += (TimeCosts[i] - AveTimeCosts)*(TimeCosts[i] - AveTimeCosts);
				SDNumItersMove += ((double)(NumItersMove[i]) - AveNumItersMove)*((double)(NumItersMove[i]) - AveNumItersMove);
				SDTimeCostsMove += (TimeCostsMove[i] - AveTimeCostsMove)*(TimeCostsMove[i] - AveTimeCostsMove);
			}
			SDNumIters = sqrt(SDNumIters / num_of_iter);
			SDTimeCosts = sqrt(SDTimeCosts / num_of_iter);
			SDNumItersMove = sqrt(SDNumItersMove / num_of_iter);
			SDTimeCostsMove = sqrt(SDTimeCostsMove / num_of_iter);

			cout << "Runing tests for objects " << obj_i << " and " << obj_j << " ..." << endl;
			if (obj_i == 0) {
				cout << "Object 0 is an ellipse with three semi axes of 1.5, 1.2, and 0.7, respectively." << endl;
			}

			cout << endl;
			cout << "Stationary cases:" << endl;
			cout << "Number of iteration:  Average, standard deviation: " << endl;
			cout << AveNumIters << " +/- " << SDNumIters << endl;
			cout << "Time cost:  Average, standard deviation: " << endl;
			cout << AveTimeCosts << " +/- " << SDTimeCosts << endl;
			cout << "Number of failed: " << num_of_failed << endl;
			cout << endl << endl;

			cout << "Moving cases:" << endl;
			cout << "Number of iteration:  Average, standard deviation: " << endl;
			cout << AveNumItersMove << " +/- " << SDNumItersMove << endl;
			cout << "Time cost:  Average, standard deviation: " << endl;
			cout << AveTimeCostsMove << " +/- " << SDTimeCostsMove << endl;
			cout << "Number of failed: " << num_of_failed_move << endl;
			
			cout << endl;
			cout <<"-----------------------------------------------------------"<< endl;
			cout << endl;

			/*
			for (int k = 0;k < NumIter.size();k++) {
				cout << "i" << NumRange[k + 1] << ": \t" << NumIter[k] << " \t" << 100 * (double)(NumIter[k]) / num_of_iter << "%" << endl;
			}

			std::cout << endl;
			cout << "Time cost: " << sum_time_cost / num_of_iter << endl;
			cout << endl;
			cout << endl;
			*/

		}
	}

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
